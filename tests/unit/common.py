"""
Common code for unit testing of rules generated with Snakemake 6.15.1.
"""

import re
import shutil
import sys
from pathlib import Path, PurePosixPath
import subprocess as sp
import os
from tempfile import TemporaryDirectory

import geopandas
import pandas


def run_test(target_name, command):
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(f"tests/unit/{target_name}/data")
        expected_path = PurePosixPath(f"tests/unit/{target_name}/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree('tests/config', f"{workdir}/config")
        shutil.copytree('tests/external_files', f"{workdir}/external_files")

        # dbg
        print(target_name, file=sys.stderr)

        if isinstance(command, str):
            command = command.split(" ")

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            *command,

            "-r",  # show reasons, helps with debugging
            "--configfile",
            "tests/config/config.yaml",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        OutputChecker(data_path, expected_path, workdir).check()


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            if "config" in path or "external_files" in path:
                continue
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if str(f).startswith(".snakemake"):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):
        print(f">>> Compare {generated_file} vs {expected_file}", file=sys.stderr)
        if re.search("\\.geoparquet$", str(generated_file), re.IGNORECASE):
            """
            NOTE: This test will **fail** if the geoparquet file does not contain geography data columns.
            This can happen where the convert_to_geoparquet job does not find any roads to write.
            We leave this failure in because it is usually unintentional that you're testing with a
            dataset where _none_ of the slices have road data, and these tests should be targeted at
            slices that _do_ have road data.
            """
            print(f">>> Method=geopandas.GeoDataFrame.compare", file=sys.stderr)
            read = geopandas.read_parquet
        elif re.search("\\.parquet$", str(generated_file), re.IGNORECASE):
            print(f">>> Method=pandas.GeoDataFrame.compare", file=sys.stderr)
            read = pandas.read_parquet
        if re.search("\\.(geo)?parquet$", str(generated_file), re.IGNORECASE):
            generated = read(generated_file)
            expected = read(expected_file)
            # use dataframe.equals to quickly check for complete table equality
            if generated.equals(expected):
                return
            else:
                # horribly slow but gives useful information on failures
                for r in range(len(generated)):
                    try:
                        assert str(generated[r:r+1]) == str(expected[r:r+1])
                    except AssertionError as e:
                        print(f">>> FAILURE at row {r}.")
                        print(f"{str(generated[r:r+1])} not equal to {str(expected[r:r+1])}")
                        raise e
                else:
                    raise RuntimeError("couldn't find mismatch between dataframes...")
        else:
            print(f">>> Method=cmp", file=sys.stderr)
            try:
                sp.check_output(["cmp", generated_file, expected_file])
            except sp.CalledProcessError as e:
                print(f">>> ERROR:\n>>> {e.stdout}", file=sys.stderr)
                raise e
        print(f">>> OK", file=sys.stderr)
