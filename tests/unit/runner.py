"""
Test runner and file comparison code for snakemake unit tests.
"""

import json
import os
import re
import shutil
import sys
import subprocess as sp
from pprint import pformat
from pathlib import Path
from tempfile import TemporaryDirectory

import geopandas as gpd
import pandas as pd


def printerr(s: str):
    """Print to stderr for visibility when invoked via pytest."""
    print(s, file=sys.stderr)


def run_snakemake_test(target_name: str, output: tuple[str]):
    """
    Create a temporary working directory, copy input files to it, run a
    snakemake rule and compare the outputs of that rule with some expected
    output.

    Args:
        target_name (str): Directory name holding reference test IO
        output (tuple[str]): Desired output files and directories for
            snakemake to generate
    """

    if not isinstance(output, tuple):
        raise TypeError(f"Expect desired outputs as tuple, got {type(output)=}")

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = Path(f"tests/unit/{target_name}/data")
        expected_path = Path(f"tests/unit/{target_name}/expected")

        # check necessary testing data is present
        assert data_path.exists()
        assert expected_path.exists()

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        auxilliary_dirs = ["config", "external_files", "bundled_data"]
        for folder in auxilliary_dirs:
            shutil.copytree(f"tests/{folder}", f"{workdir}/{folder}")

        # dbg
        printerr(target_name)

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                *output,
                "-r",  # show reasons, helps with debugging
                "--configfile",
                "tests/config/config.yaml",
                "-j1",  # single core
                "--keep-target-files",
                "--directory",
                workdir,
            ]
        )

        OutputChecker(data_path, expected_path, workdir, auxilliary_dirs).check()


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir, ignore_folders):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir
        self.ignore_folders = ignore_folders

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
            if any([folder for folder in self.ignore_folders if folder in path]):
                # skip comparison for these folders
                continue
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if ".snakemake" in str(f):
                    # ignore .snakemake/, foo/bar/.snakemake_timestamp, etc.
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
                "Unexpected files: {}".format(sorted(map(str, unexpected_files)))
            )

    def compare_files(self, generated_file: Path, expected_file: Path) -> None:
        """
        Compare two files to check if they are equal by some definition.

        Methods vary by filetype.
        """

        printerr(f">>> Compare files:\n{generated_file}\n{expected_file}")

        # PARQUET
        if re.search(r"\.(geo)?parquet$", str(generated_file), re.IGNORECASE):
            if re.search(r"\.geoparquet$", str(generated_file), re.IGNORECASE):
                """
                NOTE: This test will **fail** if the geoparquet file does not contain geography data columns.
                This can happen where the convert_to_geoparquet job does not find any roads to write.
                We leave this failure in because it is usually unintentional that you're testing with a
                dataset where _none_ of the slices have road data, and these tests should be targeted at
                slices that _do_ have road data.
                """
                read = gpd.read_parquet
            elif re.search(r"\.parquet$", str(generated_file), re.IGNORECASE):
                read = pd.read_parquet
            else:
                raise RuntimeError(f"couldn't identify read function for {generated_file}")

            generated = read(generated_file)
            expected = read(expected_file)

            self.compare_dataframes(generated, expected)

        # JSON
        elif re.search(r"\.(geo)?json$", str(generated_file), re.IGNORECASE):
            with open(generated_file, 'r') as fp:
                generated = json.load(fp)
            with open(expected_file, 'r') as fp:
                expected = json.load(fp)

            if json.dumps(generated, sort_keys=True) != json.dumps(expected, sort_keys=True):
                printerr(">>> Method: compare sorted JSON strings")
                printerr(f">>> generated:\n{pformat(generated)}")
                printerr(f">>> expected:\n{pformat(expected)}")
                raise AssertionError("JSON files do not match")

        # JPG, PDF, PNG & SVG images
        elif re.search(r"\.(jpg|jpeg|pdf|png|svg|tif|tiff)$", str(generated_file), re.IGNORECASE):
            try:
                sp.check_output(["tests/unit/visual_compare.sh", generated_file, expected_file])
            except sp.CalledProcessError as e:
                printerr(">>> Method: visual hash comparison (imagemagick's identify)")
                printerr(f">>> ERROR:\n>>> {e.stdout}")
                raise e

        # any other file type
        else:
            try:
                sp.check_output(["cmp", generated_file, expected_file])
            except sp.CalledProcessError as e:
                printerr(">>> Method: binary comparison (cmp)")
                printerr(f">>> ERROR:\n>>> {e.stdout}")
                raise e

        printerr(">>> Files are a match")

    @staticmethod
    def compare_dataframes(generated: pd.DataFrame, expected: pd.DataFrame) -> None:
        """
        Compare two dataframes, raise ValueError if they aren't the same.
        """
        # after sorting the columns so they're in the same order,
        # use dataframe.equals to quickly check for complete table equality
        # unfortunately there is an edge case this method doesn't catch...
        if not generated.sort_index(axis="columns").equals(expected.sort_index(axis="columns")):
            printerr(">>> Method: compare (geo)pandas dataframes")

            # do some basic shape and schema checks
            if len(generated) != len(expected):
                raise ValueError(
                    f"tables not of same length, {len(generated)=} & {len(expected)=}"
                )
            if difference := set(generated.columns) ^ set(expected.columns):
                raise ValueError(
                    f"tables do not have same schema: {difference=} cols are in one but not both"
                )

            # there is a case where df.equals(identical_df) can return False despite all elements being equal
            # this is when comparing Nones in the same position: https://github.com/pandas-dev/pandas/issues/20442
            mismatch_cols = set()
            for col in generated.columns:
                if any(generated[col] != expected[col]):
                    mismatch_cols.add(col)

            for col in mismatch_cols:

                # do the discrepancies occur only where there are null values (NaN & None)?
                unequal_only_where_null = all(expected[col].isna() == (expected[col].values != generated[col].values))
                if not unequal_only_where_null:
                    printerr(f"{col=} {unequal_only_where_null=}")

                    # let's try and find failing rows by converting to str
                    MAX_FAILURES_TO_PRINT = 5
                    failures = 0
                    for row in range(len(generated)):
                        gen_str = str(generated[col][row: row + 1].values)
                        exp_str = str(expected[col][row: row + 1].values)
                        if gen_str != exp_str:
                            failures += 1
                            if failures > MAX_FAILURES_TO_PRINT:
                                continue
                            else:
                                printerr(f">>> FAILURE at {col=}, {row=}: {gen_str} != {exp_str}")

                    if failures > 0:
                        raise ValueError(f"{failures} row mismatch(es) between tables")

                else:
                    # None != None according to pandas, and this is responsible for the apparent mismatch
                    # we can safely say that this column is equal
                    continue
