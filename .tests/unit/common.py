"""
Common code for unit testing of rules generated with Snakemake 6.8.1.
"""

from pathlib import Path
import subprocess as sp
import os

import geopandas
import pandas
from pandas.testing import assert_frame_equal

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
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if (
                    str(f).startswith(".snakemake")
                    or str(f).startswith("data/aqueduct")
                    or str(f) == "config.yaml"
                ):
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
        if ".geoparquet" == generated_file.suffix:
            actual = geopandas.read_parquet(generated_file)
            expected = geopandas.read_parquet(expected_file)
            assert_frame_equal(actual, expected)
        #elif ".parquet" == generated_file.suffix:
        #    actual = pandas.read_parquet(generated_file)
        #    expected = pandas.read_parquet(expected_file)
        #    assert_frame_equal(actual, expected)
        else:
            sp.check_output(["cp", generated_file, "/tmp/"])
            sp.check_output(["cmp", generated_file, expected_file])
