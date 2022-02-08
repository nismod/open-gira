import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_convert_to_geoparquet():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/unit/convert_to_geoparquet/data")
        expected_path = PurePosixPath("tests/unit/convert_to_geoparquet/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/geoparquet/northeast-oxford_filter-highway-core_slice-0.geoparquet", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/geoparquet/northeast-oxford_filter-highway-core_slice-0.geoparquet",
            "-F", 
            "-j1",
            "--keep-target-files",
            "--configfile",
            "tests/config.yaml",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
