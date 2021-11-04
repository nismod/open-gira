import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_slice():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/slice/data")
        expected_path = PurePosixPath(".tests/unit/slice/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(".tests/config.yaml", workdir)

        # dbg
        print("data/northeast-oxford-slice0.osm.pbf data/northeast-oxford-slice1.osm.pbf data/northeast-oxford-slice2.osm.pbf data/northeast-oxford-slice3.osm.pbf", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "data/northeast-oxford-slice0.osm.pbf data/northeast-oxford-slice1.osm.pbf data/northeast-oxford-slice2.osm.pbf data/northeast-oxford-slice3.osm.pbf",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
