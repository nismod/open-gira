import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_all():
    common.run_test('all', 'snakemake all -j1 --keep-target-files')
