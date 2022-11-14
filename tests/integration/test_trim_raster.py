from . import runner


def test_trim_raster():
    runner.run_snakemake_test(
        "trim_raster",
        ("results/input/hazard-aqueduct-river/djibouti-latest/inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp01000.tif",)
    )
