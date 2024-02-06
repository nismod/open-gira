"""
'Target' rules that describe final outputs.

N.B. For the --batch CLI argument to work with a rule, inputs to a rule must be
explicit (i.e. not require calculation of the DAG). Given n input files, there may be up to n batches.
"""


def STORM_exposure_by_country_by_climate_model(wildcards):
    return expand(
        "{OUTPUT_DIR}/power/by_storm_set/{STORM_SET}/exposure/admin-level-0.geoparquet",
        OUTPUT_DIR = wildcards.OUTPUT_DIR,
        STORM_SET = [
            "STORM-constant",
            "STORM-CMCC-CM2-VHR4",
            "STORM-CNRM-CM6-1-HR",
            "STORM-EC-Earth3P-HR",
            "STORM-HadGEM3-GC31-HM",
        ],
    )


rule STORM_exposure_by_storm_set:
    """
    STORM target rule, including present day ('constant') and all climate models.
    """
    input:
        STORM_exposure_by_country_by_climate_model
    output:
        "{OUTPUT_DIR}/power/by_storm_set/STORM_exposure.txt"
    shell:
        """
        # one output file per line
        echo {input} | tr ' ' '\n' > {output}
        """

"""
Test with:
snakemake --cores all --batch STORM_exposure_by_storm_set=1/5 -- results/power/by_storm_set/STORM_exposure.txt
"""
