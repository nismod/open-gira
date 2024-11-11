"""Download Harwood et al. (2024) NbS opportunities
"""

rule download_nbs_suitability:
    """
    Fetch our preprocessed data from zenodo.
    """
    output:
        treelim_current="{OUTPUT_DIR}/input/nbs-suitability/raw/treelim_current.tif",
        treelim_2050s="{OUTPUT_DIR}/input/nbs-suitability/raw/treelim_2050s.tif",
        ls_nbs_current="{OUTPUT_DIR}/input/nbs-suitability/raw/ls_nbs_current.tif",
        ls_nbs_2050="{OUTPUT_DIR}/input/nbs-suitability/raw/ls_nbs_2050.tif",
        ls_range_current="{OUTPUT_DIR}/input/nbs-suitability/raw/ls_range_current.tif",
        ls_range_2050="{OUTPUT_DIR}/input/nbs-suitability/raw/ls_range_2050.tif",
    shell:
        """
        pushd $(dirname {output.zip})
            zenodo_get -w links.txt --record=XXXXXX # TODO
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """

rule extract_nbs_suitability:
    output:
