rule download_IBTrACS:
    output:
        "{OUTPUT_DIR}/input/IBTrACS/raw/v4.csv"
    shell:
        """
        wget --output-document {output} \
            https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.ALL.list.v04r00.csv
        """


"""
Test with:
snakemake -c1 results/input/IBTrACS/v4.csv
"""
