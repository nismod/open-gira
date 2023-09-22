"""
Download the Imperial tropical cyclone storm tracks
"""

rule download_IRIS:
    """
    As of 20230626, this data is not publically available. You will need
    appropriate keys to access the files on the OUCE file store.
    """
    output:
        zip_file = "{OUTPUT_DIR}/input/IRIS/archive.zip"
    shell:
        """
        mkdir -p $(dirname {output.zip_file})
        scp /ouce-home/projects/mistral/iris/iris-data.zip {output.zip_file}
        """

"""
Test with:
snakemake -n -c1 -- results/input/IRIS/archive.zip
"""


rule extract_IRIS:
    input:
        zip_file = rules.download_IRIS.output.zip_file
    output:
        unzipped_dir = directory("{OUTPUT_DIR}/input/IRIS/iris-data/")
    shell:
        """
        unzip {input.zip_file} -d $(dirname {output.unzipped_dir})
        """
