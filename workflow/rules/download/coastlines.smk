"""Download coastline data

Reference
---------
https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip
"""


rule download_coastlines:
    output:
        directory("{OUTPUT_DIR}/input/coastlines/ne_10m_ocean/"),
    shell:
        """
        wget https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip \
            --output-document {output}.zip
        unzip -o {output}.zip -d {output}
        rm -f {output}.zip
        """


"""
Test with:
snakemake --cores 1 results/input/coastlines/ne_10m_ocean/
"""
