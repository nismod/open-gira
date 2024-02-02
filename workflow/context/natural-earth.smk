"""
Rules for getting naturalearthdata.com datasets
"""

# Not very DRY, but okay for now
rule download_ne_50m:
    output:
        "{OUTPUT_DIR}/input/admin-boundaries/zip/ne_50m.zip",
    shell:
        """
        wget -O {output} https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/ne_50m_admin_0_countries.zip
        """


"""
Test with:
snakemake --cores 1 results/input/admin-boundaries/zip/ne_50.zip
"""


rule unzip_ne_50m:
    input:
        "{OUTPUT_DIR}/input/admin-boundaries/zip/ne_50m.zip",
    output:
        directory("{OUTPUT_DIR}/input/admin-boundaries/ne_50m/"),
    shell:
        """
        unzip -o {input} -d {output}
        """


"""
Test with:
snakemake --cores 1 results/input/admin-boundaries/ne_50m/
"""
