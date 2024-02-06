"""
Download JRC Global Human Settlement Layer population

Reference
---------
https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

Dataset:

> Schiavina, Marcello; Freire, Sergio; MacManus, Kytt (2022): GHS-POP R2022A - GHS population
> grid multitemporal (1975-2030). European Commission, Joint Research Centre (JRC) [Dataset]
> DOI: 10.2905/D6D86A90-4351-4508-99C1-CB074B022C4A PID:
> http://data.europa.eu/89h/d6d86a90-4351-4508-99c1-cb074b022c4a

Concept & Methodology:

> Freire, Sergio; MacManus, Kytt; Pesaresi, Martino; Doxsey-Whitfield, Erin; Mills, Jane
> (2016): Development of new open and free multi-temporal global population grids at 250 m
> resolution. Geospatial Data in a Changing World; Association of Geographic Information
> Laboratories in Europe (AGILE). AGILE 2016.

"""

rule download_ghsl:
    output:
        "{OUTPUT_DIR}/input/ghsl/GHS_POP_E{YEAR}_GLOBE_{RELEASE}_54009_{RESOLUTION}_V1_0.tif"
    wildcard_constraints:
        YEAR=range(1975, 2031, 5),
        RESOLUTION="100|1000"
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_{wildcards.RELEASE}/GHS_POP_E{wildcards.YEAR}_GLOBE_{wildcards.RELEASE}_54009_{wildcards.RESOLUTION}/V1-0/GHS_POP_E{wildcards.YEAR}_GLOBE_{wildcards.RELEASE}_54009_{wildcards.RESOLUTION}_V1_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_POP_E{wildcards.YEAR}_GLOBE_{wildcards.RELEASE}_54009_{wildcards.RESOLUTION}_V1_0.zip \
            -d $output_dir
        """

rule download_ghsl_all:
    input:
        expand(
            os.path.join(
                "{{OUTPUT_DIR}}",
                "input",
                "ghsl",
                "GHS_POP_E{year}_GLOBE_{release}_54009_{resolution}_V1_0.tif",
            ),
            resolution=(100, 1000),
            year=(2020, ),
            release="R2022A" # TODO bump to R2023A
        )
