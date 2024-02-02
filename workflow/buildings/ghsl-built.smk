"""
Download JRC Global Human Settlement Layer Built-up Surface

GHS-BUILT-S R2023A - GHS built-up surface grid, derived from Sentinel2 composite
and Landsat, multitemporal (1975-2030)

Reference
---------
https://ghsl.jrc.ec.europa.eu/download.php?ds=bu

Dataset:

> Pesaresi M., Politis P. (2023): GHS-BUILT-S R2023A - GHS built-up surface
> grid, derived from Sentinel2 composite and Landsat, multitemporal
> (1975-2030)European Commission, Joint Research Centre (JRC) PID:
> http://data.europa.eu/89h/9f06f36f-4b11-47ec-abb0-4f8b7b1d72ea,
> doi:10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA

Concept & Methodology:

> European Commission GHSL Data Package 2023, Publications Office of the
> European Union, Luxembourg, 2023, JRC133256, ISBN 978-92-68-02341-9
> doi:10.2760/098587

"""

rule download_ghsl_built_s:
    output:
        "{OUTPUT_DIR}/input/ghsl/GHS_{RES_NRES}_E{YEAR}_GLOBE_R2023A_4326_3ss_V1_0.tif"
    wildcard_constraints:
        YEAR=range(1975, 2031, 5),
        RES_NRES="BUILT_S|BUILT_S_NRES"
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_{wildcards.RES_NRES}_GLOBE_R2023A/GHS_{wildcards.RES_NRES}_E{wildcards.YEAR}_GLOBE_R2023A_4326_3ss/V1-0/GHS_{wildcards.RES_NRES}_E{wildcards.YEAR}_GLOBE_R2023A_4326_3ss_V1_0.zip
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_{wildcards.RES_NRES}_E{wildcards.YEAR}_GLOBE_R2023A_4326_3ss_V1_0.zip \
            -d $output_dir
        """

rule download_ghsl_built_s_all:
    input:
        expand(
            os.path.join(
                "{{OUTPUT_DIR}}",
                "input",
                "ghsl",
                "GHS_{RES_NRES}_E{YEAR}_GLOBE_R2023A_4326_3ss_V1_0.tif",
            ),
            RES_NRES=("BUILT_S", "BUILT_S_NRES"),
            YEAR=(2020, )
        )
