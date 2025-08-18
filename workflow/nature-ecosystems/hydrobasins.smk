"""Hydrobasins

https://www.hydrosheds.org/products/hydrobasins

HydroBASINS represents a series of vectorized polygon layers that depict
sub-basin boundaries at a global scale. The goal of this product is to provide a
seamless global coverage of consistently sized and hierarchically nested
sub-basins at different scales (from tens to millions of square kilometers),
supported by a coding scheme that allows for analysis of catchment topology such
as up- and downstream connectivity. HydroBASINS has been extracted from the
gridded HydroSHEDS core layers at 15 arc-second resolution.

## Technical Documentation

For more information on HydroBASINS please refer to the [HydroBASINS Technical
Documentation](https://data.hydrosheds.org/file/technical-documentation/HydroBASINS_TechDoc_v1c.pdf).

## License

The HydroBASINS database is freely available for scientific, educational and
commercial use. The data are distributed under the same license agreement as the
HydroSHEDS core products, which is included in the [HydroSHEDS Technical
Documentation](https://data.hydrosheds.org/file/technical-documentation/HydroSHEDS_TechDoc_v1_4.pdf).
For all regulations regarding license grants, copyright, redistribution
restrictions, required attributions, disclaimer of warranty, indemnification,
liability, and waiver of damages, please refer to the license agreement.

By downloading and using the data the user agrees to the terms and conditions of
this license.

See `hydrobasins_license.txt` for Appendix A of the technical documentation
linked above.

### REQUIRED ATTRIBUTIONS

The following copyright statement must be displayed with, attached to or
embodied in (in a reasonably prominent manner) the documentation or metadata of
any Licensee Product or Program provided to an End User when utilizing the
Licensed Materials:

This product [insert Licensee Derivative Product name] incorporates data from
the HydroSHEDS version 1 database which is © World Wildlife Fund, Inc.
(2006-2022) and has been used herein under license. WWF has not evaluated the
data as altered and incorporated within [insert Licensee Derivative Product
name], and therefore gives no warranty regarding its accuracy, completeness,
currency or suitability for any particular purpose. Portions of the HydroSHEDS
v1 database incorporate data which are the intellectual property rights of ©
USGS (2006-2008), NASA (2000-2005), ESRI (1992-1998), CIAT (2004-2006),
UNEP-WCMC (1993), WWF (2004), Commonwealth of Australia (2007), and Her Royal
Majesty and the British Crown and are used under license. The HydroSHEDS v1
database and more information are available at https://www.hydrosheds.org.

The scientific citation for the HydroSHEDS version 1 database is:

Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from
spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.

## References

Lehner, B., Grill G. (2013). Global river hydrography and network routing:
baseline data and new approaches to study the world's large river systems.
Hydrological Processes, 27(15): 2171-2186. DOI:
[10.1002/hyp.9740](https://doi.org/10.1002/hyp.9740)
"""


rule download_hydrobasins:
    output:
        zip="{OUTPUT_DIR}/input/hydrobasins/raw/hybas_{HYDROBASINS_REGION}_lev01-12_v1c.zip",
    run:
        import os
        import requests
        fname = output.zip
        url = f"https://data.hydrosheds.org/file/hydrobasins/standard/{os.path.basename(fname)}"
        r = requests.get(url)
        with open(fname, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

rule extract_hydrobasins:
    input:
        zip="{OUTPUT_DIR}/input/hydrobasins/raw/hybas_{HYDROBASINS_REGION}_lev01-12_v1c.zip"
    output:
        shp=expand(
            "{{OUTPUT_DIR}}/input/hydrobasins/raw/hybas_{{HYDROBASINS_REGION}}_lev{{HYDROBASINS_LEVEL}}_v1c.{EXT}",
            EXT=["dbf", "prj", "sbn", "sbx", "shp", "shp.xml", "shx"]
        )
    run:
        import zipfile
        with zipfile.ZipFile(input.zip, 'r') as zip_ref:
            for output_fname in output.shp:
                fname = os.path.basename(output_fname)
                dirname = os.path.dirname(output_fname)
                zip_ref.extract(fname, dirname)

rule extract_hydrobasins_level:
    input:
        shps=expand(
            "{{OUTPUT_DIR}}/input/hydrobasins/raw/hybas_{HYDROBASINS_REGION}_lev{{HYDROBASINS_LEVEL}}_v1c.shp",
            HYDROBASINS_REGION=["af", "ar", "as", "au", "eu", "gr", "na", "sa", "si"],
        )
    output:
        parquet="{OUTPUT_DIR}/input/hydrobasins/hybas_lev{HYDROBASINS_LEVEL}_v1c.geoparquet",
    run:
        import geopandas
        import pandas
        dfs = []
        for shp in input.shps:
            df = geopandas.read_file(shp)
            dfs.append(df)
        df = pandas.concat(dfs)
        df.to_parquet(output.parquet)

rule extract_hydrobasins_all:
    input:
        parquets=expand(
            "results/input/hydrobasins/hybas_lev{HYDROBASINS_LEVEL}_v1c.geoparquet",
            HYDROBASINS_LEVEL=[f"{n:02}" for n in range(1,13)]
        )

rule hydrobasins_add_admin_codes:
    input:
        hydrobasins="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c.geoparquet",
        adm0="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
        adm1="{OUTPUT_DIR}/input/admin-boundaries/admin-level-1.geoparquet",
        adm2="{OUTPUT_DIR}/input/admin-boundaries/admin-level-2.geoparquet",
    output:
        hydrobasins_adm="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
    run:
        import geopandas
        # Read inputs
        hydrobasins = geopandas.read_parquet(input.hydrobasins, columns=["HYBAS_ID", "geometry"])
        adm0 = geopandas.read_parquet(input.adm0, columns=['GID_0', 'NAME_0', 'geometry'])
        adm1 = geopandas.read_parquet(input.adm1, columns=['GID_1', 'NAME_1', 'geometry'])
        adm2 = geopandas.read_parquet(input.adm2, columns=['GID_2', 'NAME_2', 'geometry'])
        # Use centroids for join
        hydrobasins.geometry = hydrobasins.geometry.centroid
        hb_adm = (
            hydrobasins
            .sjoin(adm0, how="left", predicate="within", lsuffix="", rsuffix="0")
            .sjoin(adm1, how="left", predicate="within", lsuffix="", rsuffix="1")
            .sjoin(adm2, how="left", predicate="within", lsuffix="", rsuffix="2")
            .drop(columns=["geometry", "index_0", "index_1"])
            .drop_duplicates(subset="HYBAS_ID", keep="first")
            .set_index("HYBAS_ID")
        )
        # Set geoms back, save
        hb_geoms = geopandas.read_parquet(input.hydrobasins, columns=["HYBAS_ID", "geometry"]) \
            .set_index("HYBAS_ID")
        hb_adm = hb_geoms.join(hb_adm)
        hb_adm.to_parquet(output.hydrobasins_adm)

rule hydrobasins_add_population:
    input:
        hydrobasins="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
        population="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif"
    output:
        hydrobasins="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes_pop.geoparquet",
        csv="{OUTPUT_DIR}/input/hydrobasins/hybas_lev12_v1c_with_gadm_codes_pop.csv",
    run:
        import geopandas
        from rasterstats import gen_zonal_stats
        from tqdm import tqdm

        hydrobasins = geopandas.read_parquet(input.hydrobasins, columns=["HYBAS_ID", "geometry"]).to_crs("+proj=moll")
        zone_pops = [
            p['sum']
            for p in gen_zonal_stats(
                hydrobasins.geometry,
                input.population,
                stats=['sum']
            )
        ]
        del hydrobasins
        hydrobasins = geopandas.read_parquet(input.hydrobasins)
        hydrobasins['population'] = zone_pops
        hydrobasins.to_parquet(output.hydrobasins)

        hydrobasins_table = hydrobasins[["GID_0", "population"]]
        hydrobasins_table.to_csv(output.csv)
