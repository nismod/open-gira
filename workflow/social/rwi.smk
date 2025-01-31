rule download_rwi:
    """
    Meta Relative Wealth Index

    Description:
    https://dataforgood.facebook.com/dfg/tools/relative-wealth-index

    Catalogue page: https://data.humdata.org/dataset/relative-wealth-index


    Citation:

    Microestimates of wealth for all low- and middle-income countries Guanghua
    Chi, Han Fang, Sourav Chatterjee, Joshua E. Blumenstock Proceedings of the
    National Academy of Sciences Jan 2022, 119 (3) e2113658119; DOI:
    10.1073/pnas.2113658119

    License: CC-BY-NC 4.0
    """
    output:
        archive = "{OUTPUT_DIR}/input/rwi/country_data/relative-wealth-index-april-2021.zip",
        csv = "{OUTPUT_DIR}/input/rwi/country_data/relative-wealth-index-april-2021.csv",
    shell:
        """
        output_dir=$(dirname {output.csv})
        mkdir -p $output_dir

        # Combined April 2021 release
        wget -nc \
            https://data.humdata.org/dataset/76f2a2ea-ba50-40f5-b79c-db95d668b843/resource/de2f953e-940c-43bb-b1f8-4d02d28124b5/download/relative-wealth-index-april-2021.zip \
            --directory-prefix=$output_dir

        # Turkey released later
        wget -nc \
            https://data.humdata.org/dataset/76f2a2ea-ba50-40f5-b79c-db95d668b843/resource/ca1dea6b-9cc1-4bed-9465-8fd8fd7e0941/download/tur_relative_wealth_index.csv \
            --directory-prefix=$output_dir

        unzip -j -n {output.archive} -d $output_dir

        pushd $output_dir
            # fix mismatched columns
            cat tur_relative_wealth_index.csv | cut -d ',' -f 2,3,4,5 > tur_relative_wealth_index.noquadkey.csv
            mv tur_relative_wealth_index.noquadkey.csv tur_relative_wealth_index.csv
            # concatenate all
            awk '(NR == 1) || (FNR > 1)' *_relative_wealth_index.csv > $(basename {output.csv})
        popd
        """

rule process_rwi_points:
    input:
        csv = "{OUTPUT_DIR}/input/rwi/country_data/relative-wealth-index-april-2021.csv",
    output:
        gpkg = "{OUTPUT_DIR}/input/rwi/relative-wealth-index-april-2021.gpkg",
    run:
        import pandas as pd
        import geopandas as gpd
        df = pd.read_csv(input.csv)
        gdf = gpd.GeoDataFrame(
            data=df[["rwi", "error"]],
            geometry=gpd.points_from_xy(df.longitude, df.latitude),
            crs="EPSG:4326"
        )
        # reproject to 3857 (point coordinates are provided in lat/lon but the
        # underlying regular grid seems to be based on Web Mercator via Quadkeys)
        gdf.to_crs("EPSG:3857", inplace=True)
        gdf.to_file(output.gpkg, engine="pyogrio")

rule process_rwi_grid:
    """Rasterise points to a grid data format

    -tr specifies x/y resolution, guessed from most common difference between point
    coordinates:

    def guess_resolution(coords):
        coords.sort()
        diffs = np.diff(coords)
        diffs.sort()
        # vals, counts = np.unique(diffs, return_counts=True)
        counts, vals = np.histogram(diffs)
        return vals[np.argmax(counts)]
    print("x res:", guess_resolution(gdf.geometry.x.unique())
    print("y res:", guess_resolution(gdf.geometry.y.unique())
    """
    input:
        gpkg = "{OUTPUT_DIR}/input/rwi/relative-wealth-index-april-2021.gpkg",
    output:
        tiff = "{OUTPUT_DIR}/input/rwi/rwi.tif",
    shell:
        """
        gdal_rasterize \
            -a rwi \
            -init -999 \
            -a_nodata -999 \
            -tr 2445.9786434 2445.96770335 \
            -ot Float64 \
            -of GTiff \
            {input.gpkg} \
            {output.tiff}
        """



