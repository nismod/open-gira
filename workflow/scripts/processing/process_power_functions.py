"""
functions and data required to perform preprocessing
"""


from importing_modules import *
import shapely.wkt as sw


def changedir():
    """For use on personal pc"""
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

def timer(s):
    print("timer : ",round((time.time() - s)/60, 2)," mins\n")

def get_powerplants(code):
    """returns all powerplants (processed data)"""
    code = code.upper()
    powerplant_file = os.path.join("data","powerplants", "global_power_plant_database.csv")
    powerplants = pd.read_csv(powerplant_file, dtype={'other_fuel3':object})  # dtype added to supress error

    powerplants = powerplants[powerplants.country == code].copy().reset_index(drop=True)

    if len(powerplants) == 0:  # no powerplants (eg vatican)
        powerplants['geometry'] = []
        return powerplants
    else:
        powerplants['geometry'] = powerplants.apply(
            lambda row: shape({'type': 'Point', 'coordinates':[row.longitude, row.latitude]}),
            axis=1
        )
        powerplants = gpd.GeoDataFrame(
            powerplants[['gppd_idnr', 'name', 'capacity_mw', 'estimated_generation_gwh_2017', 'primary_fuel', 'geometry']]
        )

        # can exclude powerplants here

        return powerplants


def get_target_areas(code):
    code = code.upper()
    # open .gpkg file and extract geometry of boundary data
    with fiona.open(os.path.join("data","adminboundaries", f"gadm36_{code}.gpkg"), "r") as src:
        shapes = [feature["geometry"] for feature in src]

    geod = Geod(ellps="WGS84")  # choice from https://gadm.org/download_country.html


    # Targets: Binary raster showing locations predicted to be connected to distribution grid.
    with rasterio.open(os.path.join("data","gridfinder","targets.tif")) as src:

        # for each country (see: code) overlay connections to distribution grid
        data, transform = rasterio.mask.mask(src, shapes, crop=True)

        geoms = []
        centroids = []
        areas = []
        # Extract feature shapes and values from the array.
        for geom, val in rasterio.features.shapes(data, transform=transform):
            if val > 0:
                feature = shape(geom)
                geoms.append(feature)
                centroids.append(feature.centroid)
                area, perimeter = geod.geometry_area_perimeter(feature)
                areas.append(area / 1e6)

    return gpd.GeoDataFrame({'area_km2':areas, 'centroid':centroids, 'geometry':geoms})


def get_population(code, targets):
    code = code.upper()

    # below: population of target areas (targets.geometry), populations is list corresponding to areas (CHECK, I think)

    if len(targets) == 0:  # if no targets
        populations = []
        population_density = []
        print("no targets for ",code)
    else:
        populations = [
            d['nansum']
            for d in zonal_stats(
                targets.geometry,
                os.path.join("data","population",f"{code}_ppp_2020_UNadj_constrained.tif"),
                stats=[],
                add_stats={'nansum': np.nansum},  # count NaN as zero for summation
                all_touched=True  # possible overestimate, but targets grid is narrower than pop
            )
        ]
        # below: find population density at centroid of target areas
        population_density = point_query(
            targets.centroid,
            os.path.join("data","population",f"{code}_ppp_2020_UNadj_constrained.tif")
        )
        #
    targets['population'] = populations
    targets['population_density_at_centroid'] = population_density
    def estimate_population_from_density(row):
        if row.population is numpy.ma.masked:
            return row.area_km2 * row.population_density_at_centroid
        else:
            return row.population
    targets['population'] = targets.apply(estimate_population_from_density, axis=1)
    return targets


def get_gdp(targets):


    # just pick centroid - GDP per capita doesn't vary at this fine granularity
    gdp_pc = []
    fn = os.path.join("data","GDP","GDP_per_capita_PPP_1990_2015_v2.nc")
    ds = nc4.Dataset(fn)
    gdp_pc_lst = []
    for ii in range(len(targets)):  # TODO vectorise
        lat_idx = round(np.interp(targets.centroid[ii].y,[-90,90],[0,2160]),0)  # convert latitude to GDP nc file
        lon_idx = round(np.interp(targets.centroid[ii].x,[-180,180],[0,4320]),0)  # convert longitude to GDP nc file

        #print(lat_idx,lon_idx)
        # Find value
        gdp_pc_lst.append(ds['GDP_per_capita_PPP'][-1,lat_idx,lon_idx])  # returns gdp of centroid of the target area (2015)

    gdp_pc_lst = [float(x) if numpy.ma.is_masked(x) == False else 0 for x in gdp_pc_lst]  # set masked to 0 (later removed)
    #print("sum", sum(gdp_pc_lst))
    targets['gdp_pc'] = gdp_pc_lst
    targets['gdp'] = targets.gdp_pc * targets.population
    return targets


def get_lines(code=None):
    """Read gridfinder lines"""

    s = time.time()

    features = []
    with fiona.open(os.path.join("data","gridfinder","grid.gpkg")) as src:



        for jj,feature in tqdm(enumerate(src), desc='grid.gpkg features', total=len(src)):
            # gridfinder GeoPackage stores an "fid" which GeoPandas ignores
            # and fiona reads as "id", not to feature['properties']
            # see https://github.com/geopandas/geopandas/issues/1035
            geom = shape(feature['geometry'])

            # if code != None:
            #     check = code_geoms_gpd.within(geom)
            #     if True in check:
            #         continue
            #     else:
            #         break



            features.append({
                'source_id': feature['id'],
                'source': feature['properties']['source'],
                'geometry': geom,
            })

    print("time for grid.gpkg processing: ",round((time.time() - s)/60, 2)," mins")

    gdf_world = gpd.GeoDataFrame(features)

    if code != None:  # preload
        with fiona.open(os.path.join("data","adminboundaries",f"gadm36_{code}.gpkg")) as src_code:
            code_geoms = []
            for feature in src_code:
                code_geoms.append(shape(feature['geometry']))
            code_geoms_gpd = gpd.GeoDataFrame({'geometry':code_geoms})
        gdf_world = gdf_world.overlay(code_geoms_gpd, how='intersection')

        print(gdf_world)
    return gdf_world


def patch_nearest_edge(point, edges):
    """Set up network

    Find nearest edge to a point
    """


    if type(point) == str:
        point = sw.loads(point)  # if point is found as string -> convert to Point(# #)
        print("changed to point")
    geom = point.buffer(1e-2)
    #print(point, " : ", type(point))




    matches_idx = edges.sindex.nearest(geom.bounds)
    nearest_geom = min(
        [edges.iloc[match_idx] for match_idx in matches_idx],
        key=lambda match: point.distance(match.geometry)
    )
    return nearest_geom

snkit.network.nearest_edge = patch_nearest_edge


def write_plants_targets(code, plantsfile, targetsfile):
    #plantsfile = f"../data/processed/{code.upper()}_plants.csv"
    #targetsfile = f"../data/processed/{code.upper()}_targets.csv"

    if os.path.exists(plantsfile) and os.path.exists(targetsfile):  # check if files already exist
        print(f"files already exist for {code}")
        pass
    else:  # do not exist
        # Power plants
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
            plants_country = get_powerplants(code) \
                .rename(columns={'gppd_idnr': 'source_id'})

            plants_country['type'] = 'source'


        # Targets
        #print("getting target areas")
        targets_country = get_target_areas(code)

        #print("getting target population")
        targets_country = get_population(code, targets_country)

        #print("getting target gdp")
        targets_country = get_gdp(targets_country)

        #print("processing targets")
        targets_country = targets_country[~targets_country.gdp.isnull()].reset_index(drop=True)  # remove the null targets
        targets_country['type'] = 'target'

        # combine
        #print("saving intermediate files")
        plants_country.to_csv(plantsfile)
        targets_country.to_csv(targetsfile)

#%% start
r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m")
COUNTRY_CODES = [row['iso3'] for row in r.json()['data']]

start = time.time()

