"""
This file downloads the plants data to geoparequet files
"""

from importing_modules import *
from process_power_functions import changedir, idx


# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)
    box_id = "box_1941"

else:  # linux
    box_id = sys.argv[1]


def combine(lsts):
    """combines lists where not masked, takes average if mulitple not masked values"""
    if len({len(i) for i in lsts}) != 1:  # must be ame length
        raise RuntimeWarning("Lists not same length")
    llen = len(lsts[0])
    num_lsts = len(lsts)
    count_overlap = 0
    fix = [None]*llen
    for j in range(llen):  # count when different list contain different non-masked values
        lst_intersect = [lst[j] for lst in lsts]
        mask_count = lst_intersect.count(numpy.ma.masked) + lst_intersect.count(None)
        if mask_count <= num_lsts - 2:  # 2 or more non-masked values
            print(j, [lst[j] for lst in lsts])
            avgval = sum([0 if np.ma.is_masked(x)==True or x==None else x for x in lst_intersect])/(num_lsts - mask_count)
            fix[j] = avgval
            count_overlap += 1

    assert count_overlap/llen < 0.4 # assert less than 40% is overlap (rough)
    print(count_overlap, count_overlap/llen)

    all = [numpy.ma.masked]*llen
    for i in range(llen):
        for lst in lsts:
            if fix[i] != None:
                all[i] = fix[i]
            else:
                val = lst[i]
                if numpy.ma.is_masked(val) == False:
                    all[i] = val
                    break
    assert len(all) == llen
    return all



def get_target_areas(box_id):
    world_boxes_path = os.path.join('data', 'processed', 'world_boxes.gpkg')
    gdf = gpd.read_file(world_boxes_path)
    box = gdf.loc[gdf.box_id==box_id]['geometry']

    geod = Geod(ellps="WGS84")

    # Targets: Binary raster showing locations predicted to be connected to distribution grid.
    with rasterio.open(os.path.join("data","gridfinder","targets.tif")) as src:

        # for each country (see: code) overlay connections to distribution grid
        try:
            data, transform = rasterio.mask.mask(src, box, crop=True)
        except:
            return gpd.GeoDataFrame()

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


def get_population(box_id, targets):

    # below: population of target areas (targets.geometry), populations is list corresponding to areas (CHECK, I think)



    with open(os.path.join('data', 'processed', 'world_boxes_metadata.txt'), 'r') as filejson:
        world_boxes_metadata = json.load(filejson)
    box_country_dict_id = world_boxes_metadata['box_country_dict'][box_id]

    pop_all = []
    pop_d_all = []
    for kk, code in enumerate(box_country_dict_id):  # run for every country in the box
        print(f"{kk+1}/{len(box_country_dict_id)} -- {code}")

        gen = gen_zonal_stats(
                targets.geometry,
                os.path.join("data","population",f"{code}_ppp_2020_UNadj_constrained.tif"),
                stats=[],
                add_stats={'nansum': np.nansum},  # count NaN as zero for summation
                all_touched=True  # possible overestimate, but targets grid is narrower than pop
            )

        populations = [
            d['nansum']
            for d in tqdm(gen, desc=f"{code} population progress", total=len(targets.geometry))]

        # populations = [
        #     d['nansum']
        #     for d in zonal_stats(
        #         targets.geometry,
        #         os.path.join("data","population",f"{code}_ppp_2020_UNadj_constrained.tif"),
        #         stats=[],
        #         add_stats={'nansum': np.nansum},  # count NaN as zero for summation
        #         all_touched=True  # possible overestimate, but targets grid is narrower than pop
        #     )
        # ]
        # below: find population density at centroid of target areas
        ss = time.time()
        population_density = point_query(
            targets.centroid,
            os.path.join("data","population",f"{code}_ppp_2020_UNadj_constrained.tif")
        )
        print(f'time for {code} pop density: ',time.time()-ss)

        pop_all.append(populations)
        pop_d_all.append(population_density)

    if len(box_country_dict_id) == 1:
        pop_all = pop_all[0]
        pop_d_all = pop_d_all[0]

    if len(box_country_dict_id) > 1:
        # join up again
        pop_all = combine(pop_all)
        pop_d_all = combine(pop_d_all)

    targets['population'] = pop_all
    targets['population_density_at_centroid'] = pop_d_all
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
    # for ii in range(len(targets)):  # TODO vectorise
    #     lat_idx = round(np.interp(targets.centroid[ii].y,[-90,90],[2160, 0]),0)  # convert latitude to GDP nc file
    #     lon_idx = round(np.interp(targets.centroid[ii].x,[-180,180],[0,4320]),0)  # convert longitude to GDP nc file



    lat_idx_arr = np.interp(targets.centroid.y,[-90,90],[2160, 0])  # convert latitude to GDP nc file
    lon_idx_arr = np.interp(targets.centroid.x,[-180,180],[0,4320])  # convert longitude to GDP nc file



    # Find value
    assert len(lat_idx_arr) == len(lon_idx_arr)
    gdp_pc_lst = [ds['GDP_per_capita_PPP'][-1,lat_idx_arr[jj],lon_idx_arr[jj]] for jj in range(len(lat_idx_arr))]  # returns gdp of centroid of the target area (2015)

    gdp_pc_lst = [float(x) if numpy.ma.is_masked(x) == False else 0 for x in gdp_pc_lst]  # set masked to 0 (later removed)
    #print("sum", sum(gdp_pc_lst))
    targets['gdp_pc'] = gdp_pc_lst
    targets['gdp'] = targets.gdp_pc * targets.population
    return targets

if __name__ == '__main__':
    # Targets
    #print("getting target areas")
    targets_box = get_target_areas(box_id)

    if len(targets_box) == 0:
        cols = ["area_km2", "centroid", "geometry", "population", "population_density_at_centroid", "gdp_pc", "gdp", "type"]
        targets_box = gpd.GeoDataFrame(columns=cols)
        targets_box.loc[0,:] = [None]*len(cols)
        targets_box['box_id'] = box_id
    else:
        #print("getting target population")
        targets_box = get_population(box_id, targets_box)

        #print("getting target gdp")
        targets_box = get_gdp(targets_box)

        #print("processing targets")
        targets_box = targets_box[~targets_box.gdp.isnull()].reset_index(drop=True)  # remove the null targets
        targets_box['type'] = 'target'
        targets_box['box_id'] = box_id

    # combine
    #print("saving intermediate files")
    targets_box.to_csv(os.path.join("data", "processed", "all_boxes", box_id, f"targets_{box_id}.csv"), index=False)

