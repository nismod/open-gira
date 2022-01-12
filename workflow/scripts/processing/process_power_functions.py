"""
functions and data required to perform preprocessing
"""


from importing_modules import *
import shapely.wkt as sw

# TODO: remove below lines once testing complete and solely on linux
def changedir():
    """For use on personal pc"""
    if "linux" not in sys.platform:
        path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
        os.chdir(path)




def idx(lat, lon):
    """Returns box index.
    boxlen - height and width of box
    lat_max - max lat value
    num_cols - number of cols (i.e. number of boxes across the equator
    lon_min - min lon value (lon_min can be negative)

    Note top left is (-,+), top right is (+,+), bottom right is (+,-) bottom left is (-,-)"""
    return (-np.round(lat/boxlen+1/2)+lat_max/boxlen)*num_cols+np.round(lon/boxlen-1/2)-lon_min/boxlen


def idxbox(lats, lons):
    """Returns a list of box_id values for the corresponding lats and lons"""
    boxes = idx(lats, lons)
    return [f"box_{int(elem)}" for elem in boxes]


def adjbox(idx):
    """Returns the indices of all boxes that touch the idx box (incl at corners). Note the top and bottom are cut offs but either left or right side is connected
    num_cols - number of cols (i.e. number of boxes across the equator
    tot_boxes - total number of boxes

    Note top left is (-,+), top right is (+,+), bottom right is (+,-) bottom left is (-,-)"""
    assert tot_boxes%num_cols == 0
    assert 0<=idx<=tot_boxes
    left, right = False, False
    adjacent = []
    if idx%num_cols==0:  # on left boundary
        adjacent+=[idx+num_cols-1, idx+1]  # adds [idx on right boundary (same row), one to the right]
        left = True
    elif (idx+1)%num_cols==0:  # on right boundary
        adjacent+=[idx-num_cols+1, idx-1]  # adds [idx on left boundary (same row), one to left]
        right = True
    else:
        adjacent+=[idx-1, idx+1]  # adds [left, right]

    if 0<=idx<=num_cols-1:  # on top row
        adjacent.append(idx+num_cols)  # adds below
        if left:
            adjacent+=[idx+num_cols+1, idx+num_cols+num_cols-1]  # adds [idx below to right one, on right boundary (same row below)]
        elif right:
            adjacent+=[idx+num_cols-1, idx+1]  # adds [idx below to left one, on left boundary (same row below)]
        else:
            adjacent+=[idx+num_cols-1, idx+num_cols+1]  # adds [below left, below right]
        assert len(adjacent) == 5
    if tot_boxes-num_cols+1<=idx<=tot_boxes:  # on bottom row
        adjacent.append(idx-num_cols)  # adds above
        if left:
            adjacent+=[idx-num_cols+1, idx-1]  # adds [idx above to right one, on right boundary (same row above)]
        if right:
            adjacent+=[idx-num_cols-1, idx-num_cols-num_cols+1]  # adds [idx above to left one, on left boundary (same row above)]
        else:
            adjacent+=[idx-num_cols-1, idx-num_cols+1]  # adds [above left, above right]
        assert len(adjacent) == 5
    if num_cols<=idx<=tot_boxes-num_cols:  # not on boundary
        if left:
            adjacent+=[idx-1, idx-num_cols, idx-num_cols+1, idx+num_cols+1, idx+num_cols+num_cols-1, idx+num_cols+1]  # add [above right (same row above), above, above right, below right (same row below), below, below right]
        elif right:
            adjacent+=[idx-num_cols-1, idx-num_cols, idx-num_cols-num_cols+1, idx+num_cols-1, idx+num_cols, idx+1]  # add [above left, above, above left (same row above), below left, below, below left (same row below)]
        else:
            adjacent+=[idx-num_cols-1, idx-num_cols, idx-num_cols+1, idx+num_cols-1, idx+num_cols, idx+num_cols+1]  # add [above left, above, above right, below left, below, below right]
        assert len(adjacent) == 8
    adjacent.sort()
    return adjacent
#
# def get_lines(code=None):
#     """Read gridfinder lines"""
#
#     s = time.time()
#
#     features = []
#     with fiona.open(os.path.join("data","gridfinder","grid.gpkg")) as src:
#
#
#
#         for jj,feature in tqdm(enumerate(src), desc='grid.gpkg features', total=len(src)):
#             # gridfinder GeoPackage stores an "fid" which GeoPandas ignores
#             # and fiona reads as "id", not to feature['properties']
#             # see https://github.com/geopandas/geopandas/issues/1035
#             geom = shape(feature['geometry'])
#
#             # if code != None:
#             #     check = code_geoms_gpd.within(geom)
#             #     if True in check:
#             #         continue
#             #     else:
#             #         break
#
#
#
#             features.append({
#                 'source_id': feature['id'],
#                 'source': feature['properties']['source'],
#                 'geometry': geom,
#             })
#
#     print("time for grid.gpkg processing: ",round((time.time() - s)/60, 2)," mins")
#
#
#     gdf_world = gpd.GeoDataFrame(features)
#
#     if code != None:  # preload
#         print(f'using {code}')
#         with fiona.open(os.path.join("data","adminboundaries",f"gadm36_{code}.gpkg"), "r", layer=3) as src_code:
#
#             code_geoms = []
#             for feature in src_code:
#                 code_geoms.append(shape(feature['geometry']))
#             print('create dataframe')
#             code_geoms_gpd = gpd.GeoDataFrame({'geometry':code_geoms})
#         print("overlay")
#         xmin, ymin, xmax, ymax = code_geoms_gpd.bounds.values[0]
#         gdf_world = gdf_world.cx[xmin:xmax,ymin:ymax]  # speed up
#
#         gdf_world = gdf_world.overlay(code_geoms_gpd, how='intersection')
#
#         print(gdf_world)
#     return gdf_world


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




#%% start
#r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m")
#COUNTRY_CODES = [row['iso3'] for row in r.json()['data']]
COUNTRY_CODES = ['PHL']  # TODO i think can remove these and two lines above

start = time.time()


changedir()  # TODO keep?
with open(os.path.join('data', 'processed', 'world_boxes_metadata.txt'), 'r') as filejson:
    world_boxes_metadata = json.load(filejson)
boxlen = world_boxes_metadata['boxlen']
lat_max = world_boxes_metadata['lat_max']
lon_min = world_boxes_metadata['lon_min']
num_cols = world_boxes_metadata['num_cols']
tot_boxes = world_boxes_metadata['tot_boxes']

