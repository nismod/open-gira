"""Common functions required to perform preprocessing
"""
import warnings

import numpy as np
import shapely.wkt as sw
import snkit

warnings.filterwarnings("ignore", category=DeprecationWarning)


def adj(idx, num_cols, tot_boxes):
    """Returns the indices of all boxes that touch the idx box (incl at corners). Note the top and bottom are cut offs but either left or right side is connected
    num_cols - number of cols (i.e. number of boxes across the equator
    tot_boxes - total number of boxes

    Note top left is (-,+), top right is (+,+), bottom right is (+,-) bottom left is (-,-)"""
    assert tot_boxes % num_cols == 0
    assert 0 <= idx <= tot_boxes
    left, right = False, False
    adjacent = []
    if idx % num_cols == 0:  # on left boundary
        adjacent += [
            idx + num_cols - 1,
            idx + 1,
        ]  # adds [idx on right boundary (same row), one to the right]
        left = True
    elif (idx + 1) % num_cols == 0:  # on right boundary
        adjacent += [
            idx - num_cols + 1,
            idx - 1,
        ]  # adds [idx on left boundary (same row), one to left]
        right = True
    else:
        adjacent += [idx - 1, idx + 1]  # adds [left, right]

    if 0 <= idx <= num_cols - 1:  # on top row
        adjacent.append(idx + num_cols)  # adds below
        if left:
            adjacent += [
                idx + num_cols + 1,
                idx + num_cols + num_cols - 1,
            ]  # adds [idx below to right one, on right boundary (same row below)]
        elif right:
            adjacent += [
                idx + num_cols - 1,
                idx + 1,
            ]  # adds [idx below to left one, on left boundary (same row below)]
        else:
            adjacent += [
                idx + num_cols - 1,
                idx + num_cols + 1,
            ]  # adds [below left, below right]
        assert len(adjacent) == 5
    if tot_boxes - num_cols + 1 <= idx <= tot_boxes:  # on bottom row
        adjacent.append(idx - num_cols)  # adds above
        if left:
            adjacent += [
                idx - num_cols + 1,
                idx - 1,
            ]  # adds [idx above to right one, on right boundary (same row above)]
        if right:
            adjacent += [
                idx - num_cols - 1,
                idx - num_cols - num_cols + 1,
            ]  # adds [idx above to left one, on left boundary (same row above)]
        else:
            adjacent += [
                idx - num_cols - 1,
                idx - num_cols + 1,
            ]  # adds [above left, above right]
        assert len(adjacent) == 5
    if num_cols <= idx <= tot_boxes - num_cols:  # not on boundary
        if left:
            adjacent += [
                idx - 1,
                idx - num_cols,
                idx - num_cols + 1,
                idx + num_cols + 1,
                idx + num_cols + num_cols - 1,
                idx + num_cols + 1,
            ]  # add [above right (same row above), above, above right, below right (same row below), below, below right]
        elif right:
            adjacent += [
                idx - num_cols - 1,
                idx - num_cols,
                idx - num_cols - num_cols + 1,
                idx + num_cols - 1,
                idx + num_cols,
                idx + 1,
            ]  # add [above left, above, above left (same row above), below left, below, below left (same row below)]
        else:
            adjacent += [
                idx - num_cols - 1,
                idx - num_cols,
                idx - num_cols + 1,
                idx + num_cols - 1,
                idx + num_cols,
                idx + num_cols + 1,
            ]  # add [above left, above, above right, below left, below, below right]
        assert len(adjacent) == 8
    adjacent.sort()
    return adjacent


def idx(lat, lon, boxlen, lat_max, num_cols, lon_min, tot_boxes):
    """Returns box index.
    input:
        lat: np array, list or float
        lon: np array, list or float
    returns
        array

    boxlen - height and width of box
    lat_max - max lat value
    num_cols - number of cols (i.e. number of boxes across the equator)
    lon_min - min lon value (lon_min can be negative)
    tot_boxes - total number of boxes

    Note top left is (-,+), top right is (+,+), bottom right is (+,-) bottom left is (-,-)"""
    if type(lat) == list:
        lat = np.array(lat)
    if type(lon) == list:
        lon = np.array(lon)
    assert type(lat) == type(lon)
    # preprocess
    eps = 1e-8
    lon = np.where(
        lon != lon_min, lon, lon + eps
    )  # ensure correct box for left boundary
    lon = np.where(
        lon != lon_min + num_cols * boxlen, lon, lon - eps
    )  # ensure correct box for right boundary
    lat = np.where(
        lat != lat_max, lat, lat - eps
    )  # ensure correct box for top boundary
    lat = np.where(
        lat != lat_max - (tot_boxes / num_cols), lat, lat + eps
    )  # ensure correct box for bottom boundary

    # compute
    ret = (
        (-np.round(lat / boxlen + 1 / 2) + lat_max / boxlen) * num_cols
        + np.round(lon / boxlen - 1 / 2)
        - lon_min / boxlen
    )  # find indices

    if ret.size == 1:
        assert 0 <= ret <= tot_boxes - 1
    else:
        assert 0 <= min(ret) and max(ret) <= tot_boxes - 1  # check all indices

    return ret


def idxbox(lats, lons, boxlen, lat_max, num_cols, lon_min, tot_boxes):
    """Returns a list of box_id values for the corresponding lats and lons"""
    boxes = idx(lats, lons, boxlen, lat_max, num_cols, lon_min, tot_boxes)
    return [f"box_{int(elem)}" for elem in boxes]


def adjbox(box_idx, num_cols, tot_boxes):
    """Returns a list of box_id values for the corresponding adjacent box"""
    boxes = adj(box_idx, num_cols, tot_boxes)
    return [f"box_{int(elem)}" for elem in boxes]


def patch_nearest_edge(point, edges):
    """Set up network

    Find nearest edge to a point
    """

    if type(point) == str:
        point = sw.loads(point)  # if point is found as string -> convert to Point(# #)
        print("changed to point")
    geom = point.buffer(1e-2)
    # print(point, " : ", type(point))

    matches_idx = edges.sindex.nearest(geom.bounds)
    nearest_geom = min(
        [edges.iloc[match_idx] for match_idx in matches_idx],
        key=lambda match: point.distance(match.geometry),
    )
    return nearest_geom


snkit.network.nearest_edge = patch_nearest_edge


#%% set variables

# start = time.time()
#
# with open(
#     os.path.join("data", "processed", "world_boxes_metadata.json"), "r"
# ) as filejson:
#     world_boxes_metadata = json.load(filejson)
# boxlen = world_boxes_metadata["boxlen"]
# lat_max = world_boxes_metadata["lat_max"]
# lon_min = world_boxes_metadata["lon_min"]
# num_cols = world_boxes_metadata["num_cols"]
# tot_boxes = world_boxes_metadata["tot_boxes"]
