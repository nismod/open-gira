"""Common functions required to perform preprocessing
"""

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
