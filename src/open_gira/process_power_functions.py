"""Common functions required to perform preprocessing
"""


def adj(idx: int, ncols: int, tot_boxes: int) -> list[int]:
    """Returns the indices of all boxes that touch the idx box (incl at corners).

    Note the assumption that index 0 (or (0,0)) is top-left.

    Parameters
    ----------
    ncols: int, number of cols (i.e. number of boxes across the equator
    tot_boxes: int, total number of boxes
    """
    assert tot_boxes % ncols == 0
    assert 0 <= idx <= tot_boxes
    nrows = tot_boxes / ncols

    x, y = xy_from_idx(idx, ncols)

    xy = [
        (x % ncols, y % nrows)  # wrap around
        for x, y in [
            (x - 1, y - 1),  # above left
            (x, y - 1),  # above
            (x + 1, y - 1),  # above right
            (x - 1, y),  # left
            (x + 1, y),  # right
            (x - 1, y + 1),  # below left
            (x, y + 1),  # below
            (x + 1, y + 1),  # below right
        ]
    ]
    return [idx_from_xy(x, y, ncols) for x, y in xy]


def xy_from_idx(idx: int, ncols: int) -> tuple[int, int]:
    """Find (x, y) tuple index from integer index, given number of columns"""
    x = idx % ncols
    y = idx // ncols
    return x, y


def idx_from_xy(x: int, y: int, ncols: int) -> int:
    """Find integer index from (x, y) tuple index, given number of columns"""
    return int(y * ncols + x)
