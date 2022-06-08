#!/usr/bin/env python
# coding: utf-8

"""
Code to determine which of the Universal Transverse Mercator projections are
most appropriate for a given longitude and latitude pair. UTM_EPSG_code returns
the EPSG indentifying code of the appropriate projection.

When run as a script, produce a labelled map of the EPSG codes.
"""

import math


def _calculate_UTM_zone_number(longitude: float, latitude: float) -> int:
    """
    Given a point specified as a latitude and longitude, return the zone number
    of the local Universal Transverse Mercator projection.

    Mostly, the earth is divided into 6 degree longitude sectors, however their
    are special cases for Norway and Svalbard (see below).

    Reference: https://pubs.usgs.gov/pp/1395/report.pdf page 62.

    Args:
        longitude (float): Decimal degrees
        latitude (float): Decimal degrees

    Returns
        int: Zone number of most appropriate UTM projection
    """

    # special cases
    # Svalbard
    if 72.0 <= latitude < 84.0:
        if 0.0 <= longitude < 9.0:
            return 31
        if 9.0 <= longitude < 21.0:
            return 33
        if 21.0 <= longitude < 33.0:
            return 35
        if 33.0 <= longitude < 42.0:
            return 37
    # Norway
    if 56.0 <= latitude < 64.0:
        if 3.0 <= longitude < 12.0:
            return 32

    # leaving aside special cases, typically to get the zone number , we count
    # up eastwards from the date line, in blocks of 6 degrees longitude
    return int(math.floor((longitude + 180) / 6) + 1)


def UTM_EPSG_code(*, longitude: float, latitude: float) -> int:
    """
    Given a point specified as a latitude and longitude, return the EPSG code
    of the local Universal Transverse Mercator projection. Mandatory kwargs.

    Args:
        longitude (float): Decimal degrees
        latitude (float): Decimal degrees

    Returns
        int: EPSG code of most appropriate UTM projection
    """

    base_epsg_code = 32600  # UTM projection EPSG codes start here

    # input checking
    if not -180.0 <= longitude <= 180:
        raise ValueError(
            "This function expects longitudes in the range [-180, 180]"
        )

    if latitude >= 84 or latitude <= -80:
        raise ValueError("UTM is inappropriate for polar regions, use UPS?")

    zone_number = _calculate_UTM_zone_number(longitude, latitude)

    epsg_code = base_epsg_code + int(zone_number)

    # southern hemisphere
    if (latitude < 0):
        epsg_code += 100

    return epsg_code


if __name__ == "__main__":
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits import basemap
    import numpy as np

    # if run as a script, show a plot of UTM codes

    fig = plt.figure(figsize=(12,8))
    m = basemap.Basemap(
        projection='merc',
        llcrnrlat=-75,
        urcrnrlat=78,
        llcrnrlon=-180,
        urcrnrlon=180,
        resolution='l'
    )
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='aqua')

    n_points = 600
    longitude = np.linspace(-180, 180, n_points)
    latitude = np.linspace(-79, 79, n_points)
    # for convenience, 'vectorise' the lookup function
    vectorised_lookup = np.vectorize(UTM_EPSG_code)

    # get the UTM codes on the grid
    xx, yy = np.meshgrid(longitude, latitude)
    UTM_field = vectorised_lookup(longitude=xx, latitude=yy)

    # make a random colormap and discretise it
    cmap = matplotlib.colors.ListedColormap(np.random.rand(256, 3))
    min_bound = np.min(UTM_field)
    max_bound = np.max(UTM_field)
    range_bound = max_bound - min_bound
    cmap_boundaries = np.linspace(min_bound, max_bound, range_bound + 1)
    norm = matplotlib.colors.BoundaryNorm(cmap_boundaries, ncolors=cmap.N)

    # plot the zones
    im = m.pcolormesh(
        longitude,
        latitude,
        UTM_field,
        latlon=True,
        cmap=cmap,
        norm=norm,
        alpha=0.3
    )

    # label all the map zones
    zone_width = 6  # typically
    label_longitudes = np.linspace(-180, 180 - zone_width, int(360 / zone_width)) + zone_width / 2
    for latitude in [15, -15]:
        for i, longitude in enumerate(label_longitudes):
            if i % 2 == 0:
                # every other label, bump latitude to ease readability
                x, y = m(longitude, latitude * 2)
            else:
                x, y = m(longitude, latitude)
            plt.text(
                x,
                y,
                UTM_EPSG_code(longitude=longitude, latitude=latitude),
                rotation=90,
                ha='center',
                va='center',
                color='white'
            )

    plt.title("EPSG codes for each Universal Transverse Mercator zone")
    plt.show()
