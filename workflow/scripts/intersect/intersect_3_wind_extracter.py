"""Process storm data and return the wind speed at each grid location.

TODO: Investigate speeding up this script. N.B. ~80% of execution time was file I/O.

"""
import os
import geopandas
import numpy
import pandas
from pyproj import CRS
from tqdm import tqdm


def main():
    edges_split_path = snakemake.input.edges_split  # type: ignore
    # TODO filter for relevance or replace with single file for (basin/model/sample)
    storm_files = snakemake.input.storm_files  # type: ignore
    output_path = snakemake.output.wind_speeds  # type: ignore

    # TODO look this up per box and/or skip basins that don't intersect
    basin = snakemake.wildcards.STORM_BASIN  # type: ignore

    # TODO check config to restrict analysis
    # snakemake.config.storm_files_sample_set
    # snakemake.config.specific_storm_analysis

    network = geopandas.read_parquet(edges_split_path)

    # Environmental pressure values (standard estimate of background pressure away from the
    # cyclone) are taken from the AIR hurricane model, table 3 in Butke (2012).
    # Available at:
    # https://www.air-worldwide.com/publications/air-currents/2012/the-pressures-on-increased-realism-in-tropical-cyclone-wind-speeds-through-attention-to-environmental-pressure/
    environmental_pressure = {
        "NI": 1006.5,
        "SA": 1014.1,
        "NA": 1014.1,
        "EP": 1008.8,
        "SI": 1010.6,
        "SP": 1008.1,
        "WP": 1008.3,
    }

    basin_environmental_pressure = environmental_pressure[basin]

    # TODO generate unique grid points from network cell indices
    points = None  # network

    # TODO check track_id in subset or do some subset read of storm files
    for storm_file in storm_files:
        storms = read_storms(storm_file)
        for track_id, track in tqdm(storms):
            print(track_id)
            track = interpolate_track(track)  # interpolate extra track points
            storm_max_ws = evaluate_speeds(track, points, basin_environmental_pressure)
            points[track_id] = storm_max_ws

    # TODO consider output - network joined with max speeds for all storms? could be large
    network_with_speeds = None
    network_with_speeds.to_parquet(output_path)



def read_storms(stormfile):
    # load in cyclone tracks for region
    tracks = pd.read_csv(
        stormfile,
        names=(
            "year",
            "month",
            "tc_number",
            "timestep",
            "basin_id",
            "lat",
            "lon",
            "min_pressure_hpa",
            "max_wind_speed_ms",
            "radius_to_max_winds_km",
            "category",
            "landfall",
            "distance_to_land_km",
        ),
    )

    # add unique number to cyclone data
    sample = os.path.basename(stormfile).replace(".txt", "")
    tracks["track_id"] = (
        sample
        + "_"
        + tracks["year"].astype(int).astype(str)
        + "_"
        + tracks["tc_number"].astype(int).astype(str)
    )

    tracks["sample"] = sample

    # change geometry from 0-360 to -180-180
    tracks.lon = numpy.where(tracks.lon > 180, tracks.lon - 360, tracks.lon)
    tracks["geometry"] = geopandas.points_from_xy(
        tracks.lon, tracks.lat, crs="EPSG:4326"
    )
    tracks = geopandas.GeoDataFrame(tracks.drop(columns=["lat", "lon"]))
    return tracks.groupby("track_id")


def interpolate_track(track, substeps=5):
    track = track.drop(columns="geometry").copy()
    track["x"] = track.geometry.x
    track["y"] = track.geometry.y

    dfs = [track]
    columns_to_interpolate = [
        "min_pressure_hpa",
        "max_wind_speed_ms",
        "radius_to_max_winds_km",
        "category",
        "landfall",
        "distance_to_land_km",
        "track_id",  # do we want to interpolate this??
        "x",
        "y",
    ]
    for increment in range(1, substeps):
        substep = increment / substeps
        tmp = track.copy()
        tmp[columns_to_interpolate] = numpy.nan
        tmp.timestep = tmp.timestep + substep
        dfs.append(tmp)

    interpolated_track = (
        pandas.concat(dfs)
        .sort_values("timestep")
        .set_index("timestep")
        .interpolate("linear")
        .fillna(0)
    )
    interpolated_track["geometry"] = geopandas.points_from_xy(
        interpolated_track.x, interpolated_track.y, crs="EPSG:4326"
    )
    return geopandas.GeoDataFrame(interpolated_track)


def holland_wind_field(
    radius_to_max_winds_km,
    wind_speed_ms,
    pressure_hpa,
    pressure_env_hpa,
    distance_m,
    lat_degrees,
):
    """Calculate wind speed at a point some distance from a cyclone track point.

    See in particular Section 3.2 in Lin and Chavas (2012).

    References
    ----------
    - Lin and Chavas (2012)
      https://www.sbafla.com/method/portals/methodology/FloodJournalArticles/Lin_Chavas_JGR12_ParametricWind.pdf
    - Holland (1980)
      https://doi.org/10.1175/1520-0493(1980)108%3C1212:AAMOTW%3E2.0.CO;2
    """
    radius_to_max_winds_m = radius_to_max_winds_km * 1000
    rho = 1.10
    f = numpy.abs(1.45842300e-4 * numpy.sin(numpy.radians(lat_degrees)))
    e = 2.71828182846
    delta_p = (pressure_env_hpa - pressure_hpa) * 100
    # case where (pressure_env_hpa == pressure_hpa) so p_drop is zero will raise ZeroDivisionError
    B = (
        numpy.power(wind_speed_ms, 2) * e * rho
        + f * wind_speed_ms * radius_to_max_winds_m * e * rho
    ) / delta_p
    Vg = (
        numpy.sqrt(
            # case where distance_m is zero will raise ZeroDivisionError
            (
                numpy.power(radius_to_max_winds_m / distance_m, B)
                * B
                * delta_p
                * numpy.exp(0 - (radius_to_max_winds_m / distance_m) ** B)
            )
            + (numpy.power(radius_to_max_winds_m, 2) * numpy.power(f, 2) / 4)
        )
        - (f * radius_to_max_winds_m) / 2
    )
    return Vg  # , B, delta_p, f


def evaluate_speeds(track, points, pressure_env_hpa):
    """Evaluate maximum wind speeds at points given a track"""
    # could calculate rolling max to keep memory down
    speeds = numpy.zeros((len(track), len(points)))
    geod_wgs84 = CRS("epsg:4326").get_geod()
    for i, track_point in enumerate(track.itertuples()):
        # distances from track to evaluation points
        _, _, distance_m = geod_wgs84.inv(
            numpy.full(len(points), track_point.geometry.x),
            numpy.full(len(points), track_point.geometry.y),
            points.geometry.x,
            points.geometry.y,
        )
        wind_speeds = holland_wind_field(
            track_point.radius_to_max_winds_km,
            track_point.max_wind_speed_ms,
            track_point.min_pressure_hpa,
            pressure_env_hpa,
            distance_m,
            track_point.geometry.y,
        )
        speeds[i] = wind_speeds
    return speeds.max(axis=0)


if __name__ == "__main__":
    main()
