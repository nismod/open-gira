"""File to calculate which wind speeds/duration etc will cause damage to that cell"""


def applythreshold(winds_df, wind_threshold):
    """Returns locations where wind is 'severe' enough to affect the power line based on some condition.
    Current condition is >= {wind_threshold} m/s"""

    winds_df = winds_df[winds_df["wind_location"] >= wind_threshold]
    return winds_df
