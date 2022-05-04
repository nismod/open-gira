"""File to calculate which wind speeds/duration etc will cause damage to that cell"""


def applythreshold(winds_df):
    """Returns locations where wind is 'severe' enough to affect the power lines. Based on some condition
    Current condition is >= {thrval} m/s"""

    thrval = 33  # m/s
    winds_df = winds_df[winds_df["wind_location"] >= thrval]
    return winds_df
