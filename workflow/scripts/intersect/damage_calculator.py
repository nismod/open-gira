"""File to calculate which wind speeds/duration etc will cause damage to that cell"""


def applythreshold(winds_df):  # TODO some threshold
    """Returns locations where wind is 'severe' enough to affect the power lines. Based on some condition [needs work]"""
    if len(winds_df) > 4:
        thrval = 12  # m/s
        winds_df = winds_df[winds_df["wind_location"] >= thrval]
    return winds_df
