"""File to calculate which wind speeds/duration etc will cause damage to that cell"""

def applythreshold(winds_df):  # TODO some threshold [for now just takes top 4]
    """Returns locations where wind is 'severe' enough to affect the power lines. Based on some condition [needs work]"""

    if len(winds_df) > 4:
        #thrval = winds_df.sort_values('wind_location')['wind_location'].iloc[-4]  # take highest 4 vals
        thrval = 12  # m/s
        winds_df = winds_df[winds_df['wind_location'] >= thrval]

    return winds_df