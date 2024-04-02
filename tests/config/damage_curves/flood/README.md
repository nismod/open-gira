# Damage curves

Each CSV file in this directory pertains to a particular `Asset`. When networks
of nodes or edges are created, each feature of interest should be tagged with
an `asset_type`. Currently accepted values may found in the open_gira.assets module.

The contents of the CSV files should be:

1) A header describing the source of the data to follow. Header lines should be prefixed with '#', e.g.
\# Source:
\# Nirandjan et al., 2024, F8.4

2) The first non-header line gives the column headings. For flooding this should be:
inundation_depth_(m),damage_fraction

3) The inundation and damage ratio data should follow, one record per line, e.g.
0.0,0.0
0.5,0.2
1.0,0.33
