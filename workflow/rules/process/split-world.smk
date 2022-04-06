"""Split the world into boxes

"""



## Process Config Inputs ##
# Preliminary Snakemake Calculations
all_boxes = [
    f"box_{int(idx)}" for idx in range(0, int((180 - -180) * (90 - -90) / config['box_width_height'] ** 2))
]  # all boxes with width and height boxlen

# all_boxes = ["box_1584", "box_1585", "box_1431", "box_1432", "box_1433", "box_1505", "box_1575", "box_1576", "box_1577", "box_1647", "box_1648", "box_1649", "box_1650", "box_1652", "box_1653", "box_1719", "box_1720", "box_1721", "box_1722", "box_1791", "box_1792", "box_1793", "box_1794", "box_1798", "box_1863", "box_1864", "box_1865", "box_1866", "box_1870", "box_1871", "box_1936", "box_1937", "box_1941", "box_1942", "box_1943", "box_2013", "box_2014"]  # AU
# all_boxes = [f"box_{num}" for num in [884, 955, 956, 957, 1028, 1029, 1030, 1031, 1102, 1103, 1104]]


if config['specific_boxes'] != 'None':
    print('Using specified boxes')
    all_boxes = [f"box_{num}" for num in config['specific_boxes']]
if len(all_boxes)==0:
    print('Specific boxes incorrectly specified')



all_box_geoms = expand(
    os.path.join('data', "processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"),
    box_id=all_boxes,
)


rule world_splitter:
    input:
        os.path.join('data', "adminboundaries", "gadm36_levels.gpkg"),
    output:
        all_box_geoms,
        os.path.join('data', "processed", "world_boxes_metadata.txt"),
    params:
        boxlen_value = config['box_width_height'],
    script:
            os.path.join("..", "..", "scripts", "process", "world_split.py"
        )