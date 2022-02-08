"""
Imports required for preprocessing
"""

import os
import glob
import requests
import fiona
import geopandas as gpd
gpd._compat.USE_PYGEOS = False
import numpy as np
import numpy.ma
import json
import pandas as pd
import rasterio
import rasterio.mask
import rasterio.features
import snkit
from pyproj import Geod
from rasterstats import zonal_stats, point_query, gen_zonal_stats
from shapely.geometry import shape, LineString, Point
from shapely import wkt
import shapely.wkt as sw
import netCDF4 as nc4  # TODO add to requirements
import time  # TODO add to requirements
import sys
import ast
from tqdm import tqdm
from collections import defaultdict
import networkx as nx
from collections import ChainMap
from shapely.errors import ShapelyDeprecationWarning
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
