# Takes a list of geoparquet files and join corresponding
# geodataframes.

# file1.geoparguet
# id obs1 obs2 geometry
# 0  a    b    geom0
# 1  c    d    geom1

# file2.geoparguet
# id obs1 obs2 geometry
# 0  A    B    GEOM0
# 1  C    D    GEOM1

# joined.geoparguet
# id obs1 obs2 geometry
# 0  a    b    geom0
# 1  c    d    geom1
# 2  A    B    GEOM0
# 3  C    D    GEOM1

# Usage: python join_data [FILE] [output]
# Example: python join_data file1.geoparguet file2.geoparguet joined.geoparguet

import sys
import geopandas as gpd
import pandas
import logging


def append_data(base, slice_files):
    slice_files.pop()
    if len(slice_files) == 0:
        return base
    base = base.append(gpd.read_parquet(slice_files[-1]))
    return append_data(base, slice_files)


def add_custom_node_references(base):
    """
    When converting to .geoparquet we added nodes at the bounding box edges.
    These nodes have no reference. We need to make it easy to identify nodes by
    ensuring that nodes in the same location have the same reference.
    We'll make it easy on ourselves by giving our inserted nodes negative reference numbers.
    """
    def fix_nodes(base, node_ref):
        node = None
        for index, r in base.iterrows():
            if pandas.isna(r['start_node_reference']):
                node = (r['start_node_longitude'], r['start_node_latitude'])
                break
            if pandas.isna(r['end_node_reference']):
                node = (r['end_node_longitude'], r['end_node_latitude'])
                break

        if node is not None:
            logging.debug(f"Fixing references for node {node} (ref {node_ref})")
            def f(x):
                if x['start_node_longitude'] == node[0] and x['start_node_latitude'] == node[1]:
                    return node_ref
                else:
                    return x['start_node_reference']
            base['start_node_reference'] = base.apply(f, axis=1)

            def f(x):
                if x['end_node_longitude'] == node[0] and x['end_node_latitude'] == node[1]:
                    return node_ref
                else:
                    return x['end_node_reference']
            base['end_node_reference'] = base.apply(f, axis=1)
            changed = True
        else:
            changed = False

        return {
            'changed': changed,
            'base': base
        }

    node_ref = -1
    max_cycles = 100000000
    while max_cycles > 0:
        max_cycles -= 1
        new_base = fix_nodes(base, node_ref)
        if not new_base['changed']:
            break
        base = new_base['base']
        node_ref -= 1

    if max_cycles == 0:
        raise RecursionError(f'Max cycles exceeded for node referencing.')
    return base


if __name__ == "__main__":
    try:
        slice_files = snakemake.input
        output_file = snakemake.output[0]
    except NameError:
        slice_files = sys.argv[1:-1]
        output_file = sys.argv[-1]

    # When getting the input files from snakemake, there is no
    # garantee that they will always in the same order. Sort them for
    # consistency. Makes testing easier.
    slice_files = sorted(slice_files)
    # We're reading the different files as a stack from the top.  Let's
    # reverse the order of files to keep the first file on top.
    slice_files = slice_files[::-1]

    base = gpd.read_parquet(slice_files[-1])
    base = append_data(base, slice_files)
    base = add_custom_node_references(base)
    base.to_parquet(output_file)
