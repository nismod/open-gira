"""
Process all files
"""


rule preprocess_all:
    input:
        out_network,
        out_connector,
