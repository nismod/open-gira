"""
Process all files
"""


rule preprocess_all:
    input:
        GENERATORS_LINES_CONSUMERS_ALL_BOXES,
        CONNECTOR_OUT,
