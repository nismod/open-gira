"""
Process all files
"""


rule preprocess_all:
    input:
        rules.process_gridfinder.output,
        GENERATORS_LINES_CONSUMERS_ALL_BOXES,
        CONNECTOR_OUT,
