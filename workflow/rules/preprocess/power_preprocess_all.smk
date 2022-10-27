"""
Process all files
"""


rule preprocess_all:
    input:
        CONNECTOR_OUT,
