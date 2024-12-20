"""
The Changing Wealth of Nations (CWON) 2021 presents the most comprehensive
accounting of global wealth. Wealth — the stock of produced, natural, and human
capital — has been firmly established as a key measure of economic prosperity.

The wealth accounting approach provides two related sets of information:
comprehensive wealth accounts (a stock measure in total and per capita values),
and adjusted net savings (a flow measure). The wealth accounts were updated in
2021, using a new methodology described in The Changing Wealth of Nations 2021.
"""

rule download_cwon:
    output:
        full="{OUTPUT_DIR}/input/capital-stocks/CWON2021_Country_Tool_Full.xlsx",
        balanced="{OUTPUT_DIR}/input/capital-stocks/CWON2021_Country_Tool_Balanced.xlsx",
        doc="{OUTPUT_DIR}/input/capital-stocks/CWON2021_Methodology_October_2021.pdf",
    shell:
        """
        wget \
            -nc \
            --output-document="{output.full}" \
            https://datacatalogfiles.worldbank.org/ddh-published/0042066/DR0084604/CWON2021%20Country%20Tool%20-%20Full%20Dataset.xlsx?versionId=2022-06-06T13:38:38.6998180Z

        wget \
            -nc \
            --output-document="{output.balanced}" \
            https://datacatalogfiles.worldbank.org/ddh-published/0042066/DR0084043/CWON2021%20Country%20Tool%20-%20Balanced%20Dataset.xlsx?versionId=2022-06-06T13:38:44.8802314Z

        wget \
            -nc \
            --output-document="{output.doc}" \
            https://datacatalogfiles.worldbank.org/ddh-published/0042066/DR0084161/CWON%202021%20Methodology%20-%20October%202021.pdf?versionId=2022-06-06T13:38:42.3087244Z
        """

rule extract_cwon_produced_capital:
    input:
        xlsx=rules.download_cwon.output.full
    output:
        csv="{OUTPUT_DIR}/input/capital-stocks/CWON2021.csv"
    run:
        import pandas
        df = pandas.read_excel(
            input.xlsx,
            sheet_name='country',
            skiprows=1,
            na_values='..'
        ).query('wb_name != 0')
        df = df[['wb_name','wb_code','year','unit','pk']]
        df.to_csv(output.csv, index=False)
