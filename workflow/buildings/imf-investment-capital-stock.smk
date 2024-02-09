"""
IMF Investment and Capital Stock Dataset, 1960-2019

(This version: May 2021)

Developed by Expenditure Policy (EP) Division, Fiscal Affairs Department (FAD),
International Monetary Fund (IMF)

This file provides comprehensive data on public investment and capital stock
(i.e. general government), private investment and capital stock, as well as
investment and capital stock arising from public-private partnerships (PPPs),
across the Fund membership countries.

The accompanying 2021 Update of the Manual "Estimating Public, Private, and PPP
Capital Stocks"
(https://infrastructuregovern.imf.org/content/dam/PIMA/Knowledge-Hub/dataset/InvestmentandCapitalStockDatabaseUserManualandFAQ_May2021.pdf)
describes in great detail the series' definitions, the investment series' data
sources, as well as the methodology in constructing the stock series. The
methodology follows the standard perpetual inventory equation and largely builds
on Gupta and others (2014) "Efficiency-Adjusted Public Capital and Growth" and
Kamps (2006) "New Estimates of Government Net Capital Stocks for 22 OECD
Countries, 1960–2001".


Please refer to it as "IMF Investment and Capital Stock Dataset, 2021” and add a
reference to the above-mentioned IMF Board Paper.


Key variables:

kgov_rppp

General government capital stock (constructed based on general government
investment flows "igov_rppp"), in billions of constant 2017 international
dollars.

kpriv_rppp

Private capital stock (constructed based on private investment flows
"igov_rppp"), in billions of constant 2017 international dollars.

kppp_rppp

Public-private partnership (PPP) capital stock (constructed based on PPP
investment flows "ippp_rppp"), in billions of constant 2017 international
dollars.

"""

rule download_imf_icsd:
    output:
        xlsx="{OUTPUT_DIR}/input/capital-stocks/icsd_dataset_2021.xlsx",
        doc="{OUTPUT_DIR}/input/capital-stocks/icsd_manual_2021.pdf",
    shell:
        """
        wget \
            -nc \
            --output-document="{output.xlsx}" \
            https://infrastructuregovern.imf.org/content/dam/PIMA/Knowledge-Hub/dataset/IMFInvestmentandCapitalStockDataset2021.xlsx

        wget \
            -nc \
            --output-document="{output.doc}" \
            https://infrastructuregovern.imf.org/content/dam/PIMA/Knowledge-Hub/dataset/InvestmentandCapitalStockDatabaseUserManualandFAQ_May2021.pdf
        """

rule extract_imf_icsd:
    input:
        xlsx=rules.download_imf_icsd.output.xlsx,
    output:
        csv="{OUTPUT_DIR}/input/capital-stocks/icsd.csv",
    run:
        import pandas
        df = pandas.read_excel(
            input.xlsx,
            sheet_name='Dataset',
        )
        # cn is "Capital stock at current PPPs (in mil. 2017US$)"
        df = df[['country','isocode','year','kgov_rppp', 'kpriv_rppp', 'kppp_rppp']]
        df.to_csv(output.csv, index=False)
