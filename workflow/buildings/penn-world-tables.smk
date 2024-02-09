"""
Penn World Table version 10.01

PWT version 10.01 is a database with information on relative levels of income,
output, input and productivity, covering 183 countries between 1950 and 2019

Feenstra, Robert C., Robert Inklaar and Marcel P. Timmer (2015), "The Next
Generation of the Penn World Table" American Economic Review, 105(10),
3150-3182, available for download at www.ggdc.net/pwt

Groningen Growth and Development Centre, 2023, "Penn World Table version 10.01",
https://doi.org/10.34894/QT5BCC, DataverseNL, V1
"""

rule download_pwt:
    output:
        zip="{OUTPUT_DIR}/input/capital-stocks/pwt.zip",
        xlsx="{OUTPUT_DIR}/input/capital-stocks/pwt/pwt1001.xlsx",
    shell:
        """
        if [ -e {output.zip} ]
        then
            echo "Skipping download"
        else
            curl -L -o {output.zip} -J  https://dataverse.nl/api/access/dataset/:persistentId/?persistentId=doi:10.34894/QT5BCC
        fi
        unzip -n -d $(dirname {output.xlsx}) {output.zip}
        """


rule extract_pwt:
    input:
        xlsx=rules.download_pwt.output.xlsx,
    output:
        csv="{OUTPUT_DIR}/input/capital-stocks/pwt.csv",
    run:
        import pandas
        df = pandas.read_excel(
            input.xlsx,
            sheet_name='Data',
        )
        # cn is "Capital stock at current PPPs (in mil. 2017US$)"
        df = df[['country','countrycode','year','cn']]
        df.to_csv(output.csv, index=False)
