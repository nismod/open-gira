"""Takes all csv files in /data/ and plots them, legend is the name of the csv file.

Takes mean of all files unless 'constant' which has separate curve"""

import os
import matplotlib.pyplot as plt
import pandas as pd


try:
    output_dir = snakemake.params['output_dir']
    inputs = snakemake.input
    output = snakemake.output
    EACA = snakemake.params['EACA'] # set True if effective pop affected, else assumes EAD
    central_threshold = snakemake.params['central_threshold']
    minimum_threshold = snakemake.params['minimum_threshold']
    maximum_threshold = snakemake.params['maximum_threshold']
except:
    raise RuntimeError("Please use snakemake to define inputs")


assert type(output) != list
output = str(output)

plot_path = os.path.join(output_dir, 'power_figures')
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


metric = 'EAD'
unit = 'USD'
if EACA:
    metric = 'Effective population affected'
    unit = 'people'




csv_files = inputs[1:]  # future

data_constant = pd.read_csv(inputs[0])


for ii, csv_file in enumerate(csv_files):
    # if ii >=1:
    #     break
    print(csv_file)
    if ii == 0:
        data = pd.read_csv(csv_file)
    else:
        data_new = pd.read_csv(csv_file)
        data = pd.merge(data, data_new, how='outer', on=['x'], suffixes=(f'_A{ii}', f'_B{ii}'))

col_lst = ['y_cen', 'y_min', 'y_max']
for col_id in col_lst:
    data[col_id] = data[[c for c in data.columns if col_id in c]].mean(axis=1)


data = data[col_lst+['x']]
if EACA:
    data = data[:-7]  # TODO why



data_dict = {'future climate': data, 'current climate': data_constant}
plt.figure(figsize=(10, 7), dpi=100)
for ii, (df_name, df) in enumerate(data_dict.items()):
    x = df['x']
    y_cen = df['y_cen']
    y_min = df['y_min']
    y_max = df['y_max']
    c3 = ((138/256/(ii+1),171/256/(ii+1), 1))
    c2 = (1/(ii+1), 0, 1-1/(ii+1))
    c1 = (0, 1-1/(ii+1), 1/(ii+1))
    plt.fill_between(x, y_min, y_max, color=c1, alpha=0.2)
    plt.plot(x, y_min, color=c1, linestyle=':', label=f'{df_name} - Minimum wind threshold ({minimum_threshold}m/s)')  # plot interpolated line
    plt.plot(x, y_max, color=c1, linestyle='--', label=f'{df_name} - Maximum wind threshold ({maximum_threshold}m/s)')  # plot interpolated line
    plt.plot(x, y_cen, color=c2, label=f'{df_name} - Central wind threshold ({central_threshold}m/s)')  # plot interpolated line
    plt.scatter(x, y_cen, s=2, color=c2)  # plot data points

plt.xlabel('Return Period [years]')
plt.ylabel(f'{metric} [{unit}]')
plt.title(f"Empirical - {metric}")

plt.grid(axis='both', which='both')
plt.xscale("log")

plt.legend()

plt.savefig(output)

plt.show()
