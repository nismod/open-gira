"""Finds EAD and EACA

must have one named constant"""



import os
import pandas as pd
import numpy as np


try:
    output_dir = snakemake.params['output_dir']
    models_future = snakemake.params['models_future']
except:
    raise RuntimeError("Please use snakemake to define inputs")




def traprule(lst, spacing):
    """Trapezium rule"""
    return 0.5*spacing*(lst[0] + lst[-1] + sum(lst[1:-1]))






metrics = ['EAD', 'EACA']
metric_unit = {'EAD': 'USD', 'EACA': 'people'}
metric_file = {'EAD': 'empirical_reconstruction cost_plotting_data.csv', 'EACA': 'empirical_effective population affected_plotting_data.csv'}


for metric in metrics:
    write_lst = []
    data_constant = pd.read_csv(os.path.join(output_dir, 'power_output-constant', 'statistics', 'empirical', 'empirical_plotting_data', metric_file[metric]))
    print(f'\n\n {metric}')
    for ii, model in enumerate(models_future):

        model_path = os.path.join(output_dir, f'power_output-{model}', 'statistics', 'empirical', 'empirical_plotting_data', metric_file[metric])
        if ii == 0:
            data = pd.read_csv(model_path)
        else:
            data_new = pd.read_csv(model_path)
            data = pd.merge(data, data_new, how='outer', on=['x'], suffixes=(f'_A{ii}', f'_B{ii}'))

    col_lst = ['y_cen', 'y_min', 'y_max']
    for col_id in col_lst:
        data[col_id] = data[[c for c in data.columns if col_id in c]].mean(axis=1)
    data = data[col_lst+['x']][:-7]  # data specific

    data_dict = {'current climate': data_constant, 'future climate': data}
    for ii, (df_name, df) in enumerate(data_dict.items()):
        x = df['x']
        y_cen = df['y_cen']
        y_min = df['y_min']
        y_max = df['y_max']


        metric_X = int(traprule(np.array(y_cen), 1/max(x)))
        metric_U = int(traprule(np.array(y_min), 1/max(x)))
        metric_L = int(traprule(np.array(y_max), 1/max(x)))


        metric_unit_indiv = metric_unit[metric]
        details_metric = f'{df_name}: {metric} = {format(metric_X, ".1E")} {metric_unit_indiv}, UR: {format(metric_L, ".1E")} {metric_unit_indiv} - {format(metric_U, ".1E")} {metric_unit_indiv}'
        print(details_metric)
        write_lst.append(details_metric)

    txt_file = os.path.join(output_dir, 'power_figures', f"{metric}_total.txt")
    with open(txt_file, 'w') as f:
        f.write(write_lst[0]+'  ||  ')
        f.write(write_lst[1])


