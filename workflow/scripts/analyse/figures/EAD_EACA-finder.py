"""Finds EAD and EACA

must have one named constant"""



import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


files_path_master = 'M:/Users/maxrob/Documents/OXFORD/ECI RA/Paper_figures/EAD_EACA/data'
plot_path = 'M:/Users/maxrob/Documents/OXFORD/ECI RA/Paper_figures/EAD_EACA'

def traprule(lst, spacing):
    """Trapezium rule"""
    return 0.5*spacing*(lst[0] + lst[-1] + sum(lst[1:-1]))






metrics = ['EAD', 'EACA']
metric_unit = {'EAD': 'USD', 'EACA': 'people'}

for metric in metrics:
    data_constant = pd.read_csv(os.path.join('data', metric, 'constant.csv'))
    print(f'\n\n {metric}')
    csv_files = [x for x in os.listdir(os.path.join('data', metric)) if 'constant' not in x]
    for ii, csv_file in enumerate(csv_files):
        # if ii >=1:
        #     break
        if ii == 0:
            data = pd.read_csv(os.path.join('data', metric, csv_file))
        else:
            data_new = pd.read_csv(os.path.join('data', metric,  csv_file))
            data = pd.merge(data, data_new, how='outer', on=['x'], suffixes=(f'_A{ii}', f'_B{ii}'))

    col_lst = ['y_cen', 'y_min', 'y_max']
    for col_id in col_lst:
        data[col_id] = data[[c for c in data.columns if col_id in c]].mean(axis=1)
    data = data[col_lst+['x']][:-7]  # TODO why

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
        print(f'{df_name}: {metric} = {format(metric_X, ".1E")} {metric_unit_indiv}, UR: {format(metric_L, ".1E")} {metric_unit_indiv} - {format(metric_U, ".1E")} {metric_unit_indiv}')

input('\n\npress enter')


