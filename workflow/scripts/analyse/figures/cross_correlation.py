"""Plots heatmap with modifications"""



import os
import pandas as pd
import matplotlib.pyplot as plt


try:
    output_dir = snakemake.params['output_dir']
    inputs = snakemake.input
    remove = snakemake.params['remove_countries']
    name_cc_constant = snakemake.params['name_cc_constant']
    name_cc_future = snakemake.params['name_cc_future']
except:
    raise RuntimeError("Please use snakemake to define inputs")


assert type(remove) == list


#
# files_path = 'M:/Users/maxrob/Documents/OXFORD/ECI RA/Paper_figures/plotter/cross-correlation/data'
# files = os.listdir(files_path)
# plot_path = 'M:/Users/maxrob/Documents/OXFORD/ECI RA/Paper_figures/plotter/cross-correlation'
# remove = ['USA', 'VEN', 'CYM', 'VCT', 'BHS']  # countries to remove
# remove = ['USA', 'VEN', 'CYM', 'VCT', 'BHS', 'ATG', 'DMA', 'LCA', 'TTO']  # countries to remove


plot_path = os.path.join(output_dir, 'power_figures')
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

def plot_relation_matrix(matrix, country_index, title, fig_num):
    """Plots and saves imshow"""
    country_index_plot = country_index.copy()
    f = plt.figure(fig_num)
    f.set_figwidth(10)
    f.set_figheight(8)

    for c in remove:  # update, remove
        if c in matrix.columns:
            matrix = matrix.drop(c, axis=1)
            matrix = matrix.drop(c, axis=0)
            del country_index_plot[c]

    country_index_plot = dict(zip(country_index_plot.keys(), range(len(country_index_plot))))  # update
    plt.imshow(matrix, cmap='viridis')
    plt.xlabel('Country B')
    plt.yticks(list(country_index_plot.values()), labels=list(country_index_plot.keys()))
    plt.xticks(list(country_index_plot.values()), labels=list(country_index_plot.keys()), rotation=90)
    plt.ylabel('Country A')

    #plt.title(title)
    cbar = plt.colorbar()
    cbar.set_label(f'JHR ({title}) [-]', rotation=90)

    plt.show()
    plt.savefig(os.path.join(plot_path, f'JHR_{title}'), bbox_inches='tight')
    print(f'Saved {title}')



#files_sorted = [[file for file in inputs if 'constant' in file], [file for file in inputs if 'constant' not in file]] # TODO fix
files_sorted = [[inputs[0]], inputs[1:]]

# files_sorted = [[file for file in files if 'constant' in file]]
# files_sorted = [[file for file in files if 'constant' not in file]]

for ii, files in enumerate(files_sorted):
    for jj, file_indiv in enumerate(files):
        if jj == 0:
            data = pd.read_csv(file_indiv, index_col=0)
        else:
            data = data + pd.read_csv(file_indiv, index_col=0)

    data = data / (jj+1)
    c_dict = dict(zip(data.columns, range(len(data.columns))))
    if ii == 0:  # constant
        name = name_cc_constant
    else:
        name = name_cc_future
    plot_relation_matrix(data, c_dict, name, ii)