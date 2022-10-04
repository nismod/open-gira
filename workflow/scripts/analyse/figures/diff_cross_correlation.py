"""Plots heatmap difference (future mean - current) with modifications"""


import os
import pandas as pd
import matplotlib.pyplot as plt

try:
    output_dir = snakemake.params["output_dir"]  # type: ignore
    inputs = snakemake.input  # type: ignore
    remove = snakemake.params["remove_countries"]  # type: ignore
    name_cc_future_perc_diff = snakemake.params["name_cc_future_perc_diff"]  # type: ignore
    name_cc_future_diff = snakemake.params["name_cc_future_diff"]  # type: ignore
except:
    raise RuntimeError("Please use snakemake to define inputs")


assert type(remove) == list


plot_path = os.path.join(output_dir, "power_figures")
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


def plot_relation_matrix(matrix, country_index, title, fig_num, perc):
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

    country_index_plot = dict(
        zip(country_index_plot.keys(), range(len(country_index_plot)))
    )  # update

    if perc == True:
        plt.imshow(matrix, cmap="RdBu_r", vmax=50, vmin=-50)
        ext = " [% of constant climate]"
    else:
        plt.imshow(matrix, cmap="RdBu_r", vmax=0.1, vmin=-0.1)
        ext = " [-]"
    plt.xlabel("Country B")
    plt.yticks(
        list(country_index_plot.values()), labels=list(country_index_plot.keys())
    )
    plt.xticks(
        list(country_index_plot.values()),
        labels=list(country_index_plot.keys()),
        rotation=90,
    )
    plt.ylabel("Country A")

    # plt.title(title)
    cbar = plt.colorbar()
    cbar.set_label(f"JHR change{ext}", rotation=90)

    plt.show()
    plt.savefig(os.path.join(plot_path, f"JHR_{title}"), bbox_inches="tight")
    print(f"Saved {title}")


files_sorted = [[inputs[0]], inputs[1:]]


datas = []
for ii, files in enumerate(files_sorted):
    for jj, file_indiv in enumerate(files):
        print(file_indiv)
        if jj == 0:
            data = pd.read_csv(file_indiv, index_col=0)
        else:
            data = data + pd.read_csv(file_indiv, index_col=0)

    data = data / (jj + 1)
    datas.append(data)
    c_dict = dict(zip(data.columns, range(len(data.columns))))

assert len(datas) == 2

# abs
data_plot = 100 * (datas[1] - datas[0]) / datas[0]  # percentage
name = name_cc_future_perc_diff
plot_relation_matrix(data_plot, c_dict, name, 0, True)

# perc
data_plot = datas[1] - datas[0]
name = name_cc_future_diff
plot_relation_matrix(data_plot, c_dict, name, 1, False)
