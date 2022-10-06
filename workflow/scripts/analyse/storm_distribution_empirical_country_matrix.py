"""Plots the empirical storm relationship matrix between (two) countries and conditional
probability relationship
"""
import itertools as it
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    output_dir = snakemake.params["output_dir"]  # type: ignore
    thrval = snakemake.params["central_threshold"]  # type: ignore
except:
    output_dir = sys.argv[1]
    thrval = sys.argv[2]
    raise RuntimeError("Please use snakemake to define inputs")


## Inputs ##


def plot_relation_matrix(matrix, title, fig_num):
    """Plots and saves imshow"""
    f = plt.figure(fig_num)
    f.set_figwidth(10)
    f.set_figheight(8)

    plt.imshow(matrix, cmap="viridis")
    plt.xlabel("Country B")
    plt.yticks(list(country_index.values()), labels=list(country_index.keys()))
    plt.xticks(
        list(country_index.values()), labels=list(country_index.keys()), rotation=90
    )
    plt.ylabel("Country A")

    plt.title(title)
    plt.colorbar()

    plt.show()
    plt.savefig(os.path.join(stat_path_empirical, title))


## load stats ##
stat_path = os.path.join(output_dir, "power_output", "statistics")
csv_path = os.path.join(stat_path, f"combined_storm_statistics_{thrval}.csv")
stats = pd.read_csv(csv_path, keep_default_na=False)
storm_count = len(stats)

stat_path_empirical = os.path.join(stat_path, "empirical")
if not os.path.exists(stat_path_empirical):
    os.makedirs(stat_path_empirical)
stat_path_empirical_data = os.path.join(stat_path_empirical, "empirical_plotting_data")
if not os.path.exists(stat_path_empirical_data):
    os.makedirs(stat_path_empirical_data)

country_storm_count = dict()  # dictionary {country1: number_of_storms, country2: ... }
countries_overlap_master = (
    dict()
)  # dictionary {country1_country2: total_overlap_storm_count, country1_country3: ... }  note that the values represent the intersection
all_countries = set()  # set of all countries investigated
for jj, stats_indiv in tqdm(
    enumerate(stats["affected countries"]), desc="Iterating targets", total=len(stats)
):
    if type(stats_indiv) == str:  # First check for string
        if len(stats_indiv) >= 1:  # Then (if string) check not ""
            individual_countries = stats_indiv.split("_")
            all_countries.update(individual_countries)

            for country_a, country_b in it.combinations(
                individual_countries, 2
            ):  # for each unique country
                country1 = min(country_a, country_b)
                country2 = max(country_a, country_b)
                country_key = country1 + "_" + country2  # min first then max (str)
                if country_key in countries_overlap_master.keys():
                    countries_overlap_master[country_key] = (
                        countries_overlap_master[country_key] + 1
                    )  # one more storm which has both countries (intersection)
                else:
                    countries_overlap_master[country_key] = 1

            for country in individual_countries:
                if country in country_storm_count.keys():
                    country_storm_count[country] = (
                        country_storm_count[country] + 1
                    )  # one more storm
                else:
                    country_storm_count[country] = 1  # first storm

# set country index
all_countries_list = list(all_countries)
all_countries_list.sort()  # ensure each run same output
country_index = dict(
    zip(all_countries_list, list(range(len(all_countries))))
)  # {country1: 0, country2: 1, country3: 2, ...}

# create matrix nore that the row number corresponds to the country in the country index (same for column)
country_matrix_unint = -1 * np.ones(
    (len(country_index), len(country_index))
)  # union intersection matrix
country_matrix_condprob = country_matrix_unint.copy()  # conditional probability matrix


# note any non-overlapping
for country_a, country_b in it.combinations(
    all_countries, 2
):  # for each unique country
    country1 = min(country_a, country_b)
    country2 = max(country_a, country_b)
    key1 = country1 + "_" + country2
    key2 = country2 + "_" + country1
    if key1 not in countries_overlap_master.keys():
        countries_overlap_master[key1] = 0  # note no intersection

    if country1 not in country_storm_count.keys():
        country_storm_count[country1] = 0  # note no storms

    if country2 not in country_storm_count.keys():
        country_storm_count[country2] = 0  # note no storms

# assign data
for country_pair, count in countries_overlap_master.items():
    country1, country2 = country_pair.split("_")
    tot_storms = (
        country_storm_count[country1] + country_storm_count[country2] - count
    )  # total number of storms for both countries (union)
    idx1 = max(country_index[country1], country_index[country2])
    idx2 = min(country_index[country1], country_index[country2])
    country_matrix_unint[idx1, idx2] = (
        count / tot_storms
    )  # add the count (intersection) for the countries using the country_index index for the respecive countries as a fraction: intersction/union

    country_matrix_condprob[country_index[country1], country_index[country2]] = (
        count / country_storm_count[country1]
    )  # given country 1 is hit, what is likelyhood of country 2 being hit
    country_matrix_condprob[country_index[country2], country_index[country1]] = (
        count / country_storm_count[country2]
    )  # given country 2 is hit, what is likelyhood of country 1 being hit

# remove 0
country_matrix_unint[country_matrix_unint == -1] = np.nan
country_matrix_condprob[country_matrix_condprob == -1] = np.nan


dfA = pd.DataFrame(country_matrix_unint)
dfA.columns = all_countries_list
dfA.index = all_countries_list
dfA.to_csv(os.path.join(stat_path_empirical_data, "country_matrix_other.csv"))

dfB = pd.DataFrame(country_matrix_condprob)
dfB.columns = all_countries_list
dfB.index = all_countries_list
dfB.to_csv(os.path.join(stat_path_empirical_data, "country_matrix_both.csv"))


x = list(range(len(country_index))) * len(country_index)
y = []
_ = [y.append(len(country_index) * [idx]) for idx in range(len(country_index))]

title_unint = (
    "Empirical fraction of a storm hitting both countries given it hits one of them"
)
plot_relation_matrix(country_matrix_unint, title_unint, 0)

title_condprob = "Given country A is hit, what is the likelihood also B is hit"
plot_relation_matrix(country_matrix_condprob, title_condprob, 1)
#  plt.close('all')
