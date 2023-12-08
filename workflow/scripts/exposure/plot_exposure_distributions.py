"""
Plot exposure distributions (over event) for a given country.
"""

import logging
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_event_distributions(
    thresholds: list[str],
    exposure_by_event: pd.DataFrame,
    plot_dir: str,
    storm_set: str,
    geography_name: str
) -> None:
    """
    Given exposure data (lengths of edge exposed to wind speed in excess of a
    given threshold), plot the distribution of events.

    Args:
        thresholds: Wind speed thresholds for exposure data, used for indexing
            column from `exposure_by_event`.
        exposure_by_event: Lengths exposed in meters, with event index, threshold columns.
        plot_dir: Location to save histograms.
        storm_set: Name of storm set (for title)
        geography_name: Name of geography (for title)
    """
    # find extrema and number of bins based on all thresholds, so x-axis is
    # constant between thresholds
    x_min_log10 = np.log10(max([1000, 0.1 * min(exposure_by_event.min())]))
    x_max_log10 = np.log10(max(exposure_by_event.max()))
    n_bins = max([10, int(np.ceil(np.sqrt(len(exposure_by_event))))])
    logging.info(f"Using {n_bins} for {len(exposure_by_event)} events")
    bins = np.logspace(x_min_log10, x_max_log10, n_bins)

    # find y (frequency) maxima across thresholds
    y_max = 0
    for threshold in thresholds:
        frequency, bin_edges = np.histogram(exposure_by_event.loc[:, threshold], bins=bins)
        y_max = max([y_max, max(frequency)])
    y_max *= 5

    for threshold in thresholds:

        # filter out zero values
        data = exposure_by_event[threshold]
        non_zero_data = data[data > 0]

        f, ax = plt.subplots()

        ax.hist(
            non_zero_data,
            bins=bins,
            color="green",
            alpha=0.5,
            label="Distribution"
        )
        p90 = np.quantile(non_zero_data, 0.9)
        p95 = np.quantile(non_zero_data, 0.95)
        p99 = np.quantile(non_zero_data, 0.99)
        ax.axvline(p90, label=r"$p_{90}$ = " + f"{p90:.2E}", ls="--", color="pink", alpha=0.7)
        ax.axvline(p95, label=r"$p_{95}$ = " + f"{p95:.2E}", ls="--", color="red", alpha=0.7)
        ax.axvline(p99, label=r"$p_{99}$ = " + f"{p99:.2E}", ls="--", color="purple", alpha=0.7)
        ax.set_ylim(0.1, y_max)
        ax.set_xlim(10**x_min_log10, 5 * 10**x_max_log10)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Exposed edge length [m]")
        ax.set_ylabel("Frequency")
        ax.set_title(
            f"{geography_name} with {storm_set}\n"
            f"Exposure distribution, n={len(non_zero_data):,d} events\n"
            f"Exposed edge length @ {threshold} [m s-1] threshold"
        )
        ax.grid(alpha=0.1)
        ax.legend()
        plt.subplots_adjust(top=0.82, bottom=0.15)

        f.savefig(os.path.join(plot_dir, threshold.replace(".", "p") + ".png"))
        plt.close(f)

    return


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Reading exposure data")
    exposure_by_event = pd.read_parquet(snakemake.input.exposure_by_event)

    logging.info("Plotting event distributions")
    os.makedirs(snakemake.output.country_event_distributions)
    plot_event_distributions(
        list(exposure_by_event.columns),
        exposure_by_event,
        snakemake.output.country_event_distributions,
        snakemake.wildcards.STORM_SET,
        snakemake.wildcards.COUNTRY_ISO_A3
    )
