"""
This file can unpickle interactive matplotlib plots.
"""

import logging
import pickle
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


matplotlib.use("TkAgg")


def show_figure(fig: Figure) -> None:
    """
    When serialising a matplotlib figure with pickle, its original canvas is lost.
    This function creates a dummy figure and use its manager to display `fig`.
    """

    # TODO: prevent the dummy figure from also showing
    dummy = plt.figure()
    new_manager = dummy.canvas.manager
    new_manager.canvas.figure = fig
    fig.set_canvas(new_manager.canvas)


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    if len(sys.argv) == 1:
        logging.info(f"Usage:\n{sys.argv[0]} path/to/plot.pickle")

    for path in sys.argv[1:]:

        if path.lower().endswith(".pickle"):
            with open(path, "rb") as fp:
                fig: Figure = pickle.load(fp)

            show_figure(fig)

        else:
            logging.info(f"Do not recognise {path} as pickle file.")

    # start the event loop
    plt.show()
