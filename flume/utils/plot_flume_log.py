import matplotlib.pyplot as plt
import os
from flume.utils.read_flume_log import read_flume_log
import argparse
import math
from icecream import ic
from matplotlib.gridspec import GridSpec


def main():
    # Construct the argument parser
    parser = argparse.ArgumentParser(
        description="Plots the data from a Flume log file against the iteration number."
    )

    # Add an argument that specifies the directory
    parser.add_argument(
        "--directory",
        "-d",
        type=str,
        default=".",
        nargs="?",
        help="Directory where the log file is located, which defaults to the current directory.",
    )

    # Add an argument that specifies the name of the log file
    parser.add_argument(
        "--filename",
        "-f",
        type=str,
        default="flume.log",
        nargs="?",
        help="Filename for the lot file, which defaults to 'flume.log'.",
    )

    # Add an argument that optionally specifies the headers to plot.
    parser.add_argument(
        "--headers",
        "-hd",
        nargs="*",
        default=None,
        help="Optional argument that, when provided, will only plot the provided header information. These headers must correspond exactly to the headers contained within the log file, including the descriptors (i.e. 'obj:', 'con:', and/or 'other':).",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Read the log file and extract the data
    data = read_flume_log(directory=args.directory, filename=args.filename)

    # Plot either all of the data or just the select data based on the input argument
    if args.headers is not None:
        # Get the number of headers
        num_headers = len(args.headers)

        # Set the keys to loop over to exclude the "iter" header
        filtered_keys = {k: v for k, v in data.items() if k in args.headers}

    else:
        # Get the number of headers
        num_headers = len(data.keys())

        # Set the keys to loop over to exclude the "iter" header
        filtered_keys = {k: v for k, v in data.items() if k != "iter"}

    # Set the number of rows/columsn for the figure
    ncols = 3
    nrows = math.ceil(num_headers / ncols)

    # Construct the figure and gridspec
    fig = plt.figure(figsize=(12, 4 * nrows))

    gs = GridSpec(nrows=nrows, ncols=ncols, figure=fig, wspace=0.35, hspace=0.25)

    # Get the color cycle
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # Loop over all of the headers and plot the data
    for i, key in enumerate(filtered_keys):

        # Add the subplot to the figure
        ax = fig.add_subplot(gs[math.floor(i / ncols), i % ncols])

        # Select the color for the data
        color = colors[i % len(colors)]

        # Plot the data
        ax.plot(data["iter"], data[key], color=color)
        ax.set_xlabel("Iteration")
        ax.set_ylabel(f"{key}")

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
