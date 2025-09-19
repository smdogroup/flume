import os
import re


def read_flume_log(directory, filename="flume.log"):
    """
    Utility script that reads the Flume log file named filename within the specified directory.

    Parameters
    ----------
    directory : str
        String that specifies the directory where the log file is located
    filename : str
        String that specifies the name of the Flume log file, which is nominally 'flume.log'

    Returns
    -------
    data : dict
        A dictionary that stores the data for each column in the log file, corresponding to the iterations from an optimization problem.
    """

    # Construct the filepath
    filepath = os.path.join(directory, filename)

    # Open the file
    with open(filepath, "r") as f:
        # Parse the lines
        lines = [line.strip() for line in f if line.strip()]

    # Set the pattern for the regular expression (extracts all of the parts from the header row that follow obj:, con:, or other: formats)
    pattern = r"(obj:\s*.*?)(?=obj:|con:|other:|$)|(con:\s*.*?)(?=obj:|con:|other:|$)|(other:\s*.*?)(?=obj:|con:|other:|$)"

    # Initialize the header list
    headers = ["iter"]

    # Find all matches with the pattern
    matches = re.findall(pattern, lines[0])
    for match in matches:
        headers.append(next(filter(None, match)))

    # Strip the extra white space from the end
    headers = [header.strip() for header in headers]

    # headers += re.findall(pattern, lines[0])
    # headers = fixed_line.split()

    data = {h: [] for h in headers}
    # Loop through the rest of the lines
    for line in lines[1:]:

        # Skip the current line if it is a header line
        if line.startswith("iter"):
            continue

        # Extract the values
        parts = line.split()

        # Otherwise, extract the data from the line
        for h, v in zip(headers, parts):
            data[h].append(float(v))

    return data
