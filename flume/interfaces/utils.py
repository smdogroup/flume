import os


class Logger:

    def __init__(self, log_path, log_name="flume.log"):

        self.log_name = log_name
        self.log_path = log_path

        # Set the log file path
        self.log_filepath = self.set_log_file_path()

    def set_log_file_path(self):

        # Set the filepath
        filepath = os.path.join(self.log_path, self.log_name)

        # Remove if it exists already
        if os.path.exists(filepath):
            os.remove(filepath)

        return filepath

    def log(self, txt="", end="\n"):
        # Open the filepath
        with open(self.log_filepath, "a") as f:
            # Write the information
            f.write("%s%s" % (txt, end))

        return
