import os
from icecream import ic


class Logger:
    log_name = "stdout.log"

    @staticmethod
    def set_log_path(log_path):
        if os.path.exists(log_path):
            os.remove(log_path)
        Logger.log_name = log_path

    @staticmethod
    def log(txt="", end="\n"):
        with open(Logger.log_name, "a") as f:
            f.write("%s%s" % (txt, end))
        return
