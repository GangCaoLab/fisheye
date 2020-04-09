import sys
import logging


def get_logger(name):
    log = logging.getLogger(name)
    handler = logging.StreamHandler(sys.stderr)
    LOGGING_FMT = "%(name)-20s %(levelname)-7s @ %(asctime)s: %(message)s"
    LOGGING_DATE_FMT = "%m/%d/%y %H:%M:%S"
    handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)
    return log