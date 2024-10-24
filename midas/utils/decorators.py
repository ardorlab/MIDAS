## Import Block ##
import logging
import time
from functools import wraps
import cProfile
import pstats
import io

"""
These are useful decorators for the MIDAS framework.
"""

logger = logging.getLogger("MIDAS_logger")

## Exception handler decorator ##
def error_handler(func):
    def wrapper(*args, **kwargs):
        try:
            exitcode = func(*args, **kwargs)
        except Exception as e:
            msg="An error in the MIDAS execution has occurred with the following message: \n\n"
            msg+="- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
            logger.error(msg, exc_info=True)
            msg="- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
            logger.info(msg)
            exitcode = 1
        return exitcode
    return wrapper

## Timer decorator ##
def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Execution time of {func.__name__}: {execution_time:.4f} seconds")
        return result
    return wrapper

## Profiler decorator ##
def profiler(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        result = func(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        logger.debug(f"Profile stats for {func.__name__}:\n{s.getvalue()}")
        return result
    return wrapper