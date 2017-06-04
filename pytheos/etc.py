import numpy as np
import uncertainties as uct


def isuncertainties(arg_list):
    """
    check if the input list contains any elements with uncertainties class

    :param arg_list: list of arguments
    :return: True/False
    """
    for arg in arg_list:
        if isinstance(arg, (list, tuple)) and isinstance(arg[0], uct.UFloat):
            return True
        elif isinstance(arg, np.ndarray) and isinstance(
                np.atleast_1d(arg)[0], uct.UFloat):
            return True
        elif isinstance(arg, (float, uct.UFloat)) and \
                isinstance(arg, uct.UFloat):
            return True
    return False
