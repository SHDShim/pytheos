import numpy as np
import uncertainties as uct
from .etc import isuncertainties


def debye_E_single(x):
    """
    calculate Debye energy using old fortran routine

    :params x: Debye x value
    :return: Debye energy
    """
    # make the function handles both scalar and array
    if ((x > 0.0) & (x <= 0.1)):
        result = 1. - 0.375 * x + x * x * \
            (0.05 - (5.952380953e-4) * x * x)
    # for 0.1 < x <= 7.25
    if ((x > 0.1) & (x <= 7.25)):
        result = ((((.0946173 * x - 4.432582) * x +
                    85.07724) * x - 800.6087) * x +
                  3953.632) / ((((x + 15.121491) * x +
                                 143.155337) * x + 682.0012) *
                               x + 3953.632)
    # for x > 7.25
    # it appears there might be error for this part, but never been exposed
    # because of rarity of such high x value.
    if (x > 7.25):
        exx = np.exp(-x)
        nn = np.round(25. / x)
        n = nn.astype(np.int64)
        temp = 0.
        if (n > 0):
            temp2 = 1.
            end = n + 1
            for i in range(1, end):
                temps = i * 1.
                temp2 = temp2 * exx
                x3 = temps * x
                temp = temp + temp2 * \
                    (6. + x3 * (6. + x3 * (3. + x3))) / \
                    (temps * temps * temps * temps)
        result = 3.0 * (6.493939402 - temp) / (x * x * x)
    return result


def debye_E(x):
    """
    calculate Debye energy using old fortran routine

    :params x: Debye x value
    :return: Debye energy
    :note: this is a wraper function (vectorization) for debye_E_single
    """
    if isuncertainties([x]):
        f_u = uct.wrap(debye_E_single)
        f_v = np.vectorize(f_u)
    else:
        f_v = np.vectorize(debye_E_single)
    return f_v(x)
