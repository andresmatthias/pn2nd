import fenics as fe 
import numpy as np


def getSystemMatricesDomain(par):

#-------- parameter renaming --------
    sigma_a = par['sigma_a']
    sigma_s = par['sigma_s']
    sigma_t = fe.Expression('sigma_a + sigma_s', degree=1, sigma_a=sigma_a, sigma_s=sigma_s)

#-------- system matrices --------

    systemMatricesDomain = dict()
    systemMatricesDomain['Kzz'] = [[0.33333333/(0.33333333*sigma_s - 1.0*sigma_t), 0.2981424/(0.33333333*sigma_s - 1.0*sigma_t)],
        [0.2981424/(0.33333333*sigma_s - 1.0*sigma_t), -(1.0*(0.085714286*sigma_s - 0.52380952*sigma_t))/(0.33333333*sigma_s*sigma_t - 1.0*sigma_t**2)]]

    systemMatricesDomain['S'] = [[-1.0*sigma_a, 0.0],
        [0.0, - 1.0*sigma_a - 1.0*sigma_s]]

    return systemMatricesDomain