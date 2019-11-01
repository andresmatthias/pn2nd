import fenics as fe 
import numpy as np


def getSystemMatricesDomain(par):

#-------- parameter renaming --------
    sigma_a = par['sigma_a']
    sigma_s = par['sigma_s']
    sigma_t = fe.Expression('sigma_a + sigma_s', degree=1, sigma_a=sigma_a, sigma_s=sigma_s)

#-------- system matrices --------

    systemMatricesDomain = dict()
    systemMatricesDomain['Kxx'] = [[0.33333333/(0.5*sigma_s - 1.0*sigma_t)]]

    systemMatricesDomain['Kxy'] = [[0.0]]

    systemMatricesDomain['Kyx'] = [[0.0]]

    systemMatricesDomain['Kyy'] = [[0.33333333/(0.5*sigma_s - 1.0*sigma_t)]]

    systemMatricesDomain['S'] = [[-1.0*sigma_a]]

    return systemMatricesDomain