import fenics as fe 
import numpy as np


def getSystemMatricesDomain(par):

#-------- parameter renaming --------
    sigma_a = par['sigma_a']
    sigma_s = par['sigma_s']
    sigma_t = fe.Expression('sigma_a + sigma_s', degree=1, sigma_a=sigma_a, sigma_s=sigma_s)

#-------- system matrices --------

    systemMatricesDomain = dict()
    systemMatricesDomain['Kzz'] = [[-0.33333333/sigma_t, -0.2981424/sigma_t, 0.0, 0.0, 0.0],
        [-0.2981424/sigma_t, -0.52380952/sigma_t, -0.25555063/sigma_t, 0.0, 0.0],
        [0.0, -0.25555063/sigma_t, -0.50649351/sigma_t, -0.25213645/sigma_t, 0.0],
        [0.0, 0.0, -0.25213645/sigma_t, -0.5030303/sigma_t, -0.25113118/sigma_t],
        [0.0, 0.0, 0.0, -0.25113118/sigma_t, -0.50175439/sigma_t]]

    systemMatricesDomain['S'] = [[-1.0*sigma_a, 0.0, 0.0, 0.0, 0.0],
        [0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0, 0.0],
        [0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0],
        [0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0],
        [0.0, 0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s]]

    return systemMatricesDomain