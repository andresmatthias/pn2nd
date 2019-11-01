import fenics as fe 
import numpy as np


def getSystemMatricesDomain(par):

#-------- parameter renaming --------
    sigma_a = par['sigma_a']
    sigma_s = par['sigma_s']
    sigma_t = fe.Expression('sigma_a + sigma_s', degree=1, sigma_a=sigma_a, sigma_s=sigma_s)

#-------- system matrices --------

    systemMatricesDomain = dict()
    systemMatricesDomain['Kxx'] = [[0.33333333/(0.5*sigma_s - 1.0*sigma_t), 0.0, -0.1490712/(0.5*sigma_s - 1.0*sigma_t), 0.25819889/(0.5*sigma_s - 1.0*sigma_t)],
        [0.0, (0.13928571*sigma_s - 0.42857143*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), 0.0, 0.0],
        [-0.1490712/(0.5*sigma_s - 1.0*sigma_t), 0.0, (0.094047619*sigma_s - 0.23809524*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), -(1.0*(0.03917734*sigma_s - 0.16495722*sigma_t))/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2)],
        [0.25819889/(0.5*sigma_s - 1.0*sigma_t), 0.0, -(1.0*(0.03917734*sigma_s - 0.16495722*sigma_t))/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), (0.13928571*sigma_s - 0.42857143*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2)]]

    systemMatricesDomain['Kxy'] = [[0.0, 0.25819889/(0.5*sigma_s - 1.0*sigma_t), 0.0, 0.0],
        [0.25819889/(0.5*sigma_s - 1.0*sigma_t), 0.0, -(1.0*(0.03917734*sigma_s - 0.16495722*sigma_t))/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), (0.075*sigma_s)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2)],
        [0.0, -(1.0*(0.03917734*sigma_s - 0.16495722*sigma_t))/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), 0.0, 0.0],
        [0.0, -(0.075*sigma_s)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), 0.0, 0.0]]

    systemMatricesDomain['Kyx'] = [[0.0, 0.25819889/(0.5*sigma_s - 1.0*sigma_t), 0.0, 0.0],
        [0.25819889/(0.5*sigma_s - 1.0*sigma_t), 0.0, -(1.0*(0.03917734*sigma_s - 0.16495722*sigma_t))/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), -(0.075*sigma_s)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2)],
        [0.0, -(1.0*(0.03917734*sigma_s - 0.16495722*sigma_t))/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), 0.0, 0.0],
        [0.0, (0.075*sigma_s)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), 0.0, 0.0]]

    systemMatricesDomain['Kyy'] = [[0.33333333/(0.5*sigma_s - 1.0*sigma_t), 0.0, -0.1490712/(0.5*sigma_s - 1.0*sigma_t), -0.25819889/(0.5*sigma_s - 1.0*sigma_t)],
        [0.0, (0.13928571*sigma_s - 0.42857143*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), 0.0, 0.0],
        [-0.1490712/(0.5*sigma_s - 1.0*sigma_t), 0.0, (0.094047619*sigma_s - 0.23809524*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), (0.03917734*sigma_s - 0.16495722*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2)],
        [-0.25819889/(0.5*sigma_s - 1.0*sigma_t), 0.0, (0.03917734*sigma_s - 0.16495722*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2), (0.13928571*sigma_s - 0.42857143*sigma_t)/(0.0625*sigma_s**2 - 0.625*sigma_s*sigma_t + 1.0*sigma_t**2)]]

    systemMatricesDomain['S'] = [[-1.0*sigma_a, 0.0, 0.0, 0.0],
        [0.0, - 1.0*sigma_a - 0.75*sigma_s, 0.0, 0.0],
        [0.0, 0.0, - 1.0*sigma_a - 0.75*sigma_s, 0.0],
        [0.0, 0.0, 0.0, - 1.0*sigma_a - 0.75*sigma_s]]

    return systemMatricesDomain