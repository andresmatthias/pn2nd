import fenics as fe 
import numpy as np


def getSystemMatricesDomain(par):

#-------- parameter renaming --------
    sigma_a = par['sigma_a']
    sigma_s = par['sigma_s']
    sigma_t = fe.Expression('sigma_a + sigma_s', degree=1, sigma_a=sigma_a, sigma_s=sigma_s)

#-------- system matrices --------

    systemMatricesDomain = dict()
    systemMatricesDomain['Kxx'] = [[-0.33333333/sigma_t, 0.0, 0.1490712/sigma_t, -0.25819889/sigma_t, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, -0.42857143/sigma_t, 0.0, 0.0, -0.21821789/sigma_t, 0.08247861/sigma_t, 0.0, 0.0, 0.0],
        [0.1490712/sigma_t, 0.0, -0.23809524/sigma_t, 0.16495722/sigma_t, 0.0, 0.0, 0.12777531/sigma_t, -0.14285714/sigma_t, 0.0],
        [-0.25819889/sigma_t, 0.0, 0.16495722/sigma_t, -0.42857143/sigma_t, 0.0, 0.0, -0.036885556/sigma_t, 0.08247861/sigma_t, -0.21821789/sigma_t],
        [0.0, -0.21821789/sigma_t, 0.0, 0.0, -0.45454545/sigma_t, 0.068720813/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.08247861/sigma_t, 0.0, 0.0, 0.068720813/sigma_t, -0.2987013/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.12777531/sigma_t, -0.036885556/sigma_t, 0.0, 0.0, -0.24675325/sigma_t, 0.17423906/sigma_t, 0.0],
        [0.0, 0.0, -0.14285714/sigma_t, 0.08247861/sigma_t, 0.0, 0.0, 0.17423906/sigma_t, -0.2987013/sigma_t, 0.068720813/sigma_t],
        [0.0, 0.0, 0.0, -0.21821789/sigma_t, 0.0, 0.0, 0.0, 0.068720813/sigma_t, -0.45454545/sigma_t]]

    systemMatricesDomain['Kxy'] = [[0.0, -0.25819889/sigma_t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.25819889/sigma_t, 0.0, 0.16495722/sigma_t, 0.0, 0.0, 0.0, -0.036885556/sigma_t, 0.0, 0.21821789/sigma_t],
        [0.0, 0.16495722/sigma_t, 0.0, 0.0, 0.0, -0.14285714/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -0.21821789/sigma_t, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -0.21821789/sigma_t, 0.0, 0.0, 0.0, 0.068720813/sigma_t, 0.0],
        [0.0, 0.0, -0.14285714/sigma_t, 0.0, 0.0, 0.0, 0.17423906/sigma_t, 0.0, -0.068720813/sigma_t],
        [0.0, -0.036885556/sigma_t, 0.0, 0.0, 0.0, 0.17423906/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.068720813/sigma_t, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.21821789/sigma_t, 0.0, 0.0, 0.0, -0.068720813/sigma_t, 0.0, 0.0, 0.0]]

    systemMatricesDomain['Kyx'] = [[0.0, -0.25819889/sigma_t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.25819889/sigma_t, 0.0, 0.16495722/sigma_t, 0.0, 0.0, 0.0, -0.036885556/sigma_t, 0.0, 0.21821789/sigma_t],
        [0.0, 0.16495722/sigma_t, 0.0, 0.0, 0.0, -0.14285714/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -0.21821789/sigma_t, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -0.21821789/sigma_t, 0.0, 0.0, 0.0, 0.068720813/sigma_t, 0.0],
        [0.0, 0.0, -0.14285714/sigma_t, 0.0, 0.0, 0.0, 0.17423906/sigma_t, 0.0, -0.068720813/sigma_t],
        [0.0, -0.036885556/sigma_t, 0.0, 0.0, 0.0, 0.17423906/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.068720813/sigma_t, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.21821789/sigma_t, 0.0, 0.0, 0.0, -0.068720813/sigma_t, 0.0, 0.0, 0.0]]

    systemMatricesDomain['Kyy'] = [[-0.33333333/sigma_t, 0.0, 0.1490712/sigma_t, 0.25819889/sigma_t, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, -0.42857143/sigma_t, 0.0, 0.0, 0.21821789/sigma_t, 0.08247861/sigma_t, 0.0, 0.0, 0.0],
        [0.1490712/sigma_t, 0.0, -0.23809524/sigma_t, -0.16495722/sigma_t, 0.0, 0.0, 0.12777531/sigma_t, 0.14285714/sigma_t, 0.0],
        [0.25819889/sigma_t, 0.0, -0.16495722/sigma_t, -0.42857143/sigma_t, 0.0, 0.0, 0.036885556/sigma_t, 0.08247861/sigma_t, 0.21821789/sigma_t],
        [0.0, 0.21821789/sigma_t, 0.0, 0.0, -0.45454545/sigma_t, -0.068720813/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.08247861/sigma_t, 0.0, 0.0, -0.068720813/sigma_t, -0.2987013/sigma_t, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.12777531/sigma_t, 0.036885556/sigma_t, 0.0, 0.0, -0.24675325/sigma_t, -0.17423906/sigma_t, 0.0],
        [0.0, 0.0, 0.14285714/sigma_t, 0.08247861/sigma_t, 0.0, 0.0, -0.17423906/sigma_t, -0.2987013/sigma_t, -0.068720813/sigma_t],
        [0.0, 0.0, 0.0, 0.21821789/sigma_t, 0.0, 0.0, 0.0, -0.068720813/sigma_t, -0.45454545/sigma_t]]

    systemMatricesDomain['S'] = [[-1.0*sigma_a, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, - 1.0*sigma_a - 1.0*sigma_s]]

    return systemMatricesDomain