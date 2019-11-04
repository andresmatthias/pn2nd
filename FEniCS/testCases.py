"""Define various test cases for numerical experiments.

%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
"""

import fenics as fe
import numpy as np
import scipy.io as sio

def boundarySource0(vx, vy, vz):
    if type(vx) is np.ndarray:
        return np.zeros(vx.shape)
    else:
        return 0.0

def boundarySource1(vx, vy, vz):
    if type(vx) is np.ndarray:
        return 1.0 / 4.0 / np.pi * np.ones(vx.shape)
    else:
        return 1.0 / 4.0 / np.pi

def boundarySource2(vx, vy, vz):
    return vz ** 2.0

def boundarySource3(vx, vy, vz):
    return vz + 2.0

def boundarySource4(vx, vy, vz):
    return vz + 1.0

def boundarySource5(vx, vy, vz):
    return vx

def boundarySource6(vx, vy, vz):
    if type(vx) is np.ndarray:
        return 1.0 * np.ones(vx.shape)
    else:
        return 1.0
    

# -----------------------------------------------------------------------------
    
def loadTestCase1():
    """Define test case 1.
    Boundary information (rho, S) ordered according to left (0) / right (1) boundary of mesh.
    """
        
    p = dict()
    p['spatialDimension'] = 1
    p['kernelName'] = 'isotropic'
    p['domain'] = {'zMin': 0.0, 'zMax': 1.0}
    p['rho'] = [0.0, 0.0]
    p['sigma_a'] = fe.Expression('1.0', degree=1)
    p['sigma_s'] = fe.Expression('0.0', degree=1)
    p['boundarySource'] = [boundarySource1,
                           boundarySource0]
    
    return p    
  
# -----------------------------------------------------------------------------

def loadTestCase2():
    """Define test case 2.
    Boundary information (rho, S) ordered according to left (0) / right (1) boundary of mesh.
    """
        
    p = dict()
    p['spatialDimension'] = 1
    p['kernelName'] = 'isotropic'
    p['domain'] = {'zMin': 0.0, 'zMax': 1.0}
    p['rho'] = [0.5, 0.5]
    p['sigma_a'] = fe.Expression('(2 + sin(2 * PI * x[0])) / 10', degree=1, PI=np.pi)
    p['sigma_s'] = fe.Expression('(3 - pow(x[0], 2)) / 10', degree=1)
    p['boundarySource'] = [boundarySource2,
                           boundarySource0]
    
    return p   

# -----------------------------------------------------------------------------

def loadTestCase3():
    """Define test case 3.
    Boundary information (rho, S) ordered according to left (0) / right (1) boundary of mesh.
    """
        
    p = dict()
    p['spatialDimension'] = 1
    p['kernelName'] = 'bilinear'
    p['domain'] = {'zMin': 0.0, 'zMax': 1.0}
    p['rho'] = [0.0, 0.0]
    p['sigma_a'] = fe.Expression('0.0', degree=1)
    p['sigma_s'] = fe.Expression('1 + x[0]', degree=1)
    p['boundarySource'] = [boundarySource3,
                           boundarySource4]
    
    return p

# ----------------------------------------------------------------------------- 

def loadTestCase4():
    """Define test case 4.
    Boundary information (rho, S) ordered according to boundary ID in 
    meshTestCase4 (rectangle):
        1: right
        2: top
        3: left
        4: bottom
    """
        
    p = dict()
    p['spatialDimension'] = 2
    p['kernelName'] = 'isotropic'
    p['meshPrefix'] = 'meshTestCase4'
    p['rho'] = [0.0, 0.5, 0.0, 0.5]
    f1 = fe.Expression('( x[0] <= 0.6 ) ? exp(-pow(x[0] - 0.6, 2) / 0.01) : 1.0', degree=1)
    f2 = fe.Expression('( x[0] >= 0.7 ) ? exp(-pow(x[0] - 0.7, 2) / 0.01) : 1.0', degree=1)
    f3 = fe.Expression('( x[1] >= 0.4 ) ? exp(-pow(x[1] - 0.4, 2) / 0.01) : 1.0', degree=1)
    
    p['sigma_a'] = fe.Expression('100.0 * f1 * f2 * f3', degree=1, f1=f1, f2=f2, f3=f3)
    p['sigma_s'] = fe.Expression('0.01', degree=1)
    p['boundarySource'] = [boundarySource0,
                           boundarySource0,
                           boundarySource1,
                           boundarySource0]
    
    return p  

# -----------------------------------------------------------------------------

def loadTestCase5(g):
    """Define test case 5.
    Boundary information (rho, S) ordered according to boundary ID in 
    meshTestCase5 (rectangle):
        1: right
        2: top
        3: left
        4: bottom
    """
        
    p = dict()
    p['spatialDimension'] = 2
    tmp = 'henyeyGreenstein_g%1.3f'%g
    p['kernelName'] = tmp.replace('.', ',')
    p['meshPrefix'] = 'meshTestCase5'
    p['rho'] = [0.0, 0.99, 0.0, 0.99]
    p['sigma_a'] = fe.Expression('0.0', degree=1)
    p['sigma_s'] = fe.Expression('1.0', degree=1)
    p['boundarySource'] = [boundarySource0,
                           boundarySource0,
                           boundarySource5,
                           boundarySource0]
    
    return p  

# -----------------------------------------------------------------------------

def loadTestCase6():
    """Define test case 6.
    Boundary information (rho, S) ordered according to boundary ID in 
    meshTestCase6 (quarter circle)
    """
    
    buildFolder = '../build/testCase6/'
    aux = sio.loadmat(buildFolder + 'meshTestCase6Normals.mat')
    nLeft = int(aux['geoGenParams']['nLeft']) - 1  # check .geo file for this number
    nTop = int(aux['geoGenParams']['nTop']) - 1 
    nInnerOuter = int(aux['geoGenParams']['nInnerOuter']) - 1 
        
    p = dict()
    p['spatialDimension'] = 2
    p['kernelName'] = 'isotropic'
    p['meshPrefix'] = 'meshTestCase6'
    p['sigma_a'] = fe.Expression('0.0', degree=1)
    p['sigma_s'] = fe.Expression('0.1', degree=1)
    
    # reflectivity
    p['rho'] = np.hstack((0.5 * np.ones(nInnerOuter), np.zeros(nTop), 0.5 * np.ones(nInnerOuter), np.zeros(nLeft)))
    
    # boundarySource
    p['boundarySource'] = []
    
    # outer curve
    for k in range(0, nInnerOuter):
        p['boundarySource'].append(boundarySource0)
    # top
    for k in range(0, nTop):
        p['boundarySource'].append(boundarySource0)
    # inner curve
    for k in range(0, nInnerOuter):
        p['boundarySource'].append(boundarySource0)
    # left
    for k in range(0, nLeft):
        p['boundarySource'].append(boundarySource6)
    
    
    return p  
