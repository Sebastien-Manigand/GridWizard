#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 09:01:50 2020

@author: sebastien
"""

import numpy as np
import matplotlib.pyplot as plt
from GridWizard import GridWizard
#exec(open('GridWizard.py').read())

def precondCGA(A, b, tol=1e-1, KMAX = 10000):
    
    x0 = np.zeros_like(b)
    r_j = b - A.dot(x0)
    p_k = np.copy(r_j)
    w = A.dot(p_k)
    a_k = (r_j.T.dot(r_j)) /(p_k.T.dot(w))
    x_k = x0 + a_k*p_k
    r_k = r_j - a_k*w
    k = 1
    r_k_norm_ = (np.linalg.norm(r_k))
    dr_norm = r_k_norm_
    while(dr_norm > tol and k < KMAX):
        
        r_k_norm = r_k_norm_
        b_k = (r_k.T.dot(r_k)) / (r_j.T.dot(r_j))
        p_l = r_k + b_k*p_k
        w = A.dot(p_l)
        a_l = (r_k.T.dot(r_k)) / (p_l.T.dot(w))
        x_l = x_k + a_l*p_l
        r_l = r_k - a_l*w
        k += 1
        dr_norm = abs(np.linalg.norm(x_l) - np.linalg.norm(x_k))
#        print(dr_norm)
        r_j = np.copy(r_k)
#        p_j = numpy.copy(p_k)
        r_k = np.copy(r_l)
        p_k = np.copy(p_l)
        x_k = np.copy(x_l)
        r_k_norm_ = (np.linalg.norm(r_k))
#        dr_norm = abs(r_k_norm_ - r_k_norm)
        
    if k < KMAX:
        print('Solution found after {0} iterations !'.format(k))
    else:
        print('Solution not found after {0} iterations...'.format(KMAX))
        
    return x_k
        
    


XMIN, XMAX = -5, 5
YMIN, YMAX = -5, 5
LEFT_T, RIGHT_T, UP_T, BOTTOM_T = 50, 50, 50, 50


grid = GridWizard(domain=[-5, 5,  -5, 5], NX = 15, NY = 15, thickness=1.0)
#grid.addObstacle(nodetype='NACA4', nodes='4415', size=4.0, rotate=0.0, refining=4, diffuse=2, faceType=grid.FACE_DIRICHLET, faceParam=250 + 273.15)
grid.addObstacle(nodetype='NACA4', nodes='0015', size=6.0, rotate=-18.0, refining=3, diffuse=2, faceType=grid.FACE_DIRICHLET, faceParam=200 + 273.15)
for i in range(len(grid.cells)):
    for j in range(len(grid.cells[i]['faces'])):
        if grid.cells[i]['faces'][j][1][0] == XMIN:
            grid.cells[i]['faces'][j][4] = grid.FACE_DIRICHLET
            grid.cells[i]['faces'][j][6] = LEFT_T + 273.15 # K
        elif grid.cells[i]['faces'][j][1][0] == XMAX:
            grid.cells[i]['faces'][j][4] = grid.FACE_DIRICHLET
            grid.cells[i]['faces'][j][6] = RIGHT_T + 273.15 # K
        elif grid.cells[i]['faces'][j][1][1] == YMIN:
            grid.cells[i]['faces'][j][4] = grid.FACE_DIRICHLET
            grid.cells[i]['faces'][j][6] = BOTTOM_T + 273.15 # K
        elif grid.cells[i]['faces'][j][1][1] == XMAX:
            grid.cells[i]['faces'][j][4] = grid.FACE_DIRICHLET
            grid.cells[i]['faces'][j][6] = UP_T + 273.15 # K
            
            
#grid.plotMPL(fsize=(8,8), plotObs = True, plotVectors=False, plotCellNum=False)


# Heat flux from the wall
qwall = 0

# Heat source per volume
Svol = 0.0 # W/m3

# Thermal Conductivity  (50 for stell)
k_cond = 500  #W/mk




print("- - - - - - - - - - - - - - -")
print("  let's the magic happen!")
print("- - - - - - - - - - - - - - -")
A = np.zeros((len(grid.cells), len(grid.cells)))
B = np.zeros(len(grid.cells))
T = np.zeros(len(grid.cells))

print("Fill A matrix and B vector...")
for i in range(len(grid.cells)):
    for j in range(len(grid.cells[i]['faces'])):
        if grid.cells[i]['faces'][j][4] == grid.FACE_INTERN and not grid.cells[i]['faces'][j][5] is None:
            ii = grid.cells[i]['faces'][j][5]
            d = np.sqrt((grid.cells[i]['center'][0] - grid.cells[ii]['center'][0])**2 + (grid.cells[i]['center'][1] - grid.cells[ii]['center'][1])**2)#(cell[i][j][7][0] - cell[ii][jj][7][0])**2 + (cell[i][j][7][1] - cell[ii][jj][7][1])**2)
            area = grid.cells[i]['faces'][j][2] * grid.cells[i]['thickness']
            A[i][i] += k_cond/d * area
            A[i][ii] += -k_cond/d * area
        elif grid.cells[i]['faces'][j][4] == grid.FACE_DIRICHLET:
            d = np.sqrt((grid.cells[i]['center'][0] - grid.cells[i]['faces'][j][1][0])**2 + (grid.cells[i]['center'][1] - grid.cells[i]['faces'][j][1][1])**2)#(cell[i][j][7][0] - cell[ii][jj][7][0])**2 + (cell[i][j][7][1] - cell[ii][jj][7][1])**2)
            area = grid.cells[i]['faces'][j][2] * grid.cells[i]['thickness']
            A[i][i] += 2*k_cond/d * area
            B[i] += 2*k_cond/d * area * grid.cells[i]['faces'][j][6]


print("Solving...")
try:
    Tvector = np.linalg.solve(A, B)
except:
    print("LinAlg failed, trying precondCGA...")
    Tvector = precondCGA(A, B, KMAX = int(1e5))



print('Plotting the solution...')
grid.plotMPL(fsize=(12,12), plotVectors=False, plotCellNum=False, plotObs=False, plotCellColor=Tvector,
             figName='simple-diffusion.pdf')




