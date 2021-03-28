#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 20:36:41 2020

@author: sebastien

This class define a Griding tool to build a refining grid which adapt to
an obstacle in the middle. 
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.path import Path
from matplotlib.patches import PathPatch
matplotlib.rcParams['font.family'] = 'serif'


class GridWizard():
    
    def __init__(self, domain = [-5, 5, -5, 5], NX = 5, NY = 5, thickness = 1.0):
        self.FACE_INTERN =      0x0000
        self.FACE_DIRICHLET =   0x0001
        self.FACE_NEUMANN =     0x0002
        self.FACE_NULL =        0x0003
        
        self.REFINE_LVL_MAX =   5
        self.xmin, self.xmax, self.NXinit = domain[0], domain[1], NX
        self.ymin, self.ymax, self.NYinit = domain[2], domain[3], NY
        X, Y = np.meshgrid(np.linspace(self.xmin, self.xmax, num=self.NXinit+1), np.linspace(self.ymin, self.ymax, num=self.NYinit+1))
        
        self.cells = []
        for i in range(NX):
            for j in range(NY):
                nodes = [[X[i][j], Y[i][j]],
                         [X[i][j+1], Y[i][j+1]],
                         [X[i+1][j+1], Y[i+1][j+1]],
                         [X[i+1][j], Y[i+1][j]]
                         ]
                self.addCell(thickness = thickness, nodes = nodes)
        self.checkNeighbors()
        self.obs = None
        
#        self.refineCell(24, refineLimit=1)
        return
    
    
    def addCell(self, surface = 1.0, thickness = 1.0, 
                nodes=[[-0.5, -0.5], [-0.5, 0.5], [0.5, 0.5], [0.5, -0.5]]):
        self.cells.append({
                'nodes': nodes.copy(),
                'center': self.centers_nodes(nodes).tolist(), #[ (nodes[0][0] + nodes[1][0] + nodes[2][0] + nodes[3][0])/4, (nodes[0][1] + nodes[1][1] + nodes[2][1] + nodes[3][1])/4 ],
                'surface': surface,
                'thickness': thickness,
                'volume': surface * thickness,
                'faces': [ ],
                'path': Path(nodes + [nodes[0]], closed=True),
                'refine.lvl': 0
                })
        for i in range(len(nodes)):
            c = [(nodes[i][0] + nodes[(i+1)%len(nodes)][0])/2, (nodes[i][1] + nodes[(i+1)%len(nodes)][1])/2]
            d = np.sqrt((nodes[(i+1)%len(nodes)][0] - nodes[i][0])**2 + (nodes[(i+1)%len(nodes)][1] - nodes[i][1])**2)
            n = [(nodes[(i+1)%len(nodes)][1] - nodes[i][1]) / d, -(nodes[(i+1)%len(nodes)][0] - nodes[i][0]) / d]
            self.cells[-1]['faces'].append([[[nodes[i][0], nodes[(i+1)%len(nodes)][0]], [nodes[i][1], nodes[(i+1)%len(nodes)][1]]], # face-nodes [[x0, x1], [y0, y1]]
                                            c.copy(), # center coordinate [xc, yc]
                                            d, # face length
                                            n.copy(), # facenormal vector [u, v]
                                            self.FACE_INTERN, # face type
                                            None, # indice of the neighbor cell
                                            0, # Dirichlet wall parameter
                                            ])
    
    
    def refineCells(self, indices, diffuse=0, refineLimit=9):
        for i in range(len(indices)):
            if self.cells[indices[0]]['refine.lvl'] != self.cells[indices[i]]['refine.lvl']:
                print('ERROR: The cells to refine are not at the same refine.lvl... ')
                exit()
        newRefineLvl = self.cells[indices[0]]['refine.lvl'] + 1
        indicesList = indices.copy()
        for d in range(diffuse):
            toAdd = []
            for i in indicesList:
                for j in range(len(self.cells[i]['faces'])):
                    if not self.cells[i]['faces'][j][5] is None:
                        if (not( self.cells[i]['faces'][j][5] in indicesList) and
                            not( self.cells[i]['faces'][j][5] in toAdd) and
                            self.cells[self.cells[i]['faces'][j][5]]['refine.lvl'] == self.cells[i]['refine.lvl']):
                            toAdd.append(self.cells[i]['faces'][j][5])
            indicesList = indicesList + toAdd
        
        toFace = []
        for i in indicesList:
            for j in range(len(self.cells[i]['faces'])):
                if not self.cells[i]['faces'][j][5] is None:
                    if (not( self.cells[i]['faces'][j][5] in indicesList) and
                        not( self.cells[i]['faces'][j][5] in toFace) and
                        self.cells[self.cells[i]['faces'][j][5]]['refine.lvl'] == self.cells[i]['refine.lvl']):
                        toFace.append(self.cells[i]['faces'][j][5])
        print(toFace)
        
        for i in indicesList:
            self.refineCell(i, refineLimit = refineLimit)
            
        for i in range(len(toFace)):
            toRem = []
            for j in range(len(self.cells[toFace[i]]['faces'])):
                if self.cells[toFace[i]]['faces'][j][5] is None:
                    pass
                elif self.cells[toFace[i]]['faces'][j][5] in indicesList:
                    x0, x1, x2 = self.cells[toFace[i]]['faces'][j][0][0][0], self.cells[toFace[i]]['faces'][j][1][0], self.cells[toFace[i]]['faces'][j][0][0][1]
                    y0, y1, y2 = self.cells[toFace[i]]['faces'][j][0][1][0], self.cells[toFace[i]]['faces'][j][1][1], self.cells[toFace[i]]['faces'][j][0][1][1]
                    self.cells[toFace[i]]['faces'].append([
                                            [[x0, x1], [y0, y1]], # face-nodes [[x0, x1], [y0, y1]]
                                            [(x0 + x1)/2, (y0 + y1)/2], # center coordinate [xc, yc]
                                            self.cells[toFace[i]]['faces'][j][2]/2, # face length
                                            self.cells[toFace[i]]['faces'][j][3].copy(), # facenormal vector [u, v]
                                            self.FACE_INTERN, # face type
                                            None, # indice of the neighbor cell
                                            0, # Dirichlet wall parameter
                                            ])
                    self.cells[toFace[i]]['faces'].append([
                                            [[x1, x2], [y1, y2]], # face-nodes [[x0, x1], [y0, y1]]
                                            [(x1 + x2)/2, (y1 + y2)/2], # center coordinate [xc, yc]
                                            self.cells[toFace[i]]['faces'][j][2]/2, # face length
                                            self.cells[toFace[i]]['faces'][j][3].copy(), # facenormal vector [u, v]
                                            self.FACE_INTERN, # face type
                                            None, # indice of the neighbor cell
                                            0, # Dirichlet wall parameter
                                            ])
                    toRem.append(j)
            for j_ in range(len(toRem)):
                for i_ in range(len(toRem)):
                    if toRem[i_] == max(toRem):
                        self.cells[toFace[i]]['faces'].pop(toRem[i_])
                        toRem.pop(i_)
                        break
            
        for j in range(len(indicesList)):
            for i in range(len(indicesList)):
                if indicesList[i] == max(indicesList):
                    self.cells.pop(indicesList[i])
                    indicesList.pop(i)
                    break
        
        toNei = []
        for i in range(len(self.cells)):
            if self.cells[i]['refine.lvl'] == newRefineLvl:
                toNei.append(i)
        self.checkNeighbors(cellList=toNei, printout=False)
    
    
    def refineCell(self, index, refineLimit = 9):
        if self.cells[index]['refine.lvl'] < min(refineLimit, self.REFINE_LVL_MAX):
            idNodes = [0, 1, 2, 3]
            for i in idNodes:
                x0, x1 = min(self.cells[index]['nodes'][i][0], self.cells[index]['center'][0]), max(self.cells[index]['nodes'][i][0], self.cells[index]['center'][0])
                y0, y1 = min(self.cells[index]['nodes'][i][1], self.cells[index]['center'][1]), max(self.cells[index]['nodes'][i][1], self.cells[index]['center'][1])
                self.addCell(surface = self.cells[index]['surface']/4, thickness = self.cells[index]['thickness'],
                             nodes = [[x0, y0],
                                      [x1, y0],
                                      [x1, y1],
                                      [x0, y1]
                                      ])
                self.cells[-1]['refine.lvl'] = self.cells[index]['refine.lvl'] + 1
            
            
    
    
    
    def checkNeighbors(self, cellList='all', printout=False):
        if cellList is 'reset':
            print("Re-setting the neighbors...")
            for i in range(len(self.cells)):
                for j in range(len(self.cells[i]['faces'])):
                    self.cells[i]['faces'][j][5] = None
            cellList = 'all'
            
        if cellList is 'all':
            for i in range(len(self.cells)-1):
                for k in range(len(self.cells[i]['faces'])):
                    for j in range(i, len(self.cells)):
                        if self.cells[i]['faces'][k][5] is None:
                            if self.cells[j]['path'].contains_point([self.cells[i]['faces'][k][1][0] + 0.0001*self.cells[i]['faces'][k][3][0], self.cells[i]['faces'][k][1][1] + 0.0001*self.cells[i]['faces'][k][3][1]]):
                                self.cells[i]['faces'][k][5] = j
                                for l in range(len(self.cells[j]['faces'])):
                                    if (self.cells[i]['faces'][k][1][0] == self.cells[j]['faces'][l][1][0] and 
                                        self.cells[i]['faces'][k][1][1] == self.cells[j]['faces'][l][1][1]):
                                        self.cells[j]['faces'][l][5] = i
                        else:
                            pass
            
        else:
            for i in cellList:
                for k in range(len(self.cells[i]['faces'])):
                    for j in range(len(self.cells)):
                        if i != j:
                            if self.cells[i]['faces'][k][5] is None:
                                if self.cells[j]['path'].contains_point([self.cells[i]['faces'][k][1][0] + 0.0001*self.cells[i]['faces'][k][3][0], self.cells[i]['faces'][k][1][1] + 0.0001*self.cells[i]['faces'][k][3][1]]):
                                    self.cells[i]['faces'][k][5] = j
                                    for l in range(len(self.cells[j]['faces'])):
                                        if (self.cells[i]['faces'][k][1][0] == self.cells[j]['faces'][l][1][0] and 
                                            self.cells[i]['faces'][k][1][1] == self.cells[j]['faces'][l][1][1]):
                                            self.cells[j]['faces'][l][5] = i
                            else:
                                pass
        if printout:
            for i in range(len(self.cells)):
                s = 'cell {0:3}: ['.format(i)
                for j in range(len(self.cells[i]['faces'])):
                    if not self.cells[i]['faces'][j][5] is None:
                        s += '{0:5}'.format(self.cells[i]['faces'][j][5])
                s += ']'
                print(s)
                
                
    
    def naca4_profile(self, naca, N = 100, center = [0, 0], length = 1.0, rotate=0.0):
        # definition of the edges
        M = int(naca[0])
        P = int(naca[1])
        XX = int(naca[2:])
        t = XX/100.0
        m = M/100.0
        p = P/10.0
        lX = np.linspace(0, 1, num=int(N/2))
        lX = lX**2 * length
        lY = t*length/0.2*(0.2969*np.sqrt(lX/length) - 0.1260*(lX/length) - 0.3516*(lX/length)**2 + 0.2843*(lX/length)**3 - 0.1015*(lX/length)**4)
        
        def cambrure(x, l):
            if x < p*l:
                yc = l * m*x/l/p**2 * (2*p - x/l)
            else:
                yc = m*(l-x)/(1-p)**2 * (1 + x/l - 2*p)
            return yc
        cambrure = np.vectorize(cambrure)
        cY = cambrure(lX, length)
        
        def dcambrure(x, l):
            if x < p*l:
                yc = l * 2*m/p**2 * (p - x/l)
            else:
                yc = l * 2*m/(1-p)**2 * (p - x/l)
            return yc
        dcambrure = np.vectorize(dcambrure)
        dcY = dcambrure(lX, length)
        theta = np.arctan(dcY)
            
        lX = lX - length/2*np.ones_like(lX)
              
        nodes = []
        for i in range(len(lX)-1): # extrados
            nodes.append([lX[i] - lY[i]*np.sin(theta[i]), cY[i] + lY[i]*np.cos(theta[i])])
        nodes.append([length/2, 0])
        for i in range(1,len(lX)-1): # intrados
            nodes.append([lX[len(lX)-1 - i] + lY[len(lX)-1 - i]*np.sin(theta[len(lX)-1 - i]), cY[len(lX)-1 - i] - lY[len(lX)-1 - i]*np.cos(theta[len(lX)-1 - i])])
           
        for i in range(len(nodes)):
            x_ = (nodes[i][0])*np.cos(rotate*np.pi/180) - (nodes[i][1])*np.sin(rotate*np.pi/180) + center[0]
            y_ = (nodes[i][0])*np.sin(rotate*np.pi/180) + (nodes[i][1])*np.cos(rotate*np.pi/180) + center[1]
            nodes[i] = [x_, y_]
        return nodes
                
                
    
    def addObstacle(self, nodetype='nodes', nodes = None, N = 100, center=[0, 0], size=1.0, rotate=0.0, refining=2, diffuse=2, faceType=0x0003, faceParam=0, reshaping=True):
        if nodetype == 'nodes' and not(nodes is None):
            self.obs = {'nodes': nodes.copy(),
                        'path': Path(nodes + [nodes[0]], closed=True),
                        'faceType': self.FACE_DIRICHLET}
        elif nodetype == 'test':
            R = size/2
            t = np.linspace(0, 2*np.pi, num=N)
            x = R * np.cos(t - rotate/180.0*np.pi)
            y = R * np.sin(t - rotate/180.0*np.pi)
            P = np.array([x, y]).T
            self.obs = {'nodes': P[:-1].tolist(),
                        'path': Path(P.tolist(), closed=True),
                        'faceType': self.FACE_DIRICHLET}
        elif nodetype == 'NACA4' and not(nodes is None):
            P = self.naca4_profile(naca = nodes, N = N, center=center, length=size, rotate=rotate)
#            print(np.array(P).shape)
            self.obs = {'nodes': P,
                        'path': Path(P + [P[0]], closed=True),
                        'faceType': self.FACE_DIRICHLET}
        
        for r in range(refining):
            toRefine = []
            for i in range(len(self.cells)):
                if self.obs['path'].intersects_path(self.cells[i]['path'], filled=False):
                    toRefine.append(i)
            self.refineCells(toRefine, diffuse=diffuse, refineLimit=r+1)
        
        toRem = []
        for i in range(len(self.cells)):
            if( not self.obs['path'].intersects_path(self.cells[i]['path'], filled=False) and
               self.obs['path'].contains_path(self.cells[i]['path'])):
                toRem.append(i)
        for j in range(len(toRem)):
            for i in range(len(toRem)):
                if toRem[i] == max(toRem):
                    self.cells.pop(toRem[i])
                    toRem.pop(i)
                    break
        
        if reshaping:
            self.checkNeighbors(cellList='reset', printout=False)
            
            self.inters = [[], []]
            toCut = []
            for i in range(len(self.cells)):
                if(self.obs['path'].intersects_path(self.cells[i]['path'], filled=False)):
                    toCut.append(i)
            for i in toCut:
                print('Reshapping the cell {0}...'.format(i))
                faceToRem = []
                newNodes = []
                for j in range(len(self.cells[i]['faces'])):
                    faceNodes = [[self.cells[i]['faces'][j][0][0][0], self.cells[i]['faces'][j][0][1][0]], [self.cells[i]['faces'][j][0][0][1], self.cells[i]['faces'][j][0][1][1]]]
                    if Path(faceNodes, closed=False).intersects_path(self.obs['path']):
                        faceToRem.append(j)
                        xI, yI = 0, 0
                        for k in range(len(self.obs['nodes'])):
                            obsNodes = [self.obs['nodes'][k], self.obs['nodes'][(k+1)%len(self.obs['nodes'])]]
                            if(Path(faceNodes + [faceNodes[0]], closed=True).intersects_path(Path(obsNodes, closed=False)) or
                               Path([[faceNodes[0][0] - 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[0][1] - 0.0001*self.cells[i]['faces'][j][3][1]],
                                     [faceNodes[1][0] - 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[1][1] - 0.0001*self.cells[i]['faces'][j][3][1]],
                                     [faceNodes[1][0] + 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[1][1] + 0.0001*self.cells[i]['faces'][j][3][1]],
                                     [faceNodes[0][0] + 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[0][1] + 0.0001*self.cells[i]['faces'][j][3][1]],
                                     [faceNodes[0][0] - 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[0][1] - 0.0001*self.cells[i]['faces'][j][3][1]]], closed=True).contains_point(obsNodes[0], transform=None, radius=0.0001)):# or
    #                           Path([[faceNodes[0][0] - 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[0][1] - 0.0001*self.cells[i]['faces'][j][3][1]],
    #                                 [faceNodes[1][0] - 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[1][1] - 0.0001*self.cells[i]['faces'][j][3][1]],
    #                                 [faceNodes[1][0] + 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[1][1] + 0.0001*self.cells[i]['faces'][j][3][1]],
    #                                 [faceNodes[0][0] + 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[0][1] + 0.0001*self.cells[i]['faces'][j][3][1]],
    #                                 [faceNodes[0][0] - 0.0001*self.cells[i]['faces'][j][3][0], faceNodes[0][1] - 0.0001*self.cells[i]['faces'][j][3][1]]], closed=True).contains_point(obsNodes[1], transform=None, radius=0.0001)):
                                xI, yI = self.intersect_segments(faceNodes, obsNodes)
                                self.inters[0].append(xI)
                                self.inters[1].append(yI)
                                newNodes.append([xI, yI])
                                
    #                            if self.obs['path'].contains_point([self.cells[i]['faces'][j][0][0][0], self.cells[i]['faces'][j][0][1][0]]):
                                if(self.obs['path'].contains_point([xI + 0.0001*(self.cells[i]['faces'][j][0][0][0] - xI), yI + 0.0001*(self.cells[i]['faces'][j][0][1][0] - yI)])):
                                    x0 = self.cells[i]['faces'][j][0][0][1]
                                    y0 = self.cells[i]['faces'][j][0][1][1]
                                else:
                                    x0 = self.cells[i]['faces'][j][0][0][0]
                                    y0 = self.cells[i]['faces'][j][0][1][0]
                                self.cells[i]['faces'].append([
                                                        [[x0, xI], [y0, yI]], # face-nodes [[x0, x1], [y0, y1]]
                                                        [(x0 + xI)/2, (y0 + yI)/2], # center coordinate [xc, yc]
                                                        np.sqrt((xI-x0)**2 + (yI-y0)**2), # face length
                                                        self.cells[i]['faces'][j][3].copy(), # facenormal vector [u, v]
                                                        self.cells[i]['faces'][j][4], # face type
                                                        self.cells[i]['faces'][j][5], # indice of the neighbor cell
                                                        self.cells[i]['faces'][j][6], # Dirichlet wall parameter
                                                        ])
    #                            break
                    elif self.obs['path'].contains_path(Path(faceNodes, closed=False)):
                        faceToRem.append(j)
                        
                if len(newNodes) == 2:
                    c = [(newNodes[0][0] + newNodes[1][0])/2, (newNodes[0][1] + newNodes[1][1])/2]
                    d = np.sqrt((newNodes[1][0] - newNodes[0][0])**2 + (newNodes[1][1] - newNodes[0][1])**2)
                    n = [(newNodes[1][1] - newNodes[0][1]) / d, -(newNodes[1][0] - newNodes[0][0]) / d]
                    self.cells[i]['faces'].append([[[newNodes[0][0], newNodes[1][0]], [newNodes[0][1], newNodes[1][1]]], # face-nodes [[x0, x1], [y0, y1]]
                                                    c.copy(), # center coordinate [xc, yc]
                                                    d, # face length
                                                    n.copy(), # facenormal vector [u, v]
                                                    faceType, # face type
                                                    None, # indice of the neighbor cell
                                                    faceParam, # Dirichlet wall parameter
                                                    ])
                    cellNodes = newNodes.copy()
                    for j in range(len(self.cells[i]['nodes'])):
                        if not self.obs['path'].contains_point(self.cells[i]['nodes'][j]):
                            cellNodes.append(self.cells[i]['nodes'][j])
                            
                elif len(newNodes) == 3:
                    if newNodes[0] == newNodes[1]:
                        c = [(newNodes[0][0] + newNodes[2][0])/2, (newNodes[0][1] + newNodes[2][1])/2]
                        d = np.sqrt((newNodes[2][0] - newNodes[0][0])**2 + (newNodes[2][1] - newNodes[0][1])**2)
                        n = [(newNodes[2][1] - newNodes[0][1]) / d, -(newNodes[2][0] - newNodes[0][0]) / d]
                        self.cells[i]['faces'].append([[[newNodes[0][0], newNodes[2][0]], [newNodes[0][1], newNodes[2][1]]], # face-nodes [[x0, x1], [y0, y1]]
                                                        c.copy(), # center coordinate [xc, yc]
                                                        d, # face length
                                                        n.copy(), # facenormal vector [u, v]
                                                        faceType, # face type
                                                        None, # indice of the neighbor cell
                                                        faceParam, # Dirichlet wall parameter
                                                        ])
                        cellNodes = [newNodes[0], newNodes[2]]
                    else:
                        c = [(newNodes[0][0] + newNodes[1][0])/2, (newNodes[0][1] + newNodes[1][1])/2]
                        d = np.sqrt((newNodes[1][0] - newNodes[0][0])**2 + (newNodes[1][1] - newNodes[0][1])**2)
                        n = [(newNodes[1][1] - newNodes[0][1]) / d, -(newNodes[1][0] - newNodes[0][0]) / d]
                        self.cells[i]['faces'].append([[[newNodes[0][0], newNodes[1][0]], [newNodes[0][1], newNodes[1][1]]], # face-nodes [[x0, x1], [y0, y1]]
                                                        c.copy(), # center coordinate [xc, yc]
                                                        d, # face length
                                                        n.copy(), # facenormal vector [u, v]
                                                        faceType, # face type
                                                        None, # indice of the neighbor cell
                                                        faceParam, # Dirichlet wall parameter
                                                        ])
                        cellNodes = [newNodes[0], newNodes[1]]
                    for j in range(len(self.cells[i]['nodes'])):
                        if not self.obs['path'].contains_point(self.cells[i]['nodes'][j]):
                            cellNodes.append(self.cells[i]['nodes'][j])
                            
                elif len(newNodes) == 4:
                    extNodes = []
                    cellNodes = []
                    cellNodes2 = []
                    newCouple = []
                    newCouple2 = []
                    whichCell = True # True for the current mitosing cell, False for the product cell
                    for j in range(len(self.cells[i]['nodes'])):
                        if not self.obs['path'].contains_point(self.cells[i]['nodes'][j]):
                            extNodes.append(j) 
                    if len(extNodes) < 2:
                        print('There is an issue with cell #{0} mitose...'.format(i))
                    cellNodes.append(self.cells[i]['nodes'][extNodes[0]])
                    for j in range(len(extNodes)-1):
                        if self.obs['path'].intersects_path(Path([self.cells[i]['nodes'][extNodes[j]], self.cells[i]['nodes'][extNodes[j+1]]], closed=False)):
                            whichCell = not whichCell
                        if whichCell:
                            cellNodes.append(self.cells[i]['nodes'][extNodes[j+1]])
                        else:
                            cellNodes2.append(self.cells[i]['nodes'][extNodes[j+1]])
    #                print(cellNodes)
    #                print(newNodes)
    #                print(cellNodes2)
                    
    #                surfCouple = [self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]),
    #                              self.surface_nodes(cellNodes + [newNodes[0] + newNodes[2]]),
    #                              self.surface_nodes(cellNodes + [newNodes[0] + newNodes[3]]),
    #                              self.surface_nodes(cellNodes + [newNodes[1] + newNodes[2]]),
    #                              self.surface_nodes(cellNodes + [newNodes[1] + newNodes[3]]),
    #                              self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]])
    #                              ]
    #                
    #                if surfCouple.index(max(surfCouple)) == 0: #self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    newCouple = [0, 1]
    #                elif surfCouple.index(max(surfCouple)) == 1: #self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    newCouple = [0, 2]
    #                elif surfCouple.index(max(surfCouple)) == 2: #self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    newCouple = [0, 3]
    #                elif surfCouple.index(max(surfCouple)) == 3: #self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    newCouple = [1, 2]
    #                elif surfCouple.index(max(surfCouple)) == 4: #self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    newCouple = [1, 3]
    #                else: #self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    newCouple = [2, 3]
    #                print(newCouple)
                    
                    for j in range(len(self.cells[i]['faces'])):
                        for k_ in range(len(cellNodes)):
                            if((self.cells[i]['faces'][j][0][0][0] == cellNodes[k_][0] and self.cells[i]['faces'][j][0][1][0] == cellNodes[k_][1]) or
                               (self.cells[i]['faces'][j][0][0][1] == cellNodes[k_][0] and self.cells[i]['faces'][j][0][1][1] == cellNodes[k_][1])):
                                for k in range(len(newNodes)):
                                    if((self.cells[i]['faces'][j][0][0][0] == newNodes[k][0] and self.cells[i]['faces'][j][0][1][0] == newNodes[k][1]) or
                                       (self.cells[i]['faces'][j][0][0][1] == newNodes[k][0] and self.cells[i]['faces'][j][0][1][1] == newNodes[k][1])):
                                        newCouple.append(k)
                                
    #                print(newCouple)
                    newCouple2 = [i for i in range(4) if i not in newCouple]
    #                print(newCouple2)
                            
    #                        if len(np.intersect1d([np.array(self.cells[i]['faces'][j][0]).T[0].tolist()], newNodes)) > 0:
    #                            newCouple.append(newNodes.index(np.array(self.cells[i]['faces'][j][0]).T[0]))
    #                        elif len(np.intersect1d([np.array(self.cells[i]['faces'][j][0]).T[1].tolist()], newNodes)) > 0:
    #                            newCouple.append(newNodes.index(np.array(self.cells[i]['faces'][j][0]).T[1]))
                                
                    
                    cellNodes = cellNodes + [newNodes[newCouple[0]], newNodes[newCouple[1]]]
                    c = [(newNodes[newCouple[0]][0] + newNodes[newCouple[1]][0])/2, (newNodes[newCouple[0]][1] + newNodes[newCouple[1]][1])/2]
                    d = np.sqrt((newNodes[newCouple[1]][0] - newNodes[newCouple[0]][0])**2 + (newNodes[newCouple[1]][1] - newNodes[newCouple[0]][1])**2)
                    n = [(newNodes[newCouple[1]][1] - newNodes[newCouple[0]][1]) / d, -(newNodes[newCouple[1]][0] - newNodes[newCouple[0]][0]) / d]
                    self.cells[i]['faces'].append([[[newNodes[newCouple[0]][0], newNodes[newCouple[1]][0]], [newNodes[newCouple[0]][1], newNodes[newCouple[1]][1]]], # face-nodes [[x0, x1], [y0, y1]]
                                                    c.copy(), # center coordinate [xc, yc]
                                                    d, # face length
                                                    n.copy(), # facenormal vector [u, v]
                                                    faceType, # face type
                                                    None, # indice of the neighbor cell
                                                    faceParam, # Dirichlet wall parameter
                                                    ])
                    
                    cellNodes2 = cellNodes2 + [newNodes[newCouple2[0]], newNodes[newCouple2[1]]]
                    c = [(newNodes[newCouple2[0]][0] + newNodes[newCouple2[1]][0])/2, (newNodes[newCouple2[0]][1] + newNodes[newCouple2[1]][1])/2]
                    d = np.sqrt((newNodes[newCouple2[1]][0] - newNodes[newCouple2[0]][0])**2 + (newNodes[newCouple2[1]][1] - newNodes[newCouple2[0]][1])**2)
                    n = [(newNodes[newCouple2[1]][1] - newNodes[newCouple2[0]][1]) / d, -(newNodes[newCouple2[1]][0] - newNodes[newCouple2[0]][0]) / d]
                    self.cells[i]['faces'].append([[[newNodes[newCouple2[0]][0], newNodes[newCouple2[1]][0]], [newNodes[newCouple2[0]][1], newNodes[newCouple2[1]][1]]], # face-nodes [[x0, x1], [y0, y1]]
                                                    c.copy(), # center coordinate [xc, yc]
                                                    d, # face length
                                                    n.copy(), # facenormal vector [u, v]
                                                    faceType, # face type
                                                    None, # indice of the neighbor cell
                                                    faceParam, # Dirichlet wall parameter
                                                    ])
                    
    #                if self.surface_nodes(cellNodes + [newNodes[0] + newNodes[1]]) < self.surface_nodes(cellNodes + [newNodes[2] + newNodes[3]]):
    #                    cellNodes = cellNodes + [newNodes[0], newNodes[1]]
    #                    c = [(newNodes[0][0] + newNodes[1][0])/2, (newNodes[0][1] + newNodes[1][1])/2]
    #                    d = np.sqrt((newNodes[1][0] - newNodes[0][0])**2 + (newNodes[1][1] - newNodes[0][1])**2)
    #                    n = [(newNodes[1][1] - newNodes[0][1]) / d, -(newNodes[1][0] - newNodes[0][0]) / d]
    #                    self.cells[i]['faces'].append([[[newNodes[0][0], newNodes[1][0]], [newNodes[0][1], newNodes[1][1]]], # face-nodes [[x0, x1], [y0, y1]]
    #                                                    c.copy(), # center coordinate [xc, yc]
    #                                                    d, # face length
    #                                                    n.copy(), # facenormal vector [u, v]
    #                                                    faceType, # face type
    #                                                    None, # indice of the neighbor cell
    #                                                    faceParam, # Dirichlet wall parameter
    #                                                    ])
    #                else:
    #                    cellNodes = cellNodes + [newNodes[2], newNodes[3]]
    #                    c = [(newNodes[2][0] + newNodes[3][0])/2, (newNodes[2][1] + newNodes[3][1])/2]
    #                    d = np.sqrt((newNodes[3][0] - newNodes[2][0])**2 + (newNodes[3][1] - newNodes[2][1])**2)
    #                    n = [(newNodes[3][1] - newNodes[2][1]) / d, -(newNodes[3][0] - newNodes[2][0]) / d]
    #                    self.cells[i]['faces'].append([[[newNodes[2][0], newNodes[3][0]], [newNodes[2][1], newNodes[3][1]]], # face-nodes [[x0, x1], [y0, y1]]
    #                                                    c.copy(), # center coordinate [xc, yc]
    #                                                    d, # face length
    #                                                    n.copy(), # facenormal vector [u, v]
    #                                                    faceType, # face type
    #                                                    None, # indice of the neighbor cell
    #                                                    faceParam, # Dirichlet wall parameter
    #                                                    ])
                        
                        
                        
                for j in range(len(faceToRem)):
                    for k in range(len(faceToRem)):
                        if faceToRem[k] == max(faceToRem):
                            self.cells[i]['faces'].pop(faceToRem[k])
                            faceToRem.pop(k)
                            break
                
                # Update the parameters of the cell
    #            cellNodes = newNodes.copy()
    #            for j in range(len(self.cells[i]['nodes'])):
    #                if not self.obs['path'].contains_point(self.cells[i]['nodes'][j]):
    #                    cellNodes.append(self.cells[i]['nodes'][j])
                angNodes = []
                for j in range(len(cellNodes)):
                    scalr = cellNodes[j][0] - self.cells[i]['center'][0]
                    vectl = self.cells[i]['center'][1] - cellNodes[j][1]
                    d = np.sqrt((cellNodes[j][0] - self.cells[i]['center'][0])**2 + (cellNodes[j][1] - self.cells[i]['center'][1])**2)
                    a = np.arccos(scalr/d)
                    if (vectl > 0): a = -a
                    angNodes.append(a)
                
                for r in range(len(angNodes)):
                    for s in range(len(angNodes)-1):
                        if angNodes[s] > angNodes[s+1]:
                            c = cellNodes[s+1]
                            a = angNodes[s+1]
                            cellNodes[s+1] = cellNodes[s]
                            angNodes[s+1] = angNodes[s]
                            cellNodes[s] = c
                            angNodes[s] = a
                            
                if len(newNodes) == 4:
                    angNodes2 = []
                    for j in range(len(cellNodes2)):
                        scalr = cellNodes2[j][0] - self.cells[i]['center'][0]
                        vectl = self.cells[i]['center'][1] - cellNodes2[j][1]
                        d = np.sqrt((cellNodes2[j][0] - self.cells[i]['center'][0])**2 + (cellNodes2[j][1] - self.cells[i]['center'][1])**2)
                        a = np.arccos(scalr/d)
                        if (vectl > 0): a = -a
                        angNodes2.append(a)
                    for r in range(len(angNodes2)):
                        for s in range(len(angNodes2)-1):
                            if angNodes2[s] > angNodes2[s+1]:
                                c = cellNodes2[s+1]
                                a = angNodes2[s+1]
                                cellNodes2[s+1] = cellNodes2[s]
                                angNodes2[s+1] = angNodes2[s]
                                cellNodes2[s] = c
                                angNodes2[s] = a
                    
                
                self.cells[i]['center'] = self.centers_nodes(cellNodes) #np.average(np.array(cellNodes), axis=0).tolist()
                self.cells[i]['nodes'] = cellNodes.copy()
                self.cells[i]['path'] = Path(cellNodes + [cellNodes[0]], closed=True)
                if self.cells[i]['path'].contains_point([self.cells[i]['faces'][-1][1][0] + 0.0001*self.cells[i]['faces'][-1][3][0], self.cells[i]['faces'][-1][1][1] + 0.0001*self.cells[i]['faces'][-1][3][1]]):
                    self.cells[i]['faces'][-1][3][0] *= -1
                    self.cells[i]['faces'][-1][3][1] *= -1
                self.cells[i]['surface'] = self.surface_nodes(self.cells[i]['nodes'])
                self.cells[i]['volume'] = self.cells[i]['surface'] * self.cells[i]['thickness']
                
                if len(newNodes) == 4:
                    faceToRem = []
                    for j in range(len(self.cells[i]['faces'])):
                        if not self.cells[i]['path'].intersects_path(Path([[self.cells[i]['faces'][j][0][0][0], self.cells[i]['faces'][j][0][1][0]],
                                                                       [self.cells[i]['faces'][j][0][0][1], self.cells[i]['faces'][j][0][1][1]]], closed=False)):
                            faceToRem.append(j)
                    self.addCell(surface = self.surface_nodes(cellNodes2), 
                                 thickness = self.cells[i]['thickness'], 
                                 nodes=cellNodes2)
                    self.cells[-1]['refine.lvl'] = self.cells[i]['refine.lvl']
                    self.cells[-1]['surface'] = self.surface_nodes(self.cells[-1]['nodes'])
                    self.cells[-1]['volume'] = self.cells[-1]['surface'] * self.cells[-1]['thickness']
                    self.cells[-1]['faces'] = []
                    for j in faceToRem:
                        self.cells[-1]['faces'].append(self.cells[i]['faces'][j].copy())
                    if self.cells[-1]['path'].contains_point([self.cells[-1]['faces'][-1][1][0] + 0.0001*self.cells[-1]['faces'][-1][3][0], self.cells[-1]['faces'][-1][1][1] + 0.0001*self.cells[-1]['faces'][-1][3][1]]):
                        self.cells[-1]['faces'][-1][3][0] *= -1
                        self.cells[-1]['faces'][-1][3][1] *= -1
                
                    for j in range(len(faceToRem)):
                        for k in range(len(faceToRem)):
                            if faceToRem[k] == max(faceToRem):
                                self.cells[i]['faces'].pop(faceToRem[k])
                                faceToRem.pop(k)
                                break
                    if self.cells[i]['path'].contains_point([self.cells[i]['faces'][-1][1][0] + 0.0001*self.cells[i]['faces'][-1][3][0], self.cells[i]['faces'][-1][1][1] + 0.0001*self.cells[i]['faces'][-1][3][1]]):
                        self.cells[i]['faces'][-1][3][0] *= -1
                        self.cells[i]['faces'][-1][3][1] *= -1
            
            

    def intersect_segments(self, seg1, seg2):
        xA, yA = seg1[0][0], seg1[0][1]
        xB, yB = seg1[1][0], seg1[1][1]
        xC, yC = seg2[0][0], seg2[0][1]
        xD, yD = seg2[1][0], seg2[1][1]
        xAB, yAB = xB - xA, yB - yA
        AB = np.sqrt(xAB**2 + yAB**2)
        AC = np.sqrt((xC-xA)**2 + (yC-yA)**2)
        AD = np.sqrt((xD-xA)**2 + (yD-yA)**2)
#        print('AB = [{0}, {1}]'.format(xAB, yAB))
        
        prod_ACAB = ((xC-xA)*xAB + (yC-yA)*yAB) / (AB**2)
        prod_ADAB = ((xD-xA)*xAB + (yD-yA)*yAB) / (AB**2)
        
        xAP, yAP = prod_ACAB*xAB, prod_ACAB*yAB
        xAQ, yAQ = prod_ADAB*xAB, prod_ADAB*yAB
        xP, yP = xAP + xA, yAP + yA
        xQ, yQ = xAQ + xA, yAQ + yA
        CP = np.sqrt((xP - xC)**2 + (yP - yC)**2)
        CDort = CP + np.sqrt((xQ - xD)**2 + (yQ - yD)**2)
        xI, yI = (xQ - xP)*(CP/CDort) + xP, (yQ - yP)*(CP/CDort) + yP
        return xI, yI
    
    def surface_3nodes(self, nodes):
        xAB, yAB = nodes[1][0]-nodes[0][0], nodes[1][1]-nodes[0][1]
        xAC, yAC = nodes[2][0]-nodes[0][0], nodes[2][1]-nodes[0][1]
        s = abs(xAB*yAC - yAB*xAC)/2
        return s
    
    def surface_nodes(self, nodes):
        s = 0
        for i in range(len(nodes)-2):
            s += self.surface_3nodes([nodes[0], nodes[i+1], nodes[i+2]])
        return s
    
    def centers_nodes(self, nodes):
        CG = np.zeros(2)
        w = 0
        for i in range(len(nodes)-2):
            s = self.surface_3nodes([nodes[0], nodes[i+1], nodes[i+2]])
            CG += np.average([nodes[0], nodes[i+1], nodes[i+2]], axis = 0) * s
            w += s
        return CG / w
    
    
    def saveGrid(self, filename):
        f = open(filename, 'w')
        f.write("# GridWizard\n")
        f.write("# - - - - - - - - - - - -\n")
        f.write("{0}".format(len(self.cells)))
        for i in range(len(self.cells)):
            f.write("\n{0:4} ".format(int(i)))
            f.write("{0:4} ".format(self.cells[i]['refine.lvl']))
            f.write("{0:8.4f} {1:8.4f} ".format(self.cells[i]['center'][0], self.cells[i]['center'][1]))
            f.write("{0:8.4f} ".format(self.cells[i]['surface']))
            f.write("{0:8.4f} ".format(self.cells[i]['thickness']))
            f.write("{0:8.4f} ".format(self.cells[i]['volume']))
            f.write("{0:3} ".format(len(self.cells[i]['faces'])))
            for j in range(len(self.cells[i]['faces'])):
                f.write("{0:8.4f} {1:8.4f} {2:8.4f} {3:8.4f} ".format(self.cells[i]['faces'][j][0][0][0], self.cells[i]['faces'][j][0][0][1], self.cells[i]['faces'][j][0][1][0], self.cells[i]['faces'][j][0][1][1]))
                f.write("{0:8.4f} {1:8.4f} ".format(self.cells[i]['faces'][j][1][0], self.cells[i]['faces'][j][1][1]))
                f.write("{0:8.4f} ".format(self.cells[i]['faces'][j][2]))
                f.write("{0:8.4f} {1:8.4f} ".format(self.cells[i]['faces'][j][3][0], self.cells[i]['faces'][j][3][1]))
                f.write("{0:04x} ".format(self.cells[i]['faces'][j][4]))
                if self.cells[i]['faces'][j][5] is None:
                    f.write("None ")
                else:
                    f.write("{0:4} ".format(self.cells[i]['faces'][j][5]))
                if type(self.cells[i]['faces'][j][6]) == int:
                    f.write("I {0:4}".format(self.cells[i]['faces'][j][6]))
                elif type(self.cells[i]['faces'][j][6]) == float:
                    f.write("F {0:8.4f}".format(self.cells[i]['faces'][j][6]))
                
        f.close()                
        
    
    def colorScale(self, xlist):
        r, g, b = 0,0,0
        color = []
        emin = min(xlist)
        emax = max(xlist)
        D = float(emax-emin)
        for l in range(len(xlist)):
            e = xlist[l]-emin
            if e < 0: 
                r=0
                g=0
                b=255
            elif e >= 0 and e < D/4.0:
                r=0
                g=int(4*e/D*255)
                b=255
            elif e >= D/4.0 and e < 2*D/4.0:
                r=0
                g=255
                b=int((1-(4*e-D)/D)*255)
            elif e >= 2*D/4.0 and e < 3*D/4.0:
                r=int((4*e-2*D)/D*255)
                g=255
                b=0
            elif e >= 3*D/4.0 and e < D:
                r=255
                g=int((1-(4*e-3*D)/D)*255)
                b=0
            else: 
                r=255
                g=0
                b=0
            color.append("#" + ("{:0>2X}{:0>2X}{:0>2X}".format(r, g, b))) 
        return color
    
    
    def plotMPL(self, fsize=(14,14), plotVectors = True, plotObs = True, plotCellNum = True, plotCellColor = None, figName = None):
        ARROW_len = 0.2
        ARROW_head = 0.15
        fig, ax = plt.subplots(figsize = fsize)
        ax.set_xlim(self.xmin - 0.5, self.xmax + 0.5)
        ax.set_ylim(self.ymin - 0.5, self.ymax + 0.5)
        
        for i in range(len(self.cells)):
#            ax.scatter(*self.cells[i]['center'], marker='+', s=64, linewidth=0.5, color='black')
            arrow_len = ARROW_len / np.power(2, self.cells[i]['refine.lvl'])
            arrow_head = ARROW_head / np.power(2, self.cells[i]['refine.lvl'])
            if not plotCellColor is None:
                if len(plotCellColor) == len(self.cells):
                    clr = self.colorScale(plotCellColor)
                    ax.add_patch(PathPatch(self.cells[i]['path'], edgecolor='black', facecolor=clr[i], fill=True, linewidth=0))
            if plotCellNum:
                ax.text(*self.cells[i]['center'], '{0}'.format(i), fontsize=max(3, 12 - int(3*self.cells[i]['refine.lvl'])), color='black')
            for j in range(len(self.cells[i]['faces'])):
                ax.plot(self.cells[i]['faces'][j][0][0],
                        self.cells[i]['faces'][j][0][1],
                        linewidth=0.5, linestyle='dashed', color='black')
                if plotVectors:
                    ax.arrow(*self.cells[i]['faces'][j][1], arrow_len*self.cells[i]['faces'][j][3][0], arrow_len*self.cells[i]['faces'][j][3][1], 
                             linewidth=0.5, head_width=0.8*arrow_head, head_length=arrow_head, facecolor='black')
        
        if (not(self.obs is None) and plotObs):
            P = np.array(self.obs['nodes'] + [self.obs['nodes'][0]]).T
            ax.plot(P[0], P[1], linewidth=1.0, color='black')
#            ax.scatter(self.inters[0], self.inters[1], marker='o', color='red', s=32)
            
        
        fig.tight_layout()
        if not figName is None:
            fig.savefig(figName)
    




if __name__ == '__main__':
    grid = GridWizard(domain = [-5, 5, -5, 5], NX = 10, NY = 10, thickness = 1.0)
#    grid.addObstacle(nodetype='test', nodes=None, size=1.0, refining=2, diffuse=1)
    grid.addObstacle(nodetype='NACA4', nodes='0025', size=6.0, rotate=-25.0, refining=3, diffuse=1, reshaping = False, faceType=grid.FACE_DIRICHLET, faceParam=273.15)
    
    grid.saveGrid('testGridWizard.dat')
    
    grid.plotMPL(plotVectors = False,
                 plotObs = False,
                 plotCellNum = True,
                 figName = 'testGridWizard.pdf')




