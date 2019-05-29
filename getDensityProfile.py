#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 09:36:41 2019
currently written for lammps trajectory file where box is sliced 
in the z axis
@author: nvthaomy
"""
import numpy as np
import argparse, os
import matplotlib.pyplot as plt
import matplotlib

showPlots = True
try:
    os.environ["DISPLAY"] #Detects if display is available
except KeyError:
    showPlots = False
    matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",required=True,help="traj file")
parser.add_argument("-zmin",required=True, type=float,help="min value of z")
parser.add_argument("-zmax",required=True, type=float, help="max value of z")     
parser.add_argument("-ns",required=True, type = int, help="number of slabs")
parser.add_argument("-N", required=True, type = float, help="number of atoms in molecule of interest")
parser.add_argument("-at", nargs='+',required=True, help="atom types in molecule of interest")
args = parser.parse_args()

def getindices(traj):
    trjFile = open(traj,'r')
    line =trjFile.readline()
    gotIndices = False
    while len(line):
        if "ITEM: ATOMS" in line:
            atomtype_ind = (line.split()).index('type') -2
            z_ind = (line.split()).index('z') -2
            gotIndices = True
        if gotIndices:
            line = ""
        line =trjFile.readline()
    return atomtype_ind,z_ind

z = np.linspace(args.zmin, args.zmax, args.ns+1,endpoint = True)
dz = (args.zmax -args.zmin)/float(args.ns)
atomtype_ind,z_ind = getindices(args.input)
rho = []
for i in range(len(z)-1):
    print "\nGetting density for {} < z < {}".format(z[i],z[i+1])
    local_rho = []
    frame = 0
    density = 0
    trjFile = open(args.input,'r')
    line =trjFile.readline()
    while len(line): 
        oldframe = frame
        if "ITEM: TIMESTEP" in line:
            frame  += 1
            readBox = False
            readingTraj = False
            #print "Frame: {}".format(frame)
        if "ITEM: BOX" in line:
            readBox = True
            counter = 0
        if readBox and not "ITEM: BOX" in line: #get box size in x and y axis
            if counter == 0:
                x = float(line.split()[1])-float(line.split()[0])
            elif counter == 1:
                y = float(line.split()[1])-float(line.split()[0])
            counter += 1
            if counter > 1:
                slabVol = dz * x * y
                readBox = False
        if oldframe != frame and frame !=1: #check if moving to the next frame
            density = float(density)/slabVol/args.N #get number density of molecule in each slab
            local_rho.append(density)
            density = 0
        else:
            if "ITEM: ATOMS" in line:
                readingTraj = True
            if readingTraj and not "ITEM: ATOMS" in line:
                atomtype = line.split()[atomtype_ind] 
                current_z = float(line.split()[z_ind])
                if atomtype in args.at and current_z >= z[i] and current_z <= z[i+1]:
                    density += 1
        line =trjFile.readline()
    local_rho = np.mean(np.array(local_rho)) #averaging slab density over trj frames
    rho.append(local_rho)
    
z_plot = []
for i in range(len(z)-1): #making z array for plotting 
    zavg = (z[i]+z[i+1])/2
    z_plot.append(zavg)
plt.plot(z_plot,rho)
plt.ylabel('$\\rho$')
plt.xlabel('z')
title = 'Density Profile'
plt.title(title,loc='center')
plt.savefig('DensityProfile.png',dpi=500,transparent=True)
plt.show()                    
            
                    
                
            
    

            