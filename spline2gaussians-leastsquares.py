#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:47:10 2019

@author: nvthaomy
"""
import numpy as np
from scipy.optimize import least_squares
import spline, sys, argparse, re
import matplotlib.pyplot as plt
import matplotlib
"""Fitting Gaussians to spline using least squares 
   obj: objective function, calculates the Boltzmann weighted residuals 
   constrains: 
       Gaussians with even index to be repulsive, and odd index to be attractive
       upper bound of repulsive Gaussian is maximum potential from spline, can modify in getBounds
   by default, optimize with incremental number of Gaussians. Initial guess for 1st Gaussian opt is B = max value of spline potential, K = 1,
       initial values for the remaining opt are optimized parameters from previous opt + [B=0,K=0] for the newly added Gaussian
   can also optimize in one stage by providind intial values and using -nostage flag

   Outputs:
       Gaussian parameters,energy scale and kappa, in the form: [B1, K1, B2, K2,...,Bn,Kn]   
   """
#test command
#python spline2gaussians-leastsquares.py  -k "2.7835e+02 , 3.3541e+00 , -5.8015e-01, 1.6469e-01 ,-1.1965e-01, 5.2720e-02 , -2.3451e-02, 2.6243e-03" -cut 11 -n 2

parser = argparse.ArgumentParser(description="decomposing spline into Gaussians using least squares")
parser.add_argument("-k",required = True ,type = str, help="cubic spline knots, e.g. '1,2,3' or '1 2 3'")
parser.add_argument("-cut", required = True, type = float, help = "cut off distance")
parser.add_argument("-n", default = 2, type = int, help="number of Gaussians")
parser.add_argument("-N", default = 1000, type = int, help="number of points used for fitting")
parser.add_argument("-nostage", action = 'store_true')
parser.add_argument("-x0", type = str,help="initial values for Gaussian parameters, format '1 0.5 -10  0.1'")
args = parser.parse_args() 

knots = [float(i) for i in re.split(' |,',args.k) if len(i)>0]
rcut = args.cut
n = args.n
N = args.N

def obj(x,w,rs,u_spline): 
    """Calculate Boltzmann weighted residuals"""
    n = int(len(x)/2) #number of Gaussians
    u_gauss = getUgauss(x,rs,n)
    return w*(u_gauss-u_spline)

def getUgauss(x,rs,n):
    u_gauss = np.zeros(len(rs))
    for i in range(n):
        B = x[i*2]
        K = x[i*2+1]
        u_gauss += B*np.exp(-K*rs**2)
    return u_gauss

def getUspline(knots, rcut, rs, MaxPairEnekBT = 20, kB = 1,TempSet = 1):
    """calculate spline potential and adjust the hard core region
    to match the tabulated pair potentials generated by sim 
    (PotentialTablePair class in sim/export/lammps.py)"""
    myspline = spline.Spline(rcut,knots)
    u_spline = []
    du_spline = []
    for r in rs:
        u_spline.append(myspline.Val(r))
        du_spline.append(myspline.DVal(r)) 
    u_spline = np.array(u_spline)
    #get the maximum pair energy
    MaxPairEne = kB * TempSet * MaxPairEnekBT
    #indices where energy is greater
    ind = np.where(u_spline> MaxPairEne)[0]
    if len(ind):
        #find the first index where energy is valid
        i = ind[-1] + 1
        #do a linear extrapolation in the hard core region
        u_spline[:i] = (rs[i] - rs[:i]) * -du_spline[i] + u_spline[i]
        for j in ind:
            du_spline[j] = du_spline[i]
    return u_spline,du_spline
    
def weight(rs,u_spline):
    w = np.exp(-u_spline)
    w = w/np.sum(w)
    return w

    
rs = np.linspace(0,rcut,N)
u_spline, du_spline = getUspline(knots,rcut,rs)
u_max = np.max(u_spline)

w = weight(rs,u_spline)

def getBounds(n):
    bounds = ([],[]) 
    lower_energybound = -np.inf
    upper_energybound = u_max
    for i in range(n):
        if i % 2 == 0: #bounds of repulsive Gaussian
            bounds[0].extend([0,0]) #lower bound of B and K
            bounds[1].extend([upper_energybound,np.inf]) #upper bound of B and K
        else: #bounds of attractive Gaussian
            bounds[0].extend([lower_energybound,0])
            bounds[1].extend([0,np.inf])
    return bounds
def plot(xopt,rs,n_G,u_spline):
    u_gauss = getUgauss(xopt,rs,n_G)
    if n_G == 1:
        plt.figure(figsize=[10,4])
        plt.subplot(1,2,1)
        plt.plot(rs,u_spline,label="spline",linewidth = 2)
        plt.scatter(np.linspace(0,rcut,len(knots)),knots,label = "spline knots",c='r')
        plt.ylim(min(np.min(u_spline),np.min(u_gauss))*2)
    plt.plot(rs,u_gauss,label="{}-Gaussian".format(n_G),linewidth = 1)
    plt.xlim(0,rcut)
    plt.xlabel('r')
    plt.ylabel('u(r)')
    plt.legend(loc='best')

#    plt.show()

LSQ = []
outfile = open("Spline_to_{}G.txt".format(n),'w')    
if not args.nostage:   
    for i in range(n):
        bounds = getBounds(i+1)
        if i == 0:
            x0 = np.array([u_max,1.]) #initial vals for B and kappa of first Gaussian
            sys.stdout.write('\nInitial guess for 1st Gaussian:')
            sys.stdout.write('\nB: {}, K: {}'.format(x0[0],x0[1]))
            sys.stdout.write('\nParameters from optimizing {} Gaussian:'.format(i+1))            
            outfile.write('Spline knots:\n{}\n'.format(args.k))
            outfile.write('Cut:  {}'.format(rcut))
            outfile.write('\n\nDecomposing spline into Gaussians, parameters are in the form: [B1, K1, B2, K2,...,Bn,Kn]\n')
            outfile.write('\nParameters from optimizing {} Gaussian:\n'.format(i+1))
        else:
            x0 = [p for p in xopt]
            x0.extend([0,0])
            sys.stdout.write('\nInitial guess: {}'.format(x0))
            sys.stdout.write('\nParameters from optimizing {} Gaussians:'.format(i+1))
            outfile.write('Parameters from optimizing {} Gaussians:\n'.format(i+1))
        gauss = least_squares(obj,x0, args = (w,rs,u_spline),bounds=bounds)
        xopt = gauss.x
        sys.stdout.write('\n{}'.format(xopt))
        outfile.write('{}\n'.format(xopt))
        sys.stdout.write('\nLSQ: {}\n'.format(gauss.cost))
        outfile.write('LSQ: {}\n\n'.format(gauss.cost))
        LSQ.append(gauss.cost)
        plot(xopt,rs,i+1,u_spline)
    outfile.close()
else:
    if len(args.x0) == 0:
        raise Exception('Need initial values of Gaussian parameters')
    else:
        x0 = [float(i) for i in re.split(' |,',args.x0) if len(i)>0]
        if len(x0) != 2*n:
            raise Exception('Wrong number of initial values')
    bounds = getBounds(n)
    sys.stdout.write('\nInitial guess:')
    sys.stdout.write('\n{}'.format(x0))
    gauss = least_squares(obj,x0, args = (w,rs,u_spline),bounds=bounds)
    xopt = gauss.x
    sys.stdout.write('\nParameters from optimizing {} Gaussians:'.format(n))
    sys.stdout.write('\n{}'.format(xopt))
    sys.stdout.write('\nLSQ: {}\n'.format(gauss.cost))
    LSQ.append(gauss.cost)
    plot(xopt,rs,n,u_spline)

u_gauss = getUgauss(xopt,rs,n)
#plt.figure()
plt.subplot(1,2,2)
plt.plot(rs,u_spline,label="spline",linewidth = 2)
plt.plot(rs,u_gauss,label="{}-Gaussian".format(n),linewidth = 1)
plt.scatter(np.linspace(0,rcut,len(knots)),knots,label = "spline knots",c='r')
plt.ylim(min(np.min(u_spline),np.min(u_gauss))*1.1,1)
plt.xlim(0,rcut)
plt.xlabel('r')
plt.ylabel('u(r)')
plt.legend(loc='best')
plt.subplots_adjust(wspace=0.2)
plt.savefig("Spline_to_{}G.png".format(n),dpi=500)

plt.figure()
plt.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
plt.scatter(range(1,n+1),LSQ)
plt.xlim(0,n+1)
plt.xticks(range(0,n+2,1))
plt.ylim(np.min(LSQ)*0.9,np.max(LSQ)*1.1)
plt.xlabel('Number of Gaussians')
plt.ylabel('Objective')

plt.savefig("Spline_to_{}G_LSQ.png".format(n),dpi=500)
plt.show()


