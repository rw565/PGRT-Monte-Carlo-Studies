#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import collections as mc
import scipy.stats as ss
import scipy
import numpy as np
import os
from pathlib import Path
import uproot
import re
import click
from matplotlib.patches import Rectangle


def hello():
    print('Hello Ryan')

def get_stat_value(s, v):
    """
    This function retrive a value from a stat file
    """
    g = r''+v+'\w+'
    a = re.search(g, s)
    if a == None:
        return -1
    a = a.group(0)[len(v):]
    return float(a)  

def get_energy(filename):
    data = open(filename,'r')
    lines = data.readlines()[9:]
    #print("lines:",lines)
    data.close()
    energy = []
    bin_width = []
    frequency = []

    for line in lines:
        p = line.split()
        if ("!" not in line) and ("#" not in line):
            energy.append(float(p[0]))
            bin_width.append(float(p[1]))
            frequency.append(float(p[2]))
    return energy, bin_width, frequency
    
def tget(t, array_name):
    """ 
    it is faster to access to root array like this dont know exactly why
    """
    return t.arrays([array_name])[array_name]


def plot_transaxial_position(ax, coinc, slice_time):
    times = tget(coinc, 'time1')
    runID = tget(coinc, 'runID')
    gpx1 = tget(coinc, 'globalPosX1')
    gpx2 = tget(coinc, 'globalPosX2')
    gpy1 = tget(coinc, 'globalPosY1')
    gpy2 = tget(coinc, 'globalPosY2')
    # only consider coincidences  with time lower than time_slice
    # (assuming 2 time slices only)
    mask = (times < slice_time)
    n = 1000 # restrict to the n first values
    r0_gpx1 = gpx1[mask][:n]
    r0_gpx2 = gpx2[mask][:n]
    r0_gpy1 = gpy1[mask][:n]
    r0_gpy2 = gpy2[mask][:n]
    r0x = np.concatenate((r0_gpx1,r0_gpx2, r0_gpx1))
    r0y = np.concatenate((r0_gpy1,r0_gpy2, r0_gpy1))
    ax.scatter(r0x, r0y, s=1)
    mask = (times > slice_time)
    r1_gpx1 = gpx1[mask][:n]
    r1_gpx2 = gpx2[mask][:n]
    r1_gpy1 = gpy1[mask][:n]
    r1_gpy2 = gpy2[mask][:n]
    r1x = np.concatenate((r1_gpx1,r1_gpx2, r1_gpx1))
    r1y = np.concatenate((r1_gpy1,r1_gpy2, r1_gpy1))
    ax.scatter(r1x, r1y, s=1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
    ax.set_title('Transaxial detection position ({} first events only)'.format(n))


def plot_axial_detection(ax, coinc):
    # Axial Detection
    ad1 = tget(coinc, 'globalPosZ1')
    ad2 = tget(coinc, 'globalPosZ2')
    ad = np.concatenate((ad1, ad2))
    ax.hist(ad, histtype='step', bins=100)
    ax.set_xlabel('mm',fontsize = 13)
    ax.set_ylabel('counts',fontsize =13)
    ax.set_title('Average Z-coordinate of All Coincidences',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
def get_counts(coinc):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    ad1 = tget(coinc, 'globalPosZ1')
    ad2 = tget(coinc, 'globalPosZ2')
    z = (ad1+ad2)/2
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    trues = z[mask]
    scatters = z[~mask]
    # Randoms
    eventID1 = tget(coinc, 'eventID1')
    eventID2 = tget(coinc, 'eventID2')
    time = tget(coinc, 'time1')
    randoms = time[eventID1 != eventID2]
    Ctot = len(trues) + len(scatters) + len(randoms)
    return trues, scatters, randoms, Ctot

def get_annihilation(hits):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    proccess = tget(hits, 'processName')
    x = tget(hits, 'posX')
    y = tget(hits, 'posY')
    z = tget(hits, 'posZ')
    t = tget(hits, 'time')
    e = tget(hits, 'edep')
    mask =  ((proccess=='annihil'))
    annx = x[mask]
    anny = y[mask]
    annz = z[mask]
    annt = t[mask]
    anne = e[mask]
    return annx, anny, annz, annt, anne
  
def get_process(hits):
    process = tget(hits, 'processName')
    return process
    
def get_hits(hits):
    #return the number of hits
    hitsX = tget(hits, 'posX')
    hitsY = tget(hits, 'posY')
    hitsZ = tget(hits, 'posZ')
    hits = len(hitsX)
    return hits, hitsX, hitsY, hitsZ

def get_hit_time(hits):
    time = tget(hits, 'time')
    return time

def get_sin(coinc,Emin,Emax):
    theta = tget(coinc, 'sinogramTheta')*180/np.pi
    s = tget(coinc, 'sinogramS')
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    all1 = tget(coinc, 'energy1')
    all2 = tget(coinc, 'energy2')
    totE = all1 + all2
    mask =  ((totE/0.001 < Emax) & (totE/0.001 > Emin))
    theta = theta[mask]
    s = s[mask]
    return theta, s, totE/0.001
    
def get_singles(singles):
    xpos = tget(singles, 'globalPosX')
    ypos = tget(singles, 'globalPosY')
    zpos = tget(singles, 'globalPosZ')
    d = np.sqrt(xpos**2 + ypos**2 + zpos**2)
    return xpos,ypos,zpos,d

def get_scatt_singles(singles):
    xpos = tget(singles, 'globalPosX')
    ypos = tget(singles, 'globalPosY')
    zpos = tget(singles, 'globalPosZ')
    compt = tget(singles, 'comptonPhantom')
    rayl = tget(singles, 'RayleighPhantom')
    print(compt)
    print(rayl)
    mask = ((compt != 0) | (rayl != 0))
    xpos = xpos[mask]
    ypos = ypos[mask]
    zpos = zpos[mask]
    d = np.sqrt(xpos**2 + ypos**2 + zpos**2)
    return xpos,ypos,zpos,d

def get_filter_singles(singles, energy, erange):
    x1 = tget(singles, 'globalPosX')
    y1 = tget(singles, 'globalPosY')
    z1 = tget(singles, 'globalPosZ')
    all1 = tget(singles, 'energy')
    compt1 = tget(singles, 'comptonPhantom')
    rayl1 = tget(singles, 'RayleighPhantom')
    mask = (all1 >= (energy - erange/2)) & (all1 <= (energy + erange/2)) & ((compt1==0) & (rayl1==0))
    fx1 = x1[mask]
    fy1 = y1[mask]
    fz1 = z1[mask]
    fe1 = all1[mask]
    return fx1, fy1, fz1, fe1

def get_noise_singles_filter(singles,coinc,energy,erange):
    #create an array of True/False determining whether the single ID matches a coincidence ID
    #to mask the time of singles 
    es = tget(singles,'energy')
    maskSingle = (es >= (energy - erange/2)) & (es <= (energy + erange/2)) 
    ts = tget(singles,'time')
    ts = ts[maskSingle]
    all1 = tget(coinc, 'energy1')
    all2 = tget(coinc, 'energy2')
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    t1 = tget(coinc,'time1')
    t2 = tget(coinc,'time2')
    mask = (all1 >= (energy - erange/2)) & (all1 <= (energy + erange/2)) & ((compt1==0) & (rayl1==0)) & (all2 >= (energy - erange/2)) & (all2 <= (energy + erange/2)) &  ((compt2==0) & (rayl2==0))
    tcoinc1 = t1[mask]
    tcoinc2 = t2[mask] 
    tcoinc = np.concatenate((tcoinc1,tcoinc2))
    #find non-similar times between
    noisySingles = np.setxor1d(ts,tcoinc)
    return noisySingles,ts
        
def get_filter_singles_time(singles, energy, erange):
    x1 = tget(singles, 'globalPosX')
    y1 = tget(singles, 'globalPosY')
    z1 = tget(singles, 'globalPosZ')
    all1 = tget(singles, 'energy')
    t1 = tget(singles, 'time')
    mask = (all1 >= (energy - erange/2)) & (all1 <= (energy + erange/2))
    fx1 = x1[mask]
    fy1 = y1[mask]
    fz1 = z1[mask]
    fe1 = all1[mask]
    ft1 = t1[mask]
    return ft1

def get_unscatt_singles(singles):
    xpos = tget(singles, 'globalPosX')
    ypos = tget(singles, 'globalPosY')
    zpos = tget(singles, 'globalPosZ')
    compt = tget(singles, 'comptonPhantom')
    rayl = tget(singles, 'RayleighPhantom')
    print(compt)
    print(rayl)
    mask = ((compt == 0) & (rayl == 0))
    xpos = xpos[mask]
    ypos = ypos[mask]
    zpos = zpos[mask]
    d = np.sqrt(xpos**2 + ypos**2 + zpos**2)
    return xpos,ypos,zpos,d

def get_singles_eventID(singles):
    runID = tget(singles, 'eventID')
    return runID
    
def get_singles_time(singles):
    time = tget(singles, 'time')
    return time
    
def get_time(coinc):
    time1 = tget(coinc, 'time1')
    time2 = tget(coinc, 'time2')
    return time1, time2
    
def plot_axial_sensitivity_detection(ax, trues):
    ax.hist(trues, bins=100)
    ax.set_xlabel('mm',fontsize = 13)
    ax.set_ylabel('counts',fontsize =13)
    ax.set_title('Average Z-coordinate of True Coincidences',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)

def plot_axial_scatter_fraction(ax, coinc, scatters): 
    #% of coincidences whose average signal location results from scattered coincidence
    ad1 = tget(coinc, 'globalPosZ1')
    ad2 = tget(coinc, 'globalPosZ2')
    z = (ad1+ad2)/2
    countsa, binsa = np.histogram(scatters, bins=100)
    #print(countsa)
    countsr, binsr = np.histogram(z, bins=100)
    #print(countsr)
    ax.hist(binsa[:-1], bins=100, weights=countsa/countsr*100)
    #print(countsa/countsr)
    ax.set_xlabel('mm')
    ax.set_ylabel('%')
    ax.set_title('Coincidence Scatter fraction')
    
def get_decays(coinc):
    time = tget(coinc, 'time1')
    decayF18 = time
    return decayF18

def plot_transaxial_hits(ax, hits):
    runID = tget(hits, 'runID')
    x = tget(hits, 'posX')
    y = tget(hits, 'posY')
    ax.scatter(x, y, s=1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x position (mm)')
    ax.set_ylabel('y position (mm)')
    ax.set_title('Transaxial detection position')
    
def plot_XZ_hits(ax, hits):
    runID = tget(hits, 'runID')
    x = tget(hits, 'posX')
    z = tget(hits, 'posZ')
    ax.scatter(x, z, s=1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x position (mm)')
    ax.set_ylabel('z position (mm)')
    ax.set_title('X-Z detection position')
    
def plot_axial_hits(ax, hits):
    # Axial Detection
    z = tget(hits, 'posZ')
    y = tget(hits, 'posY')
    ax.scatter(z, y, s=1)
    ax.set_xlabel('z position (mm)')
    ax.set_ylabel('y positoin (mm)')
    ax.set_title('Axial detection position')

def plot_rad_decay(ax, end_time, decayF18):
    # histogram of decayO15
    bin_heights, bin_borders = np.histogram(np.array(decayF18), bins=100, density=False)
    bin_widths = np.diff(bin_borders)
    bin_centers = bin_borders[:-1] + bin_widths / 2
    # exponential fit
    # (ignore the warning for overflow error)
    np.seterr(all='warn', over='ignore')
    def exponenial_func(x, a, b):
        return a*np.exp(-b*x)
    popt, pcov = scipy.optimize.curve_fit(exponenial_func, bin_centers, bin_heights)
    xx = np.linspace(0, int(end_time), int(end_time))
    yy = exponenial_func(xx, *popt)
    hl = np.log(2)/popt[1]
    # plot
    ax.hist(decayF18, bins=100, label='F18 HL = 6586.2 sec', alpha=0.5, density=False)
    ax.plot(xx, yy, label=f'F18 fit HL = {hl:.2f} sec')
    ax.legend()
    ax.set_xlabel('time (s)')
    ax.set_ylabel('decay')
    ax.set_title('Rad decays')

def plot_randoms_delays(ax, randoms, delays):
    t1 = tget(delays, 'time1')
    ax.hist(randoms, bins=100, histtype='stepfilled', alpha=0.6, label=f'Real randoms = {len(randoms)}')
    ax.hist(t1, bins=100, histtype='step', label=f'Delays (estimated randoms) = {len(delays)}')
    ax.legend()
    ax.set_xlabel('time (s)')
    ax.set_ylabel('events')
    ax.set_title('Randoms')
    
def plot_sinogram(ax, theta, s):
    plt.scatter(s, theta)
    ax.set_xlabel('s (mm)',fontsize=14)
    ax.set_ylabel('theta (degrees)',fontsize = 14)
    ax.set_title('Sinogram',fontsize=16)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)

    
def plot_LOR(ax, coinc, nb, center):    
    x1 = coinc.arrays()['globalPosX1']
    y1 = coinc.arrays()['globalPosY1']
    x2 = coinc.arrays()['globalPosX2']
    y2 = coinc.arrays()['globalPosY2']
    x1 = x1[0:nb]
    y1 = y1[0:nb]
    x2 = x2[0:nb]
    y2 = y2[0:nb]
    ax.plot([x1,x2],[y1,y2],zorder = 1)    
    ax.autoscale()
    ax.set_xlabel('X-position (mm)',fontsize=14)
    ax.set_ylabel('Y-poisition (mm)',fontsize=14)
    ax.set_title('Lines of response (LOR); All Coincidences',fontsize=16)
    ax.scatter(center[0],center[1],s=100,color="black",zorder=2)
    phan = plt.Circle((0,0), 50, fill = False,zorder = 3,color="blue")
    ax.set_aspect( 1 )
    ax.add_artist(phan)
    detin = plt.Circle((0,0), 150, fill = False, zorder = 4,color="orange")
    ax.set_aspect( 1 )
    ax.add_artist(detin)
    detout = plt.Circle((0,0), 170, fill = False, zorder = 5,color="orange")
    ax.set_aspect( 1 )
    ax.add_artist(detout)
    #LSO = plt.Circle((0,0), 450, fill = False, zorder = 6,color="yellow")
    #ax.set_aspect( 1 )
    #ax.add_artist(LSO)
    #BGO = plt.Circle((0,0), 465, fill = False, zorder = 7,color="green")
    #ax.set_aspect( 1 )
    ##ax.add_artist(BGO)
    #shield = plt.Circle((0,0), 480, fill = False, zorder = 8,color="black")
    #ax.set_aspect( 1 )
    #ax.add_artist(shield)
    
def plot_true_LOR_filter(ax, coinc, nb, center,energy,erange,angle = True):    
    x1 = coinc.arrays()['globalPosX1']
    y1 = coinc.arrays()['globalPosY1']
    x2 = coinc.arrays()['globalPosX2']
    y2 = coinc.arrays()['globalPosY2']
    all1 = tget(coinc, 'energy1')
    all2 = tget(coinc, 'energy2')
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    t1 = tget(coinc,'time1')
    t2 = tget(coinc,'time2')
    mask = (all1 >= (energy - erange/2)) & (all1 <= (energy + erange/2)) & ((compt1==0) & (rayl1==0)) & (all2 >= (energy - erange/2)) & (all2 <= (energy + erange/2)) &  ((compt2==0) & (rayl2==0))
    x1 = x1[mask]
    x2 = x2[mask]
    y1 = y1[mask]
    y2 = y2[mask]
    t1 = t1[mask]
    t2 = t2[mask]
    print(len(x1))
    if angle == True:
        #create temporary lists
        x1t = []
        x2t = []
        y1t = []
        y2t = []
        t1t = []
        t2t = []
        for i in range(0,len(x1)):
            #take the x any y coordinate of the single to calculate azimuthal angle
            phi = np.arctan(y2[i]/abs(x2[i]))*180/np.pi

            #calculate the angle relative to beam (negative y-direction)
            thetai = 90 + phi

            if thetai >= 60 and thetai <= 120:
                x1t.append(x1[i])
                x2t.append(x2[i])
                y1t.append(y1[i])
                y2t.append(y2[i])
                t1t.append(t1[i])
                t2t.append(t2[i])

        x1 = x1t[0:nb]
        x2 = x2t[0:nb]
        y1 = y1t[0:nb]
        y2 = y2t[0:nb]
        t1 = t1t[0:nb]
        t2 = t2t[0:nb]
        print(len(x1))
    elif angle != True:
        print("ok")
    ax.plot([x1,x2],[y1,y2],zorder = 1)    
    plt.scatter(x2,y2,zorder=4,color="black")
    plt.scatter(x1,y1,zorder=3,color="red")
    ax.autoscale()
    ax.set_xlabel('X-position (mm)',fontsize=14)
    ax.set_ylabel('Y-poisition (mm)',fontsize=14)
    ax.set_title('Lines of response (LOR); All Coincidences',fontsize=16)
    ax.scatter(center[0],center[1],s=100,color="black",zorder=2)
    phan = plt.Circle((0,0), 50, fill = False,zorder = 3,color="blue")
    ax.set_aspect( 1 )
    ax.add_artist(phan)
    detin = plt.Circle((0,0), 150, fill = False, zorder = 4,color="orange")
    ax.set_aspect( 1 )
    ax.add_artist(detin)
    detout = plt.Circle((0,0), 170, fill = False, zorder = 5,color="orange")
    ax.set_aspect( 1 )
    ax.add_artist(detout)
    #LSO = plt.Circle((0,0), 450, fill = False, zorder = 6,color="yellow")
    #ax.set_aspect( 1 )
    #ax.add_artist(LSO)
    #BGO = plt.Circle((0,0), 465, fill = False, zorder = 7,color="green")
    #ax.set_aspect( 1 )
    ##ax.add_artist(BGO)
    #shield = plt.Circle((0,0), 480, fill = False, zorder = 8,color="black")
    #ax.set_aspect( 1 )
    #ax.add_artist(shield)
    return x1, x2, y1, y2, t1, t2
    
def plot_trans_LOR(ax, coinc, nb,center):        
    z1 = coinc.arrays()['globalPosZ1']
    y1 = coinc.arrays()['globalPosY1']
    z2 = coinc.arrays()['globalPosZ2']
    y2 = coinc.arrays()['globalPosY2']
    z1 = z1[0:nb]
    y1 = y1[0:nb]
    z2 = z2[0:nb]
    y2 = y2[0:nb]
    ax.plot([z1,z2],[y1,y2],zorder = 1)    
    ax.autoscale()
    ax.set_xlabel('Z-position in mm')
    ax.set_ylabel('Y-Poisition in mm')
    ax.set_title('Lines of response (LOR); Y-Z Plane')
    ax.scatter(center[0],center[1],s=100,color="black",zorder=2)
    ax.add_patch(Rectangle((-200,-400),400,800,edgecolor="black",fill=False,lw = 5,zorder=2))
    ax.add_patch(Rectangle((-200,-520),400,1040,edgecolor="black",fill=False,lw = 5,zorder=3))
    ax.add_patch(Rectangle((-100,-100),200,200,edgecolor="blue",fill=False,lw = 5,zorder=4))
    
def plot_true_LOR(ax, coinc, nb,center):  
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    
    x1 = coinc.arrays()['globalPosX1']
    y1 = coinc.arrays()['globalPosY1']
    x2 = coinc.arrays()['globalPosX2']
    y2 = coinc.arrays()['globalPosY2']
    x1all = x1[mask]
    y1all = y1[mask]
    x2all = x2[mask]
    y2all = y2[mask]
    x1 = x1all[0:nb]
    x2 = x2all[0:nb]
    y1 = y1all[0:nb]
    y2 = y2all[0:nb]
    ax.plot([x1,x2],[y1,y2],zorder = 1)    
    #d1 = np.sqrt(x1all**2 + y1all**2)
    #d2 = np.sqrt(x2all**2 + y2all**2)
    #print(f'Min first coincidence distance:,{min(d1)}')
    #print(f'Max first coincidence distance:,{max(d1)}')
    #print(f'Min second coincidence distance:,{min(d2)}')
    #print(f'Max second coincidence distance:,{max(d2)}')    
    ax.autoscale()
    ax.set_xlabel('X-position (mm)',fontsize=14)
    ax.set_ylabel('Y-poisition (mm)',fontsize=14)
    ax.set_title('Lines of response (LOR); True Coincidences',fontsize=16)
    ax.scatter(center[0],center[1],s=100,color="black",zorder=2)
    phan = plt.Circle((0,0), 200, fill = False,zorder = 3,color="blue")
    ax.set_aspect( 1 )
    ax.add_artist(phan)
    detin = plt.Circle((0,0), 400, fill = False, zorder = 4,color="black")
    ax.set_aspect( 1 )
    ax.add_artist(detin)
    detout = plt.Circle((0,0), 520, fill = False, zorder = 5,color="black")
    ax.set_aspect( 1 )
    ax.add_artist(detout)
    
    #cant approximate rectangular modules with a circle without some leakage outside of circle
    #LSO = plt.Circle((0,0), 450, fill = False, zorder = 6,color="yellow")
    #ax.set_aspect( 1 )
    #ax.add_artist(LSO)
    #BGO = plt.Circle((0,0), 465, fill = False, zorder = 7,color="green")
    #ax.set_aspect( 1 )
    #ax.add_artist(BGO)
    #shield = plt.Circle((0,0), 480, fill = False, zorder = 8,color="black")
    #ax.set_aspect( 1 )
    #ax.add_artist(shield)

def plot_true_trans_LOR(ax, coinc, nb,center):  
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    
    z1 = coinc.arrays()['globalPosZ1']
    y1 = coinc.arrays()['globalPosY1']
    z2 = coinc.arrays()['globalPosZ2']
    y2 = coinc.arrays()['globalPosY2']
    z1all = z1[mask]
    y1all = y1[mask]
    z2all = z2[mask]
    y2all = y2[mask]
    z1 = z1all[0:nb]
    z2 = z2all[0:nb]
    y1 = y1all[0:nb]
    y2 = y2all[0:nb]
    ax.plot([z1,z2],[y1,y2],zorder = 1)    
    ax.autoscale()
    ax.set_xlabel('Position in mm')
    ax.set_ylabel('Poisition in mm')
    ax.set_title('Lines of response (LOR)')
    ax.scatter(center[0],center[1],s=100,color="black",zorder=2)
    ax.add_patch(Rectangle((-200,-400),400,800,edgecolor="black",fill=False,lw = 5,zorder=2))
    ax.add_patch(Rectangle((-200,-520),400,1040,edgecolor="black",fill=False,lw = 5,zorder=3))
    ax.add_patch(Rectangle((-100,-100),200,200,edgecolor="blue",fill=False,lw = 5,zorder=4))
    
def plot_z_detection(ax, coinc):
    # Axial Detection
    ad1 = tget(coinc, 'globalPosZ1')
    ad2 = tget(coinc, 'globalPosZ2')
    ad = np.concatenate((ad1, ad2))
    ax.hist(ad, histtype='step', bins=100)
    ax.set_xlabel('mm',fontsize = 13)
    ax.set_ylabel('counts',fontsize =13)
    ax.set_title('Z-coordinate of All Coincidences',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
def plot_x_detection(ax, coinc):
    # Axial Detection
    ad1 = tget(coinc, 'globalPosX1')
    ad2 = tget(coinc, 'globalPosX2')
    ad = np.concatenate((ad1, ad2))
    ax.hist(ad, histtype='step', bins=100)
    ax.set_xlabel('mm',fontsize = 13)
    ax.set_ylabel('counts',fontsize =13)
    ax.set_title('X-coordinate of All Coincidences',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
def plot_y_detection(ax, coinc):
    # Axial Detection
    ad1 = tget(coinc, 'globalPosY1')
    ad2 = tget(coinc, 'globalPosY2')
    ad = np.concatenate((ad1, ad2))
    ax.hist(ad, histtype='step', bins=100)
    ax.set_xlabel('mm',fontsize = 13)
    ax.set_ylabel('counts',fontsize =13)
    ax.set_title('Y-coordinate of All Coincidences',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)

def get_coinc_energy(coinc):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    all1 = tget(coinc, 'energy1')
    all2 = tget(coinc, 'energy2')
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    trues1 = all1[mask]
    scatters1 = all1[~mask]
    trues2 = all2[mask]
    scatters2 = all2[~mask]
    return trues1, scatters1, trues2, scatters2

def get_coinc_pos(coinc,true):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    
    x1 = tget(coinc, 'globalPosX1')
    x2 = tget(coinc, 'globalPosX2')
    y1 = tget(coinc, 'globalPosY1')
    y2 = tget(coinc, 'globalPosY2')
    z1 = tget(coinc, 'globalPosZ1')
    z2 = tget(coinc, 'globalPosZ2')
    if true == False or true == "None":
        tx1 = x1
        tx2 = x2
        ty1 = y1
        ty2 = y2
        tz1 = z1
        tz2 = z2
    elif true == True:
        compt1 = tget(coinc, 'comptonPhantom1')
        compt2 = tget(coinc, 'comptonPhantom2')
        rayl1 = tget(coinc, 'RayleighPhantom1')
        rayl2 = tget(coinc, 'RayleighPhantom2')
        mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
        tx1 = x1[mask]
        tx2 = x2[mask]
        ty1 = y1[mask]
        ty2 = y2[mask]
        tz1 = z1[mask]
        tz2 = z2[mask]
    return tx1, tx2, ty1, ty2, tz1, tz2

def filter_coinc(coinc, energy, erange):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    #erange is range around energy that is used for filtering
    x1 = tget(coinc, 'globalPosX1')
    x2 = tget(coinc, 'globalPosX2')
    y1 = tget(coinc, 'globalPosY1')
    y2 = tget(coinc, 'globalPosY2')
    z1 = tget(coinc, 'globalPosZ1')
    z2 = tget(coinc, 'globalPosZ2')
    all1 = tget(coinc, 'energy1')
    all2 = tget(coinc, 'energy2')
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    comptcryst1 = tget(coinc,'comptonCrystal1')
    comptcryst2 = tget(coinc,'comptonCrystal2')
    raycryst1 = tget(coinc,'RayleighCrystal1')
    raycryst2 = tget(coinc,'RayleighCrystal2')
    t1 = tget(coinc,'time1')
    t2 = tget(coinc,'time2')
    mask = (all1 >= (energy - erange/2)) & (all1 <= (energy + erange/2)) & ((compt1==0) & (rayl1==0)) & (all2 >= (energy - erange/2)) & (all2 <= (energy + erange/2)) &  ((compt2==0) & (rayl2==0))
    fx1 = x1[mask]
    fx2 = x2[mask]
    fy1 = y1[mask]
    fy2 = y2[mask]
    fz1 = z1[mask]
    fz2 = z2[mask]
    fe1 = all1[mask]
    fe2 = all2[mask]
    t1 = t1[mask]
    t2 = t2[mask]
    return fx1, fx2, fy1, fy2, fz1, fz2, fe1, fe2, t1, t2

def get_true_singles(singles):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    energy = tget(singles, 'energy')
    compt = tget(singles, 'comptonPhantom')
    rayl = tget(singles, 'RayleighPhantom')
    mask =  ((compt==0) & (rayl==0) )
    trues = energy[mask]/0.001
    scatters = energy[~mask]/0.001
    return trues, scatters

def plot_basic_xy_reconstruction(ax, coinc, nb, cmin):
    # Axial Detection
    x1 = tget(coinc, 'globalPosX1')
    x2 = tget(coinc, 'globalPosX2')
    y1 = tget(coinc, 'globalPosY1')
    y2 = tget(coinc, 'globalPosY2')
    x = (x1+x2)/2
    y = (y1+y2)/2
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    truex = x[mask]
    truey = y[mask]
    #ax.hist2d(x, y, bins=nb)
    #ax.set_xlabel('x (mm)',fontsize = 13)
    #ax.set_ylabel('y (mm)',fontsize =13)
    #ax.set_title('Axial Reconstruction; Trues',fontsize=13)
    #plt.yticks(fontsize=12)
    #plt.xticks(fontsize=12)
    return x1, x2, y1, y2
    
def plot_basic_yz_reconstruction(ax, coinc, nb, cmin):
    # Axial Detection
    z1 = tget(coinc, 'globalPosZ1')
    z2 = tget(coinc, 'globalPosZ2')
    y1 = tget(coinc, 'globalPosY1')
    y2 = tget(coinc, 'globalPosY2')
    z = (z1+z2)/2
    y = (y1+y2)/2
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    truez = z[mask]
    truey = y[mask]
    ax.hist2d(truez, truey, bins=nb, cmin = cmin)
    ax.set_xlabel('z (mm)',fontsize = 13)
    ax.set_ylabel('y (mm)',fontsize =13)
    ax.set_title('YZ-Transaxial Reconstruction; Trues',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
def plot_basic_xz_reconstruction(ax, coinc, nb, cmin):
    # Axial Detection
    z1 = tget(coinc, 'globalPosZ1')
    z2 = tget(coinc, 'globalPosZ2')
    x1 = tget(coinc, 'globalPosX1')
    x2 = tget(coinc, 'globalPosX2')
    z = (z1+z2)/2
    x = (x1+x2)/2
    compt1 = tget(coinc, 'comptonPhantom1')
    compt2 = tget(coinc, 'comptonPhantom2')
    rayl1 = tget(coinc, 'RayleighPhantom1')
    rayl2 = tget(coinc, 'RayleighPhantom2')
    mask =  ((compt1==0) & (compt2==0) & (rayl1==0) & (rayl2==0))
    truez = z[mask]
    truex = x[mask]
    ax.hist2d(truez, truex, bins=nb, cmin = cmin)
    ax.set_xlabel('z (mm)',fontsize = 13)
    ax.set_ylabel('x (mm)',fontsize =13)
    ax.set_title('XZ-Transaxial Reconstruction; Trues',fontsize=13)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
def plot_sinogram(ax, theta, s):
    plt.scatter(s, theta)
    ax.set_xlabel('s (mm)',fontsize=14)
    ax.set_ylabel('theta (degrees)',fontsize = 14)
    ax.set_title('Sinogram',fontsize=16)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    
def get_number_coinc(coinc):
    z1 = tget(coinc, 'globalPosZ1')
    print("There are ",len(z1), " coincidences")
    
def get_TOF_data(coinc):
    # trues, scatters, randoms, Ctot = p.get_counts(coinc)
    # 
    # Cscat : scatter counts is the number of falsely located coincidence events resulting from gamma rays scattering inside the phantom
    # Ctrue : is the number of true coincidences
    # Crnd  : the number of random (accidental) coincidences
    # Ctot  : Ctot = Cscat + Ctrue + Crnd is the total number of detected coincidences, sometimes called 'prompts'
    # 
    
    x1 = tget(coinc, 'globalPosX1')
    x2 = tget(coinc, 'globalPosX2')
    y1 = tget(coinc, 'globalPosY1')
    y2 = tget(coinc, 'globalPosY2')
    z1 = tget(coinc, 'globalPosZ1')
    z2 = tget(coinc, 'globalPosZ2')
    t1 = tget(coinc, 'time1')
    t2 = tget(coinc, 'time2')
    
    return x1, x2, y1, y2, z1, z2, t1, t2