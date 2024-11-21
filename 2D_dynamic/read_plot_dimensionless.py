#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:44:03 2024
@author: lcbt87
"""
import numpy as np
import os as os
import pandas as pd
import regex as re
from scipy.special import lambertw
import matplotlib as mpl
from matplotlib import pyplot as plt
#%%
''' MODEL 1 '''
# READ input file for parameters
os.chdir('/home/lcbt87/Fault2dPML/Model1')
f = open('fpar','r')
f.readline()
li=f.readline()
dx=re.split(r'\s{1,}',li)[0]
li=f.readline()
dt=re.split(r'\s{1,}',li)[0]
ntime=f.readline() 
ntime=float(re.split(r'\s{1,}',ntime)[0])
peak=f.readline() 
ssini=f.readline() 
FricType=f.readline() 
li=f.readline()
delta, vcha, vmax, tau_f, na, na = re.split(r'\s{1,}',li)
f.readline()
f.readline() 
li=f.readline()
InitAspLen,na,na=re.split(r'\s{1,}',li)
f.readline() 
f.readline() 
li = f.readline()
nxs, nys, na,na,na= re.split(r'\s{1,}',li)
f.close()
#
InitAspLen=float(InitAspLen);nxs=float(nxs);nys=float(nys)
dx=float(dx);dt=float(dt)

InitAspLen=InitAspLen+30


# #%% read and plot wavefiled snapshot 
# ########################
# itef = open('./iter','r')
# li=itef.readline()
# iter=float(re.split(r'\s{1,}',li)[3])
# itef.close()
# if iter > 1000:
#     #iter=13000
#     ite = int(int(iter/1e3)*1e3)
#     ifi=1;fig,ax=plt.subplots(num=ifi,clear=True)
#     df=pd.read_csv('./RES/fvex'+f"{ite:05}",header=None,sep='\s{1,}')
#     vx=np.array(df).T
#     ax.imshow(vx)
#     ax.plot(-vx[124,:],'-')
#%% reference velocity
rmu= 28E9
rlam=rmu
rrho=2700.
vs=np.sqrt(rmu/rrho)
v_ray=vs*0.9194
#%% read rupture position
ouf1='./RES/ruptures'
rups1 = pd.read_csv(ouf1,header=None, delimiter=r"\s+",skipinitialspace=True).to_numpy()
rup=rups1[:]
rupx=[];
for rr in rup[:]:    # append pos, time. Subtract asperity center
    if rr[0]>nxs:rupx.append([rr[0]-nxs+0.5,rr[2]]) 
rupx=np.array(rupx)
#%% dimensionless rupture location
rup_tim=rupx[:,1]*dt*v_ray/(dx*InitAspLen/2)
rup_pos=rupx[:,0]*dx/(dx*InitAspLen/2)
#%% plot reference Rayleigh wave timeline:
##########################################
v_ray=1.0
#
ifi=3;fig2,ax2=plt.subplots(1,1,num=ifi,clear=True)
c4=mpl.colors.rgb2hex((.8, .8, .8, 0.8), keep_alpha=True)
ax2.axvline(1,color=c4)
ax2.axvline(2.7,color=c4)
lmax=(4)
tmax1=(6/5)*lmax/v_ray
k=16;ax2.plot([0,lmax],[k/2,tmax1+k/2],':',color='gray', label=r'$C_{lim}$')
#%% plot numerical solution:
#############################
nlimit=14
nlimit1=225
c1=mpl.colors.rgb2hex((0.0, 0.0, 1.0, 1.0), keep_alpha=True)
c2=mpl.colors.rgb2hex((0.3, 0.3, 1.0, 0.3), keep_alpha=True)
c3=mpl.colors.rgb2hex((0.3, 0.3, 1.0, 0.3), keep_alpha=True)

ax2.plot(rup_pos[nlimit+1:nlimit1],rup_tim[nlimit+1:nlimit1],'.', color=c1, 
         label='$L_c$=1218 m',linewidth=3)
#ax2.plot(rup_pos[:nlimit],rup_tim[:nlimit],'.', color=c2)
ax2.plot(rup_pos[nlimit1:],rup_tim[nlimit1:],'.',color=c1, linewidth=3)
ax2.set_xlim([0,lmax]);
ax2.set_ylim([8,4*tmax1]);#ax.set_aspect(.1)
ax2.set_aspect(.35)

#%%
# import time
# import sys
# time.sleep(0.2)
# sys.exit()

''' MODEL 2 '''
# READ input file for parameters
os.chdir('/home/lcbt87/Fault2dPML/Model2')
f = open('fpar','r')
f.readline()
li=f.readline()
dx=re.split(r'\s{1,}',li)[0]
li=f.readline()
dt=re.split(r'\s{1,}',li)[0]
ntime=f.readline() 
ntime=float(re.split(r'\s{1,}',ntime)[0])
peak=f.readline() 
ssini=f.readline() 
FricType=f.readline() 
li=f.readline()
delta, vcha, vmax, tau_f, na, na = re.split(r'\s{1,}',li)
f.readline()
f.readline() 
li=f.readline()
InitAspLen,na,na=re.split(r'\s{1,}',li)
f.readline() 
f.readline() 
li = f.readline()
nxs, nys, na,na,na= re.split(r'\s{1,}',li)
f.close()
#
InitAspLen=float(InitAspLen);nxs=float(nxs);nys=float(nys)
dx=float(dx);dt=float(dt)

InitAspLen=InitAspLen+30


#%% read and plot wavefiled snapshot 
########################
# itef = open('./iter','r')
# li=itef.readline()
# iter=float(re.split(r'\s{1,}',li)[3])
# itef.close()
# if iter > 1000:
#     #iter=13000
#     ite = int(int(iter/1e3)*1e3)
#     ifi=1;fig,ax=plt.subplots(num=ifi,clear=True)
#     df=pd.read_csv('./RES/fvex'+f"{ite:05}",header=None,sep='\s{1,}')
#     vx=np.array(df).T
#     ax.imshow(vx)
#     ax.plot(-vx[124,:],'-')
#%% reference velocity
rmu= 28E9
rlam=rmu
rrho=2700.
vs=np.sqrt(rmu/rrho)
v_ray=vs*0.9194
#%% read rupture position
ouf1='./RES/ruptures'
rups1 = pd.read_csv(ouf1,header=None, delimiter=r"\s+",skipinitialspace=True).to_numpy()
rup=rups1[:]
rupx=[];
for rr in rup[:]:    # append pos, time. Subtract asperity center
    if rr[0]>nxs:rupx.append([rr[0]-nxs+0.5,rr[2]]) 
rupx=np.array(rupx)
#%% dimensionless rupture location
rup_tim=rupx[:,1]*dt*v_ray/(dx*InitAspLen/2) - 1.7 
rup_pos=rupx[:,0]*dx/(dx*InitAspLen/2)
v_ray=1.0

#%% plot numerical solution:
#############################
nlimit=14
nlimit1=225
c1=mpl.colors.rgb2hex((0.0, 1.0, 0.0, 1.0), keep_alpha=True)
c2=mpl.colors.rgb2hex((0.3, 1.0, 0.3, .3), keep_alpha=True)
c3=mpl.colors.rgb2hex((0.3, 1.0, 0.3, .3), keep_alpha=True)

ax2.plot(rup_pos[nlimit+1:nlimit1],rup_tim[nlimit+1:nlimit1],'-', color=c1, 
         label='$L_c$=609 m',linewidth=1)
ax2.plot(rup_pos[nlimit1:],rup_tim[nlimit1:],'-', color=c1,linewidth=1)

#ax2.plot(rup_pos[:nlimit],rup_tim[:nlimit],'.', color=c2)

#ax2.set_xlim([0,lmax]);
#ax2.set_ylim([0,4*tmax1]);#ax2.set_aspect(.1)
#%%  
#%%
########################
#### Analytic solution
# ########################
ta = np.linspace(0,20,100)
#tab = np.linspace(800*dt,12000*dt,100)
''' 1) alpha factor solution in exponential and lambertw'''
def ya1(t, t0, yc, C, alpha):
    yy = yc * (1 + (1/alpha) * lambertw( np.exp( alpha* (C/yc) * (t-t0)   )))
    return(np.real(yy))
tsh=14.4
fa1 = ya1(ta,tsh,1,v_ray,1.0)
#fa1b = ya1(tab,tsh,InitAspLen*dx/2,v_ray,1.0)

ax2.plot(fa1, ta, '--', label='analytic', color='r', linewidth=2);
#ax2.plot(fa1b, tab, '--', color='r');
ax2.set_xlabel('position (m)');ax2.set_ylabel('time (s)')
#ax2.axhline(13400*dt)
#ax2.axvline(2000)
lgr=mpl.colors.rgb2hex((0.9, .8, .9, 1.0), keep_alpha=True)
ax2.annotate('supershear transition', xy=(2.75,11.9), rotation=90)
ax2.fill_between([2.7, 5.0], 8, 20, alpha=0.5, color=lgr)
#ax2.set_xlim(0,3500)
ax2.set_xlabel(r'L/$L_c$',fontsize=14)
ax2.set_ylabel(r'$t \ C_{lim} / L_c$',fontsize=14)
ax2.legend(loc='lower right')
ax2.set_title('2D homogeneous')
fig2.savefig('../2d_6.25_0.0003125_1218_609b.pdf')

