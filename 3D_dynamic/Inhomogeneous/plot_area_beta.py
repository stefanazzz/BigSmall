# read and plot results from OctantSpace FD simulation
import pandas as pd
import os as os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
import sys
from scipy.special import lambertw


# clean trailing and extra spaces:
def clean_spaces(inf,ouf):
    awk_command = f"awk '$1=$1' {inf}"
    try:
        completed_process = subprocess.run(awk_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if completed_process.returncode == 0:
            with open(ouf, "w") as file:
                file.write(completed_process.stdout)
            print(f"Output saved to {ouf}")
        else:
            print(f"Command failed with error code {completed_process.returncode}: {completed_process.stderr}")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error code {e.returncode}: {e.output}")
#%%#############################
#### DEFINE WORKING FIRST FOLDER
##############################
os.chdir('/home/lcbt87/BigSmall/Inhomogeneous2/RES_0.2_0.05_2.0_5.0')
#os.chdir('/home/lcbt87/Big/RES')

# read parameters of simulation from parameterisation file fpar:
f=open(file='fpar',mode='r')
fpar = f.readlines();
dx = fpar[1].split(" ");dx = float(dx[0])
dt = fpar[2].split(" ");dt = float(dt[0])
rad = fpar[12].split(" ");rad = float(rad[0])
#%%#######################################
#### READ AND PLOT SNAPSHOT OF SLIP RATE
########################################
bobo = open('../iter','r')
iter= bobo.read()
print(iter)
inf = './fldz'+str(eval(iter)).zfill(5) # current or final iteration
inf = './fldz00600' # fixed iteration
ouf = inf+'s'
clean_spaces(inf,ouf)
fldz = pd.read_csv(ouf,header=None,sep=' ',skipinitialspace=True)
fz=fldz.to_numpy()
# fill symmetry
fz_lr = np.flip(fz, axis=1)
joined_arrays = np.concatenate((fz_lr, fz), axis=1)
fz_ud = np.flip(joined_arrays, axis=0)
joined_arrays = np.concatenate((fz_ud, joined_arrays), axis=0)
#################################
# plot ruptre area
#################################
fig,ax=plt.subplots(num=1,clear=True)
#ax.cla()
ax.imshow(joined_arrays)
ax.set_title(inf)
ax.set_aspect(1)

#%%########################################
#### plot limit velocity timeline for reference
#########################################
fig,ax=plt.subplots(num=2,clear=True)
''' v_lim = sqrt(C_s * C_Ray) is equivalent vel for elliptical crack'''
v_lim=np.sqrt(1.0*0.9194)
#v_lim=1.00
slow=1/v_lim
l100=int(10)
ax.plot([0,10],[0,l100],':', label=r'$C_{lim}$')
nx=200
xma=float(nx)*dx/2.0
xma=2.*xma/rad
tma=1.2*xma/v_lim
ax.set_xlim([0,xma])
ax.set_ylim([0,tma])
ax.set_xlabel(r'L/$L_c$',fontsize=14)
ax.set_ylabel(r'$t \ C_{lim} / L_c$',fontsize=14)
ax.set_aspect(.7)
c4=mpl.colors.rgb2hex((.8, .8, .8, 0.8), keep_alpha=True)
ax.axvline(1,color=c4)


#%%#####################################################
#### Read, process and plot ruptures from FIRST folder
######################################################

inf2='./ruptures'
ouf2 = inf2+'s'
clean_spaces(inf2,ouf2)
rups = pd.read_csv(ouf2,header=None,sep=' ',skipinitialspace=True)
rup = rups.to_numpy()
'''
- Append only last output from each timestep, by appending previous line
when time step changes. Aim is get complete num of rupt points for each time. 
- Only 1/4 of rupture area is computed, therefore:
    Rupture area is A = (4 * dx^2 * nbr) where nbr is number of 'broken' cells
    Equivalent radius is sqrt(A/pi) = 2* dx * sqrt(nbr/pi)
    Equivalent diameter is 4 * dx *sqrt(nbr/pi)
'''
dira=2; # 2 if L is radius, 4 if L is diameter

rupxy=[]; fi = 0; timei = 0
for i, rr in enumerate(rup[:]):    
    if rr[2] > timei: 
        fi = fi +1
        timei = rr[2]
        rupxy.append([rup[i-1,2]*dt, dira*dx*np.sqrt(float(rup[i-1,3])/np.pi) ])
        if (fi == 2):
            rad = dira*dx*np.sqrt(float(rup[i,3])/np.pi) # equiv. initial radius
            print('first',rad)
''' normalise both t and x and plot'''
rupxy=np.array(rupxy)/rad; # norm t and x for self-similar x/t''':
rupxy[:,0]=rupxy[:,0]*v_lim # scale time by limit vel
ax.plot(rupxy[:750,1],rupxy[:750,0],'.', clip_on = True,label='$L_c$=10.0 m',color='blue')
#ax.legend()
#ax.plot(rupy[:,1],rupy[:,0],'-')

#sys.exit()
# #%%#####################################################
# #### Read, process and plot ruptures from SECOND folder
# #######################################################
os.chdir('/home/lcbt87/BigSmall/Inhomogeneous2/RES_0.4_0.1_4.0_10.0')
#os.chdir('/home/lcbt87/BigSmall/Inhomogeneous2/RES_0.2_0.05_2.0_5.0')
f=open(file='fpar',mode='r')
fpar = f.readlines();
dx = fpar[1].split(" ");dx = float(dx[0])
dt = fpar[2].split(" ");dt = float(dt[0])
#rad = fpar[12].split(" ");rad = float(rad[0])
inf2='./ruptures'
ouf2 = inf2+'s'
clean_spaces(inf2,ouf2)
rups = pd.read_csv(ouf2,header=None,sep=' ',skipinitialspace=True)
rup = rups.to_numpy()
rupxy=[]; fi = 0; timei = 0
for i, rr in enumerate(rup[:]):    
    if rr[2] > timei: 
        fi = fi +1
        timei = rr[2]
        rupxy.append([rup[i-1,2]*dt, dira*dx*np.sqrt(float(rup[i-1,3])/np.pi) ]) 
        if (fi == 2):
            rad = dira*dx*np.sqrt(float(rup[i,3])/np.pi) # approx initial length L is rad.
            print('second',rad)
''' normalise and plot'''
rupxy=np.array(rupxy)/rad; # norm t and x for self-similar x/t''':
rupxy[:,0]=rupxy[:,0]*v_lim # scale time by limit vel
c1=mpl.colors.rgb2hex((0.0, 1.0, 0.0, 1.0), keep_alpha=True)
ax.plot(rupxy[:750,1],rupxy[:750,0],marker='.', clip_on = True,label='$L_c$=5.0 m',
        color=c1, markersize=.01)
ax.legend()
#ax.grid()
#ax.plot(rupy[:,1],rupy[:,0],'-')

# ########################
# #### Analytic solution
# # ########################

''' 0) original solution adjusted'''
ta = np.linspace(-5,20,100)
def ya0(t, t0, yc, C): # analytic sol.
    yy = yc * (1 + lambertw( np.exp( (C/yc) * (t-t0)   )))
    return(np.real(yy))
fa0 = ya0(ta,2.2,.5,v_lim)
fan0 = fa0 + .5
tan = ta*v_lim 
#ax.plot(fan0, tan, '--', label='fit 0');

''' 1) alpha factor solution in exponential and lambertw'''
def ya1(t, t0, yc, C, alpha):
    yy = yc * (1 + (1/alpha) * lambertw( np.exp(alpha* (C/yc) * (t-t0)   )))
    return(np.real(yy))
#fa1 = ya1(ta,10.,1.0,v_lim,1) #fpar_ok
fa1 = ya1(ta,2.45,1.0,v_lim,2)
#fa1 = ya1(ta,2.,.7,v_lim,1)
fan1= fa1
tan = ta*v_lim 
ax.plot(fan1, tan, '--', label='analytic', color='r');

''' 2) overall inverse alpha '''
def ya2(t, t0, yc, C, alpha):
    yy = (1/alpha) *yc * (1 +  lambertw( np.exp( (C/yc) * (t-t0)   )))
    return(np.real(yy))
fa2 = ya2(ta,4.4,1,v_lim,2)
fan2 = fa2 + 0.5
tan = .5*ta*v_lim # why the factor of 2 in time????
#ax.plot(fan2, tan, '--', label='fit2');

''' 3) half in exponential '''
def ya3(t, t0, yc, C):
    yy = yc * (1 + 1 * lambertw( np.exp( 0.5* (C/yc) * (t-t0)   )))
    return(np.real(yy))
fa3 = ya3(ta,4.4,1/2,v_lim)
fan3 = fa3 + 0.5
tan = 0.5*ta*v_lim # why the factor of 2 in time????
#ax.plot(fan3, tan, '--', label='fit3');

# Normalisation only for time*v_lim as Lc in ya is 1.
# Still not sure why the 0.5 factor and +1 in fan is necessary here. 
# But the scaling across the two simulation works perfectly. 
#tan = ta*v_lim
ax.legend(loc='lower right')
ax.set_title('3D inhomogeneous')
fig.savefig('../3d_inhomogeneous_two_sizes.svg')
