#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
June 24 2024
@author: stefan

Plots the outputs from Check_Areas.py

2) For different sub-fault sizes L, calculate F = mean(stress) x sqrt(L)
3) Find max of F(L) 
    Rationale is Lc > C Gamma mu' / (mean(tau))*2
    take new Gamma = C mu' Gamma' and square root, so that
    mean(tau) sqrt(L) > Gamma 
    The largest possible value of F will fail first. 
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
from operator import itemgetter
from copy import deepcopy
#%%
def read_proc(infile):
    data = np.loadtxt(infile, delimiter=',')
    ''' x = data[:, 0] # dimension x of subfault or array x**2
        y = data[:, 1] # stress mean
        z = data[:, 2] # standard deviation within each subfault 
        Sort the values based on the size of the subfault (32, 64, ...4096 squared cells)'''
    sorted_values = sorted(data, key=itemgetter(0))
    ''' Group values based on the size of the subfault (32, 64, ...4096 squared cells)'''
    grouped_values = {key: list(group) for key, group in groupby(sorted_values, key=itemgetter(0))}
    #grouped_values = {key: list(group) for key, group in groupby(data, key=itemgetter(0))}
    
    ''' Organise mean (y_values) by sub-array size. 
        NOTE: also calculates mean and std but they are NOT being used for now'''
    result = {}
    for key, group in grouped_values.items():
        print('        sub array size:',key)
        y_values = [item[1] for item in group]
        std_values = [item[2] for item in group] # not used
        var_values = [item[3] for item in group] #not used
        mean_y = np.mean(y_values) # not used - this is mean of the means so basically a const. 
        std_mean = np.std(y_values) # not used
        std_all = np.mean(std_values) #not used
        result[key] = {'mean': mean_y, 'stdmean': std_mean, 'stdall': std_all, 'y_values': y_values}
    
    x_values = list(result.keys()) # this are the sizes of sub-arrays
    #mean_values = [values['mean'] for values in result.values()]
    #std_values = [values['stdmean'] for values in result.values()]
    #stdi_values = [np.sqrt(values['stdall']) for values in result.values()]
    yvals_values = [values['y_values'] for values in result.values()]
    
    ''' Calculate F = mean(tau) * sqrt(L) for each array size and find maximum Fmax.'''
    Flist=[]; FlistL=[]; Fras=[];
    for j, xv in enumerate(x_values[:]):
       F = np.sqrt(xv) * np.array(yvals_values[j])
       indy=np.argmax(F)
       Flist.append(F[indy])
       FL = np.array(yvals_values[j])
       FlistL.append(max(FL))
       Fr = xv * np.array(yvals_values[j])
       Fras.append(max(Fr))
    return(x_values,Flist,FlistL,Fras)
#%%
res16=[];res4=[];res1=[];res2=[];res8=[]
maxall=[];
for i in range(6):
    print('_____________________________________________')
    print('reading/processing results from iteration',i,':')
    #
    print('    upper limit:',16384)
    x_values,Flist,FlistL,Fras = read_proc('Areas_stats_16384_iter_'+str(i)+'.csv')
    res16.append((x_values,Flist,FlistL,Fras))
    maxall.append([16384,x_values[np.argmax(Flist)]])
    print('    upper limit:',4096)
    x_values,Flist,FlistL,Fras = read_proc('Areas_stats_4096_iter_'+str(i)+'.csv')
    res4.append((x_values,Flist,FlistL,Fras))
    maxall.append([4096,x_values[np.argmax(Flist)]])
    print('    upper limit:',1024)
    x_values,Flist,FlistL,Fras = read_proc('Areas_stats_1024_iter_'+str(i)+'.csv')
    res1.append((x_values,Flist,FlistL,Fras))
    maxall.append([1024,x_values[np.argmax(Flist)]])
    print('    upper limit:',2048)
    x_values,Flist,FlistL,Fras = read_proc('Areas_stats_2048_iter_'+str(i)+'.csv')
    res2.append((x_values,Flist,FlistL,Fras))
    maxall.append([2048,x_values[np.argmax(Flist)]])
    print('    upper limit:',8192)
    x_values,Flist,FlistL,Fras = read_proc('Areas_stats_8192_iter_'+str(i)+'.csv')
    res8.append((x_values,Flist,FlistL,Fras))
    maxall.append([8192,x_values[np.argmax(Flist)]])
res16=np.array(res16);res4=np.array(res4);res1=np.array(res1);res2=np.array(res2);res8=np.array(res8)
len16=res16[0,0];len4=res4[0,0];len1=res1[0,0]; len2=res2[0,0]; len8=res8[0,0]

#%%
''' plot F = mean(tau) * sqrt(L)  vs. L'''
print('_____________________________________________')
print('_____________________________________________')
print('_____________________________________________')
print('_____________________________________________')
print('                PLOTTING....')
mean_max16= res16[:,1].mean(axis=0); std_max16=res16[:,1].std(axis=0)
mean_max4 =  res4[:,1].mean(axis=0);  std_max4= res4[:,1].std(axis=0)
mean_max1 =  res1[:,1].mean(axis=0);  std_max1= res1[:,1].std(axis=0)
mean_max2 =  res2[:,1].mean(axis=0);  std_max2= res2[:,1].std(axis=0)
mean_max8 =  res8[:,1].mean(axis=0);  std_max8= res8[:,1].std(axis=0)
#%%
fig,ax=plt.subplots(num=1,clear=True)
xfit1 = np.array([1e-3,1e5])
yfit1 = 1e4*xfit1**(1/3)
#yfit2 = 2e5*np.log(x_values[:-1])
for ili in range(-8,8):
    ax.plot(xfit1 , np.sqrt(2)**ili * yfit1, '-', color='lightgrey')
#
ax.errorbar(len16, mean_max16,yerr=std_max16/2, fmt='o')    
ax.errorbar(len4,  mean_max4, yerr=std_max4/2,  fmt='o')    
ax.errorbar(len1,  mean_max1, yerr=std_max1/2,  fmt='o')      
ax.errorbar(len2,  mean_max2, yerr=std_max2/2,  fmt='o')      
ax.errorbar(len8,  mean_max8, yerr=std_max8/2,  fmt='o')      
ax.set_xlabel(r'$L_{}$'); 
#ax.set_ylabel(r'max($\overline{\tau}$ $\sqrt{L}$)')
ax.set_ylabel(r'max($\zeta$)')
ax.set_yscale('log')
ax.set_xscale('log')
yallmax = max([max(mean_max16),max(mean_max4),max(mean_max1),max(mean_max2),max(mean_max8)])+1e4
yallmin = min([min(mean_max16),min(mean_max4),min(mean_max1),min(mean_max2),min(mean_max8)])-1e1
xallmax = max([max(len16),max(len4),max(len1),max(len2),max(len8)])+1e4
xallmin = min([min(len16),min(len4),min(len1),min(len2),min(len8)])-1e1
ax.set_ylim(top=yallmax, bottom=yallmin); ax.set_xlim(xallmin,xallmax)
fig.tight_layout()
fig.canvas.draw()
fig.canvas.flush_events()
fig.savefig('tau_sqrt_L.pdf')

#%%
''' plot F = mean(tau) for each array size 
    This is the case where Gamma' = Gamma x L, i.e., scalable fracture energy '''
mean_max16= res16[:,2].mean(axis=0); std_max16=res16[:,2].std(axis=0)
mean_max4 =  res4[:,2].mean(axis=0);  std_max4= res4[:,2].std(axis=0)
mean_max1 =  res1[:,2].mean(axis=0);  std_max1= res1[:,2].std(axis=0)
mean_max2 =  res2[:,2].mean(axis=0);  std_max2= res2[:,2].std(axis=0)
mean_max8 =  res8[:,2].mean(axis=0);  std_max8= res8[:,2].std(axis=0)
fig,ax=plt.subplots(num=2,clear=True)
xfit1 = np.array([1e-3,1e5])
yfit1 = 1e4*xfit1**(-1/7)
for ili in range(-8,8):
    ax.plot(xfit1 , np.sqrt(2)**ili * yfit1, '-', color='lightgrey')
ax.errorbar(len16, mean_max16,yerr=std_max16/2, fmt='o')    
ax.errorbar(len4,  mean_max4, yerr=std_max4/2,  fmt='o')    
ax.errorbar(len1,  mean_max1, yerr=std_max1/2,  fmt='o')      
ax.errorbar(len2,  mean_max2, yerr=std_max2/2,  fmt='o')      
ax.errorbar(len8,  mean_max8, yerr=std_max8/2,  fmt='o')      
ax.set_xlabel(r'$L_{}$'); ax.set_ylabel(r'max($\overline{\tau}$)')
ax.set_yscale('log')
ax.set_xscale('log')
yallmax = max([max(mean_max16),max(mean_max4),max(mean_max1),max(mean_max1),max(len8)])+1e3
yallmin = min([min(mean_max16),min(mean_max4),min(mean_max1),min(mean_max2),max(mean_max8)])-1e2
xallmax = max([max(len16),max(len4),max(len1),max(len2),max(len8)])+1e4
xallmin = min([min(len16),min(len4),min(len1),min(len2),min(len8)])-1e1
ax.set_ylim(top=yallmax, bottom=yallmin); ax.set_xlim(xallmin,xallmax)
fig.tight_layout()
fig.canvas.draw()
fig.canvas.flush_events()
#fig.savefig('tau.pdf')


#%%
''' plot F = mean(tau) * L for each array size 
    This is the case for R&S friction Lb = mu' Dc / (b sigma) assuming tau propto sigma'''
mean_max16= res16[:,3].mean(axis=0); std_max16=res16[:,3].std(axis=0)
mean_max4 =  res4[:,3].mean(axis=0);  std_max4= res4[:,3].std(axis=0)
mean_max1 =  res1[:,3].mean(axis=0);  std_max1= res1[:,3].std(axis=0)
mean_max2 =  res2[:,3].mean(axis=0);  std_max2= res2[:,3].std(axis=0)
mean_max8 =  res8[:,3].mean(axis=0);  std_max8= res8[:,3].std(axis=0)
fig,ax=plt.subplots(num=3,clear=True)
xfit1 = np.array([1e-3,1e5])
yfit1 = 1e5*xfit1**(.85)
for ili in range(-10,5):
    ax.plot(xfit1 , 2**ili * yfit1, '-', color='lightgrey')
ax.errorbar(len16, mean_max16,yerr=std_max16/2, fmt='o')    
ax.errorbar(len4,  mean_max4, yerr=std_max4/2,  fmt='o')    
ax.errorbar(len1,  mean_max1, yerr=std_max1/2,  fmt='o')  
ax.errorbar(len2,  mean_max2, yerr=std_max2/2,  fmt='o')      
ax.errorbar(len8,  mean_max8, yerr=std_max8/2,  fmt='o')          
ax.set_xlabel(r'$L_{}$'); ax.set_ylabel(r'max($\overline{\tau}$ ${L}$)')
ax.set_yscale('log')
ax.set_xscale('log')
yallmax = max([max(mean_max16),max(mean_max4),max(mean_max1)])+5e6
yallmin = min([min(mean_max16),min(mean_max4),min(mean_max1)])-1e4
xallmax = max([max(len16),max(len4),max(len1)])+1e4
xallmin = min([min(len16),min(len4),min(len1)])-1e1
ax.set_ylim(top=yallmax, bottom=yallmin); ax.set_xlim(xallmin,xallmax)
fig.tight_layout()
fig.canvas.draw()
fig.canvas.flush_events()
#fig.savefig('tau_L.pdf')

#%%%################################################
#maxall=np.array(maxall)
sortd = sorted(maxall, key=itemgetter(0))
groupd = {key: list(group) for key, group in groupby(sortd, key=itemgetter(0))}
#%%
fig,ax=plt.subplots(num=4,clear=True)

almax=[]
for k in groupd.keys():    
    valz = np.array(groupd[k])
    meaz = np.mean(valz[:,1])
    stdz = np.std(valz[:,1])
    print(k,meaz,stdz)
    #ax.plot(valz[:,0],valz[:,1],'o')
    ax.errorbar(k,meaz,yerr=stdz,fmt='o')
    for jj in valz[:,1]:
        almax.append([k,jj,jj/k])
almax=np.array(almax)
x = almax[:,0]
x=x[:,np.newaxis]
a, _, _, _ = np.linalg.lstsq(x, almax[:,1], rcond=None)
a=a[0]
n = len(almax)
#  standard error on slope
mex = np.mean(x[:,0])
SE1 = np.sqrt( np.sum( (almax[:,1]-x[:,0])**2 ) / (n-1)) 
SE2 = np.sqrt( np.sum( (mex-x[:,0])**2 ) ) 
SE = SE1 / SE2
print('slope=',a,'+/-',SE/2)
ax.plot([0,16348],[0,16348*a],'-',color='lightgrey')
ax.plot([0,16348],[0,16348*(a-SE/2)],'-',color='lightgrey')
ax.plot([0,16348],[0,16348*(a+SE/2)],'-',color='lightgrey')

#ax.plot(almax[:,0],almax[:,1],'o')
# ax.errorbar(maxall[:,0], maxall[:,1], fmt='o')    
ax.set_xlabel(r'$L_{f}$'); ax.set_ylabel(r'$L_c$')
ax.set_yscale('log')
ax.set_xscale('log')
fig.tight_layout()
fig.canvas.draw()
fig.canvas.flush_events()
fig.savefig('Lf_Lc_linearfit.pdf')
#%%####################
fig,ax=plt.subplots(num=5,clear=True)
ax.hist(almax[:,2], bins=6, density=False,
         histtype ='bar',
         ec='black')
fig.tight_layout()
fig.canvas.draw()
fig.canvas.flush_events()
#%%
sortalmax=np.sort(almax[:,2])
p = 1. * np.arange(len(sortalmax)) / (len(sortalmax) - 1)
fig,ax=plt.subplots(num=6,clear=True)
ax.plot(sortalmax,p,'o')
fig.tight_layout()
fig.canvas.draw()
fig.canvas.flush_events()
#%%
import collections
counter = collections.Counter(sortalmax)