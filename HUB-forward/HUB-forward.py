import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import matplotlib.ticker 
from scipy.optimize import minimize
from random import seed
from random import random
import random
from scipy.stats import gumbel_l
import math
from scipy.stats import poisson
import os
seed(123434)  
plt.rc('font', family = 'Sans', serif = 'Computer Modern Roman', size = 50)
plt.rc('xtick', labelsize=40) 
plt.rc('ytick', labelsize=40)
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
modes=[]
################################################################################
print("In the HUB-forward code we propose an underlying distribution (probability distribution function), and get the fraction of frozen droplets (fice) and the cumulative freezing spectrum (Nm).")
print("Enter below the input parameters: ")
check=True
while(check):
    ndroplets=input("Number of droplets (10-10000): "); ndroplets=int(ndroplets)
    if ndroplets >=10 and ndroplets<=10000:
        check=False
while(check):
    distchoice=input("Do you want to change the type of distribution? Default is Gaussian. Yes or No? ")
    if distchoice == 'No':
        check=False; disttype=1
    if distchoice == 'Yes':
        check=False; 
        disttype=input("Type 2 for Log-normal or 3 for Gumbel with left tail: "); disttype=int(disttype) 
check=True
while(check):
    nconc=input("Number of concentrations (5-9): "); nconc=int(nconc)
    if nconc>=5 and nconc<=9:
        check=False
check=True
while(check):
    nsubpop=input("Number of subpopulations (1, 2 or 3): "); nsubpop=int(nsubpop)
    if nsubpop>=1 and nsubpop<=3:
        check=False
if nsubpop==1:
   Tmode1=input("Most likely freezing temperature #1 (in Celsius): "); Tmode1=float(Tmode1)
   s1=input("Spread #1: "); s1=float(s1)
   param=np.array([Tmode1,s1]) 
   modes.append(Tmode1)
if nsubpop==2:
   Tmode1=input("Most likely freezing temperature #1 (in Celsius): "); Tmode1=float(Tmode1)
   s1=input("Spread #1: "); s1=float(s1)
   Tmode2=input("Most likely freezing temperature #2 (in Celsius): "); Tmode2=float(Tmode2)
   s2=input("Spread #2: "); s2=float(s2)
   c2=input("Weight for subpopulation #2: "); c2=float(c2)
   param=np.array([Tmode1,s1,Tmode2,s2,c2]) 
   modes.append(Tmode1)
   modes.append(Tmode2)
if nsubpop==3:
   Tmode1=input("Most likely freezing temperature #1 (in Celsius): "); Tmode1=float(Tmode1)
   s1=input("Spread #1: "); s1=float(s1)
   Tmode2=input("Most likely freezing temperature #2 (in Celsius): "); Tmode2=float(Tmode2)
   s2=input("Spread #2: "); s2=float(s2)
   c2=input("Weight for subpopulation #2: "); c2=float(c2)
   Tmode3=input("Most likely freezing temperature #3 (in Celsius): "); Tmode3=float(Tmode3)
   s3=input("Spread #3: "); s3=float(s3)
   c3=input("Weight for subpopulation #3: "); c3=float(c3)
   param=np.array([Tmode1,s1,Tmode2,s2,c2,Tmode3,s3,c3])
   modes.append(Tmode1)
   modes.append(Tmode2)
   modes.append(Tmode3)
print("\n################### Summary ")
print("Sample "+str(nconc)+" concentrations with "+str(ndroplets)+" droplets each.")
if nsubpop==1:
    print("Underlying distribution composed of 1 subpopulation of IN with Tmode1= "+str(Tmode1)+"^oC and spread1="+str(s1))
if nsubpop==2:
    print("Underlying distribution composed of 2 subpopulations of INs with Tmode1= "+str(Tmode1)+"^oC, spread1="+str(s1)+", weight1="+str(1-c2)+", Tmode2="+str(Tmode2)+"^oC, spread2="+str(s2)+", weight2="+str(c2))
if nsubpop==3:
    print("Underlying distribution composed of 3 subpopulations of INs with Tmode1= "+str(Tmode1)+"^oC, spread1="+str(s1)+", weight1="+str(1-c2-c3)+", Tmode2="+str(Tmode2)+"^oC, spread2="+str(s2)+", weight2="+str(c2)+", Tmode3="+str(Tmode3)+"^oC, spread3="+str(s3)+", weight3="+str(c3))
################################################################################
#The underlying probability distribution function can be modified here
def normalized_PDF(x,mode,scale,disttype):
    if disttype==1:
       return (1/(scale*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mode)/scale)**2)                  # Gaussian distribution
    if disttype==2:
       return (1/((x-mode)*scale*np.sqrt(2*np.pi)))*(np.exp(-0.5*(np.log(x-mode)/scale)**2)) # Log-normal distribution
    if disttype==3:
       return (1/scale)*np.exp((x-mode)/scale - np.exp((x-mode)/scale))                      # Gumbel with left-tail
################################################################################
x = np.linspace(-30, 0, num=10000)
if nsubpop==1:
   #pdf=normalized_PDF(x,param[0],param[1])
   pdf=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf); pdf[ii] = 0
if nsubpop==2:
   #pdf=(1-param[4])*normalized_PDF(x,param[0],param[1])+param[4]*normalized_PDF(x,param[2],param[3])
   pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
   pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
   pdf=(1-param[4])*pdf1+param[4]*pdf2
if nsubpop==3:
   #pdf=(1-param[4]-param[5])*normalized_PDF(x,param[0],param[1])+param[4]*normalized_PDF(x,param[2],param[3])+param[5]*normalized_PDF(x,param[6],param[7])
   pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
   pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
   pdf3=normalized_PDF(temp_fit,param[5],param[6],disttype); ii = np.isnan(pdf3); pdf3[ii] = 0
   pdf=(1-param[4]-param[7])*pdf1+param[4]*pdf2+param[7]*pdf3
norm=np.sum(pdf)
pdf=pdf/norm
######################################################################################################################
def N_bacteria_varied(N_bacteria_set_mean,N):
    Nbacteria_random = poisson.rvs(mu=N_bacteria_set_mean, size=N)  
    sets = []
    for i in range(N): 
        # for each droplet, sample the number of INs in each droplet
        if Nbacteria_random[i]!=0:
           Temp=np.random.choice(a=x, p=pdf,size=Nbacteria_random[i])
           sets.append(max(Temp))
    return sets
######################################################################################################################
Thet1e0= []; Thet1e1 = []; Thet1e2 = []; Thet1e3 = []; Thet1e4 = []; 
Thet1e5= []; Thet1e6 = []; Thet1e7 = []; Thet1e8 = []; Thet1e9 = []; 
print("\n################### Generating the data")
if nconc==5:
   Thet1e0=N_bacteria_varied(N_bacteria_set_mean=1,N=ndroplets) 
   print("10^0 INs per droplet ")   
   Thet1e1=N_bacteria_varied(N_bacteria_set_mean=10,N=ndroplets) 
   print("10^1 INs per droplet ")   
   Thet1e2=N_bacteria_varied(N_bacteria_set_mean=100,N=ndroplets) 
   print("10^2 INs per droplet ") 
   Thet1e3=N_bacteria_varied(N_bacteria_set_mean=1000,N=ndroplets) 
   print("10^3 INs per droplet ") 
   Thet1e4=N_bacteria_varied(N_bacteria_set_mean=10000,N=ndroplets) 
   print("10^4 INs per droplet ")
if nconc==6:
   Thet1e0=N_bacteria_varied(N_bacteria_set_mean=1,N=ndroplets) 
   print("10^0 INs per droplet ")   
   Thet1e1=N_bacteria_varied(N_bacteria_set_mean=10,N=ndroplets) 
   print("10^1 INs per droplet ")   
   Thet1e2=N_bacteria_varied(N_bacteria_set_mean=100,N=ndroplets) 
   print("10^2 INs per droplet ") 
   Thet1e3=N_bacteria_varied(N_bacteria_set_mean=1000,N=ndroplets) 
   print("10^3 INs per droplet ") 
   Thet1e4=N_bacteria_varied(N_bacteria_set_mean=10000,N=ndroplets) 
   print("10^4 INs per droplet ")
   Thet1e5=N_bacteria_varied(N_bacteria_set_mean=100000,N=ndroplets) 
   print("10^5 INs per droplet ")
if nconc==7:
   Thet1e0=N_bacteria_varied(N_bacteria_set_mean=1,N=ndroplets) 
   print("10^0 INs per droplet ")   
   Thet1e1=N_bacteria_varied(N_bacteria_set_mean=10,N=ndroplets) 
   print("10^1 INs per droplet ")   
   Thet1e2=N_bacteria_varied(N_bacteria_set_mean=100,N=ndroplets) 
   print("10^2 INs per droplet ") 
   Thet1e3=N_bacteria_varied(N_bacteria_set_mean=1000,N=ndroplets) 
   print("10^3 INs per droplet ") 
   Thet1e4=N_bacteria_varied(N_bacteria_set_mean=10000,N=ndroplets) 
   print("10^4 INs per droplet ")
   Thet1e5=N_bacteria_varied(N_bacteria_set_mean=100000,N=ndroplets) 
   print("10^5 INs per droplet ")    
   Thet1e6=N_bacteria_varied(N_bacteria_set_mean=1000000,N=ndroplets) 
   print("10^6 INs per droplet ")  
if nconc==8:
   Thet1e0=N_bacteria_varied(N_bacteria_set_mean=1,N=ndroplets) 
   print("10^0 INs per droplet ")   
   Thet1e1=N_bacteria_varied(N_bacteria_set_mean=10,N=ndroplets) 
   print("10^1 INs per droplet ")   
   Thet1e2=N_bacteria_varied(N_bacteria_set_mean=100,N=ndroplets) 
   print("10^2 INs per droplet ") 
   Thet1e3=N_bacteria_varied(N_bacteria_set_mean=1000,N=ndroplets) 
   print("10^3 INs per droplet ") 
   Thet1e4=N_bacteria_varied(N_bacteria_set_mean=10000,N=ndroplets) 
   print("10^4 INs per droplet ")
   Thet1e5=N_bacteria_varied(N_bacteria_set_mean=100000,N=ndroplets) 
   print("10^5 INs per droplet ")    
   Thet1e6=N_bacteria_varied(N_bacteria_set_mean=1000000,N=ndroplets) 
   print("10^6 INs per droplet ") 
   Thet1e7=N_bacteria_varied(N_bacteria_set_mean=10000000,N=ndroplets) 
   print("10^7 INs per droplet ") 
if nconc==9:
   Thet1e0=N_bacteria_varied(N_bacteria_set_mean=1,N=ndroplets) 
   print("10^0 INs per droplet ")   
   Thet1e1=N_bacteria_varied(N_bacteria_set_mean=10,N=ndroplets) 
   print("10^1 INs per droplet ")   
   Thet1e2=N_bacteria_varied(N_bacteria_set_mean=100,N=ndroplets) 
   print("10^2 INs per droplet ") 
   Thet1e3=N_bacteria_varied(N_bacteria_set_mean=1000,N=ndroplets) 
   print("10^3 INs per droplet ") 
   Thet1e4=N_bacteria_varied(N_bacteria_set_mean=10000,N=ndroplets) 
   print("10^4 INs per droplet ")
   Thet1e5=N_bacteria_varied(N_bacteria_set_mean=100000,N=ndroplets) 
   print("10^5 INs per droplet ")    
   Thet1e6=N_bacteria_varied(N_bacteria_set_mean=1000000,N=ndroplets) 
   print("10^6 INs per droplet ") 
   Thet1e7=N_bacteria_varied(N_bacteria_set_mean=10000000,N=ndroplets) 
   print("10^7 INs per droplet ") 
   Thet1e8=N_bacteria_varied(N_bacteria_set_mean=100000000,N=ndroplets) 
   print("10^8 INs per droplet ") 
print("Done!")
######################################################################################################################
fig1, ax = plt.subplots()
initial=min(modes)-8
#initial=-30
final=0
binwidth=0.1
bins=np.arange(initial, final + binwidth, binwidth)
######################################################################################################################
def histogram_PDF(Thet):
    p,edges = np.histogram(Thet, bins=bins, density=False)
    T = 0.5*(edges[1:]+ edges[:-1])
    norm=np.trapz(p,T)
    Nfrozen=np.sum(p)
    p=p/norm
    return T, p, Nfrozen
######################################################################################################################
T1e0,p1e0,Nfrozen1e0=histogram_PDF(Thet1e0)
ax.plot(T1e0,p1e0,'-s',color='blue',linewidth=3,label='$10^{0}$',fillstyle='none',markersize=8,mew=2)
if nsubpop==1:
   pdf=normalized_PDF(T1e0,param[0],param[1])
if nsubpop==2:
   pdf=(1-param[4])*normalized_PDF(T1e0,param[0],param[1])+param[4]*normalized_PDF(T1e0,param[2],param[3])
if nsubpop==3:
   pdf=(1-param[4]-param[5])*normalized_PDF(T1e0,param[0],param[1])+param[4]*normalized_PDF(T1e0,param[2],param[3])+param[5]*normalized_PDF(T1e0,param[6],param[7])
if nconc==5:
   T1e1,p1e1,Nfrozen1e1=histogram_PDF(Thet1e1)
   ax.plot(T1e1,p1e1,'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   T1e2,p1e2,Nfrozen1e2=histogram_PDF(Thet1e2)
   ax.plot(T1e2,p1e2,'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   T1e3,p1e3,Nfrozen1e3=histogram_PDF(Thet1e3)
   ax.plot(T1e3,p1e3,'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   T1e4,p1e4,Nfrozen1e4=histogram_PDF(Thet1e4)
   ax.plot(T1e4,p1e4,'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("Pu_and_Pmax.txt", np.column_stack((T1e0,pdf,p1e0,p1e1,p1e2,p1e3,p1e4)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4") 
if nconc==6:
   T1e1,p1e1,Nfrozen1e1=histogram_PDF(Thet1e1)
   ax.plot(T1e1,p1e1,'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   T1e2,p1e2,Nfrozen1e2=histogram_PDF(Thet1e2)
   ax.plot(T1e2,p1e2,'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   T1e3,p1e3,Nfrozen1e3=histogram_PDF(Thet1e3)
   ax.plot(T1e3,p1e3,'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   T1e4,p1e4,Nfrozen1e4=histogram_PDF(Thet1e4)
   ax.plot(T1e4,p1e4,'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   T1e5,p1e5,Nfrozen1e5=histogram_PDF(Thet1e5)
   ax.plot(T1e5,p1e5,'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("Pu_and_Pmax.txt", np.column_stack((T1e0,pdf,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5") 
if nconc==7:
   T1e1,p1e1,Nfrozen1e1=histogram_PDF(Thet1e1)
   ax.plot(T1e1,p1e1,'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   T1e2,p1e2,Nfrozen1e2=histogram_PDF(Thet1e2)
   ax.plot(T1e2,p1e2,'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   T1e3,p1e3,Nfrozen1e3=histogram_PDF(Thet1e3)
   ax.plot(T1e3,p1e3,'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   T1e4,p1e4,Nfrozen1e4=histogram_PDF(Thet1e4)
   ax.plot(T1e4,p1e4,'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   T1e5,p1e5,Nfrozen1e5=histogram_PDF(Thet1e5)
   ax.plot(T1e5,p1e5,'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   T1e6,p1e6,Nfrozen1e6=histogram_PDF(Thet1e6)
   ax.plot(T1e6,p1e6,'-^',color='magenta',linewidth=3,label='$10^{6}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("Pu_and_Pmax.txt", np.column_stack((T1e0,pdf,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5,p1e6)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6") 
if nconc==8:
   T1e1,p1e1,Nfrozen1e1=histogram_PDF(Thet1e1)
   ax.plot(T1e1,p1e1,'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   T1e2,p1e2,Nfrozen1e2=histogram_PDF(Thet1e2)
   ax.plot(T1e2,p1e2,'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   T1e3,p1e3,Nfrozen1e3=histogram_PDF(Thet1e3)
   ax.plot(T1e3,p1e3,'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   T1e4,p1e4,Nfrozen1e4=histogram_PDF(Thet1e4)
   ax.plot(T1e4,p1e4,'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   T1e5,p1e5,Nfrozen1e5=histogram_PDF(Thet1e5)
   ax.plot(T1e5,p1e5,'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   T1e6,p1e6,Nfrozen1e6=histogram_PDF(Thet1e6)
   ax.plot(T1e6,p1e6,'-^',color='magenta',linewidth=3,label='$10^{6}$',fillstyle='none',markersize=8,mew=2)
   T1e7,p1e7,Nfrozen1e7=histogram_PDF(Thet1e7)
   ax.plot(T1e7,p1e7,'-o',color='gray',linewidth=3,label='$10^{7}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("Pu_and_Pmax.txt", np.column_stack((T1e0,pdf,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5,p1e7)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7") 
if nconc==9:
   T1e1,p1e1,Nfrozen1e1=histogram_PDF(Thet1e1)
   ax.plot(T1e1,p1e1,'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   T1e2,p1e2,Nfrozen1e2=histogram_PDF(Thet1e2)
   ax.plot(T1e2,p1e2,'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   T1e3,p1e3,Nfrozen1e3=histogram_PDF(Thet1e3)
   ax.plot(T1e3,p1e3,'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   T1e4,p1e4,Nfrozen1e4=histogram_PDF(Thet1e4)
   ax.plot(T1e4,p1e4,'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   T1e5,p1e5,Nfrozen1e5=histogram_PDF(Thet1e5)
   ax.plot(T1e5,p1e5,'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   T1e6,p1e6,Nfrozen1e6=histogram_PDF(Thet1e6)
   ax.plot(T1e6,p1e6,'-^',color='magenta',linewidth=3,label='$10^{6}$',fillstyle='none',markersize=8,mew=2)
   T1e7,p1e7,Nfrozen1e7=histogram_PDF(Thet1e7)
   ax.plot(T1e7,p1e7,'-o',color='gray',linewidth=3,label='$10^{7}$',fillstyle='none',markersize=8,mew=2)
   T1e8,p1e8,Nfrozen1e8=histogram_PDF(Thet1e8)
   ax.plot(T1e8,p1e8,'-s',color='orange',linewidth=3,label='$10^{8}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("Pu_and_Pmax.txt", np.column_stack((T1e0,pdf,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5,p1e7,p1e8)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8") 
ax.plot(T1e0,pdf,'--',color='black',linewidth=4, label='underlying',ms=6)
ax.set_ylabel('$P_u$, $P_{max}$')
ax.set_xlabel('T ($^{o}C$)')
#if nsubpop>=2:
#  if c2 < 0.1:
#    ax.set_yscale('log')
#    ax.set_ylim([c2*10**(-2),2])
ax.legend(loc='best',fontsize=20,frameon=False)
fig1.set_size_inches(10.0, 10.0, forward=True)
fig1.tight_layout()
fig1.savefig("Figure1_Pu_and_Pmax.pdf")
######################################################################################################################
######################################################################################################################
######################################################################################################################
print("\n################### Plotting the fraction of frozen droplets...")
fig2, ax = plt.subplots()
def ice_fraction(PDF,Nfrozen):
    cdf = np.cumsum(PDF)*0.1
    fice=(1-cdf)*(Nfrozen/ndroplets) #*(1-np.exp(-lambdaa))
    fice=np.around(fice,8)
    return fice
if nconc==5:
   ax.plot(T1e0,ice_fraction(p1e0,Nfrozen1e0),'-s',color='blue',linewidth=3,label='$10^{0}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e1,ice_fraction(p1e1,Nfrozen1e1),'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e2,ice_fraction(p1e2,Nfrozen1e2),'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e3,ice_fraction(p1e3,Nfrozen1e3),'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e4,ice_fraction(p1e4,Nfrozen1e4),'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("f_ice.txt", np.column_stack((T1e0,p1e0,p1e1,p1e2,p1e3,p1e4)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4") 
if nconc==6:
   ax.plot(T1e0,ice_fraction(p1e0,Nfrozen1e0),'-s',color='blue',linewidth=3,label='$10^{0}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e1,ice_fraction(p1e1,Nfrozen1e1),'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e2,ice_fraction(p1e2,Nfrozen1e2),'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e3,ice_fraction(p1e3,Nfrozen1e3),'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e4,ice_fraction(p1e4,Nfrozen1e4),'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e5,ice_fraction(p1e5,Nfrozen1e5),'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("f_ice.txt", np.column_stack((T1e0,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5") 
if nconc==7:
   ax.plot(T1e0,ice_fraction(p1e0,Nfrozen1e0),'-s',color='blue',linewidth=3,label='$10^{0}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e1,ice_fraction(p1e1,Nfrozen1e1),'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e2,ice_fraction(p1e2,Nfrozen1e2),'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e3,ice_fraction(p1e3,Nfrozen1e3),'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e4,ice_fraction(p1e4,Nfrozen1e4),'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e5,ice_fraction(p1e5,Nfrozen1e5),'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e6, ice_fraction(p1e6,Nfrozen1e6),'-^',color='magenta',linewidth=3,label='$10^{6}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("f_ice.txt", np.column_stack((T1e0,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5,p1e6)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6") 
if nconc==8:
   ax.plot(T1e0,ice_fraction(p1e0,Nfrozen1e0),'-s',color='blue',linewidth=3,label='$10^{0}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e1,ice_fraction(p1e1,Nfrozen1e1),'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e2,ice_fraction(p1e2,Nfrozen1e2),'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e3,ice_fraction(p1e3,Nfrozen1e3),'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e4,ice_fraction(p1e4,Nfrozen1e4),'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e5,ice_fraction(p1e5,Nfrozen1e5),'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e6,ice_fraction(p1e6,Nfrozen1e6),'-^',color='magenta',linewidth=3,label='$10^{6}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e7,ice_fraction(p1e7,Nfrozen1e7),'-o',color='gray',linewidth=3,label='$10^{7}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("f_ice.txt", np.column_stack((T1e0,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5,p1e7)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7") 
if nconc==9:
   ax.plot(T1e0,ice_fraction(p1e0,Nfrozen1e0),'-s',color='blue',linewidth=3,label='$10^{0}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e1,ice_fraction(p1e1,Nfrozen1e1),'-o',color='cyan',linewidth=3,label='$10^{1}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e2,ice_fraction(p1e2,Nfrozen1e2),'-D',color='green',linewidth=3,label='$10^{2}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e3,ice_fraction(p1e3,Nfrozen1e3),'-x',color='#e5bf05',linewidth=3,label='$10^{3}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e4,ice_fraction(p1e4,Nfrozen1e4),'-^',color='red',linewidth=3,label='$10^{4}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e5,ice_fraction(p1e5,Nfrozen1e5),'-^',color='purple',linewidth=3,label='$10^{5}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e6,ice_fraction(p1e6,Nfrozen1e6),'-^',color='magenta',linewidth=3,label='$10^{6}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e7,ice_fraction(p1e7,Nfrozen1e7),'-o',color='gray',linewidth=3,label='$10^{7}$',fillstyle='none',markersize=8,mew=2)
   ax.plot(T1e8,ice_fraction(p1e8,Nfrozen1e8),'-s',color='orange',linewidth=3,label='$10^{8}$',fillstyle='none',markersize=8,mew=2)
   np.savetxt("f_ice.txt", np.column_stack((T1e0,p1e0,p1e1,p1e2,p1e3,p1e4,p1e5,p1e7,p1e8)), fmt='%10.8f', newline='\n', header="Temperature, underlying, INs per droplet: 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8") 
ax.legend(loc='best',fontsize=20,frameon=False)
ax.set_ylabel('Fraction of ice')
ax.set_xlabel('T ($^{o}C$)')
fig2.set_size_inches(10.0, 10.0, forward=True)
fig2.tight_layout()
fig2.savefig("Figure2_fractionofice.pdf")
######################################################################################################################
######################################################################################################################
######################################################################################################################
fig3, ax = plt.subplots()
np.seterr(divide='ignore')
print("\n################### Plotting the cumulative freezing spectrum...")
print("Density of the initial suspension is 1e-3 grams/1e-3 liters and the droplet volume is 0.1*1e-6 liters.")
normalize=input("Do you want to change the mass of the particles and volume in the initial suspension and droplet volume (use the format 1e3, 1e-3)? Yes or No?: ")
if normalize =='Yes':
  m=input("Mass of the particles in the initial suspension in grams: "); m=float(m)
  Vwash=input("Volume of the initial suspension in liters: "); Vwash=float(Vwash)
  Vdrop=input("Droplet volume in liters: "); Vdrop=float(Vdrop)
else:
  m=1*10**(-3)                              # mass of the particles in the initial suspension in grams
  Vwash=1*10**(-3)                          # volume of the initial suspension in liters 
  Vdrop=0.1*10**(-6)                        # droplet volume in micro-liters

def remove_inf(arr1,arr2):
    new1=np.array([]); new2=np.array([])
    for i in range (0, len(arr2)):
     if arr2[i] != np.inf and arr2[i]>0:
      new1=np.append(new1,arr1[i])
      new2=np.append(new2,arr2[i]) 
    return new1, new2

if nconc==5:
    d=np.array([10**4,10**3,10**2,10**1,10**0]) # dilution factor
    Nm1e0 = -np.log(1 - ice_fraction(p1e0,Nfrozen1e0))*(Vwash/Vdrop)*(d[0]/m)
    T1e0,Nm1e0=remove_inf(T1e0,Nm1e0)
    ax.plot(T1e0,Nm1e0,'s',color='blue',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{0}$')
    Nm1e1 = -np.log(1 - ice_fraction(p1e1,Nfrozen1e1))*(Vwash/Vdrop)*(d[1]/m)
    T1e1,Nm1e1=remove_inf(T1e1,Nm1e1)
    ax.plot(T1e1,Nm1e1,'o',color='cyan',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{1}$')
    Nm1e2 = -np.log(1 - ice_fraction(p1e2,Nfrozen1e2))*(Vwash/Vdrop)*(d[2]/m)
    T1e2,Nm1e2=remove_inf(T1e2,Nm1e2)
    ax.plot(T1e2,Nm1e2,'D',color='green',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{2}$')
    Nm1e3 = -np.log(1 - ice_fraction(p1e3,Nfrozen1e3))*(Vwash/Vdrop)*(d[3]/m)
    T1e3,Nm1e3=remove_inf(T1e3,Nm1e3)
    ax.plot(T1e3,Nm1e3,'x',color='#e5bf05',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{3}$')
    Nm1e4 = -np.log(1 - ice_fraction(p1e4,Nfrozen1e4))*(Vwash/Vdrop)*(d[4]/m)
    T1e4,Nm1e4=remove_inf(T1e4,Nm1e4)
    ax.plot(T1e4,Nm1e4,'^',color='red',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{4}$')
    Nm_all = np.concatenate((Nm1e0,Nm1e1,Nm1e2,Nm1e3,Nm1e4))
    T_all = np.concatenate((T1e0,T1e1,T1e2,T1e3,T1e4))    
    np.savetxt("Nm_1e0_INperDroplet.txt", np.column_stack((T1e0, Nm1e0)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1 IN per droplet') 
    np.savetxt("Nm_1e1_INperDroplet.txt", np.column_stack((T1e1, Nm1e1)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10 INs per droplet') 
    np.savetxt("Nm_1e2_INperDroplet.txt", np.column_stack((T1e2, Nm1e2)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100 INs per droplet') 
    np.savetxt("Nm_1e3_INperDroplet.txt", np.column_stack((T1e3, Nm1e3)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000 INs per droplet') 
    np.savetxt("Nm_1e4_INperDroplet.txt", np.column_stack((T1e4, Nm1e4)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000 INs per droplet') 
if nconc==6:
    d=np.array([10**5,10**4,10**3,10**2,10**1,10**0]) # dilution factor
    Nm1e0 = -np.log(1 - ice_fraction(p1e0,Nfrozen1e0))*(Vwash/Vdrop)*(d[0]/m)
    T1e0,Nm1e0=remove_inf(T1e0,Nm1e0)
    ax.plot(T1e0,Nm1e0,'s',color='blue',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{0}$')
    Nm1e1 = -np.log(1 - ice_fraction(p1e1,Nfrozen1e1))*(Vwash/Vdrop)*(d[1]/m)
    T1e1,Nm1e1=remove_inf(T1e1,Nm1e1)
    ax.plot(T1e1,Nm1e1,'o',color='cyan',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{1}$')
    Nm1e2 = -np.log(1 - ice_fraction(p1e2,Nfrozen1e2))*(Vwash/Vdrop)*(d[2]/m)
    T1e2,Nm1e2=remove_inf(T1e2,Nm1e2)
    ax.plot(T1e2,Nm1e2,'D',color='green',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{2}$')
    Nm1e3 = -np.log(1 - ice_fraction(p1e3,Nfrozen1e3))*(Vwash/Vdrop)*(d[3]/m)
    T1e3,Nm1e3=remove_inf(T1e3,Nm1e3)
    ax.plot(T1e3,Nm1e3,'x',color='#e5bf05',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{3}$')
    Nm1e4 = -np.log(1 - ice_fraction(p1e4,Nfrozen1e4))*(Vwash/Vdrop)*(d[4]/m)
    T1e4,Nm1e4=remove_inf(T1e4,Nm1e4)
    ax.plot(T1e4,Nm1e4,'^',color='red',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{4}$')
    Nm1e5 = -np.log(1 - ice_fraction(p1e5,Nfrozen1e5))*(Vwash/Vdrop)*(d[5]/m)
    T1e5,Nm1e5=remove_inf(T1e5,Nm1e5)
    ax.plot(T1e5,Nm1e5,'^',color='purple',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{5}$')
    Nm_all = np.concatenate((Nm1e0,Nm1e1,Nm1e2,Nm1e3,Nm1e4))
    T_all = np.concatenate((T1e0,T1e1,T1e2,T1e3,T1e4)) 
    np.savetxt("Nm_1e0_INperDroplet.txt", np.column_stack((T1e0, Nm1e0)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1 IN per droplet') 
    np.savetxt("Nm_1e1_INperDroplet.txt", np.column_stack((T1e1, Nm1e1)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10 INs per droplet') 
    np.savetxt("Nm_1e2_INperDroplet.txt", np.column_stack((T1e2, Nm1e2)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100 INs per droplet') 
    np.savetxt("Nm_1e3_INperDroplet.txt", np.column_stack((T1e3, Nm1e3)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000 INs per droplet') 
    np.savetxt("Nm_1e4_INperDroplet.txt", np.column_stack((T1e4, Nm1e4)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000 INs per droplet') 
    np.savetxt("Nm_1e5_INperDroplet.txt", np.column_stack((T1e5, Nm1e5)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100000 INs per droplet') 
if nconc==7:
    d=np.array([10**6,10**5,10**4,10**3,10**2,10**1,10**0]) # dilution factor
    Nm1e0 = -np.log(1 - ice_fraction(p1e0,Nfrozen1e0))*(Vwash/Vdrop)*(d[0]/m)
    T1e0,Nm1e0=remove_inf(T1e0,Nm1e0)
    ax.plot(T1e0,Nm1e0,'s',color='blue',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{0}$')
    Nm1e1 = -np.log(1 - ice_fraction(p1e1,Nfrozen1e1))*(Vwash/Vdrop)*(d[1]/m)
    T1e1,Nm1e1=remove_inf(T1e1,Nm1e1)
    ax.plot(T1e1,Nm1e1,'o',color='cyan',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{1}$')
    Nm1e2 = -np.log(1 - ice_fraction(p1e2,Nfrozen1e2))*(Vwash/Vdrop)*(d[2]/m)
    T1e2,Nm1e2=remove_inf(T1e2,Nm1e2)
    ax.plot(T1e2,Nm1e2,'D',color='green',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{2}$')
    Nm1e3 = -np.log(1 - ice_fraction(p1e3,Nfrozen1e3))*(Vwash/Vdrop)*(d[3]/m)
    T1e3,Nm1e3=remove_inf(T1e3,Nm1e3)
    ax.plot(T1e3,Nm1e3,'x',color='#e5bf05',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{3}$')
    Nm1e4 = -np.log(1 - ice_fraction(p1e4,Nfrozen1e4))*(Vwash/Vdrop)*(d[4]/m)
    T1e4,Nm1e4=remove_inf(T1e4,Nm1e4)
    ax.plot(T1e4,Nm1e4,'^',color='red',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{4}$')
    Nm1e5 = -np.log(1 - ice_fraction(p1e5,Nfrozen1e5))*(Vwash/Vdrop)*(d[5]/m)
    T1e5,Nm1e5=remove_inf(T1e5,Nm1e5)
    ax.plot(T1e5,Nm1e5,'^',color='purple',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{5}$')
    Nm1e6 = -np.log(1 - ice_fraction(p1e6,Nfrozen1e6))*(Vwash/Vdrop)*(d[6]/m)
    T1e6,Nm1e6=remove_inf(T1e6,Nm1e6)
    ax.plot(T1e6,Nm1e6,'^',color='magenta',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{6}$')
    Nm_all = np.concatenate((Nm1e0,Nm1e1,Nm1e2,Nm1e3,Nm1e4,Nm1e5,Nm1e6))
    T_all = np.concatenate((T1e0,T1e1,T1e2,T1e3,T1e4,T1e5,T1e6))    
    np.savetxt("Nm_1e0_INperDroplet.txt", np.column_stack((T1e0, Nm1e0)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1 IN per droplet') 
    np.savetxt("Nm_1e1_INperDroplet.txt", np.column_stack((T1e1, Nm1e1)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10 INs per droplet') 
    np.savetxt("Nm_1e2_INperDroplet.txt", np.column_stack((T1e2, Nm1e2)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100 INs per droplet') 
    np.savetxt("Nm_1e3_INperDroplet.txt", np.column_stack((T1e3, Nm1e3)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000 INs per droplet') 
    np.savetxt("Nm_1e4_INperDroplet.txt", np.column_stack((T1e4, Nm1e4)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000 INs per droplet') 
    np.savetxt("Nm_1e5_INperDroplet.txt", np.column_stack((T1e5, Nm1e5)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100000 INs per droplet') 
    np.savetxt("Nm_1e6_INperDroplet.txt", np.column_stack((T1e6, Nm1e6)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000000 INs per droplet') 
if nconc==8:
    d=np.array([10**7,10**6,10**5,10**4,10**3,10**2,10**1,10**0]) # dilution factor
    Nm1e0 = -np.log(1 - ice_fraction(p1e0,Nfrozen1e0))*(Vwash/Vdrop)*(d[0]/m)
    T1e0,Nm1e0=remove_inf(T1e0,Nm1e0)
    ax.plot(T1e0,Nm1e0,'s',color='blue',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{0}$')
    Nm1e1 = -np.log(1 - ice_fraction(p1e1,Nfrozen1e1))*(Vwash/Vdrop)*(d[1]/m)
    T1e1,Nm1e1=remove_inf(T1e1,Nm1e1)
    ax.plot(T1e1,Nm1e1,'o',color='cyan',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{1}$')
    Nm1e2 = -np.log(1 - ice_fraction(p1e2,Nfrozen1e2))*(Vwash/Vdrop)*(d[2]/m)
    T1e2,Nm1e2=remove_inf(T1e2,Nm1e2)
    ax.plot(T1e2,Nm1e2,'D',color='green',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{2}$')
    Nm1e3 = -np.log(1 - ice_fraction(p1e3,Nfrozen1e3))*(Vwash/Vdrop)*(d[3]/m)
    T1e3,Nm1e3=remove_inf(T1e3,Nm1e3)
    ax.plot(T1e3,Nm1e3,'x',color='#e5bf05',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{3}$')
    Nm1e4 = -np.log(1 - ice_fraction(p1e4,Nfrozen1e4))*(Vwash/Vdrop)*(d[4]/m)
    T1e4,Nm1e4=remove_inf(T1e4,Nm1e4)
    ax.plot(T1e4,Nm1e4,'^',color='red',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{4}$')
    Nm1e5 = -np.log(1 - ice_fraction(p1e5,Nfrozen1e5))*(Vwash/Vdrop)*(d[5]/m)
    T1e5,Nm1e5=remove_inf(T1e5,Nm1e5)
    ax.plot(T1e5,Nm1e5,'^',color='purple',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{5}$')
    Nm1e6 = -np.log(1 - ice_fraction(p1e6,Nfrozen1e6))*(Vwash/Vdrop)*(d[6]/m)
    T1e6,Nm1e6=remove_inf(T1e6,Nm1e6)
    ax.plot(T1e6,Nm1e6,'^',color='magenta',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{6}$')
    Nm1e7 = -np.log(1 - ice_fraction(p1e7,Nfrozen1e7))*(Vwash/Vdrop)*(d[7]/m)
    T1e7,Nm1e7=remove_inf(T1e7,Nm1e7)
    ax.plot(T1e7,Nm1e7,'o',color='gray',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{7}$')
    Nm_all = np.concatenate((Nm1e0,Nm1e1,Nm1e2,Nm1e3,Nm1e4,Nm1e5,Nm1e6,Nm1e7))
    T_all = np.concatenate((T1e0,T1e1,T1e2,T1e3,T1e4,T1e5,T1e6,T1e7))  
    np.savetxt("Nm_1e0_INperDroplet.txt", np.column_stack((T1e0, Nm1e0)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1 IN per droplet') 
    np.savetxt("Nm_1e1_INperDroplet.txt", np.column_stack((T1e1, Nm1e1)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10 INs per droplet') 
    np.savetxt("Nm_1e2_INperDroplet.txt", np.column_stack((T1e2, Nm1e2)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100 INs per droplet') 
    np.savetxt("Nm_1e3_INperDroplet.txt", np.column_stack((T1e3, Nm1e3)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000 INs per droplet') 
    np.savetxt("Nm_1e4_INperDroplet.txt", np.column_stack((T1e4, Nm1e4)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000 INs per droplet') 
    np.savetxt("Nm_1e5_INperDroplet.txt", np.column_stack((T1e5, Nm1e5)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100000 INs per droplet') 
    np.savetxt("Nm_1e6_INperDroplet.txt", np.column_stack((T1e6, Nm1e6)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000000 INs per droplet') 
    np.savetxt("Nm_1e7_INperDroplet.txt", np.column_stack((T1e7, Nm1e7)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000000 INs per droplet') 
if nconc==9:
    d=np.array([10**8,10**7,10**6,10**5,10**4,10**3,10**2,10**1,10**0]) # dilution factor
    Nm1e0 = -np.log(1 - ice_fraction(p1e0,Nfrozen1e0))*(Vwash/Vdrop)*(d[0]/m)
    T1e0,Nm1e0=remove_inf(T1e0,Nm1e0)
    ax.plot(T1e0,Nm1e0,'s',color='blue',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{0}$')
    Nm1e1 = -np.log(1 - ice_fraction(p1e1,Nfrozen1e1))*(Vwash/Vdrop)*(d[1]/m)
    T1e1,Nm1e1=remove_inf(T1e1,Nm1e1)
    ax.plot(T1e1,Nm1e1,'o',color='cyan',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{1}$')
    Nm1e2 = -np.log(1 - ice_fraction(p1e2,Nfrozen1e2))*(Vwash/Vdrop)*(d[2]/m)
    T1e2,Nm1e2=remove_inf(T1e2,Nm1e2)
    ax.plot(T1e2,Nm1e2,'D',color='green',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{2}$')
    Nm1e3 = -np.log(1 - ice_fraction(p1e3,Nfrozen1e3))*(Vwash/Vdrop)*(d[3]/m)
    T1e3,Nm1e3=remove_inf(T1e3,Nm1e3)
    ax.plot(T1e3,Nm1e3,'x',color='#e5bf05',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{3}$')
    Nm1e4 = -np.log(1 - ice_fraction(p1e4,Nfrozen1e4))*(Vwash/Vdrop)*(d[4]/m)
    T1e4,Nm1e4=remove_inf(T1e4,Nm1e4)
    ax.plot(T1e4,Nm1e4,'^',color='red',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{4}$')
    Nm1e5 = -np.log(1 - ice_fraction(p1e5,Nfrozen1e5))*(Vwash/Vdrop)*(d[5]/m)
    T1e5,Nm1e5=remove_inf(T1e5,Nm1e5)
    ax.plot(T1e5,Nm1e5,'^',color='purple',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{5}$')
    Nm1e6 = -np.log(1 - ice_fraction(p1e6,Nfrozen1e6))*(Vwash/Vdrop)*(d[6]/m)
    T1e6,Nm1e6=remove_inf(T1e6,Nm1e6)
    ax.plot(T1e6,Nm1e6,'^',color='magenta',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{6}$')
    Nm1e7 = -np.log(1 - ice_fraction(p1e7,Nfrozen1e7))*(Vwash/Vdrop)*(d[7]/m)
    T1e7,Nm1e7=remove_inf(T1e7,Nm1e7)
    ax.plot(T1e7,Nm1e7,'o',color='gray',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{7}$')
    Nm1e8 = -np.log(1 - ice_fraction(p1e8,Nfrozen1e8))*(Vwash/Vdrop)*(d[8]/m)
    T1e8,Nm1e8=remove_inf(T1e8,Nm1e8)
    ax.plot(T1e8,Nm1e8,'o',color='orange',linewidth=3,fillstyle='none',markersize=10,mew=2,label='$10^{8}$')
    Nm_all = np.concatenate((Nm1e0,Nm1e1,Nm1e2,Nm1e3,Nm1e4,Nm1e5,Nm1e6,Nm1e7,Nm1e8))
    T_all = np.concatenate((T1e0,T1e1,T1e2,T1e3,T1e4,T1e5,T1e6,T1e7,T1e8))  
    np.savetxt("Nm_1e0_INperDroplet.txt", np.column_stack((T1e0, Nm1e0)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1 IN per droplet') 
    np.savetxt("Nm_1e1_INperDroplet.txt", np.column_stack((T1e1, Nm1e1)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10 INs per droplet') 
    np.savetxt("Nm_1e2_INperDroplet.txt", np.column_stack((T1e2, Nm1e2)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100 INs per droplet') 
    np.savetxt("Nm_1e3_INperDroplet.txt", np.column_stack((T1e3, Nm1e3)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000 INs per droplet') 
    np.savetxt("Nm_1e4_INperDroplet.txt", np.column_stack((T1e4, Nm1e4)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000 INs per droplet') 
    np.savetxt("Nm_1e5_INperDroplet.txt", np.column_stack((T1e5, Nm1e5)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100000 INs per droplet') 
    np.savetxt("Nm_1e6_INperDroplet.txt", np.column_stack((T1e6, Nm1e6)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 1000000 INs per droplet') 
    np.savetxt("Nm_1e7_INperDroplet.txt", np.column_stack((T1e7, Nm1e7)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 10000000 INs per droplet') 
    np.savetxt("Nm_1e8_INperDroplet.txt", np.column_stack((T1e8, Nm1e8)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1), 100000000 INs per droplet') 
   
ax.set_yscale('log')
ax.legend(loc='best',fontsize=20,ncol=1,frameon=False)
locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,numticks=100)
ax.yaxis.set_major_locator(locmaj)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylabel('N$_m$ ($g^{-1}$)')
ax.set_xlabel('T ($^{o}C$)')
fig3.set_size_inches(10.0, 10.0, forward=True)
fig3.tight_layout()
fig3.savefig("Figure3_Nm_spectrum.pdf")
np.savetxt("Nm_all.txt", np.column_stack((T_all, Nm_all)), fmt='%10.8f', newline='\n',header='Temperature, Nm (g^-1)') 
print("\n################### Done!")
print("Check for files in your directory.")
plt.show()
