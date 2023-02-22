import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import inf
import numpy as np
from scipy.optimize import minimize
from random import seed
from random import random
import random
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import matplotlib.ticker 
from operator import itemgetter
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import dual_annealing
import os
import warnings
warnings.filterwarnings('ignore')
seed(12334) 
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
def normalized_PDF(x,mode,scale,disttype):
    if disttype==1:
       return (1/(scale*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mode)/scale)**2)                  # Gaussian distribution
    if disttype==2:
       return (1/((x-mode)*scale*np.sqrt(2*np.pi)))*(np.exp(-0.5*(np.log(x-mode)/scale)**2)) # Log-normal distribution
    if disttype==3:
       return (1/scale)*np.exp((x-mode)/scale - np.exp((x-mode)/scale))                      # Gumbel with left-tail
######################################################################################################################
######################################################################################################################
######################################################################################################################
print("In the HUB-backward code we obtain the differential freezing spectrum (n_m) that better fits the target data.")
target=0
check=True
while(check):
   target=input("What is the target data? \n1 - cumulative freezing spectrum (N_m)\n2 - fraction of frozen droplets (f_ice)\n"); target=int(target)
   if target==1:
      fice=False; check=False
   if target==2:
      fice=True; check=False
dir_list = os.listdir('./input'); check=False
while(check is False):
  file_name=input("Enter the file's name (file format is .txt, the 1st column is temperature and the second is Nm or fice): ")
  check=file_name in dir_list
  if check is False:
     print("\nCould not find the file. Try again.")
     print("List of files: "+str(dir_list))
file = file_name
print("\nThe data is interpolated to have equally spaced points. A spline fit and a Savitzky–Golay filter are used.")
check=True
while(check):
    numberofpoints=input("Do you want to change the number of points (default is 100)? Yes or No? ")
    if numberofpoints == 'No':
        check=False; npoints=100
    if numberofpoints == 'Yes':
        check=False; 
        npoints=input("How many? "); npoints=int(npoints)
check=True
while(check):
    windowlength=input("Do you want to change the window length (default is 3)? Yes or No? ")
    if windowlength == 'No':
        check=False; window_length=3
    if windowlength == 'Yes':
        check=False; 
        window_length=input("What is the new window length? "); window_length=int(window_length)
check=True
while(check):
    poly=input("Do you want to change the polynomial order (default is 1)? Yes or No? ")
    if poly == 'No':
        check=False; polyorder=1
    if poly == 'Yes':
        check=False; 
        polyorder=input("What is the new polynomial order? "); polyorder=int(polyorder)
check=True
while(check):
    distchoice=input("Do you want to change the type of distribution? Default is Gaussian. Yes or No? ")
    if distchoice == 'No':
        check=False; disttype=1
    if distchoice == 'Yes':
        check=False; 
        disttype=input("Type 2 for Log-normal or 3 for Gumbel with left tail: "); disttype=int(disttype) 
nsubpop=0
while nsubpop<1 or nsubpop>3:
   nsubpop=input("Number of subpopulations (1, 2 or 3): "); nsubpop=int(nsubpop)
print("\nDo you want to change the standard bounds for parameters?")
if nsubpop==1: 
    print("The default is Tmode #1=[-30,0], spread #1=[0.1,10]")
if nsubpop==2: 
    print("The default is Tmode #1=[-30,0], spread #1=[0.1,10], Tmode #2=[-30,0], spread #2=[0.1,10], weight #2=[1e-10,1]")
if nsubpop==3: 
    if file_name=='fractionofice_cholesterol_fig10.txt':
       print("The default is Tmode #1=[-17,-15], spread #1=[0.1,1], Tmode #2=[-12,-9], spread #2=[0.1,1], weight #2=[1e-3,1], Tmode #3=[-9,6], spread #3=[0.1,1], weight #3=[1e-3,1]")
    elif file_name=='fractionofice_cholesterol_fig13.txt':
       print("The default is Tmode #1=[-30,-12], spread #1=[0.1,0.8], Tmode #2=[-12,-9], spread #2=[0.1,1], weight #2=[1e-3,1], Tmode #3=[-9,6], spread #3=[0.1,0.8], weight #3=[1e-3,1]") 
    elif file_name=='Nm_pollen_thesis.txt':
          print("The default is Tmode #1=[-20,-15], spread #1=[0.1,5], Tmode #2=[-16,-11], spread #2=[0.1,2], weight #2=[1e-8,1e-3], Tmode #3=[-11,-6.0], spread #3=[0.1,2], weight #3=[1e-8,1e-3]")
    else: 
       print("The default is Tmode #1=[-30,0], spread #1=[0.1,10], Tmode #2=[-30,0], spread #2=[0.1,10], weight #2=[1e-10,1], Tmode #3=[-30,0], spread #3=[0.1,10], weight #3=[1e-10,1e-2]")
check=True
while(check):
    changebnds=input("Yes or No? ")
    if changebnds == 'No':
        check=False
    if changebnds == 'Yes':
        check=False
if changebnds == 'No':
    if nsubpop==1:
       bnds = [[-30,0],[0.1,10],[1e-3,1]]
    if nsubpop==2:
       bnds = [[-30,0],[0.1,10],[-30,0],[0.1,10],[1e-10,1],[1e-3,1]]
    if nsubpop==3:
       bnds = [[-30,0],[0.1,10],[-30,0],[0.1,10],[1e-10,1],[-30,0],[0.1,10],[1e-10,1e-2],[1e-3,1]]
       if file_name=='fractionofice_cholesterol_fig10.txt':
          bnds =[[-17,-15],[0.1,1],[-12,-9],[0.1,1],[1e-3,1],[-9,6],[0.1,1],[1e-3,1]]
       if file_name=='fractionofice_cholesterol_fig13.txt':
          bnds =[[-30,-12],[0.1,0.8],[-12,-9],[0.1,1],[1e-3,1],[-9,6],[0.1,0.8],[1e-3,1]]
       if file_name=='Nm_pollen_thesis.txt':
          bnds =[[-20,-15],[0.1,5],[-16,-11],[0.1,2],[1e-8,0.001],[-11,-6.0],[0.1,2],[1e-8,1e-3],[1e-3,1]]
else:
    if nsubpop==1:
       Tmode1_l=input("Tmode #1 lower bound: "); Tmode1_l=float(Tmode1_l)
       Tmode1_u=input("Tmode #1 upper bound: "); Tmode1_u=float(Tmode1_u)
       spread1_l=input("Spread #1 lower bound: "); spread1_l=float(spread1_l)
       spread1_u=input("Spread #1 upper bound: "); spread1_u=float(spread1_u)
       bnds=[[Tmode1_l,Tmode1_u],[spread1_l,spread1_u],[1e-3,1]]
    if nsubpop==2:
       Tmode1_l=input("Tmode #1 lower bound: "); Tmode1_l=float(Tmode1_l)
       Tmode1_u=input("Tmode #1 upper bound: "); Tmode1_u=float(Tmode1_u)
       spread1_l=input("Spread #1 lower bound: "); spread1_l=float(spread1_l)
       spread1_u=input("Spread #1 upper bound: "); spread1_u=float(spread1_u)
       Tmode2_l=input("Tmode #2 lower bound: "); Tmode2_l=float(Tmode2_l)
       Tmode2_u=input("Tmode #2 upper bound: "); Tmode2_u=float(Tmode2_u)
       spread2_l=input("Spread #2 lower bound: "); spread2_l=float(spread2_l)
       spread2_u=input("Spread #2 upper bound: "); spread2_u=float(spread2_u)
       weight2_l=input("Weight #2 lower bound: "); weight2_l=float(weight2_l)
       weight2_u=input("Weight #2 upper bound: "); weight2_u=float(weight2_u)
       bnds=[[Tmode1_l,Tmode1_u],[spread1_l,spread1_u],[Tmode2_l,Tmode2_u],[spread2_l,spread2_u],[weight2_l,weight2_u],[1e-3,1]]
    if nsubpop==3:
       Tmode1_l=input("Tmode #1 lower bound: "); Tmode1_l=float(Tmode1_l)
       Tmode1_u=input("Tmode #1 upper bound: "); Tmode1_u=float(Tmode1_u)
       spread1_l=input("Spread #1 lower bound: "); spread1_l=float(spread1_l)
       spread1_u=input("Spread #1 upper bound: "); spread1_u=float(spread1_u)
       Tmode2_l=input("Tmode #2 lower bound: "); Tmode2_l=float(Tmode2_l)
       Tmode2_u=input("Tmode #2 upper bound: "); Tmode2_u=float(Tmode2_u)
       spread2_l=input("Spread #2 lower bound: "); spread2_l=float(spread2_l)
       spread2_u=input("Spread #2 upper bound: "); spread2_u=float(spread2_u)
       weight2_l=input("Weight #2 lower bound: "); weight2_l=float(weight2_l)
       weight2_u=input("Weight #2 upper bound: "); weight2_u=float(weight2_u)
       Tmode3_l=input("Tmode #3 lower bound: "); Tmode3_l=float(Tmode3_l)
       Tmode3_u=input("Tmode #3 upper bound: "); Tmode3_u=float(Tmode3_u)
       spread3_l=input("Spread #3 lower bound: "); spread3_l=float(spread3_l)
       spread3_u=input("Spread #3 upper bound: "); spread3_u=float(spread3_u)
       weight3_l=input("Weight #3 lower bound: "); weight3_l=float(weight3_l)
       weight3_u=input("Weight #3 upper bound: "); weight3_u=float(weight3_u)
       bnds=[[Tmode1_l,Tmode1_u],[spread1_l,spread1_u],[Tmode2_l,Tmode2_u],[spread2_l,spread2_u],[weight2_l,weight2_u],[Tmode3_l,Tmode3_u],[spread3_l,spread3_u],[weight3_l,weight3_u],[1e-3,1]]
print("\nDo you want to give an initial guess? ")
check=True
while(check):
    initialguess=input("Yes or No? ")
    if initialguess == 'No':
        check=False; initial_guess=False
    if initialguess == 'Yes':
        check=False; initial_guess=True
if initial_guess is True:
    if nsubpop==1:
       Tmode1_0=input("Tmode #1: "); Tmode1_0=float(Tmode1_0)
       spread1_0=input("Spread #1: "); spread1_0=float(spread1_0)
       p0=[Tmode1_0, spread1_0, 0.6]
    if nsubpop==2:
       Tmode1_0=input("Tmode #1: "); Tmode1_0=float(Tmode1_0)
       spread1_0=input("Spread #1: "); spread1_0=float(spread1_0)
       Tmode2_0=input("Tmode #2: "); Tmode2_0=float(Tmode2_0)
       spread2_0=input("Spread #2: "); spread2_0=float(spread2_0)
       weight2_0=input("Weight #2: "); weight2_0=float(weight2_0)
       p0=[Tmode1_0, spread1_0, Tmode2_0, spread2_0, weight2_0, 0.6]
    if nsubpop==3:
       Tmode1_0=input("Tmode #1: "); Tmode1_0=float(Tmode1_0)
       spread1_0=input("Spread #1: "); spread1_0=float(spread1_0)
       Tmode2_0=input("Tmode #2: "); Tmode2_0=float(Tmode2_0)
       spread2_0=input("Spread #2: "); spread2_0=float(spread2_0)
       weight2_0=input("Weight #2: "); weight2_0=float(weight2_0)
       Tmode3_0=input("Tmode #3: "); Tmode3_0=float(Tmode2_0)
       spread3_0=input("Spread #3: "); spread3_0=float(spread3_0)
       weight3_0=input("Weight #3: "); weight3_0=float(weight3_0)
       p0=[Tmode1_0, spread1_0, Tmode2_0, spread2_0, weight2_0, Tmode3_0, spread3_0, weight3_0, 0.6]
######################################################################################################################
######################################################################################################################
######################################################################################################################
# reading the input data
initial=0
data = np.loadtxt("input/"+str(file),usecols=[0,1])
nsteps=data.shape[0]
Temps_exp = np.copy(data[initial:nsteps,0])#-273.15
data_exp = np.copy(data[initial:nsteps,1])
bothlists = list(zip(Temps_exp, data_exp))
bothlists.sort(key = itemgetter(0)) # sort according to the first list info
Temps_exp,data_exp =list(zip(*bothlists))
######################################################################################################################
######################################################################################################################
######################################################################################################################
fig1, ax = plt.subplots()
ax.plot(Temps_exp,data_exp,'o',color='black',linewidth=3,label='input data',fillstyle='none',markersize=8,mew=2)
##### interpolating the input data to produce equally spaced data points
f = interp1d(Temps_exp, data_exp)
temp_fit = np.linspace(min(Temps_exp), max(Temps_exp), num=npoints)
delta_temp = (temp_fit[1]-temp_fit[0])
data_fit=f(temp_fit)
data_fit_smooth = savgol_filter(data_fit, window_length, polyorder) 
ax.plot(temp_fit,data_fit_smooth, '-s',color='magenta',label='spline',lw=2,ms=4)
# generating the data for the initial guess
if initial_guess is True:
   if nsubpop==3:
       param=p0
       if fice is False:
          beta=param[8]
       if fice is True:
          beta=1
       #Prob = (1-param[4]-param[7])*  normalized_PDF(temp_fit,param[0],param[1])+param[4]*  normalized_PDF(temp_fit,param[2],param[3])+param[7]*  normalized_PDF(temp_fit,param[5],param[6])
       pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
       pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
       pdf3=normalized_PDF(temp_fit,param[5],param[6],disttype); ii = np.isnan(pdf3); pdf3[ii] = 0
       Prob=(1-param[4]-param[7])*pdf1+param[4]*pdf2+param[7]*pdf3
   if nsubpop==2:
       param=p0
       if fice is False:
          beta=param[5]
       if fice is True:
          beta=1
       #Prob = (1-param[4])*  normalized_PDF(temp_fit,param[0],param[1])+param[4]*  normalized_PDF(temp_fit,param[2],param[3]) 
       pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
       pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
       Prob=(1-param[4])*pdf1+param[4]*pdf2
   if nsubpop==1:
       param=p0
       if fice is False:
          beta=param[2]
       if fice is True:
          beta=1
       Prob =   normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(Prob); Prob[ii] = 0
   integral_normalized=np.cumsum(Prob)*delta_temp
   if fice is False:
      #g=beta*(1-integral_normalized)
      g=-np.log(1-beta*(1-integral_normalized))
      Nm_trial=g*max(data_exp)
      ax.plot(temp_fit,Nm_trial,'-D',color='blue',linewidth=3, label='initial guess')
   if fice is True:
      f=beta*(1-integral_normalized/max(integral_normalized))
      ax.plot(temp_fit,f,'-D',color='blue',linewidth=3, label='initial guess')
######################################################################################################################
if fice is False:
   ax.set_yscale('log')
   locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
   locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,numticks=100)
   ax.yaxis.set_major_locator(locmaj)
   ax.yaxis.set_minor_locator(locmin)
   ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())   
   ax.set_ylabel('N$_m$ ($g^{-1}$)')
if fice is True:
   ax.set_ylabel('Fraction of ice')
ax.legend(loc='best',fontsize=20,ncol=1)
ax.set_xlabel('T ($^{o}C$)')
ax.grid(True,which="both")
fig1.set_size_inches(10.0, 10.0, forward=True)
fig1.tight_layout()
fig1.savefig("Figure1_target_spline.pdf")
######################################################################################################################
######################################################################################################################
######################################################################################################################
# defining the objetive function 
scores=[]
def objective_function(param):
   if nsubpop==3:
       if fice is False:
          beta=param[8]
       if fice is True:
          beta=1
       #Prob = (1-param[4]-param[7])*normalized_PDF(temp_fit,param[0],param[1])+param[4]*  normalized_PDF(temp_fit,param[2],param[3])+param[7]*  normalized_PDF(temp_fit,param[5],param[6])
       pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
       pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
       pdf3=normalized_PDF(temp_fit,param[5],param[6],disttype); ii = np.isnan(pdf3); pdf3[ii] = 0
       Prob=(1-param[4]-param[7])*pdf1+param[4]*pdf2+param[7]*pdf3
   if nsubpop==2:
       if fice is False:
          beta=param[5]
       if fice is True:
          beta=1
       #Prob = (1-param[4])*normalized_PDF(temp_fit,param[0],param[1])+param[4]*  normalized_PDF(temp_fit,param[2],param[3]) 
       pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
       pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
       Prob=(1-param[4])*pdf1+param[4]*pdf2
   if nsubpop==1:
       if fice is False:
          beta=param[2]
       if fice is True:
          beta=1
       Prob = normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(Prob); Prob[ii] = 0
            
   integral_normalized=np.cumsum(Prob)*(temp_fit[1]-temp_fit[0])
   if fice is False:
      #g=beta*(1-integral_normalized)
      g=-np.log(1-beta*(1-integral_normalized))
      Nm_trial=g*max(data_exp)
      deltaNm=abs(np.log10(data_fit_smooth)-np.log10(Nm_trial))
      new_deltaNm2=[]
      for i in range(0,len(Nm_trial)):
         if deltaNm[i]!=np.inf and int(np.isnan(deltaNm[i])) != 1 and Nm_trial[i]>min(data_fit_smooth):
            new_deltaNm2.append(deltaNm[i]**2)
      if len(new_deltaNm2)!=0:
         score=(1/len(new_deltaNm2))*np.sum(new_deltaNm2) # Mean Squared Error (MSE) 
         scores.append(score) 
      else:
         score=100
      if len(new_deltaNm2)<len(temp_fit) and nsubpop==3: ####
         score=100
   if fice is True:
      f=beta*(1-integral_normalized/max(integral_normalized))
      deltaf=abs(data_fit_smooth-f)
      new_deltaf2=[]
      for i in range(0,len(f)):
         if deltaf[i]!=np.inf and int(np.isnan(deltaf[i])) != 1:
            new_deltaf2.append(deltaf[i]**2)
      if len(new_deltaf2)!=0:
         score=(1/len(new_deltaf2))*np.sum(new_deltaf2) # Mean Squared Error (MSE) 
         scores.append(score)
      else:
         score=100
   return score  
######################################################################################################################
print("\n################### Summary ")
print('Input file: '+str(file))
print('Defined bounds: '+str(bnds))
print('Initial guess? '+str(initialguess))
random_seed = np.random.randint(low=0,high=2**32-1,dtype=np.int64)
print('Seed for the dual annealing search: '+str(random_seed))
if fice is False:
   print('Fitting the cumulative freezing spectrum...')
if fice is True:
   print('Fitting the fraction of ice...')
######################################################################################################################
if initial_guess is True:
   min_res = dual_annealing(objective_function, bounds=bnds,maxfun=100000000,seed=random_seed,x0=p0) 
if initial_guess is False:
   min_res = dual_annealing(objective_function, bounds=bnds,maxfun=100000000,seed=random_seed) 
print("\n################### Best set of parameters ")
best=min_res.x
MSE=round(min_res.fun,4)
print('The mean squared error (MSE) to the target data is '+str(MSE))   
######################################################################################################################
######################################################################################################################
######################################################################################################################
temp_fit = np.linspace(min(Temps_exp), max(Temps_exp), num=npoints)
def function_generate_artificial_Nm(param,factor_Nm_underlying):
   if nsubpop==3:
       if fice is False:
          beta=param[8]
       if fice is True:
          beta=1
       #Prob = (1-param[4]-param[7])*normalized_PDF(temp_fit,param[0],param[1])+param[4]*normalized_PDF(temp_fit,param[2],param[3])+param[7]*normalized_PDF(temp_fit,param[5],param[6])
       pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
       pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
       pdf3=normalized_PDF(temp_fit,param[5],param[6],disttype); ii = np.isnan(pdf3); pdf3[ii] = 0
       Prob=(1-param[4]-param[7])*pdf1+param[4]*pdf2+param[7]*pdf3
   if nsubpop==2:
       if fice is False:
          beta=param[5]
       if fice is True:
          beta=1
       #Prob = (1-param[4])*normalized_PDF(temp_fit,param[0],param[1])+param[4]*normalized_PDF(temp_fit,param[2],param[3]) 
       pdf1=normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(pdf1); pdf1[ii] = 0
       pdf2=normalized_PDF(temp_fit,param[2],param[3],disttype); ii = np.isnan(pdf2); pdf2[ii] = 0
       Prob=(1-param[4])*pdf1+param[4]*pdf2
   if nsubpop==1:
       if fice is False:
          beta=param[2]
       if fice is True:
          beta=1
       Prob = normalized_PDF(temp_fit,param[0],param[1],disttype); ii = np.isnan(Prob); Prob[ii] = 0
   integral_normalized=np.cumsum(Prob)*(temp_fit[1]-temp_fit[0])
   if fice is False:
      #g=beta*(1-integral_normalized)
      g=-np.log(1-beta*(1-integral_normalized))
      Nm_trial=g*max(data_exp)
      return temp_fit,Nm_trial
   if fice is True:
      f=beta*(1-integral_normalized/max(integral_normalized))
      return temp_fit,f
######################################################################################################################
fig3, ax = plt.subplots()
ax.plot(Temps_exp,data_exp,'o',color='black',linewidth=3,label='input data',fillstyle='none',markersize=10,mew=2)
######################################################################################################################
temp_fit,Nm_underlying=function_generate_artificial_Nm(best,max(data_fit_smooth))
ax.plot(temp_fit,Nm_underlying,'-',color='red',linewidth=4,ms=6,mew=2,fillstyle='full',zorder=20,label='optimized')
if fice is False:
   ax.set_yscale('log')
   ax.set_ylabel('N$_m$ ($g^{-1}$)')
if fice is True:
   ax.set_ylabel('Fraction of ice')

ax.set_xlabel('T ($^{o}C$)')
ax.legend(loc='best',fontsize=20,frameon=False)
fig3.set_size_inches(10.0, 10.0, forward=True)
fig3.tight_layout()
np.savetxt("Target_optimized.txt", np.column_stack((temp_fit,Nm_underlying)), fmt='%10.8f', newline='\n',header='Temperature, Nm or fice') 
fig3.savefig('Figure2_target_optimized.pdf')
######################################################################################################################
######################################################################################################################
######################################################################################################################
fig4, ax = plt.subplots()
temp_fit2 = np.linspace(min(Temps_exp)-2, max(Temps_exp)+2, num=npoints*2)
def function_generate_pdf(param):
   if nsubpop==3:
       Prob = (1-param[4]-param[7])*normalized_PDF(temp_fit2,param[0],param[1])+param[4]*normalized_PDF(temp_fit2,param[2],param[3])+param[7]*normalized_PDF(temp_fit2,param[5],param[6])
   if nsubpop==2:
       Prob = (1-param[4])*normalized_PDF(temp_fit2,param[0],param[1])+param[4]*normalized_PDF(temp_fit2,param[2],param[3]) 
   if nsubpop==1:
       Prob = normalized_PDF(temp_fit2,param[0],param[1]) 
   return Prob
######################################################################################################################
prob1=function_generate_pdf(best)
print("The parameters of the underlying distribution are ")
ax.plot(temp_fit2,prob1,'-',color='red',zorder=2,linewidth=4,fillstyle='none',ms=8,mew=2)
if nsubpop==1:
   ax.axvline(x=best[0],color='black',lw=2,linestyle='--',label='{:.1f}'.format(best[0]))
   print("Tmode #1="+str(round(best[0],2)))
   print("spread #1="+str(round(best[1],2)))
if nsubpop==2:
   ax.axvline(x=best[0],color='black',lw=2,linestyle='--',label='{:.1f}'.format(best[0]))
   ax.axvline(x=best[2],color='black',lw=2,linestyle='--',label='{:.1f}'.format(best[2]))
   print("Tmode #1="+str(round(best[0],2)))
   print("spread #1="+str(round(best[1],2)))
   print("weight #1="+str(round(1-best[4],8)))
   print("Tmode #2="+str(round(best[2],2)))
   print("spread #2="+str(round(best[3],2)))
   print("weight #2="+str(round(best[4],8)))
   if fice is False:
      print("beta="+str(round(best[5],2)))
   print("\nThe probability of freezing for the sub population at "+str(round(best[0],2))+" is "+str(round((1-best[4])*100,8))+"%.")
   print("The probability of freezing for the sub population at "+str(round(best[2],2))+" is "+str(round(best[4]*100,8))+"%.\n\n")
if nsubpop==3:
   ax.axvline(x=best[0],color='black',lw=2,linestyle='--',label='{:.1f}'.format(best[0]))
   ax.axvline(x=best[2],color='black',lw=2,linestyle='--',label='{:.1f}'.format(best[2]))
   ax.axvline(x=best[5],color='black',lw=2,linestyle='--',label='{:.1f}'.format(best[5]))
   print("Tmode #1="+str(round(best[0],2)))
   print("spread #1="+str(round(best[1],2)))
   print("weight #1="+str(round(1-best[4]-best[7],8)))
   print("Tmode #2="+str(round(best[2],2)))
   print("spread #2="+str(round(best[3],2)))
   print("weight #2="+str(round(best[4],8)))
   print("Tmode #3="+str(round(best[5],2)))
   print("spread #3="+str(round(best[6],2)))
   print("weight #3="+str(round(best[7],2)))
   if fice is False:
      print("beta="+str(round(best[8],8)))
   print("\nThe probability of freezing for the sub population at "+str(round(best[0],2))+" is "+str(round((1-best[4])*100,8))+"%.")
   print("The probability of freezing for the sub population at "+str(round(best[2],2))+" is "+str(round(best[4]*100,8))+"%.")
   print("The probability of freezing for the sub population at "+str(round(best[5],2))+" is "+str(round(best[7]*100,8))+"%.\n\n")
ax.set_yscale('log')
ax.set_ylabel('n$_m$')
ax.set_xlabel('T ($^{o}C$)')
ax.set_ylim([min(abs(best))/1e5,2])
ax.legend(loc='best',fontsize=25,ncol=1)
fig4.set_size_inches(10.0, 10.0, forward=True)
fig4.tight_layout()
fig4.savefig('Figure3_nm_differential_freezing_spectrum.pdf', dpi=300)
np.savetxt("nm_differential_freezing_spectrum.txt", np.column_stack((temp_fit2,prob1)), fmt='%10.8f', newline='\n',header='Temperature, nm') 
plt.show()