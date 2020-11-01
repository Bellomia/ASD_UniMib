""" Standard Libraries """
import numpy as np
from pylab import *
from math import *

""" Pearson's Test Library """
from scipy.stats import chisquare 

## Fraunhofer Theory for Diffraction [N: #{apertures}]

def f(x):
    
        N = 5
        I0 = 1
        a = 1
        b = 0.5
        return I0*(np.sin(b*x)/(b*x))**2*np.sin(N*a*x)**2/np.sin(a*x)**2

## Smearing Parameter Choice [0 < c < 10]

c = 1
          
## First MC-simulation ['hit or miss'] to compute N_true

x1 = np.linspace(-3.7,-0.01,200)  # Have to decompose the domain to
x2 = np.linspace(0.01,3.7,200)    # eliminate the x = 0 singularity!
v = max(f(x1))

N_true = np.array([0])
while(sum(N_true==0)>=1): # Have to assure no bin results empty

 u = np.random.uniform(-3.7,3.7,500000)
 t = f(u)
 y = np.random.uniform(0,v,500000)
 boolean = y<t
 u = u[boolean]
 N_true,bin_edges1 = np.histogram(u,bins=100) # Pure Theoretical MC-Simulation

## Plotting 'Theory' and 'Theory Vs MC-Simulation'

half_bin = 0.5*(bin_edges1[len(bin_edges1)-1]-bin_edges1[len(bin_edges1)-2])
bin_values = bin_edges1[:-1]+half_bin
half_bin = 0.5*(bin_values[len(bin_values)-1]-bin_values[len(bin_values)-2])
bin_edges1 = np.append(bin_values-half_bin,bin_values[len(bin_values)-1]+half_bin)

plt.figure('Fraunhofer Theory')
plt.plot(x1,f(x1), color='red', linewidth=2)
plt.plot(x2,f(x2), color='red', linewidth=2)
plt.xlabel('ϑ-angle')
plt.ylabel('Intensity')

plt.figure('TheoryVsMontecarlo')
plt.plot(x1,f(x1)*(2*3.7/half_bin), color='red', linewidth=2,label='Theory')
plt.plot(x2,f(x2)*(2*3.7/half_bin), color='red', linewidth=2)
plt.bar(bin_edges1[:-1],N_true,width=half_bin*2,color='black',alpha=0.5,edgecolor='black',label='Pure MC')
plt.xlabel('ϑ-angle')
plt.ylabel('Counts')
plt.legend()

## Computing Smearing Matrix from Normal-Random-Generator [...from numpy library]

n = len(bin_values)
delta_x = half_bin*2
sigma = c*delta_x

M = np.zeros((n,n))
for j in range(0,n):
  if N_true[j] == 0:
   print('error')
  p = np.random.normal(bin_values[j],sigma,N_true[j])
  for k in range(0,len(p)):
   if p[k] < -3.7:                      #
    p[k] = p[k]+2*(bin_values[j]-p[k])  # Folding the tails in...
   if p[k] > 3.7:                       # -----------------------
    p[k] = p[k]-2*(p[k]-bin_values[j])  # 
  for i in range(0,n):
   M[i,j] = sum((p>bin_values[i]-half_bin) & (p<bin_values[i]+half_bin))/(N_true[j]*1.0)
   
N_smeared = np.matmul(M,N_true) # ROWxCOL product by numpy package

T = N_true      # For subsequent use
R = N_smeared   # ------------------

## Plotting Smearing Matrix [colormap]

fig, ax = subplots()
im = imshow(M, cmap=cm.tab20, vmin=M.min(), vmax=M.max(), extent=[0, 100, 100, 0])
cb = fig.colorbar(im)
ax.xaxis.tick_top()

## Second MC-simulation ['hit or miss'] to compute a new, independent, N_smeared [D != R]

N_true = np.array([0])
while(sum(N_true==0)>=1): # Have to assure no bin result empty

 u = np.random.uniform(-3.7,3.7,500000)
 t = f(u)
 y = np.random.uniform(0,v,500000)
 boolean = y<t
 u = u[boolean]
 N_true,bin_edges1 = np.histogram(u,bins=100) # Pure MC

n = len(bin_values)
delta_x = half_bin*2
sigma = c*delta_x

N_smeared = np.matmul(M,N_true) # Represents a 'Simulated Experiment' to unfold...

## Plotting 'True Signal Vs Detected Signal'

plt.figure('TrueVsSmeared')
plt.bar(bin_edges1[:-1],N_true,width=half_bin*2,color='black',alpha=0.5,edgecolor='black',label='Pure MC')
plt.bar(bin_edges1[:-1],N_smeared,width=half_bin*2,color='orange',alpha=0.6,edgecolor='black',label='Smeared')
plt.xlabel('ϑ-angle')
plt.ylabel('Counts')
plt.legend()

D = N_smeared # For subsequent use

## 'Bin-By-Bin Correction Factors' and Unfolding

C = T/R
U = C*D
sigma_U = C*D**0.5
dev = abs(U-T)

for i in range(U.shape[0]):         	#
    if dev[i] > 2.3548200*sigma_U[i]:   # Checking error bars...
        print('Unfolding FAILED!')  	#

N_unfold = U

## Chi-square test [Pearson] with respect to the last pure simulation

chi,p = chisquare(N_unfold,N_true)
print('chi: ', chi)
print('p-value: ', p)

## Plotting 'Theory Vs Experimental Vs Unfolded Signal'

plt.figure('Unfolded')
plt.plot(x1,f(x1)*(2*3.7/half_bin), color='red', linewidth=1.5,linestyle='--',label='Theory')
plt.plot(x2,f(x2)*(2*3.7/half_bin), color='red', linewidth=1.5, linestyle='--')
plt.bar(bin_edges1[:-1],N_smeared,width=half_bin*2,color='orange',alpha=0.6,edgecolor='black',label='Experimental')
plt.bar(bin_edges1[:-1],N_unfold,width=half_bin*2,color='green',alpha=0.5,edgecolor='black',label='Unfolded')
plt.xlabel('ϑ-angle')
plt.ylabel('Counts')
plt.legend()

## Showing all plots

plt.show()
