import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.mlab as mlab
import matplotlib.lines as mlines
import scipy.stats as stats

## Background Class H0 Definition

mean0=[0,1] 
sigmax0=1.5
sigmay0=1.5
rho=0.7
cov0=[[sigmax0**2,rho*sigmax0*sigmay0],[rho*sigmax0*sigmay0,sigmay0**2]]
n0=1000

## Signal Class H1 Definition

mean1=[4,4] 
sigmax1=1
sigmay1=1
cov1=[[sigmax1**2,0],[0,sigmay1**2]]
n1=1000

## Pseudo-Random generation for H0 e H1 gaussians with defined parameters

x0,y0=np.random.multivariate_normal(mean0,cov0,n0).T
T0=np.column_stack((x0,y0))*1.0 ## Data on two columns
m0=np.sum(T0,axis=0)/n0 ## Mean over rows

x1,y1=np.random.multivariate_normal(mean1,cov1,n1).T
T1=np.column_stack((x1,y1))*1.0 ## Data on two columns
m1=np.sum(T1,axis=0)/n1 ## Mean over rows

## Covariance matrices calculation

S_w=(n1)*np.cov(T1,rowvar=False,bias=1)+(n0)*np.cov(T0,rowvar=False,bias=1) ## S_w=S0+S1
S_inv_w=np.linalg.inv(S_w) ## S_w inversion with linear algebra numpy library
S_b=np.outer(m1-m0,m1-m0)

## Fisher vector definition

v=np.matmul(S_inv_w, m1-m0) ## ROWxCOL product by numpy package

## Data projection on \vec{v}

x11=np.matmul(v,T1.transpose()) ## Signal data projection
x00=np.matmul(v,T0.transpose()) ## Background data projection

## Plotting Data, \vec{v} and projections

plt.figure('Data and v')

plt.scatter(x1,y1,color='red', edgecolor='black')
plt.scatter(x0,y0, color='yellow', edgecolor='black')
plt.xlabel('x')
plt.ylabel('y')

def newline(p1, p2): ## Straight line from two given points...
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if(p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
        ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

    l = mlines.Line2D([xmin,xmax], [ymin,ymax], color='black', linestyle='-')
    ax.add_line(l)
    return l
newline([0,0], v)

plt.figure('Data projections') ## (bins choosen for graphic-only purposes)

plt.hist(x11, bins=33, color='red', edgecolor='black', alpha=0.8)
plt.hist(x00, bins=66, color='yellow', edgecolor='black', alpha=0.8)
plt.xlabel('Projected Data')
plt.ylabel('Counts')

## Showing plots

plt.show()

## Significance and Specificity optimization and consequent t_cut determination

t=np.concatenate((x11,x00),axis=0)  ## Definition of Fisher Statistic
q=np.linspace(min(t),max(t),len(t)) ## Homogeneuous domain for optimum-value search
alpha=np.zeros(len(t))              ## alpha ~ 1 - E_1 -> #{False Negatives}
beta=np.zeros(len(t))               ## beta  ~ E_0 -----> #{False Positives}

for i in range(0,len(q)):
 alpha[i]=(x11<q[i]).sum()/n1 ## Rejection of H1 when H1 is true
 beta[i]=(x00>q[i]).sum()/n0  ## Acceptance of H1 when H1 is false

ToMinimize=alpha+beta
best=min(ToMinimize)

t_cut=np.mean(q[ToMinimize==best]) ## There can be more than one best-match...
print (round(t_cut,6) , ': t_cut optimum value')

plt.axvline(t_cut) 
plt.xlabel('Projected Data')
plt.ylabel('Counts')

## Plotting alpha+beta as a function of trial t_cut (i.e. q)

plt.figure('Minimization Functional')
plt.plot(q,ToMinimize, linewidth = 5)
plt.xlabel('Trial Boundary [t]')
plt.ylabel('(False Positives + False Negatives)')
plt.axvline(np.mean(q[ToMinimize==best]))
plt.show()

## Identified signal purity computation

bestalpha=np.mean(alpha[np.where(((q-t)**2)**0.5==np.min((q-t)**2)**0.5)])
bestbeta=np.mean(beta[np.where(((q-t)**2)**0.5==np.min((q-t)**2)**0.5)])

print (round(bestalpha,6), ': alpha optimum value')
print (round(bestbeta,6), ': beta optimum value')
print (round(1-bestalpha,6)/(round(bestbeta,6)+round(1-bestalpha,6)), ': obtained signal purity')

## ROC curve determination

plt.figure('ROC curve')
plt.plot(beta,1-alpha, linewidth = 5, color='green')
plt.xlabel('False Positive Rate [Background Efficiency]')
plt.ylabel('True Positive Rate [Signal Efficiency]')
plt.title('ROC curve')
plt.show()

## Area under the ROC curve
    
SamplingPts = beta.shape[0]
    
# Integrating in [0.1] with "trapezoid-rule"
    
I = 0.5*alpha[0] + 0.5*alpha[beta.shape[0]-1]
for j in range(1,SamplingPts):
    I += alpha[j]
    dx = beta[j] - beta[j-1]
    I *= dx
I += 1
print('ROC area: ', I)