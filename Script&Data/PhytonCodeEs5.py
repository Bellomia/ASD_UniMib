""" Standard Libraries """
import numpy as np
import pylab as plt
from math import *

""" I/O Libraries """
import os
import sys

""" Gaussian-Smoothing Library """
from scipy.ndimage import gaussian_filter as smoothing

""" Best-Fit Library """
from scipy.optimize import curve_fit

## Choosing the source

SourceName = str(input('Please insert the desired source [Cs137|Co57|Na22|Ba133|Eu152]: '))

## Opening the ASCII file related to the chosen source

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
    
path = get_script_path()
data=np.loadtxt(path+'/Script&Data/CeBr3_'+SourceName+'.dat')

## Defining Effective-Range

i = len(data)
while(data[i-1] <= 0.01*np.max(data)/2): 
    i = i - 1
    
x_max = i
data = data[range(1,x_max)]

## Plotting retrieved experimental data

x = np.array(range(1,x_max))
y = data

plt.figure('Raw'+SourceName)
plt.plot(x,y,'ko',markersize=0.5)
plt.xlim((1, x_max))
plt.xlabel('Channels')
bottom, top = plt.ylim()
plt.ylabel('Counts')
plt.title('RAW Data ['+SourceName+']')
plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0),useMathText=True)
plt.show()

## Peak-Recognition Routine [RAW]

# Derivatives of data [RAW]
y1 = np.gradient(y)
y2 = np.gradient(y1)

X_Peaks1 = np.array([])
Y_Peaks1 = np.array([])
X_Peaks2 = np.array([])
Y_Peaks2 = np.array([])
Y1 = np.array([])
Y2 = np.array([])

steps = 20
for i in range(1,steps+1):
    
    # Minimization only in a partial domain
    X = x[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y = y[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y1 = y1[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y2 = y2[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    X_CandidatePeak = X[np.argmax(Y)]
    Y_CandidatePeak = np.max(Y)
    Y1_CandidatePeak = Y1[np.argmax(Y)]
    Y2_CandidatePeak = Y2[np.argmax(Y)]
    
    a = 1*np.max(Y)/2           # Threshold on Background
    b = 1*np.max(abs(Y1))/2     # Threshold on Slope
    c = 0.1*np.max(abs(Y2))/2   # Threshold on Concavity
    
    if (Y_CandidatePeak > a)and(abs(Y1_CandidatePeak)< b)and(Y2_CandidatePeak<-c)and(Y_CandidatePeak > 0.13*np.max(y)):
        
        X_Peaks1 = np.append(X_Peaks1, X_CandidatePeak)
        Y_Peaks1 = np.append(Y_Peaks1, Y_CandidatePeak)
        
steps = 8
for i in range(1,steps+1):
    
    # Minimization only in a partial domain
    X = x[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y = y[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y1 = y1[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y2 = y2[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    X_CandidatePeak = X[np.argmax(Y)]
    Y_CandidatePeak = np.max(Y)
    Y1_CandidatePeak = Y1[np.argmax(Y)]
    Y2_CandidatePeak = Y2[np.argmax(Y)]
    
    a = 1*np.max(Y)/2           # Threshold on Background
    b = 1*np.max(abs(Y1))/2     # Threshold on Slope
    c = 0.1*np.max(abs(Y2))/2   # Threshold on Concavity
    
    if (Y_CandidatePeak > a)and(abs(Y1_CandidatePeak)< b)and(Y2_CandidatePeak<-c)and(Y_CandidatePeak > 0.13*np.max(y)):
        
        for j in range(len(X_Peaks1)):
            if X_CandidatePeak == X_Peaks1[j]:
                X_Peaks2 = np.append(X_Peaks2, X_CandidatePeak)
                Y_Peaks2 = np.append(Y_Peaks2, Y_CandidatePeak)
        
""" Showing determined peaks """
plt.figure('Raw'+SourceName+' with determined peaks')
plt.plot(x,y,'ko',markersize=0.5)
plt.xlim((1, x_max))
plt.xlabel('Channels')
bottom, top = plt.ylim()
plt.ylabel('Counts')
plt.title('RAW Data with determined peaks ['+SourceName+']')
plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0),useMathText=True)
plt.show()
plt.plot(X_Peaks2,Y_Peaks2,'X',c='orange',markeredgecolor='k',markersize=10,alpha=0.75)

## Smoothing Data

y=smoothing(data,sigma=5)
plt.figure('Smoothing'+SourceName)
plt.plot(x,data,'ko',markersize=0.5,label='RAW Data   ')
plt.plot(x,y,'ro',markersize=0.5,label='Smoothed   ')
plt.xlim((1, x_max))
plt.xlabel('Channels')
plt.ylim((bottom, top))
plt.ylabel('Counts')
plt.title('Smoothing ['+SourceName+']')
plt.legend(markerscale=10)
plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0),useMathText=True)
plt.figure('Smoothed'+SourceName)
plt.plot(x,y,'ro',markersize=0.5)
plt.xlim((1, x_max))
plt.xlabel('Channels')
plt.ylim((bottom, top))
plt.ylabel('Counts')
plt.title('Smoothed Data with determined peaks ['+SourceName+']')
plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0),useMathText=True)
plt.show()

## Peak-Recognition Routine [Smoothed]

# Derivatives of data [Smoothed]
y1 = np.gradient(y)
y2 = np.gradient(y1)

X_Peaks1 = np.array([])
Y_Peaks1 = np.array([])
X_Peaks2 = np.array([])
Y_Peaks2 = np.array([])
Y1 = np.array([])
Y2 = np.array([])

steps = 20
for i in range(1,steps+1):
    
    # Minimization only in a partial domain
    X = x[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y = y[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y1 = y1[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y2 = y2[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    X_CandidatePeak = X[np.argmax(Y)]
    Y_CandidatePeak = np.max(Y)
    Y1_CandidatePeak = Y1[np.argmax(Y)]
    Y2_CandidatePeak = Y2[np.argmax(Y)]
    
    a = 1*np.max(Y)/2           # Threshold on Background
    b = 1*np.max(abs(Y1))/2     # Threshold on Slope
    c = 0.1*np.max(abs(Y2))/2   # Threshold on Concavity
    
    if (Y_CandidatePeak > a)and(abs(Y1_CandidatePeak)< b)and(Y2_CandidatePeak<-c)and(Y_CandidatePeak > 0.13*np.max(y)):
        
        X_Peaks1 = np.append(X_Peaks1, X_CandidatePeak)
        Y_Peaks1 = np.append(Y_Peaks1, Y_CandidatePeak)
        
steps = 8
for i in range(1,steps+1):
    
    # Minimization only in a partial domain
    X = x[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y = y[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y1 = y1[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    Y2 = y2[(x > (x_max/steps)*(i-1)) & (x <= (x_max/steps)*i)]
    X_CandidatePeak = X[np.argmax(Y)]
    Y_CandidatePeak = np.max(Y)
    Y1_CandidatePeak = Y1[np.argmax(Y)]
    Y2_CandidatePeak = Y2[np.argmax(Y)]
    
    a = 1*np.max(Y)/2           # Threshold on Background
    b = 1*np.max(abs(Y1))/2     # Threshold on Slope
    c = 0.1*np.max(abs(Y2))/2   # Threshold on Concavity
    
    if (Y_CandidatePeak > a)and(abs(Y1_CandidatePeak)< b)and(Y2_CandidatePeak<-c)and(Y_CandidatePeak > 0.13*np.max(y)):
        
        for j in range(len(X_Peaks1)):
            if X_CandidatePeak == X_Peaks1[j]:
                X_Peaks2 = np.append(X_Peaks2, X_CandidatePeak)
                Y_Peaks2 = np.append(Y_Peaks2, Y_CandidatePeak)
        
""" Showing determined peaks """
plt.plot(X_Peaks2,Y_Peaks2,'X',c='orange',markeredgecolor='k',markersize=10,alpha=0.75)

## Principal-Peak [P.P] Analysis

xmax_peak = X_Peaks2[np.argmax(Y_Peaks2)]   # Principal-Peak position
ymax_peak = np.max(Y_Peaks2)                # Principal-Peak value
index = np.where(x==xmax_peak)[0]           # Finding P.P in our data array

# Definition of a ROI: we must have |yROI-yPeak| <= d and |yROI'| < t, on both sides of P.P.
############################################################################################
d = ymax_peak*0.75 # ----------------------
t = 50             # BEST VALUES [HEURISTIC]
#################### ----------------------
x1=x[(abs(y1)<t) & (x<xmax_peak) & (abs(y-ymax_peak)>d)]
index1=np.where(x==x1[abs(x1-xmax_peak)==min(abs(x1-xmax_peak))])[0]
i1=np.array(index1[0])
x2=x[(abs(y1)<t) & (x>xmax_peak) & (abs(y-ymax_peak)>d)]
index2=np.where(x==x2[abs(x2-xmax_peak)==min(abs(x2-xmax_peak))])[0]
i2=np.array(index2[0])
xnew=x[i1:i2]
ynew=y[i1:i2]

#Plotting the ROI-confined Peak
plt.figure('Peak-Fitting_'+SourceName)
plt.plot(x[range(i1-1,i2+1)],y[range(i1-1,i2+1)],'ko',markersize=1.5,label='ROI-data')

# Background Esimate: Linear-Fit beetwen ROI edges... 
a=xnew[0]
b=xnew[len(xnew)-1]
x_to_fit = np.array([xnew[0],xnew[len(xnew)-1]])
y_to_fit = np.array([ynew[0],ynew[len(xnew)-1]])
coefficients = np.polyfit(x_to_fit, y_to_fit, deg=1) # Linear-Fit
polynomial = np.poly1d(coefficients)                 # ----------
x_pol = np.linspace(xnew[0],xnew[len(xnew)-1],100)
y_pol = polynomial(x_pol)
plt.plot(x_pol,y_pol,'m-.',linewidth=1.5,alpha=0.5,label='Background')

# Background-correcting the Peak [and plotting]
ynew=ynew-polynomial(xnew)
plt.plot(xnew,ynew,'o',c='orange',markersize=1.5,label='Corrected-Peak')
plt.ylim(0,ymax_peak*2-1.25*d)

## Gaussian Best-Fit of the corrected P.P. data

def Gauss(x,A,x0,sigma):
    return A*np.exp(-(x-x0)**2/(2*sigma**2))

A_init = np.max(ynew)       # -----------
x0_init = np.mean(xnew)     # First guess
sigma_init = np.std(xnew)   # -----------

# Least-Squares Fitting...
popt,pcov = curve_fit(Gauss,xnew,ynew,p0=[A_init,x0_init,sigma_init])

A = popt[0]                                  # ----------------
x0 = popt[1]                                 # Best-Fit results
sigma = popt[2]                              # ----------------
StDev_vec = 2.3548200*np.sqrt(np.diag(pcov)) # ----------------

W = 2.3548200*popt[2] # Peak-Width ~ FWHM
S = 0.5*A*W # Area with Triangle-Approx.
R = W/x0 # Resolution ~ FWHM/x0

dA = StDev_vec[0]
dx = StDev_vec[1]
dW = 2.3548200*StDev_vec[2]
dS = (dA/A + dW/W)*S
dR = (dW/W + dx/x0)*R

# Printing Results...
print('PRINCIPAL-PEAK POSITION', x0,'±', dx)
print('PRINCIPAL-PEAK HEIGHT: ', A,'±', dA)
print('PRINCIPAL-PEAK WIDTH: ', W,'±', dW)
print('PRINCIPAL-PEAK AREA :', S,'±', dS)
print('CRYSTAL RESOLUTION: ', R,'±', dR)

# Plotting & Showing All
yFIT = Gauss(xnew,A,x0,sigma)
plt.plot(xnew,yFIT,'g--',linewidth=1.5,label='Gaussian-Fit')
plt.title('Principal Peak Best-Fitting ['+SourceName+']')
plt.ticklabel_format(axis='both',style='sci',scilimits=(0,0),useMathText=True)
plt.legend(markerscale=3)
plt.show()