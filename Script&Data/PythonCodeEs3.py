from pylab import *
import numpy as np
import operator

## Search-Routine Parameters

n = 2 
mu_exp = 2.0
mu_rifl = 1.0
mu_contr_ex = 1.0/2
mu_contr_int = -1.0/2
mu_red = 1.0/2

## Random Trial Simplex [may need a careful range-definition]

def vertici_iniziali():

 import random 
 random.seed()
 vertici=[]
 
 for i in range(n+1):

  x=random.uniform(-12,-8);
  y=random.uniform(0,3);
  vertici.append([x,y])
 
 return vertici 

vertex=vertici_iniziali() 
x1=vertex[0]
x2=vertex[1]
x3=vertex[2]

## Target Functions [uncomment the desired one]

def f(x,y):
    
 #return -1.0*np.cos(x)*np.cos(y)*np.exp(-((x-np.pi)**2+(y-np.pi)**2))          # EASOM
 #return np.exp(0.5*(x**2+y**2-25)**2)+(np.sin(4*x-3*y))**4+0.5*(2*x+y-10)**2   # GOLDSTEIN-PRICE
 return 100*abs(y-0.01*x**2)+0.01*abs(x+10)                                    # BUKIN 6th
 #return  (x+2*y-7)**2+(2*x+y-5)**2                                             # BOOTH

data_x=[x1,x2,x3]

data_f=[f(x1[0],x1[1]),f(x2[0],x2[1]),f(x3[0],x3[1])]
data=[[f(x1[0],x1[1]),x1],[f(x2[0],x2[1]),x2],[f(x3[0],x3[1]),x3]]
print (data)

## Ordering [increasing f(x)]

data=sorted(data,key=operator.itemgetter(0)) 
data_f= [item[0] for item in data]
data_x= [item[1] for item in data]
print (data_f, 'f(trial simplex)')
print ( data_x ,'trial simplex')

## Plotting the Target Function and the Trial Simplex

xvec = np.linspace(-12, -8, 1000)
yvec = np.linspace(-0, 3, 1000)
X,Y = np.meshgrid(xvec, yvec)
Z = f(X, Y).T
fig, ax = subplots()
im = imshow(Z, cmap=cm.magma, vmin=Z.min(), vmax=Z.max(), extent=[-12, -8, 0, 3])
im.set_interpolation('bilinear')
cb = fig.colorbar(im)

Xvertex = np.array([])
Yvertex = np.array([])

for i in range(n+1):
    
    Xvertex = np.append(Xvertex, data_x[i][0])
    Yvertex = np.append(Yvertex, data_x[i][1])

plt.scatter(Xvertex,Yvertex, color='green', edgecolor='black', s=200)

coord = data_x
coord.append(coord[0]) # have to repeat the first point to create a 'closed loop'

xs, ys = zip(*coord) # creates lists of x and y values
plt.plot(xs,ys, color='white', alpha=0.3, ls='--') # Polytope draws up 

## Evolving the Simplex

epsilon = 10**(-5) # See...
loop=1

while data_f[n]- data_f[0] > epsilon: # ...this!

    loop=loop+1

    ## Centroid
    
    a=np.array(data_x[0:n])
    a=a/n
    xc=a.sum(axis=0)
    xr=(1+mu_rifl)*xc-mu_rifl*np.array(data_x[n])
    fr=f(xr[0],xr[1])
    print ( fr, 'fr')
    
    ## Reflection step
    
    if data_f[0]<=fr<data_f[n-1]:   
    
        data[n][1]=xr
        data[n][0]=f(xr[0],xr[1])
        data=sorted(data,key=operator.itemgetter(0)) 
        print('loop',loop,'riflessione')
        data_f= [item[0] for item in data]
        data_x= [item[1] for item in data]
        print (data_f, 'valori funzione ')
        print ( data_x ,'vertici ')
        for i in range(n+1):
            
            Xvertex = np.append(Xvertex, data_x[i][0])
            Yvertex = np.append(Yvertex, data_x[i][1])
        
        plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
        
        coord = data_x
        coord.append(coord[0]) # repeat the first point to create a 'closed loop'
        
        xs, ys = zip(*coord) # create lists of x and y values
        plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
        continue
    
    ## Expansion step
    
    if fr<data_f[0]: 
    
        a=np.array(data_x[0:n])
        a=a/n
        xc=a.sum(axis=0)
        xe=(1+mu_exp)*xc-mu_exp*np.array(data_x[n])
        fe=f(xe[0],xe[1])
        
        if fe<fr: 
        
            data[n][1]=xe  
            data[n][0]=f(xe[0],xe[1])
            data=sorted(data,key=operator.itemgetter(0)) 
            print('loop' ,loop,'espansione')
            
            data_f= [item[0] for item in data]
            data_x= [item[1] for item in data]
            print (data_f, 'valori funzione ')
            print ( data_x ,'vertici ')
            for i in range(n+1):
                
                Xvertex = np.append(Xvertex, data_x[i][0])
                Yvertex = np.append(Yvertex, data_x[i][1])
            
            plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
            
            coord = data_x
            coord.append(coord[0]) # repeat the first point to create a 'closed loop'
            
            xs, ys = zip(*coord) # create lists of x and y values
            plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
            continue 
        
        else: 
        
            data[n][1]=xr
            data[n][0]=f(xr[0],xr[1])
            data=sorted(data,key=operator.itemgetter(0))
            print('loop',loop,'riflessione')
            
            data_f= [item[0] for item in data]
            data_x= [item[1] for item in data]
            print (data_f, 'valori funzione ')
            print ( data_x ,'vertici ')
            for i in range(n+1):
                
                Xvertex = np.append(Xvertex, data_x[i][0])
                Yvertex = np.append(Yvertex, data_x[i][1])
            
            plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
            
            coord = data_x
            coord.append(coord[0]) # repeat the first point to create a 'closed loop'
            
            xs, ys = zip(*coord) # create lists of x and y values
            plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
            continue
    
    ## External-Contraction step
       
    if data_f[n-1]<=fr<data_f[n]:  
    
    
        a=np.array(data_x[0:n])
        a=a/n
        xc=a.sum(axis=0)
        xoc=(1+mu_contr_ex)*xc-mu_contr_ex*np.array(data_x[n])
        foc=f(xoc[0],xoc[1])
        
        if foc<fr:
            
            data[n][1]=xoc
            data[n][0]=f(xoc[0],xoc[1])
            data=sorted(data,key=operator.itemgetter(0)) 
            print( 'loop' ,loop,'contrazione esterna')
            
            data_f= [item[0] for item in data]
            data_x= [item[1] for item in data]
            print (data_f, 'valori funzione ')
            print ( data_x ,'vertici ')
            for i in range(n+1):
                
                Xvertex = np.append(Xvertex, data_x[i][0])
                Yvertex = np.append(Yvertex, data_x[i][1])
            
            plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
            
            coord = data_x
            coord.append(coord[0]) # repeat the first point to create a 'closed loop'
            
            xs, ys = zip(*coord) # create lists of x and y values
            plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
            continue
        
        ## Reduction step [!]
        
        else: 
        
            a=np.array(data_x)
            
            for i in range(1,n+1):
                data[i][1]=a[0]+mu_red*(a[i]-a[0])
                data[i][0]=f(data[i][1][0],data[i][1][1])
            
            
            
            data=sorted(data,key=operator.itemgetter(0)) 
            print('loop',loop,'riduzione')
            
            data_f= [item[0] for item in data]
            data_x= [item[1] for item in data]
            print (data_f, 'valori funzione ')
            print ( data_x ,'vertici ')
            for i in range(n+1):
                
                Xvertex = np.append(Xvertex, data_x[i][0])
                Yvertex = np.append(Yvertex, data_x[i][1])
            
            plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
            
            coord = data_x
            coord.append(coord[0]) # repeat the first point to create a 'closed loop'
            
            xs, ys = zip(*coord) # create lists of x and y values
            plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
            continue
    
    ## Internal-Contraction step
    
    if fr>=data_f[n]:  
    
        a=np.array(data_x[0:n])
        a=a/n
        xc=a.sum(axis=0)
        xic=(1+mu_contr_int)*xc-mu_contr_int*np.array(data_x[n])
        fic=f(xic[0],xic[1])
        
        if fic<data_f[n]: 
        
            data[n][1]=xic
            data[n][0]=f(xic[0],xic[1])
            data=sorted(data,key=operator.itemgetter(0)) 
            print('loop',loop ,'contrazione interna') 
            
            data_f= [item[0] for item in data]
            data_x= [item[1] for item in data]
            print (data_f, 'valori funzione ')
            print ( data_x ,'vertici ')
            for i in range(n+1):
                
                Xvertex = np.append(Xvertex, data_x[i][0])
                Yvertex = np.append(Yvertex, data_x[i][1])
            
            plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
            
            coord = data_x
            coord.append(coord[0]) # repeat the first point to create a 'closed loop'
            
            xs, ys = zip(*coord) # create lists of x and y values
            plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
            continue 
            
        ## Reduction step [!]    
            
        else: 
        
            a=np.array(data_x)
            for i in range(1,n+1):
                data[i][1]=a[0]+mu_red*(a[i]-a[0])
                data[i][0]=f(data[i][1][0],data[i][1][1])
                
            
            data=sorted(data,key=operator.itemgetter(0)) 
            print('loop',loop, 'riduzione')
            
            data_f= [item[0] for item in data]
            data_x= [item[1] for item in data]
            print (data_f, 'valori funzione ')
            print ( data_x ,'vertici ')
            for i in range(n+1):
                
                Xvertex = np.append(Xvertex, data_x[i][0])
                Yvertex = np.append(Yvertex, data_x[i][1])
            
            plt.scatter(Xvertex,Yvertex, color='white', edgecolor='black')
            
            coord = data_x
            coord.append(coord[0]) # repeat the first point to create a 'closed loop'
            
            xs, ys = zip(*coord) # create lists of x and y values
            plt.plot(xs,ys, color='white', alpha=0.3, ls='--') 
        
            continue
 
## Plotting the Final Simplex [a single point if the routine has converged!]
            
Xvertex = np.array([])
Yvertex = np.array([])

for i in range(n+1):
    
    Xvertex = np.append(Xvertex, data_x[i][0])
    Yvertex = np.append(Yvertex, data_x[i][1])

plt.scatter(Xvertex,Yvertex, color='red', edgecolor='black', s=100)

## Showing all the plots...

plt.show()
