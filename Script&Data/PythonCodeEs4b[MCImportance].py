import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import triang
from scipy.stats import beta

## Integrand function
def f(x): return x**5+x**4
 
x=np.linspace(-1,1,endpoint=True,dtype=float)
plt.plot(x,f(x),'.')

## 1st Interval: Triangular-approximation
s=np.linspace(-1,-0.25,50)
plt.plot(s,0.04*triang.pdf(s,0.2/0.7,a,b-a), '--')
a=-1
b=-0.25
n1=50
x1=np.random.triangular(a,-0.8,b,n1)  
eval_funct1=f(x1) 
y1=triang.pdf(x1,0.2/0.7,a,b-a)
n1=n1*1.0

i1=(1/n1)*np.sum(eval_funct1/y1)

mu_camp1=(np.sum(eval_funct1/y1))/n1
var_camp1=1/(n1-1)*np.sum((eval_funct1/y1-mu_camp1)**2)
var1=(1/n1)*var_camp1

## 2nd Interval: Beta-approximation
k=np.linspace(-0.25,1,100)
plt.plot(k,0.4*beta.pdf(k,5,1), '-.')

n2=9950
x2=beta.rvs(5,1,size=n2)  
eval_funct2=f(x2) 
y2=beta.pdf(x2,5,1)                     

n2=n2*1.0
i2=(1/n2)*np.sum(eval_funct2/y2)

mu_camp2=(np.sum(eval_funct2/y2))/n2
var_camp2=1/(n2-1)*np.sum((eval_funct2/y2-mu_camp2)**2)
var2=(1/n2)*var_camp2

## Results
print(i1+i2,'Integral Mean with Importance Sampling')
print(var1+var2,'Total Variance with Importance Sampling')
print(sqrt(var1+var2), 'Total Standard Deviation with Importance Sampling')
plt.xlabel('x')
plt.ylabel('y')
plt.show()