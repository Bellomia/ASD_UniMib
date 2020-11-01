import numpy as np
import matplotlib.pyplot as plt
import math

## Integrand function
def f(x): return x**5+x**4

## Integration interval
a =-1.0
b = 1.0

## Number of random number generations
n = 10000

## Standard MC implementation
h=[]
for k in range(1,1500):  

 x=np.random.uniform(a,b,n)  # [a,b]=[-1.0,1.0]
 eval_funct=f(x) 
 h.append((b-a)*np.sum(eval_funct)/(n))

S=(b-a)*(np.sum(eval_funct))/n
n=10000.0
mu_camp=(np.sum(eval_funct))/n
var_camp=1/(n-1)*np.sum((eval_funct-mu_camp)**2)
var=(b-a)**2*(1/n)*var_camp

print (S,'Integral Mean with Standard MC')
print (var,'Variance with Standard MC')
print(sqrt(var), 'Standard Deviation with Standard MC')

## Plotting a histogram of the generated gaussian
hist,bin_edges=np.histogram(h,bins=100)
plt.figure()
plt.hist(h,bin_edges)
plt.xlabel('Bins')
plt.ylabel('Counts')
plt.title('MC evaluation of Integral value [Standard Sampling]')
axes=plt.gca()

plt.show()