import numpy as np
import math

## Integrand function
def f(x): return x**5+x**4

## 1st Interval 
a=-1.0
b=-0.4
n1=100
x1=np.random.uniform(a,b,n1)  # [a,b]=[-1.0,-0.4]
k=np.linspace(0.35,1,50)
eval_funct1=f(x1)
n1=n1*1.0 
i1=(b-a)*(np.sum(eval_funct1)/n1) # mean
mu_camp1=(np.sum(eval_funct1))/n1
var_camp1=(1/(n1-1))*np.sum((eval_funct1-mu_camp1)**2)  # variance
var1=(b-a)**2*(1/n1)*var_camp1

## 2nd Interval 
a=-0.4
b=0.4
n2=100
x2=np.random.uniform(a,b,n2)  # [a,b]=[-0.4,0.4]
eval_funct2=f(x2) 
n2=n2*1.0
i2=(b-a)*(np.sum(eval_funct2)/n2) # mean
mu_camp2=(np.sum(eval_funct2))/n2
var_camp2=1/(n2-1)*np.sum((eval_funct2-mu_camp2)**2)  # variance
var2=(b-a)**2*(1/n2)*var_camp2

## 3rd Interval 
a=0.4
b=0.6
n3=800
x3=np.random.uniform(a,b,n3)  # [a,b]=[0.4,0.6]
eval_funct3=f(x3) 
n3=n3*1.0
i3=(b-a)*(np.sum(eval_funct3)/n3) # mean
mu_camp3=(np.sum(eval_funct3))/n3
var_camp3=1/(n3-1)*np.sum((eval_funct3-mu_camp3)**2)  # variance
var3=(b-a)**2*(1/n3)*var_camp3

## 4th Interval 
a=0.6
b=0.7
n4=1000
x4=np.random.uniform(a,b,n4)  # [a,b]=[0.6,0.7]
eval_funct4=f(x4) 
n4=n4*1.0
i4=(b-a)*(np.sum(eval_funct4)/n4) # mean
mu_camp4=(np.sum(eval_funct4))/n4
var_camp4=1/(n4-1)*np.sum((eval_funct4-mu_camp4)**2)  # variance
var4=(b-a)**2*(1/n4)*var_camp4

## 5th Interval 
a=0.7
b=0.8
n5=1500
x5=np.random.uniform(a,b,n5)  # [a,b]=[0.7,0.8]
eval_funct5=f(x5) 
n5=n5*1.0
i5=(b-a)*(np.sum(eval_funct5)/n5) # mean
mu_camp5=(np.sum(eval_funct5))/n5
var_camp5=1/(n5-1)*np.sum((eval_funct5-mu_camp5)**2)  # variance
var5=(b-a)**2*(1/n5)*var_camp5

## 6th Interval 
a=0.8
b=0.9
n6=3000
x6=np.random.uniform(a,b,n6)  # [a,b]=[0.8,0.9]
eval_funct6=f(x6) 
n6=n6*1.0
i6=(b-a)*(np.sum(eval_funct6)/n6) # mean
mu_camp6=(np.sum(eval_funct6))/n6
var_camp6=1/(n6-1)*np.sum((eval_funct6-mu_camp6)**2)  # variance
var6=(b-a)**2*(1/n6)*var_camp6

## 7th Interval 
a=0.9
b=1.0
n7=3500
x7=np.random.uniform(a,b,n7)  # [a,b]=[0.9,1]
eval_funct7=f(x7) 
n7=n7*1.0
i7=(b-a)*(np.sum(eval_funct7)/n7) # mean
mu_camp7=(np.sum(eval_funct7))/n7
var_camp7=1/(n7-1)*np.sum((eval_funct7-mu_camp7)**2)  # variance
var7=(b-a)**2*(1/n7)*var_camp7

## Results
print(n1+n2+n3+n4+n5+n6+n7,'Check on total number of generations [has to be 10^4]')

S=i1+i2+i3+i4+i5+i6+i7
print(S,'Total mean for S, with Stratified Sampling')

var=var1+var2+var3+var4+var5+var6+var7
print(var ,'Total variance for S, with Stratified Sampling')

dev = sqrt(var)
print(dev ,'Total standard deviation for S, with Stratified Sampling')