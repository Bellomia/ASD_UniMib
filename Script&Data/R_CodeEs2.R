# Entering the desired order for the harmonic sum
n = 50

# Direct: literal transposition of sum in for cycle
H = 0
for (k in seq(1,n, by = 1))
{ H = H+ 1/k	} 
print(H, digits = 16)

# Trick: invert sum order (avoids round-offs!)
H = 0
k = n
while (k > 0)
{ H = H + 1/k; k = k - 1 }
print(H, digits = 16)

# Setting Eulero-Mascheroni constant (correct up to 16 digits)
EM = 0.5772156649015329

# Entering the desired truncation-parameter for the series
n = 10^5

# First attempt: direct algorithm from definition
em = -log(n)
for (k in seq(1,n, by = 1))
{ em = em + 1/k	} 
print(EM, digits = 16)
print(em, digits = 16)
print(abs(EM-em), digits = 1)

# Second attempt: inverted sum order
em = 0
k = n
while (k > 0)
{ em = em + 1/k; k = k - 1 }
em = em-log(n)
print(EM, digits = 16)
print(em, digits = 16)
print(abs(EM-em), digits = 1)

# Third attempt: Macys formula
em = 0
k = n^2
while (k > n)
	{ em = em - 1/k; k = k - 1 }
while (k > 0)
{ em = em + 1/k; k = k -1 }
print(EM, digits = 16)
print(em, digits = 16)
print(abs(EM-em), digits = 1)
	
# Final refinement: some magic
em = 0
k = n^2+n # 1st Trick
while (k > n)
{ em = em - 1/k; k = k - 1 }
while (k > 0)
{ em = em + 1/k; k = k - 1 }
em = em + 1/(6*n^2) - 1/(6*n^3)# 2nd Trick
print(EM, digits = 16)
print(em, digits = 16)
print(abs(EM-em), digits = 1)
