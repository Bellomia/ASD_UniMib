# Entering here the desired Breit-Wigner parameters:
x0 = 0
HWHM = 1

# Setting here an appropriate x domain (as a linear spaced vector)
x = seq(-10, 10, by = 0.04) # Note that "by" fixes the number of bins

# Defining the Breit-Wigner function
f_x = HWHM / (pi*HWHM^2 + pi*(x-x0)^2)

# Plotting and saving a graph of the function
cairo_pdf('Breit-Wigner.pdf', width = 16, height = 9, pointsize = 18)
mycolor = rgb(0.3, 0.7, 0)
plot(x, f_x, xlab='x', ylab='f(x)', col=mycolor)
dev.off()

# Extracting the required N uniform random numbers (N is a parameter of choice)
N = 10^5
uniform_set = runif(N, 0, 1) 

# Using inverse-cumulant equation to get Breit-Wigner dataset
BW_set = x0 + HWHM*tan(pi*(uniform_set-1/2))
	
# Defining histogram range and bins from x vector
hist_range = range(x)
bins = x

# Retrieving a proper subset of Breit-Wigner data 
BW_subset = subset(BW_set, BW_set <= max(hist_range) & BW_set >= min(hist_range))

# Plotting and saving the histogram of Breit-Wigner dataset (with analitical check)
myhistogram = hist(BW_subset, bins)
A = (myhistogram$counts / myhistogram$density)[1] # Normalization factor
cairo_pdf('BW_histogram.pdf', width = 16, height = 9, pointsize = 18)
plot(myhistogram)
curve(A*dcauchy(x,location=x0,scale=HWHM,log=FALSE), lwd=5, add=TRUE,  col=mycolor)
dev.off()

# Defining Minimal Standard parameters for LCG
a = 7^5
b = 0
c = 2^31-1

# Entering user's choice parameters (length of sequence, interval radius and seed)
n = 10^5
r = 10
s = 1
	
# Inizializing recursive vectors
x = numeric(length = n)
u = numeric(length = n)
x[1] = s

# Recursive cycle
t = 1
while(t < n)
{		u[t] = x[t] * r/c; 						#Normalizing the output
		x[t+1] = (a*x[t] + b) %% c;		#Recursive action (%% is R syntax for mod)
		t = t+1
};  u[t] = x[t] * r/c 
	
# Saving a histogram plot of {u}
cairo_pdf('LCG_histogram.pdf', width = 16, height = 9, pointsize = 22)
LCG_hist = hist(u, breaks = n/1000)
LCG_hist$name = 'Minimal Standard LCG Histogram'
LCG_hist$ylab = 'Frequency'
plot(LCG_hist,  ylab=LCG_hist$ylab, main=LCG_hist$name)
dev.off()

# Assuming to have a vector u[t], with t = 1,...,n and uniform distributed in [0,r]

	# Sample Mean
	sm = 0
	for (t in seq(1, n, by = 1))
		{ sm = sm + u[t] / n }
	
	# Sample Variance
	ss = 0
	for (t in seq(1, n, by = 1))
		{ ss = ss + (u[t] - sm)^2 / (n-1) }
		
	# Displaying results and relative deviation from 'true values'
	true_media = r/2;		true_varianza = r^2/12;
	media = sm;		varianza = ss
	e_media = abs(media - true_media)/true_media
	e_varianza = abs(varianza - true_varianza)/true_varianza
	media; varianza; e_media; e_varianza;

	# Subset of "even extrations"
	X = u[c(TRUE,FALSE)] 
	# Subset of "odd extractions"
	Y = u[c(FALSE,TRUE)]
	# Saving a scatter-plot of Y vs X
	cairo_pdf('CorrelationStudy.pdf')
	Xlab = 'Even Extractions'
	Ylab = 'Odd Extractions'
	Title = 'Correlation Study of LCG Sequence'
	plot(X , Y, pch = '.', xlab = Xlab, ylab = Ylab, main = Title)
	dev.off()
