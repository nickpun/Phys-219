## Based  on script by Evan  Kiefl 2016, Edits   by Rob Kiefl and Evan Thomas
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#LIST OF ALL INPUTS

# fname is assumed to be a .csv file (comma separated values). All rows 
# containing non-number text (example: names for each column) must begin 
# with "#". By default, the oscilloscope generates text column names 
# without "#"--it must be added manually. The .csv file must be in the same 
# folder as this fit program, otherwise the full file extension must be added
# to fname: e.g. fname = 'folder/subfolder/subsubfolder/file.csv'
fname = 'lab24p4.csv'
# the data file fname has x data (first column) and y data (second column).
# Enter the name and units of each column here so that the plots are properly
# labeled. e.g. x_name = 'Time', x_units = 'ms', y_name = 'Voltage', y_units
# = 'V'.
x_name = 'Time'
x_units = 's'
y_name = 'Voltage'
y_units = 'V'
# it is necessary to provide good initial guesses for the parameters in order 
# for the chi2 fit to be successful. A function describing the signal could have  4 parameters:
# amplitude, relaxation  time, frequency,phase  and offset. Initial guesses for these
# parameters are housed in guesses. For example, guesses = (1., 5e-5,,4e+4.  0.,0.)
# corresponds to an expoentially damped cosine   with amplitude 1.0 V at  the start 
#time defined as t0,   relaxation time of 5e-4 s  or 50 microsec , frequency of 40kHz
#and phase shift of 0.  Note if the guesses are not very good  the fitting routine 
#may not converge on the minimum in Chi^2.   
guesses = (0.05635,2.209e-4,35440,2.074)
# There is always uncertainty in measurements. For oscilloscope data, it is
# reasonable to assume each data point has the same absolute uncertainty (e.g. 
# +- 0.02V). Define y_sigma to be this uncertainty. If you would like to 
# perform the fit without uncertainties, set y_sigma = None.
y_sigma = 0.000172
#t0 is the start time in the data when the applied voltage suddenly goes to 0.
t0=0.0

###############################################################################
# loads data, plots guessed curve, then performs the fit
###############################################################################

# definition of the fit function
def fit_function(x, amplitude, tau, freq, phase):
#    """ sine function """
    return amplitude * np.exp(-(x-t0)/tau )*np.cos(2*np.pi*freq*(x-t0)+phase) + 0.0008

#load the file "fname", defined above
data = np.loadtxt(fname, delimiter=',', comments='#')
# access the data columns and assign variables xraw,yraw
xraw = data[5000:,0]
yraw = data[5000:,1]
#define packing factor
npac=50
#define a function  to pack the data
def pack(A,p):
  # A is an array, and p is the packing factor
  B = np.zeros(len(A)//p)
  i = 1
  while i-1<len(B):
    B[i-1] = np.mean(A[p*(i-1):p*i])
    i += 1
  return B
# pack the data
x=pack(xraw,npac)
y=pack(yraw,npac)    
 #define an array of frequency points for plotting the theory function
xtheory = np.linspace(min(x),max(x),500)       
#compares the guessed curve to the data for visual reference
y_guess = fit_function(xtheory,*guesses)
plt.errorbar(x,y,yerr=y_sigma,marker='.',linestyle='',label="measured data")
plt.plot(xtheory,y_guess,marker="",linestyle="-",linewidth=1,color="g",
         label="initial guess")
plt.xlabel('{} [{}]'.format(x_name,x_units))
plt.ylabel('{} [{}]'.format(y_name,y_units))
plt.title(r'Comparison between the data and the intial guess')
plt.legend(loc='best',numpoints=1)
print ('\nDisplaying plot 1')
plt.show()


    
    
# fit the data to fit_function. fit_params is the resultant parameters, fit_cov
# is the covariance matrix between all the parameters. maxfev=10**5 means the
# fitting procedure is iterated at maximum 10^5 times before giving up.
fit_params,fit_cov = curve_fit(fit_function,x,y,sigma=y_sigma*np.ones(len(x)),p0=guesses,
                               maxfev=10**5)
                               
###############################################################################
# prints the chi2 
###############################################################################

# function that  calculates the chi square value of a fit
def chi_square (fit_parameters, x, y, sigma):
    if sigma is None:
        sigma = 1.
    return np.sum((y-fit_function(x, *fit_parameters))**2/sigma**2)
    
# calculate and print chi square as well as the per degree-of-freedom value
chi2 = chi_square(fit_params,x,y,y_sigma)
dof = len(x) - len(fit_params)
print ("\nGoodness of fit - chi square measure:")
print ("Chi2 = {}, Chi2/dof = {}\n".format(chi2, chi2/dof))

###############################################################################
# prints the fit parameters (with uncertainty)
###############################################################################

# the covariance matrix is rescaled to cancel the inverse scaling that is 
# performed for numerical reasons during execution of curve_fit -- do not 
# change this line!
fit_cov = fit_cov*dof/chi2
# calculate the standard deviation as uncertainty of the fit parameters
fit_params_error = np.sqrt(np.diag(fit_cov))

# read out parameter results,
param_names = ['amplitude','relaxation time','frequency','phase']
print ("Fit parameters:")
for i in range(len(fit_params)):
    print ('{} = {:.3e} +/- {:.3e}'.format(param_names[i],
                                          fit_params[i],
                                          fit_params_error[i]))


# prints out the covariance between all variables. Uncomment this code to see
# the results of this more sophisticated statistical measure.
#print "\nCovariance between fit parameters:"
#for i,fit_covariance in enumerate(fit_cov):
#    for j in xrange(i+1,len(fit_covariance)) :
#        print "{} and {} : {:.3e}".format(param_names[i],
#                                          param_names[j],
#                                          fit_cov[i,j])
                                                                        
###############################################################################
# plots the data and the fitted curve
###############################################################################
    
# y_fit is the curve using the parameters calculated from the fit
x_fitfunc = np.linspace(min(x),max(x),len(x))
y_fitfunc = fit_function(x_fitfunc,*fit_params)    

# plots y, y_fit, and y_guess on the same plot.
# marker='.' : data points are not indicated by markers
# linestyle= '-' : a continuous line is drawn
# linewidth=2 : the line thickness is set to 2
# color='r' : the color of the line is set to red
# label=string : the string is shown in the legend
plt.errorbar(x,y,yerr=y_sigma,marker='.',linestyle='',label="measured data")
plt.plot(x_fitfunc,y_fitfunc,marker="",linestyle="-",linewidth=2,color="r",
         label=" fit")
# add axis labels and title
plt.xlabel('{} [{}]'.format(x_name,x_units))
plt.ylabel('{} [{}]'.format(y_name,y_units))
plt.title(r'Voltage vs. Time w/o first 5000 pts')
# set the x and y boundaries of your plot
#plt.xlim(lower_x,upper_x)
#plt.ylim(lower_y,upper_y)
# show a legend. loc='best' places legend where least amount of data is 
# obstructed. 
plt.legend(loc='best',numpoints=1)
print ('\nDisplaying plot 2')
# plt.show() may or may not need to be commented out depending on your python
# editor (spyder) settings.
plt.show()


###############################################################################
# plots residual and histogram of residual. Don't touch this part of the code
###############################################################################

# residual is the difference between the data and theory
y_fit=fit_function(x,*fit_params)
residual = y-y_fit
normresidual=residual/y_sigma
# creates a histogram of the residuals
hist,bins = np.histogram(normresidual,bins=30)
# this complicated code produces a figure with a plot of the residuals as well
# as a histogram of the residuals. You do not need to understand this code.
fig = plt.figure(figsize=(7,10))
ax1 = fig.add_subplot(211)
ax1.errorbar(x,residual,yerr=y_sigma,marker='.',linestyle='',
             label="residual (y-y_fit)")
ax1.hlines(0,np.min(x),np.max(x),lw=2,alpha=0.8)
ax1.set_xlabel('{} [{}]'.format(x_name,x_units))
ax1.set_ylabel('y-y_fit [{}]'.format(y_units))
plt.title('Residuals')
ax1.legend(loc='best',numpoints=1)
ax2 = fig.add_subplot(212)
ax2.bar(bins[:-1],hist,width=bins[1]-bins[0])
ax2.set_ylim(0,1.2*np.max(hist))
ax2.set_xlabel('(y-y_fit)/y_sigma')
ax2.set_ylabel('Number of occurences')
plt.title('Histogram  of normalized residuals')

if y_sigma != None:
    # within_err is the percentage of data points whose error bars overlap with 
    # the fit. If chi2/dof = 1, within_error ~ 66% (one standard deviation)
    within_err=100.*np.sum((residual<=y_sigma)&(residual>=-y_sigma))/len(residual)
    print ("\nResidual information:")
    print ('{:.1f}% of data points agree with fit to within 1 sigma'.format(within_err))
    ax2.vlines(-1.0,0,np.max(hist)*1.3,lw=2,color='r',linestyles='--')
    ax2.vlines(+1.0,0,np.max(hist)*1.3,lw=2,color='r',linestyles='--',label='+-one sigma')
    ax2.text(0.0,np.max(hist)*1.1,'{:.1f}% of data points within one sigma of best fit'.format(within_err),
             horizontalalignment='center',verticalalignment='center')
#    ax2.text(np.min(bins),np.max(hist)*1.1,'{:.1f}% of data'.format(100-within_err),
#             horizontalalignment='left',verticalalignment='center')
#    ax2.text(np.max(bins),np.max(hist)*1.1,'{:.1f}% of data'.format(100-within_err),
#             horizontalalignment='right',verticalalignment='center')
    ax2.legend(loc=(0.35,0.05))
    
print ('\nDisplaying plot 3')
plt.show()
















