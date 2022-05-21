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
fname = 'lab4pt4_2c.csv'
# the data file fname has x data (first column) and y data (second column) and 
# the uncertianty in y (third column).
# Enter the name and units of each column here so that the plots are properly
# labeled. e.g. x_name = 'Time', x_units = 'ms', y_name = 'Voltage', y_units
# = 'V'.
x_name = 'Frequency'
x_units = 'Hz'
y_name = 'Voltage Gain'
y_units = 'no units'
normresidunits='no units'
# it is often necessary to provide good initial guesses for the parameters in order 
# for the chi2 fit to be converge correctly. A Lorentizain  funtion may have  3 
#parameters: # amplitude, width and  and  a y_offset. Initial guesses for these
# parameters are housed in guesses. For example, guesses = (2., 0.004. 0.01)
# corresponds to an lorenztian with amplitude gain 100.0  at  zero 
#frequency,  the product  fcA0=1e6 s^-1 and offset 0. If the guesses are 
#not good  the fitting routine  may not converge on the global minimum in Chi^2.   
guesses = (150,1.5e6,0)

###############################################################################
# loads data, plots guessed curve
###############################################################################

# definition of the fit function
def fit_function(x, Gain,fcA0, y_offset):
    return Gain/np.sqrt(1+(x*Gain/(fcA0))**2) + y_offset
#    return amplitude/(1+(2*np.pi*x*tau)**2) + y_offset
#load the file "fname", defined above
data = np.loadtxt(fname, delimiter=',', comments='#')
# access the data columns and assign variables x,y  and  uncetainty in y y_sigma
x = data[:,0]
y = data[:,3]
y_sigma = data[:,4]
#define an array of frequency points for plotting the theory function
xtheory = np.linspace(min(x),max(x),500)
    
#compares the guessed curve to the data for visual reference
y_guess = fit_function(xtheory,*guesses)
plt.errorbar(x,y,yerr=y_sigma,marker='.',linestyle='',label="measured data")
plt.plot(xtheory,y_guess,marker="",linestyle="-",linewidth=1,color="g",
         label="initial guess")
plt.xlabel('{} [{}]'.format(x_name,x_units))
plt.ylabel('{} [{}]'.format(y_name,y_units))
plt.title('Comparison between the data and the intial guess')
plt.legend(loc='best',numpoints=1)
print ('Displaying plot 1')
plt.show()


    
    
# fit the data to fit_function. fit_params is the resultant parameters, fit_cov
# is the covariance matrix between all the parameters. maxfev=10**5 means the
# fitting procedure is iterated at maximum 10^5 times before giving up.
fit_params,fit_cov = curve_fit(fit_function,x,y,sigma=y_sigma,p0=guesses,
                               maxfev=10**5)
                               
###############################################################################
# prints the chi2 
###############################################################################

# function that  calculates the chi square value of a fit
def chi_square (fit_parameters, x, y, sigma):
#
    return np.sum((y-fit_function(x, *fit_parameters))**2/sigma**2)
    
# calculate and print chi square as well as the per degree-of-freedom value
chi2 = chi_square(fit_params,x,y,y_sigma)
dof = len(x) - len(fit_params)
print ("Goodness of fit - chi square measure:")
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

# read out parameter results
param_names = ['Gain','FcA0','V_offset']
print ("Fit parameters:")
for i in range(len(fit_params)):
    print ('{} = {:.3e} +/- {:.3e}'.format(param_names[i],
                                          fit_params[i],
                                          fit_params_error[i]))


# prints out the covariance between all variables. Uncomment this code to see
# the results of this more sophisticated statistical measure.

                  
print( "\nCovariance between fit parameters:")
for i,fit_covariance in enumerate(fit_cov):
    for j in range(i+1,len(fit_covariance)) :
        print( "{} and {} : {:.3e}".format(param_names[i],
                                          param_names[j],fit_cov[i,j]))                  
                                          
print ("\nCorrelation between fit parameters:")
for i,fit_covariance in enumerate(fit_cov):
    for j in range(i+1,len(fit_covariance)) :
        print( "{} and {} : {:.3e}".format(param_names[i],
    param_names[j],fit_cov[i,j]/(np.abs(fit_params_error[i]*fit_params_error[j]))))

                                                                        
###############################################################################
# plots the data and the fit curve
###############################################################################
    
# y_fit is the curve using the parameters calculated from the fit
y_fitfunc = fit_function(xtheory,*fit_params) 
#y_fit = fit_function(x,*fit_params)
#y_curve=fit_function(xtheory,*fit_params)

# plots y, and y_fit on the same plot.
# marker='.' : data points are not indicated by markers
# linestyle= '-' : a continuous line is drawn
# linewidth=2 : the line thickness is set to 2
# color='r' : the color of the line is set to red
# label=string : the string is shown in the legend
plt.errorbar(x,y,yerr=y_sigma,marker='.',linestyle='',label="measured data")
plt.plot(xtheory,y_fitfunc,marker="",linestyle="-",linewidth=2,color="r",
         label=" fit")
# add axis labels and title
plt.xlabel('{} [{}]'.format(x_name,x_units))
plt.ylabel('{} [{}]'.format(y_name,y_units))
plt.title(r'Voltage Gain vs. Frequency (G = 150)')
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
# plots residual and histogram of residual. Don't touch this par tof the code
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
ax1.legend(loc='best',numpoints=1)
ax2 = fig.add_subplot(212)
ax2.bar(bins[:-1],hist,width=bins[1]-bins[0])
ax2.set_ylim(0,1.2*np.max(hist))
ax2.set_xlabel('(y-y_fit)/y_sigma ')
ax2.set_ylabel('Number of occurences')

#if y_sigma != None:
# within_err is the percentage of data points whose error bars overlap with 
# the fit. If chi2/dof = 1, within_error ~ 66% (one standard deviation)
within_err=100.*np.sum((residual<=y_sigma)&(residual>=-y_sigma))/len(residual)
print ("\nResidual information:")
print ('{:.1f}% of data points agree within 1 sigma  of fit'.format(within_err))
ax2.text(0.0,np.max(hist)*1.1,'{:.1f}% of data point within 1 sigma'.format(within_err),
         horizontalalignment='center',verticalalignment='center')
ax2.vlines(-1.0,0,np.max(hist)*1.3,lw=2,color='r',linestyles='--')
ax2.vlines(+1.0,0,np.max(hist)*1.3,lw=2,color='r',linestyles='--',label='+/- one sigma')
#ax2.text(np.min(bins),np.max(hist)*1.1,'{:.1f}% of data'.format(100-within_err),
#         horizontalalignment='left',verticalalignment='center')
#ax2.text(np.max(bins),np.max(hist)*1.1,'{:.1f}% of data'.format(100-within_err),
#         horizontalalignment='right',verticalalignment='center')
#ax2.legend(loc=(0.35,0.05))

print ('\nDisplaying plot 3')
plt.show()
















