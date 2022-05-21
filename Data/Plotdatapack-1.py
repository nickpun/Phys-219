# orginal  version  Evan Kiefl; modifications  Evan Thomas  and  R Kiefl
# comments are typcially used to explain the following line of code
#
############################################################

# import the  library numpy  and rename it  np
import numpy as np
# import the library matplotlib   and rename it plot
import matplotlib.pyplot as plt
#name  the input file  with the data
fname = 'lab5x4x3.csv'

# read in data - the file is assumed to be in csv format (comma separated variables). 
#Files need to be specified with a full path OR they have to be saved in the same folder 
#as the script
data = np.loadtxt(fname, delimiter=',', comments='#',usecols=(0,1))
#data = np.loadtxt(fname, delimiter=',',comments='#' )
# access the data columns and assign variables xraw and yraw
#generate  an array  xraw  which is the first  column  of  data.  Note the first column is 
#indexed as  zero.
xraw = data[:,0]
#generate  an array  yraw  which is the second  column  of  data  (index  1)
yraw = data[:,1]
#define packing factor  npac
npac=10
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
#note the data is not copied during this process - x,y are 'pointing' to the same 
#memory as data
#define the uncertainty in y. 
sigmay=0.003
# plot the data
plt.errorbar(x, y,yerr=sigmay,marker='s',linestyle='')
## marker='o' : use markers to indicate each data point (x_1,y_1),(x_2,y_2)
## linestyle= '' : no line is drawn to connect the data points
## linestyle= '-' : a line is drawn to connect the data points

# add axis labels
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Voltage output vs. Time')
# this next command makes  sure the plot is shown. 
plt.show()
# ----- you can save a figure by right-clicking on the graph in the console
# ----- alternatively use: plt.savefig("NAMEOFFIGURE")