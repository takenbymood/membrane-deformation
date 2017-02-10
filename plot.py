import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

parser = argparse.ArgumentParser(description='Plot fitness logs from GA')
parser.add_argument('--in','-i', dest='filepath', required=False,
                  	default=None,
                    help='csv file for plotting')

args = parser.parse_args()

# example data
data = np.genfromtxt(args.filepath, delimiter='	', skip_header=2,
                     skip_footer=0, names=['GEN','N','NV','AVG', 'STD','MIN','MAX','NOV'])

x = data['GEN']
yorig = data['AVG']
print data
y = 1.0/yorig
y2 = 1.0/data['MIN']
err = data['STD']
errN = map(math.sqrt,data['NV'])
stdErr = err/errN
errPerc = (stdErr/yorig)

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Genome Fitness")    
ax1.set_xlabel('Generation Number')
ax1.set_ylabel('Fitness')

ax1.plot(x, y, color='r')
ax1.errorbar(x,y, yerr=y*errPerc, fmt='o')
ax1.plot(x,y2)


plt.ylim(np.min(y)-0.1*(np.max(y)-np.min(y)),np.max(y2)+0.2*(np.max(y)-np.min(y)))

plt.show()