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
data = np.genfromtxt(args.filepath, delimiter='	', skip_header=10,
                     skip_footer=10, names=['GEN','N','AVG', 'MIN','MAX','STD'])

x = data['GEN']
y = data['AVG']
err = data['STD']
errN = map(math.sqrt,data['N'])

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Genome Fitness")    
ax1.set_xlabel('Generation Number')
ax1.set_ylabel('Fitness')

ax1.plot(x, y, color='r')
ax1.errorbar(x,y, yerr=err/errN, fmt='o')

plt.show()