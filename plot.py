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
                     skip_footer=0, names=['GEN','N','AVG', 'STD','MIN','MAX'])

x = data['GEN']
y = 1.0/data['AVG']
y2 = 1.0/data['MIN']
err = data['STD']
errN = map(math.sqrt,data['N'])

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Genome Fitness")    
ax1.set_xlabel('Generation Number')
ax1.set_ylabel('Fitness')

ax1.plot(x, y, color='r')
ax1.errorbar(x,y, fmt='o')
ax1.plot(x,y2)

plt.show()