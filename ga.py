import random

import numpy as np
import math
import time
import sys

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

flatten = lambda l: [item for sublist in l for item in sublist]
partition = lambda l,s : [l[x:x+s] for x in xrange(0, len(l), s)]

#Genome properties


class Ligand:
	def __init__(self,rad,ang,eps,sig,mass=0.01,cutoff=2.5):
		self.rad = rad
		self.ang = ang
		self.eps = eps
		self.sig = sig
		self.mass = mass
		self.cutoff = cutoff

	def __str__(self):
		ligstr= "rad:"+str(self.rad)+", ang:"+str(self.ang)+", eps:"+str(self.eps)+", sig:"+str(self.sig)
		ligstr+= ", mass:"+str(self.mass)+", cutoff:"+str(self.cutoff)
		return ligstr

class Protein:
	ligands=[]
	def __init__(self,x=0,y=0,mass=1,eps=1,sig=4,cutoff=2**(1/6)):
		self.x = x
		self.y = y
		self.mass = mass
		self.eps = eps
		self.sig = sig
		self.cutoff = cutoff

	def addLigand(self,ligand):
		if(isinstance(ligand,Ligand)):
			self.ligands.append(ligand)

	def __str__(self):
		protstr = "x:"+str(self.x)+", y:"+str(self.y)+", m:"+str(self.mass)
		protstr += ", eps:"+str(self.eps)+", sig:"+str(self.sig)+", cut:"+str(self.cutoff)
		protstr += "\n"+str(len(self.ligands))+" ligands"
		i=0
		for l in self.ligands:
			i+=1
			protstr += "\nligand " + str(i) +" - "+ str(l)
		return protstr

class Genome:

	def __init__(self,genes=20,ljEpsPlaces=6,ljSigmaPlaces=6,ligRadPlaces=6,ligAngPlaces=6):
		self.ljEpsPlaces = ljEpsPlaces
		self.ljSigmaPlaces = ljSigmaPlaces
		self.ligRadPlaces = ligRadPlaces
		self.ligAngPlaces = ligAngPlaces
		self.geneSize = ljEpsPlaces+ljSigmaPlaces+ligRadPlaces+ligAngPlaces

		self.size = genes*self.geneSize


class Algorithm:

	def __init__(self,genome=Genome(),mutationRate=0.1,hofSize=5):

		self.genome = genome
		self.toolbox = base.Toolbox()
		self.toolbox.register("bit", random.randint,0,1)
		self.toolbox.register("individual", tools.initRepeat, creator.Individual,
						 self.toolbox.bit, n=self.genome.size)
		self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

		self.toolbox.register("mate", tools.cxTwoPoint)
		self.toolbox.register("mutate",tools.mutFlipBit,indpb=mutationRate)
		self.toolbox.register("select", tools.selTournament, tournsize=3)
		self.toolbox.register("evaluate", self.evaluate)

		self.hof = tools.HallOfFame(hofSize)

	def evaluate(self,individual):
		return sum(individual),

	def run(self,popSize=100,CXPB=0.5,MUTPB=0.2,NGEN=100,log=True):

		pop = self.toolbox.population(n=popSize)
		if(log):
			logfile = open("out/ft_"+str(int(math.floor(time.time())))+".tsv", 'w')
		# Evaluate the entire population

		
		stats = tools.Statistics(lambda ind: ind.fitness.values)
		stats.register("Avg", np.mean)
		stats.register("Std", np.std)
		stats.register("Min", np.min)
		stats.register("Max", np.max)

		if(log):
			orig_stdout = sys.stdout
			sys.stdout = logfile

		algorithms.eaSimple(pop, self.toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, stats=stats,
                        halloffame=self.hof, verbose=True)

		if(log):
			sys.stdout = orig_stdout
			logfile.close()

		return pop



class State:

	instances = []
	populations = []

	def __init__(self):
		creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
		creator.create("Individual", list, fitness=creator.FitnessMin)

	def registerInstance(self,genome = Genome(),mutationRate=0.1):
		self.instances.append(Algorithm(genome,mutationRate))

	def run(self,popSize=100,crossPb=0.5,mutPb=0.2,nGen=100):
		self.populations = []
		for i in self.instances:
			self.populations.append(i.run(popSize,crossPb,mutPb,nGen))
		return self.populations

def main():

	state = State()
	state.registerInstance(Genome(),0.1)
	p = state.run(100,0.5,0.2,100)
if __name__ == "__main__":
	main()
