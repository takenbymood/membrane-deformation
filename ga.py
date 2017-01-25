import random

import numpy as np
import math

from deap import base
from deap import creator
from deap import tools

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

	def __init__(self,genome=Genome(),mutationRate=0.1):

		self.genome = genome
		self.toolbox = base.Toolbox()
		self.toolbox.register("bit", random.randint,0,1)
		self.toolbox.register("individual", tools.initRepeat, creator.Individual,
						 self.toolbox.bit, n=self.genome.size)
		self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

		self.toolbox.register("mate", tools.cxTwoPoint)
		self.toolbox.register("mutate",tools.mutUniformInt, low = 0, up = 1,indpb=mutationRate)
		self.toolbox.register("select", tools.selTournament, tournsize=3)
		self.toolbox.register("evaluate", self.evaluate)

	def evaluate(self,individual):
		return sum(individual),

	def run(self,popSize=10,CXPB=0.5,MUTPB=0.5,NGEN=100):
		pop = self.toolbox.population(n=popSize)
		# Evaluate the entire population
		fitnesses = map(self.toolbox.evaluate, pop)
		for ind, fit in zip(pop, fitnesses):
			ind.fitness.values = fit

		for g in range(NGEN):

			# Select the next generation individuals
			offspring = tools.selBest(pop, 4)
			# Clone the selected individuals
			offspring = map(self.toolbox.clone, offspring)

			# Apply crossover and mutation on the offspring
			for child1, child2 in zip(offspring[::2], offspring[1::2]):
				if random.random() < CXPB:
					self.toolbox.mate(child1, child2)
					del child1.fitness.values
					del child2.fitness.values

			for mutant in offspring:
				if random.random() < MUTPB:
					self.toolbox.mutate(mutant)
					del mutant.fitness.values

			# Evaluate the individuals with an invalid fitness
			invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
			fitnesses = map(self.toolbox.evaluate, invalid_ind)
			for ind, fit in zip(invalid_ind, fitnesses):
				ind.fitness.values = fit

			# The population is entirely replaced by the offspring
			pop[:] = offspring

			fits = [ind.fitness.values[0] for ind in pop]
			
			length = len(pop)
			mean = sum(fits) / length
			sum2 = sum(x*x for x in fits)
			std = abs(sum2 / length - mean**2)**0.5
			
			# print("  Min %s" % min(fits))
			# print("  Max %s" % max(fits))
			# print("  Avg %s" % mean)
			# print("  Std %s" % std)
		
		return pop



class State:

	instances = []

	def __init__(self):
		creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
		creator.create("Individual", list, fitness=creator.FitnessMin)

	def registerInstance(self,genome = Genome(),mutationRate=0.1):
		self.instances.append(Algorithm(genome,mutationRate))

	def run(self,popSize=10,crossPb=0.5,mutPb=0.5,nGen=100):
		for i in self.instances:
			i.run(popSize,crossPb,mutPb,nGen)


def main():

	state = State()
	state.registerInstance(Genome(),0.1)
	state.run(10,0.5,0.5,100)

if __name__ == "__main__":
	main()
