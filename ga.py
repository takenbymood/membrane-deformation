import random

import numpy as np
import math

from deap import base
from deap import creator
from deap import tools

flatten = lambda l: [item for sublist in l for item in sublist]
partition = lambda l,s : [l[x:x+s] for x in xrange(0, len(l), s)]

#Chromosome properties

class Chromosome:

	def __init__(self,ljEpsPlaces=6,ljSigmaPlaces=6,ligRadPlaces=6,ligAngPlaces=6,ljNum=4,ligNum=6):
		self.ljEpsPlaces = ljEpsPlaces
		self.ljSigmaPlaces = ljSigmaPlaces
		self.ligRadPlaces = ligRadPlaces
		self.ligAngPlaces = ligAngPlaces
		self.ljNum = ljNum #must be a power of 2 for stable results
		self.ligNum = ligNum

		self.ljGeneSize = self.ljEpsPlaces + self.ljSigmaPlaces
		self.ligTypePlaces = int(math.floor(math.log(self.ljNum,2)))
		self.ligGeneSize = self.ligTypePlaces+self.ligRadPlaces+self.ligAngPlaces

		self.size = self.ljNum*self.ljGeneSize + self.ligNum*self.ligGeneSize


class Algorithm:

	def __init__(self,chromosome=Chromosome(),mutationRate=0.1):

		self.chromosome = chromosome
		self.toolbox = base.Toolbox()
		self.toolbox.register("bit", random.randint,0,1)
		self.toolbox.register("individual", tools.initRepeat, creator.Individual,
						 self.toolbox.bit, n=self.chromosome.size)
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

	def registerInstance(self,chromosome = Chromosome(),mutationRate=0.1):
		self.instances.append(Algorithm(chromosome,mutationRate))

	def run(self,popSize=10,crossPb=0.5,mutPb=0.5,nGen=100):
		for i in self.instances:
			i.run(popSize,crossPb,mutPb,nGen)


def main():

	state = State()
	state.registerInstance(Chromosome(6,6,6,6,4,6),0.1)
	state.run(10,0.5,0.5,100)

if __name__ == "__main__":
	main()
