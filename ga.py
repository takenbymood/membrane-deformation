import random

import numpy as np

from deap import base
from deap import creator
from deap import tools

flatten = lambda l: [item for sublist in l for item in sublist]
partition = lambda l,s : [l[x:x+s] for x in xrange(0, len(l), s)]

CHROMOSOME_SIZE = 8
GENE_SIZE = 30

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
toolbox.register("attr_int", random.randint,0,1)
toolbox.register("gene",tools.initRepeat,creator.Individual,toolbox.attr_int,n=GENE_SIZE)
toolbox.register("individual", tools.initRepeat, creator.Individual,
				 toolbox.gene, n=CHROMOSOME_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def evaluate(individual):
	return sum(flatten(individual)),

toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate",tools.mutUniformInt, low = 1, up = 1,indpb=0.1)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)


def main():

	pop = toolbox.population(n=10)
	len(pop[0])
	CXPB, MUTPB, NGEN = 0.5, 0.5, 100
	# Evaluate the entire population
	fitnesses = map(toolbox.evaluate, pop)
	for ind, fit in zip(pop, fitnesses):
		ind.fitness.values = fit

	for g in range(NGEN):

		# Select the next generation individuals
		offspring = tools.selBest(pop, 4)
		# Clone the selected individuals
		offspring = map(toolbox.clone, offspring)

		# Apply crossover and mutation on the offspring
		for child1, child2 in zip(offspring[::2], offspring[1::2]):
			if random.random() < CXPB:
				toolbox.mate(child1, child2)
				del child1.fitness.values
				del child2.fitness.values

		for mutant in offspring:
			if random.random() < MUTPB:
				toolbox.mutate(mutant[random.randint(0,CHROMOSOME_SIZE-1)])
				del mutant.fitness.values

		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses = map(toolbox.evaluate, invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit

		# The population is entirely replaced by the offspring
		pop[:] = offspring

		fits = [ind.fitness.values[0] for ind in pop]
		
		length = len(pop)
		mean = sum(fits) / length
		sum2 = sum(x*x for x in fits)
		std = abs(sum2 / length - mean**2)**0.5
		
		print("  Min %s" % min(fits))
		print("  Max %s" % max(fits))
		print("  Avg %s" % mean)
		print("  Std %s" % std)
	
	print(pop[0])
	return pop


if __name__ == "__main__":
	main()
