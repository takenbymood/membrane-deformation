import random

import numpy as np
import math
import time
import sys

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

import lammpsbuilder

flatten = lambda l: [item for sublist in l for item in sublist]
partition = lambda l,s : [l[x:x+s] for x in xrange(0, len(l), s)]

def grayToNumber(g):
	b=[g[0]]
	for i in range(len(g)-1):
		b.append(b[i]^g[i+1])
	out = 0
	for bit in b:
		out = (out << 1) | bit
	return out

class Ligand:
	def __init__(self,eps,sig,rad,ang,mass=0.01,cutoff=2.5):
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
	ligands = []
	def __init__(self,x=0,y=0,mass=1,eps=1,sig=4,cutoff=2**(1/6)):
		self.x = x
		self.y = y
		self.mass = mass
		self.eps = eps
		self.sig = sig
		self.cutoff = cutoff
		self.ligands=[]

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

	def __init__(self,genes=6,ljEpsPlaces=6,ljSigmaPlaces=6,ligRadPlaces=6,ligAngPlaces=6,maxRadius=15,maxEps=20,maxSigma=10,maxAngle=6.283185,minRadius=0,minEps=0,minSigma=0,minAngle=0):
		self.ljEpsPlaces = ljEpsPlaces
		self.ljSigmaPlaces = ljSigmaPlaces
		self.ligRadPlaces = ligRadPlaces
		self.ligAngPlaces = ligAngPlaces
		self.geneSize = ljEpsPlaces+ljSigmaPlaces+ligRadPlaces+ligAngPlaces
		self.genes = genes
		self.size = genes*self.geneSize

		self.maxRadius = maxRadius
		self.maxAngle = maxAngle
		self.maxEps = maxEps
		self.maxSigma = maxSigma

		self.minRadius = minRadius
		self.minAngle = minAngle
		self.minEps = minEps
		self.minSigma = minSigma

		self.invMaxRad = 1.0/(2**ligRadPlaces)
		self.invMaxAngle = 1.0/(2**ligAngPlaces)
		self.invMaxEps = 1.0/(2**ljEpsPlaces)
		self.invMaxSig = 1.0/(2**ljSigmaPlaces)	

	def constructProtein(self,individual):
		p = Protein()
		for i in range(self.genes):
			gPos = self.geneSize*i
			gStart,gEnd = gPos,gPos + self.ljEpsPlaces
			eps = grayToNumber(individual[gStart:gEnd])*self.invMaxEps*(self.maxEps-self.minEps)+self.minEps
			gStart,gEnd = gEnd+1,gEnd + 1 + self.ljSigmaPlaces
			sig = grayToNumber(individual[gStart:gEnd])*self.invMaxSig*(self.maxSigma-self.minSigma)+self.minSigma
			gStart,gEnd = gEnd+1,gEnd + 1 + self.ligRadPlaces
			rad = grayToNumber(individual[gStart:gEnd])*self.invMaxRad*(self.maxRadius-self.minRadius)+self.minRadius
			gStart,gEnd = gEnd+1,gEnd + 1 + self.ligAngPlaces
			ang = grayToNumber(individual[gStart:gEnd])*self.invMaxAngle*(self.maxAngle-self.minAngle)+self.minAngle
			p.addLigand(Ligand(eps,sig,rad,ang))

		return p


class Algorithm:

	def __init__(self,genome=Genome(),mutationRate=0.1,hofSize=1):

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
		p = self.genome.constructProtein(individual)
		rsum = 0
		for l in p.ligands:
			rsum += l.rad

		return rsum,

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

		for i in self.hof:
			print self.genome.constructProtein(i)

		return pop



class State:

	instances = []
	populations = []

	def __init__(self):
		creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
		creator.create("Individual", list, fitness=creator.FitnessMin)

	def registerInstance(self,genome = Genome(),mutationRate=0.1):
		self.instances.append(Algorithm(genome,mutationRate))

	def run(self,popSize=100,crossPb=0.5,mutPb=0.2,nGen=100,log=False):
		self.populations = []
		for i in self.instances:
			self.populations.append(i.run(popSize,crossPb,mutPb,nGen,log))
		return self.populations

def main():

	state = State()
	state.registerInstance(Genome(),0.1)
	p = state.run(10,0.5,0.2,100,False)
	data = lammpsbuilder.LammpsData(1,1,1)
	data.addAtom(1,5,5)
	data.addAtom(1,10,5)
	data.addAtom(1,15,5)
	data.addBond(1,1,2)
	data.addAngle(1,1,2,3)
	print (data)
if __name__ == "__main__":
	main()
