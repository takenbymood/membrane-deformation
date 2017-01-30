import random

import numpy as np
import math
import time
import sys
import subprocess
import os
import fileinput

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

import lammpsbuilder as lb

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

	def __init__(self,genes=6,ljEpsPlaces=6,ljSigmaPlaces=6,ligRadPlaces=6,ligAngPlaces=6,maxRadius=6,maxEps=5,maxSigma=2,maxAngle=6.283185,minRadius=2,minEps=0,minSigma=1,minAngle=0):
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
		num = grayToNumber(individual)
		run = "250000"
		sim = MembraneSimulation("p_"+str(num),p,run=run,dumpres=run)
		sim.saveFiles()
		dir_path = os.path.dirname(os.path.realpath(__file__))
		path = dir_path+"/"+sim.filedir
		print ("lammps -in "+path)
		try:
			proc = subprocess.Popen("cd "+ path + " && lammps -in "+sim.scriptName,shell=True)
			proc.wait()
		except: 
			return 1E10,
		sim.deleteFiles()

		outData = []

		with open(dir_path+"/out/out/"+"p_"+str(num)+"_out.xyz", 'r+') as f:
			lines = f.readlines()
			for i in range(len(lines)):
				if run in lines[i]:
					lines[i] = ""
					break
				lines[i] = ""

			for line in lines:
				if line != "":
					outData.append(line.replace("\n","").replace(" ",","))

		os.remove(dir_path+"/out/out/"+"p_"+str(num)+"_out.xyz")

		if len(outData)<100:
			return 1E10,

		outVectors = {}
		for line in outData:
			slist = line.split(",")
			if(len(slist)<3):
				return 1E10,
			if int(slist[0]) in outVectors:
				outVectors[int(slist[0])].append({'x':float(slist[1]),'y':float(slist[2])})
			else:
				outVectors[int(slist[0])] = []
				outVectors[int(slist[0])].append({'x':float(slist[1]),'y':float(slist[2])})

		magnitudes = []
		for key, value in outVectors.iteritems():
			if key > 3:
				for v in value:
					currentMin = 1E10
					for v2 in outVectors[1]:
						xd = v['x']-v2['x']
						yd = v['y']-v2['y']
						m = xd*xd + yd*yd
						if m<currentMin:
							currentMin = m
					magnitudes.append(currentMin)


		if len(magnitudes)<1:
			return 1E10,

		msum = 0
		for m in magnitudes:
			msum += m

		msum = msum/float(len(magnitudes))

		return msum,

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
			p = self.genome.constructProtein(i)
			print(p)
			num = grayToNumber(i)
			run = "250000"
			sim = MembraneSimulation("p_"+str(num),p,run=run,dumpres="500")
			sim.saveFiles()
			dir_path = os.path.dirname(os.path.realpath(__file__))
			path = dir_path+"/"+sim.filedir
			print ("lammps -in "+path)
			proc = subprocess.Popen("cd "+ path + " && lammps -in "+sim.scriptName,shell=True)
			proc.wait()
		return pop



class State:

	

	def __init__(self):
		creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
		creator.create("Individual", list, fitness=creator.FitnessMin)
		self.instances = []
		self.populations = []

	def registerInstance(self,genome = Genome(),mutationRate=0.1):
		self.instances.append(Algorithm(genome,mutationRate))

	def run(self,popSize=100,crossPb=0.5,mutPb=0.2,nGen=100,log=False):
		self.populations = []
		for i in self.instances:
			self.populations.append(i.run(popSize,crossPb,mutPb,nGen,log))
		return self.populations


class MembraneSimulation(lb.LammpsSimulation):

	def __init__(self,name,protein,mLength=100,spacing=1.5,corepos_x=0, corepos_y=6,run="250000",dumpres="100"):
		lb.LammpsSimulation.__init__(self,name,"out/",run=run)
		self.script.dump = "id all xyz "+dumpres+" out/" + name +"_out.xyz"
		self.data.atomTypes = 3+len(protein.ligands)
		self.data.bondTypes = 1
		self.data.angleTypes = 1
		self.data.addMass(1,1)
		self.data.addMass(2,1)
		self.data.addMass(3,2)
		for i in range(len(protein.ligands)):
			self.data.addMass(4+i,0.01)

		startX = -(0.5*mLength*spacing)

		self.data.addAtom(2,startX,0)

		for i in range(mLength-2):
			self.data.addAtom(1,startX+spacing*i+2,0)
			self.data.addBond(1,i+1,i+2)
			self.data.addAngle(1,i+1,i+2,i+3)

		self.data.addAtom(2,startX+spacing*mLength,0)
		self.data.addBond(1,mLength-1,mLength)

		mol = self.data.addAtom(3,corepos_x,corepos_y,0)

		self.script.addBond(1,2.0,1.3)
		self.script.addAngle(1,30,180)
		self.script.addPair("*","*",0,0,0)

		aType = 4
		for l in protein.ligands:
			self.data.addAtom(aType,corepos_x+l.rad*math.cos(l.ang),corepos_y+l.rad*math.sin(l.ang),0,mol)
			self.script.addPair("1",str(aType),l.eps,l.sig,l.cutoff)
			aType+=1
		
		self.script.addPair(1,3,1,4,5.6123)

		self.script.addGroup("move",[1])
		self.script.addGroup("anchor",[2])
		pGroup = [3]
		for i in range(len(protein.ligands)):
			pGroup.append(i+4)
		self.script.addGroup("protein",pGroup)

		self.script.addFix("move","nve")
		self.script.addFix("protein","rigid/nve molecule")
		self.script.addFix("all","langevin 1 1 1 1000")
		self.script.addFix("all","enforce2d")

def main():

	state = State()
	state.registerInstance(Genome(),0.1)
	p = state.run(10,0.5,0.2,10,False)
	
if __name__ == "__main__":
	main()
