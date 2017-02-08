import random

import numpy as np
import math
import time
import sys
import subprocess
import os
import fileinput
import time

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

from threading import Timer

import lammpsbuilder as lb

flatten = lambda l: [item for sublist in l for item in sublist]
partition = lambda l,s : [l[x:x+s] for x in xrange(0, len(l), s)]
rmcase = lambda l,n: [x for x in l if not (x==n)]

def grayToNumber(g):
	b=[g[0]]
	for i in range(len(g)-1):
		b.append(b[i]^g[i+1])
	out = 0
	for bit in b:
		out = (out << 1) | bit
	return out

def kill(p):
    try:
    	print "killing "+str(p) 
        p.kill()
    except OSError:
        pass # ignore

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

	def __init__(self,genes=6,ljEpsPlaces=4,ljSigmaPlaces=1,ligRadPlaces=1,ligAngPlaces=4,maxRadius=4,maxEps=10,maxSigma=2,maxAngle=6.283185,minRadius=4,minEps=2,minSigma=2,minAngle=0):
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

	def __init__(self,genome=Genome(),mutationRate=0.2,hofSize=10,runtime=250000):

		self.genome = genome
		self.toolbox = base.Toolbox()
		self.toolbox.register("bit", random.randint,0,1)
		self.toolbox.register("individual", tools.initRepeat, creator.Individual,
						 self.toolbox.bit, n=self.genome.size)
		self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

		self.toolbox.register("mate", self.mate)
		self.toolbox.register("mutate",tools.mutFlipBit,indpb=mutationRate)
		self.toolbox.register("crossover",self.crossover)
		self.toolbox.register("select", tools.selTournament, tournsize=4)
		self.toolbox.register("evaluate", self.evaluate)

		self.runtime = str(runtime)

		self.currentP = 0

		self.hof = tools.HallOfFame(hofSize)

		self.stamp = str(int(math.floor(time.time())))

		self.maxFit = 1E20

	def crossover(self,ind1,ind2):
		pos = random.randint(1,self.genome.genes)
		pos2 = random.randint(pos,self.genome.genes)
		a = partition(ind1,self.genome.geneSize)
		b = partition(ind2,self.genome.geneSize)


		ap1 = a[0:pos]
		ap2 = a[pos:pos2]

		bp1 = b[0:pos]
		bp2 = b[pos:pos2]

		ap3 = a[pos2:self.genome.genes]
		bp3 = b[pos2:self.genome.genes]

		ind3 = flatten(ap1+bp2+ap3)
		ind4 = flatten(bp1+ap2+bp3)

		return ind3,ind4

	def mate(self,ind1,ind2):
		child1, child2 = [self.toolbox.clone(ind) for ind in (ind1, ind2)]
		self.toolbox.crossover(child1, child2)
		del child1.fitness.values
		del child2.fitness.values
		return child1,child2



	def evaluate(self,individual):
		maxFit = self.maxFit
		p = self.genome.constructProtein(individual)
		num = self.currentP
		self.currentP+=1
		sim = MembraneSimulation("p_"+str(num),p,"xyz/","",run=self.runtime,dumpres=self.runtime)
		sim.saveFiles()
		dir_path = os.path.dirname(os.path.realpath(__file__))
		path = dir_path+"/"+sim.filedir
		try:
			proc = subprocess.Popen("cd "+ path + " && lammps -in "+sim.scriptName+" > lammps.out",shell=True)
			t = Timer(60, kill, [proc])
			t.start()
			proc.wait()
			t.cancel()
		except: 
			individual.fitness.valid = False
			#print "crashed or failed"
			return maxFit,
		

		outData = []

		time.sleep(0.5)

		sim.deleteFiles()

		outFilename = dir_path+"/out/xyz/"+"p_"+str(num)+"_out.xyz"

		if(not os.path.exists(outFilename)):
			#print "no outfile"
			return maxFit,

		with open(outFilename, 'r+') as f:
			lines = f.readlines()
			for i in range(len(lines)):
				if self.runtime in lines[i]:
					lines[i] = ""
					break
				lines[i] = ""

			for line in lines:
				if line != "":
					outData.append(line.replace("\n","").replace(" ",","))


		os.remove(outFilename)

		if len(outData)<50:
			#print str(outData)
			return maxFit,

		outVectors = {}
		for line in outData:
			slist = line.split(",")
			if(len(slist)<3):
				return maxFit,
			if int(slist[0]) in outVectors:
				outVectors[int(slist[0])].append({'x':float(slist[1]),'y':float(slist[2])})
			else:
				outVectors[int(slist[0])] = []
				outVectors[int(slist[0])].append({'x':float(slist[1]),'y':float(slist[2])})

		magnitudes = []
		boxsize = 20
		for key, value in outVectors.iteritems():
			if key == 3:
				for v in value:
					inrange = 0
					fmag = 0
					for v2 in outVectors[1]:
						xd = v['x']-v2['x']
						yd = v['y']-v2['y']
						m = math.sqrt(xd*xd+yd*yd)
						if(m<7):
							inrange+=1
					if(inrange>0):
						magnitudes.append(inrange)

		if len(magnitudes)<1:
			#print str(num) + " protein out of range"
			return maxFit,

		msum = 0
		for m in magnitudes:
			msum += m

		if(msum == 0):
			#print "no msum"
		 	return maxFit,

		return 1.0/msum,

	def clmean(self,l):
		return np.mean(rmcase(l,(self.maxFit,)))

	def clstd(self,l):
		return np.std(rmcase(l,(self.maxFit,)))

	def clmin(self,l):
		return np.min(rmcase(l,(self.maxFit,)))

	def clmax(self,l):
		return np.max(rmcase(l,(self.maxFit,)))

	def clnum(self,l):
		return len(rmcase(l,(self.maxFit,)))

	def run(self,popSize=100,CXPB=0.5,MUTPB=0.2,NGEN=100,log=True):

		pop = self.toolbox.population(n=popSize)
		
		os.mkdir("out/"+self.stamp)
		self.logfile = None
		if(log):
			self.logfile = open("out/"+self.stamp+"/ft.tsv", 'w')
		# Evaluate the entire population

		
		self.stats = tools.Statistics(lambda ind: ind.fitness.values)
		self.stats.register("NValid", self.clnum)
		self.stats.register("Avg", self.clmean)
		self.stats.register("Std", self.clstd)
		self.stats.register("Min", self.clmin)
		self.stats.register("Max", self.clmax)

		if(log):
			orig_stdout = sys.stdout
			sys.stdout = self.logfile

		algorithms.eaSimple(pop, self.toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, stats=self.stats,
                        halloffame=self.hof, verbose=True)

		if(log):
			sys.stdout = orig_stdout
			self.logfile.close()

		tag=1
		with open("out/"+self.stamp+"/hof.out", 'w') as file_:
			for i in self.hof:
				p = self.genome.constructProtein(i)
				file_.write(str(i))
				file_.write(str(p))
				num = grayToNumber(i)
				sim = MembraneSimulation("hof_"+str(tag),p,"",str(self.stamp),run=self.runtime,dumpres="100")
				sim.filedir = "out/"+self.stamp+"/"
				sim.saveFiles()
				dir_path = os.path.dirname(os.path.realpath(__file__))
				path = dir_path+"/"+sim.filedir
				proc = subprocess.Popen("cd "+ path + " && lammps -in "+sim.scriptName+" > hoflog.out",shell=True)
				proc.wait()
				tag+=1

		with open("out/"+self.stamp+"/pop.out", 'w') as file_:
			for i in pop:
				p = self.genome.constructProtein(i)
				file_.write(str(i))
				file_.write(str(p))

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

	def __init__(self,name,protein,outdir,filedir,mLength=70,spacing=2.0,corepos_x=0, corepos_y=10,run="200000",dumpres="100"):
		lb.LammpsSimulation.__init__(self,name,"out/",run=run)
		self.script.dump = "id all xyz "+dumpres+" "+outdir+name +"_out.xyz"
		self.data.atomTypes = 3+len(protein.ligands)
		self.data.bondTypes = 1
		self.data.angleTypes = 1
		self.data.addMass(1,1)
		self.data.addMass(2,1)
		self.data.addMass(3,3)
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
		self.script.addAngle(1,40,180)
		self.script.addPair("*","*",0,0,0)

		aType = 4
		for l in protein.ligands:
			self.data.addAtom(aType,corepos_x+l.rad*math.cos(l.ang),corepos_y+l.rad*math.sin(l.ang),0,mol)
			self.script.addPair("1",str(aType),l.eps,l.sig,l.sig*l.cutoff)
			aType+=1
		
		self.script.addPair(1,3,100,4,4.4)
		self.script.addPair(1,1,100,2,2.24)
		self.script.addPair(1,2,100,1,1.1)
		self.script.addPairModify("shift yes")

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
		self.script.addLine("fix 4 all wall/lj93 yhi 18 1.0 1.0 2.5")
		

def main():
	state = State()
	state.registerInstance(Genome(),0.2)
	p = state.run(40,0.5,0.2,40,False)
	
if __name__ == "__main__":
	main()
