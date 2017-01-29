class LammpsScript:
	sep = "			"
	dimension = ""
	units = ""
	atom_style = ""
	atom_modify = ""
	read_data = ""
	dump = ""
	neighbour = ""
	neigh_modify = ""

	thermo = ""
	timestep = ""
	run = ""
	velocity = ""

	bond_style=""
	angle_style=""
	pair_style=""
	
	bond_coeffs=[]
	angle_coeffs=[]
	pair_coeffs=[]
	groups = []
	fixes = []

	def addPair(self,atom1,atom2,eps=0,sig=0,cutoff=""):
		s = str(atom1)+" "+str(atom2)+" "+str(eps)+" "+str(sig)+" "+str(cutoff)
		self.pair_coeffs.append(s)

	def addBond(self,bond,K,x0):
		s = str(bond)+" "+str(K)+" "+str(x0)
		self.bond_coeffs.append(s)

	def addAngle(self,angle,K,theta0):
		s = str(angle)+" "+str(K)+" "+str(theta0)
		self.angle_coeffs.append(s)

	def addGroup(self,name,members,order="type"):
		s = str(name)+" "+str(order)+" "
		for m in members:
			s+=str(m)+" "
		self.groups.append(s)

	def addFix(self,group,action):
		s = str(group)+ " " + str(action)
		self.fixes.append(s)

	def __init__(self,read_data="",dump="id all xyz 100 twodim.xyz",thermo="300",timestep="0.004",run="150000",dimension="2",units="lj",velocity="all create 1.0 1000",atom_style="molecular",atom_modify="sort 0 1",neighbour="0.3 bin",neigh_modify="every 1 delay 1",angle_style="harmonic",bond_style="harmonic",pair_style="lj/cut 2.5"):
		self.read_data = read_data
		self.dump = dump
		self.dimension = dimension
		self.units = units
		self.atom_style = atom_style
		self.atom_modify = atom_modify
		self.neighbour = neighbour
		self.neigh_modify = neigh_modify
		self.angle_style = angle_style
		self.bond_style = bond_style
		self.pair_style = pair_style
		self.timestep = timestep
		self.thermo = thermo
		self.run = run
		self.velocity = velocity

	def __str__(self):
		d = self.sep
		s="dimension"+d+self.dimension+"\n"
		s+="units"+d+self.units+"\n"
		s+="atom_style"+d+self.atom_style+"\n"
		s+="atom_modify"+d+self.atom_modify+"\n"
		s+="\n"
		s+="read_data"+d+self.read_data+"\n"
		s+="neighbour"+d+self.neighbour+"\n"
		s+="neigh_modify"+d+self.neigh_modify+"\n"
		s+="bond_style"+d+self.bond_style+"\n"
		for b in self.bond_coeffs:
			s+="bond_coeff"+d+b+"\n"
		s+="angle_style"+d+self.angle_style+"\n"
		for a in self.angle_coeffs:
			s+="angle_coeff"+d+a+"\n"
		s+="pair_style"+d+self.pair_style+"\n"
		for p in self.pair_coeffs:
			s+="pair_coeff"+d+p+"\n"
		s+="velocity"+d+self.velocity+"\n"
		s+="dump"+d+self.dump+"\n"

		for g in self.groups:
			s+="group"+d+g+"\n"

		i=0
		for f in self.fixes:
			s+="fix"+d+str(i)+" "+f+"\n"
			i+=1

		s+="\n"
		s+="thermo"+d+self.thermo+"\n"
		s+="timestep"+d+self.timestep+"\n"
		s+="run"+d+self.run+"\n"

		return s