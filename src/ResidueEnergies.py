import re





###################################################################################################################
#Residue Energies: jede Zeile der PDB-Datei entspricht einer Instanz


#testpdb einlesen
filename='testpdb.txt'
f = open(filename)
lines = f.readlines()



scoreterms = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 'fa_elec', 'pro_close', 'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc', 'dslf_fa13', 'atom_pair_constraint', 'coordinate_constraint', 'angle_constraint', 'dihedral_constraint', 'rama', 'omega', 'fa_dun', 'p_aa_pp', 'ref', 'chainbreak', 'res_type_constraint', 'total']


#start- und endpunkt beim einlesen der datei
for i in range(len(lines)):
	if re.match('#BEGIN_POSE_ENERGIES_TABLE', lines[i]):
		start = i+4
	elif re.match('#END_POSE_ENERGIES_TABLE', lines[i]):
		end = i

#print '\n reading from line %s to line %s \n' %(start, end) 




#Klasse ResidueEnergies erstellen:
class ResidueEnergies(object):

	#legt fest wie man einzelne instances abfragen kann. muessen name und argumente enthalten
    def __init__(self, name, arguments):
        self.name = name
        self.arguments = arguments
 
 	#ermoeglicht die abfrage eines namens
    def getName(self):
        return self.name
 	
 	#ermoeglicht die abfrage von score terms (argumente) fuer bestimmten residues
    def getArguments(self):
        return self.arguments

    def getArgument(self, term):
    	args = self.arguments
    	for i in range(len(scoreterms)):
    		if term == scoreterms[i]:
    			print args[i]
    		
 	#beschreibt eine instance der klasse
    def __str__(self):
        return "%s has the following score terms \n %s" % (self.name, self.arguments)


 
#namelist enthaelt alle instances der klasse ResidueEnergies
namelist=[]

for i in range(end):
	if i >= start:
		entry = lines[i].split()
		name_entry = entry[0]
		#print name_entry
		argument_entry = entry[1:]
		name_entry = ResidueEnergies(name_entry, argument_entry)
		#print name_entry.getArguments()
		namelist.append(name_entry)


print '\n', namelist[0],'\n'


#einzelne scoreterms fuer bestimmte position abfragen
#print namelist[0].getArgument('fa_atr')


#print namelist[0].getName()
#print namelist[0].getArguments()


###################################################################################################
#Pose Energies


class PoseEnergies(object):

    #legt fest wie man einzelne instances erstellen kann. muss name und argumente enthalten.
    def __init__(self):
        self.namelist = []

    def loadFile(self, filename):
        f = open(filename)
        lines = f.readlines()
        for i in range(end):
            if i >= start:
                entry = lines[i].split()
                name_entry = entry[0]
                argument_entry = entry[1:]
                name_entry = ResidueEnergies(name_entry, argument_entry)
                #print name_entry.getArguments()
                self.namelist.append(name_entry)
#
#    def printNames(self):
#        for i in range(len(namelist)):
#            print namelist[i].getName()


	# gibt alle Argumente fuer bestimmte position im protein zurueck
    def getArguments(self, residue):
    	for i in range(len(namelist)):
    		if namelist[i].getName() == residue:
    			print argument_entry


    #gibt bestimmtes argument fuer bestimmten residue zurueck			
    def getArgument(self, residue, term):
    	for i in range(len(namelist)):
    		if namelist[i].getName() == residue:
    			print residue
    			for j in range(len(scoreterms)):
    				#print namelist[i].getArguments()[j]
    				if j == scoreterms.index(term):
    					print namelist[i].getArguments()[j]
    	


#Datei einlesen.
PoseEnergies().loadFile("testpdb.txt")

#einzelne werte abfragen
PoseEnergies().getArgument('ALA_100', 'total')
