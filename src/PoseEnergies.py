#Pose Energies
from ResidueEnergies import ResidueEnergies

class PoseEnergies(ResidueEnergies):

    #legt fest wie man einzelne instances erstellen kann. muss name und argumente enthalten.
    def __init__(self):
        self.namelist = []

    def loadfile(self, filename):
        f = open(filename)
        lines = f.readlines()
        namelist = []
        for i in range(5291):
            if i >= 4984:
                entry = lines[i].split()
                name_entry = entry[0]
                print name_entry
                argument_entry = entry[1:]
                name_entry = ResidueEnergies(name_entry, argument_entry)
                #print name_entry.getArguments()
                namelist.append(name_entry)

    def print_all_names(self):
        for i in len(namelist):
            # GetName ist funktion der Klasse ResidueEnergies
            print namelist[i].GetName()

x= PoseEnergies()
x.loadfile("testpdb.txt")
x.print_all_names()