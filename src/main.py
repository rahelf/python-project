#!/usr/bin/python
from ResidueEnergies import ResidueEnergies, PoseEnergies
from ResTypeAverageScores import ResTypeAverageScores
from ResTypesStatisticsCollector import ResTypesStatisticsCollector


#pose_energies = PoseEnergies() # creates instance of PoseEnergies
#pose_energies.loadFile('../test/testpdb.txt')

statistics_collector = ResTypesStatisticsCollector()

pdbfile1 = "../test/avetest_mock.pdb"
pdbfile2 = "../test/avetest_mock2.pdb"

pe_instance_1 = PoseEnergies()
pe_instance_1.loadFile( pdbfile1 )
print pe_instance_1.score_term_list[1:]
#print pe_instance_1.calculate_averages_and_stdevs('ALA', 'faketerm1')


pe_instance_2 = PoseEnergies()
pe_instance_2.loadFile( pdbfile2)

statistics_collector.add_pose_energies( pe_instance_1)
statistics_collector.add_pose_energies( pe_instance_2)

#print pose_energies.res_e_list[0].get_value('fa_atr') 
#print pose_energies.get_score_term_value_for_residue(1, "fa_rep")
#print pose_energies.get_number_of_residues('ALA')
#print pose_energies.get_average('ALA', 'fa_atr')
#print pose_energies.get_standard_deviation('ALA', 'fa_atr')


'''
for aa in aminoacids:
    for score_term in pose_energies.score_term_list[1:]:
        mean = str(pose_energies.calculate_averages_and_stdevs(aa, score_term)[0])
        stdev = str(pose_energies.calculate_averages_and_stdevs(aa, score_term)[1])
        print score_term.ljust(25), aa, mean.ljust(20), stdev.ljust(20)


print pose_energies.calculate_averages_and_stdevs("ALA", 'fa_atr')
print pose_energies.calculate_averages_and_stdevs("GLU", 'fa_sol')
























class Atest(object):
	def __init__(self, aval, bval):
		self.aval = aval
		self.bval = bval

	def get_aval():
		return self.aval

	def get_bval():
		return self.bval

	def get_sum():
		return self.aval + self.bval




	first_instance = Atest(3,5)

	print first_instance.get_sum()

	second_instance = Atest(4,9)
'''


















'''

FileList = []
Listfile = ''
template = ''
outfile = ""
Singlefile = ''
listmode = 0
external_template = 1

CommandArgs = sys.argv[1:]

for arg in CommandArgs:
    if arg == '-t':
        template = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-l':
        Listfile = CommandArgs[CommandArgs.index(arg)+1]
        listmode = 1
    elif arg == '-out':
        outfile = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-s':
        Singlefile = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-add_pdb_suffix':
        ADD_PDB_SUFFIX = 1

if( Listfile and not template ):
    hurz = 2

elif( (not Listfile and not Singlefile) or (not template) ):
     print 'Error, please supply name of the listfile, the resfile and the template'
     sys.exit()


if(listmode):
    inlist = open(Listfile,'r')
    FileList = inlist.readlines()
    inlist.close()
    if( not template ):
        template = FileList[0].replace("\n","")
        external_template = 0
    if outfile == ' ':
        outfile = Listfile + '.ana'
    print "Checking structures for %s structures in %s to template %s" % (len(FileList), Listfile, template)

else:
    FileList.append(Singlefile)

template_coords = {}

if template != '':
    template_coords = get_coordinates(template)

seq_prof = SequenceProfile( template_coords, external_template )

outstring = ""

for struct in FileList:

    struct_coords = get_coordinates(struct)

    seq_prof.add_struct( struct_coords )
 
    #print mutstring

outstring = seq_prof.get_outstring()
pymutstring =  seq_prof.get_pymutstring()
    #print outstring


if outfile == "":
    print outstring
    print pymutstring

else:
    outf = open(outfile,'w')
    outf.write(outstring)
    outf.close()

'''
