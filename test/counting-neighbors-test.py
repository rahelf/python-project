import sys
sys.path.append('../src')
from ResidueEnergies import ResidueEnergies, PoseEnergies





#############################################################################
#calculate number of neighbors

error4_seen = False

filename = '../../pdbdir/10mh_nohet_1_relax.pdb'
pe_instance = PoseEnergies()
pe_instance.loadFile(filename)


if pe_instance.res_e_list[209].res_type != 'LYS' or pe_instance.res_e_list[209].number_of_neighbors != 8:
	error4_seen = True
elif pe_instance.res_e_list[117].res_type != 'MET' or pe_instance.res_e_list[117].number_of_neighbors != 16:
	error4_seen = True


if not error4_seen:
	print 'No error when counting neighbors'
else:
	print 'ERROR: Error in counting neighbors'

