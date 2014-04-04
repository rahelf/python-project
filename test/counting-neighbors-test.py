import sys
sys.path.append('../src')
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
from ResidueEnergies import ResidueEnergies, PoseEnergies,critical_distance_squared
from ResTypeAverageScores import ResTypeAverageScores
from constants import *

import cPickle
import numpy as np
import sys


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




###########################################################################
#extracting subgroup depending of number of neighbors

error5_seen = False


filename1 = '../test/10gs_nohet_1_relax.pdb'
filename2 = '../test/serin-only-test-pdb.txt'

statistics_collector_from_pdb1 = ResTypesStatisticsCollector()
statistics_collector_from_pdb2 = ResTypesStatisticsCollector()


pe_instance1 = PoseEnergies()
pe_instance1.loadFile(filename1)
statistics_collector_from_pdb1.add_pose_energies(pe_instance1)

pe_instance2 = PoseEnergies()
pe_instance2.loadFile(filename2)
statistics_collector_from_pdb2.add_pose_energies(pe_instance2)


print 'SER with 10-20 neighbors:', statistics_collector_from_pdb1.calculate_averages_and_stdevs_from_subset('SER', 'rama', range(10, 21))
print 'all SER from artificial file (contains all SER with 10-20 neighbors', statistics_collector_from_pdb2.calculate_averages_and_stdevs('SER', 'rama')

    