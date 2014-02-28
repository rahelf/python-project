#!/usr/bin/python

import sys
sys.path.append('../src')
from ResidueEnergies import ResidueEnergies, PoseEnergies
from ResTypeAverageScores import ResTypeAverageScores
import numpy as np
from ResTypesStatisticsCollector import ResTypesStatisticsCollector


aminoacids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']


pose_energies = PoseEnergies() # creates instance of PoseEnergies
pose_energies.loadFile('avetest_mock.pdb')

statistics_collector = ResTypesStatisticsCollector()

pdbfile1 = "../test/avetest_mock.pdb"
pdbfile2 = "../test/avetest_mock2.pdb"

pe_instance_1 = PoseEnergies()
pe_instance_1.loadFile( pdbfile1 )
pe_instance_2 = PoseEnergies()
pe_instance_2.loadFile( pdbfile2)

statistics_collector.add_pose_energies( pe_instance_1)
statistics_collector.add_pose_energies( pe_instance_2)



########################################################################
# test: mean and stdev for single .pdb-files

error1_seen = False


if pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm1')[1] != np.std([1, 1, 1]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm1')[0] != np.mean([1, 1, 1]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm2')[1] != np.std([4, -4, 2]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm2')[0] != np.mean([4, -4, 2]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm3')[1] != np.std([6, 6, 5]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm3')[0] != np.mean([6, 6, 5]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm4')[1] != np.std([9, -10, 2]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm4')[0] != np.mean([9, -10, 2]):
	error1_seen = True

elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm1')[1] != np.std([5, 5]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm1')[0] != np.mean([5, 5]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm2')[1] != np.std([6, 3]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm2')[0] != np.mean([6, 3]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm3')[1] != np.std([7, 6]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm3')[0] != np.mean([7, 6]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm4')[1] != np.std([8, -8]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm4')[0] != np.mean([8, -8]):
	error1_seen = True

elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm1')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm1')[0] != 1:
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm2')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm2')[0] != 2:
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm3')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm3')[0] != 3:
	error1_seen = True
elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm4')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm4')[0] != 4:
	error1_seen = True


if not error1_seen:
	print "No error when calculating mean and stdev for a single .pdb-file"
else:
	print "There was an error when calculating mean and standard deviation for a single .pdb-file"





####################################################################
#test: mean and stdev out of two .pdb-files

error2_seen = False


if statistics_collector.calculate_averages_and_stdevs("GLU", 'faketerm2')[1] != np.std([0, 0, 5, 6, 6]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm2')[0] != np.mean([0, 0, 5, 6, 6]):
	error1_seen = True
elif statistics_collector.calculate_averages_and_stdevs("MET", 'faketerm1')[1] != np.std([1, 4, 4]) or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm1')[0] != np.mean([1, 4, 4]):
	error1_seen = True
elif statistics_collector.calculate_averages_and_stdevs("ASP", 'faketerm3')[1] != np.std([7, 6]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm3')[0] != np.mean([7, 6]):
	error1_seen = True
elif statistics_collector.calculate_averages_and_stdevs("TRP", 'faketerm4')[1] != 0 or pose_energies.calculate_averages_and_stdevs("TRP", 'faketerm4')[0] != 0:
	error1_seen = True

if not error2_seen:
	print "No error when calculating mean and stdev multiple .pdb-files"
else:
	print "There was an error when calculating mean and standard deviation for multiple .pdb-files"




