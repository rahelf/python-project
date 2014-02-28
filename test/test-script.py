#!/usr/bin/python


import sys
sys.path.append('../src')
from ResidueEnergies import ResidueEnergies, PoseEnergies
from ResTypeAverageScores import ResTypeAverageScores
import numpy as np


aminoacids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']


pose_energies = PoseEnergies() # creates instance of PoseEnergies
pose_energies.loadFile('avetest_mock.pdb')


'''
for aa in aminoacids:
    for score_term in pose_energies.score_term_list[1:]:
        mean = str(pose_energies.calculate_averages_and_stdevs(aa, score_term)[0])
        stdev = str(pose_energies.calculate_averages_and_stdevs(aa, score_term)[1])
        print score_term.ljust(25), aa, mean.ljust(20), stdev.ljust(20)
        '''


###### Test

error_seen = False


if pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm1')[1] != np.std([1, 1, 1]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm1')[0] != np.mean([1, 1, 1]):
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm2')[1] != np.std([4, -4, 2]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm2')[0] != np.mean([4, -4, 2]):
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm3')[1] != np.std([6, 6, 5]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm3')[0] != np.mean([6, 6, 5]):
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm4')[1] != np.std([9, -10, 2]) or pose_energies.calculate_averages_and_stdevs("GLU", 'faketerm4')[0] != np.mean([9, -10, 2]):
	error_seen = True

elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm1')[1] != np.std([5, 5]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm1')[0] != np.mean([5, 5]):
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm2')[1] != np.std([6, 3]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm2')[0] != np.mean([6, 3]):
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm3')[1] != np.std([7, 6]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm3')[0] != np.mean([7, 6]):
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm4')[1] != np.std([8, -8]) or pose_energies.calculate_averages_and_stdevs("ASP", 'faketerm4')[0] != np.mean([8, -8]):
	error_seen = True

elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm1')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm1')[0] != 1:
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm2')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm2')[0] != 2:
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm3')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm3')[0] != 3:
	error_seen = True
elif pose_energies.calculate_averages_and_stdevs("MET", 'faketerm4')[1] != 0 or pose_energies.calculate_averages_and_stdevs("MET", 'faketerm4')[0] != 4:
	error_seen = True


if not error_seen:
	print "No error seen"
else:
	print "There was an error when calculating mean and standard deviation"

