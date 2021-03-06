#!/usr/bin/python

import cPickle
import sys
sys.path.append('../src')
from ResidueEnergies import ResidueEnergies, PoseEnergies
from ResTypeAverageScores import ResTypeAverageScores
import numpy as np
from ResTypesStatisticsCollector import ResTypesStatisticsCollector


aminoacids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']


pose_energies = PoseEnergies() # creates instance of PoseEnergies
pose_energies.loadFile('stddev_and_mean/avetest_mock.pdb')

statistics_collector = ResTypesStatisticsCollector()

pdbfile1 = "stddev_and_mean/avetest_mock.pdb"
pdbfile2 = "stddev_and_mean/avetest_mock2.pdb"

pe_instance_1 = PoseEnergies()
pe_instance_1.loadFile( pdbfile1 )
pe_instance_2 = PoseEnergies()
pe_instance_2.loadFile( pdbfile2)

statistics_collector.add_pose_energies( pe_instance_1)
statistics_collector.add_pose_energies( pe_instance_2)



########################################################################
# test: mean and stddev for single .pdb-files

error1_seen = False


if pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm1')[1] != np.std([1,1]) or pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm1')[0] != np.mean([1, 1]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm2')[1] != np.std([4, 2]) or pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm2')[0] != np.mean([4, 2]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm3')[1] != np.std([6, 5]) or pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm3')[0] != np.mean([6, 5]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm4')[1] != np.std([9, -10]) or pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm4')[0] != np.mean([9, -10]):
	error1_seen = True

elif pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm1')[1] != np.std([5]) or pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm1')[0] != np.mean([5]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm2')[1] != np.std([6]) or pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm2')[0] != np.mean([6]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm3')[1] != np.std([7]) or pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm3')[0] != np.mean([7]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm4')[1] != np.std([8]) or pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm4')[0] != np.mean([8]):
	error1_seen = True

elif pose_energies.calculate_averages_and_stddevs("MET", 'faketerm1')[1] != np.std([1, 1]) or pose_energies.calculate_averages_and_stddevs("MET", 'faketerm1')[0] != np.mean([1, 1]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("MET", 'faketerm2')[1] != np.std([2, -4]) or pose_energies.calculate_averages_and_stddevs("MET", 'faketerm2')[0] != np.mean([2, -4]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("MET", 'faketerm3')[1] != np.std([3, 6]) or pose_energies.calculate_averages_and_stddevs("MET", 'faketerm3')[0] != np.mean([3, 6]):
	error1_seen = True
elif pose_energies.calculate_averages_and_stddevs("MET", 'faketerm4')[1] != np.std([4, 2]) or pose_energies.calculate_averages_and_stddevs("MET", 'faketerm4')[0] != np.mean([4, 2]):
	error1_seen = True


if not error1_seen:
	print "No error when calculating mean and stddev for a single .pdb-file"
else:
	print "There was an error when calculating mean and standard deviation for a single .pdb-file"





####################################################################
#test: mean and stddev out of two .pdb-files

error2_seen = False


if statistics_collector.calculate_averages_and_stddevs("GLU", 'faketerm2')[1] != np.std([0, 0, 5, 6, 6]) or pose_energies.calculate_averages_and_stddevs("GLU", 'faketerm2')[0] != np.mean([0, 0, 5, 6, 6]):
	error1_seen = True
elif statistics_collector.calculate_averages_and_stddevs("MET", 'faketerm1')[1] != np.std([1, 4, 4]) or pose_energies.calculate_averages_and_stddevs("MET", 'faketerm1')[0] != np.mean([1, 4, 4]):
	error1_seen = True
elif statistics_collector.calculate_averages_and_stddevs("ASP", 'faketerm3')[1] != np.std([7, 6]) or pose_energies.calculate_averages_and_stddevs("ASP", 'faketerm3')[0] != np.mean([7, 6]):
	error1_seen = True
elif statistics_collector.calculate_averages_and_stddevs("TRP", 'faketerm4')[1] != 0 or pose_energies.calculate_averages_and_stddevs("TRP", 'faketerm4')[0] != 0:
	error1_seen = True

if not error2_seen:
	print "No error when calculating mean and stddev multiple .pdb-files"
else:
	print "There was an error when calculating mean and standard deviation for multiple .pdb-files"


######################################################################
#test: serialization and deserialization:

error3_seen = False

statistics_collector_from_pdb = ResTypesStatisticsCollector()
statistics_collector_from_archive = ResTypesStatisticsCollector()


all_file_names = ["stddev_and_mean/avetest_mock.pdb", "stddev_and_mean/avetest_mock2.pdb"]

for filename in all_file_names:
    pe_instance = PoseEnergies()
    pe_instance.loadFile(filename)
    statistics_collector_from_pdb.add_pose_energies(pe_instance)


for aminoacid in aminoacids:
    statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores('pickle-test/'+aminoacid+'.txt')


for aminoacid in aminoacids:
    filename = 'pickle-test/'+ aminoacid + '.txt'
    f = file(filename, 'rb')
    loaded_object = cPickle.load(f)
    statistics_collector_from_archive.add_archived_data( loaded_object)
    #archived_res_type_average_scores[aminoacid] = loaded_object
    f.close()

if not statistics_collector_from_pdb.calculate_averages_and_stddevs('GLU', 'faketerm1') == statistics_collector_from_archive.calculate_averages_and_stddevs('GLU', 'faketerm1'):
	error3_seen = True
elif not statistics_collector_from_pdb.calculate_averages_and_stddevs('ASP', 'faketerm1') == statistics_collector_from_archive.calculate_averages_and_stddevs('ASP', 'faketerm1'):
	error3_seen = True
elif not statistics_collector_from_pdb.calculate_averages_and_stddevs('MET', 'faketerm2') == statistics_collector_from_archive.calculate_averages_and_stddevs('MET', 'faketerm2'):
	error3_seen = True
#elif not statistics_collector_from_pdb.calculate_averages_and_stddevs('ALA', 'faketerm4') == statistics_collector_from_archive.calculate_averages_and_stddevs('ALA', 'faketerm4'):
#	error3_seen = True





if not error3_seen:
	print 'No error in serialization or deserialization!'
else:
	print 'There is an error in serialization or deserialization!'
