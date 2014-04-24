#!/usr/bin/python
import sys
sys.path.append('../../src')
from ResidueEnergies import ResidueEnergies, PoseEnergies,critical_distance_squared
from ResTypeAverageScores import ResTypeAverageScores
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
from constants import *

import numpy as np


statistics_collector = ResTypesStatisticsCollector()

pdbfile1 = "../stddev_and_mean/avetest_mock.pdb"
pdbfile2 = "../stddev_and_mean/avetest_mock2.pdb"

pe_instance_1 = PoseEnergies()
pe_instance_1.loadFile( pdbfile1 )
pe_instance_2 = PoseEnergies()
pe_instance_2.loadFile( pdbfile2)

statistics_collector.add_pose_energies( pe_instance_1)
statistics_collector.add_pose_energies( pe_instance_2)

#########################################
error_combined_scored = False



#combination of score terms
score_terms_to_be_combined = ['faketerm1', 'faketerm2']
for aminoacid in aminoacids:
    statistics_collector.restype_av_scores[aminoacid].calculate_sum_of_several_score_terms(score_terms_to_be_combined)

print 'statistics for ALA and fakterm1+faketerm2: ', statistics_collector.calculate_averages_and_stddevs('ALA', 'faketerm1+faketerm2')

if statistics_collector.calculate_averages_and_stddevs('ALA', 'faketerm1+faketerm2') != (2.25, 1.75):
    error_combined_scored = True



score_terms_to_be_combined2 = ['faketerm3', 'faketerm2']
statistics_collector.restype_av_scores['GLU'].calculate_sum_of_several_score_terms(score_terms_to_be_combined2)


#print 'statistics for GLU and fakterm3+faketerm2: ', statistics_collector.calculate_averages_and_stddevs('GLU', 'faketerm3+faketerm2')

#print "hack"
#print statistics_collector.calculate_averages_and_stddevs('GLU', 'faketerm3+faketerm2')


if round(statistics_collector.calculate_averages_and_stddevs('GLU', 'faketerm3+faketerm2')[0],5) != float(4.83333):
    error_combined_scored = True

if round(statistics_collector.calculate_averages_and_stddevs('GLU', 'faketerm3+faketerm2')[1],5) != float(5.20150):
    error_combined_scored = True
    


if error_combined_scored == True:
    print 'ERROR in combination scores'
else:
    print 'No error seen in combination of scores'




