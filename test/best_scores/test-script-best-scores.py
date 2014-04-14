#test for best scores

import sys
sys.path.append('../../src')
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
from ResidueEnergies import ResidueEnergies, PoseEnergies,critical_distance_squared
from ResTypeAverageScores import ResTypeAverageScores
from constants import *
import cPickle
import numpy as np





FileList = ['12as_nohet_1_relax.pdb',  '12e8_nohet_1_relax.pdb', '12ca_nohet_1_relax.pdb',  '12gs_nohet_1_relax.pdb']
FileList_modified = ['12as_nohet_1_relax.pdb',  '12e8_nohet_1_relax.pdb', 'test_modified_12ca_nohet_1_relax.pdb',  '12gs_nohet_1_relax.pdb']



statistics_collector = ResTypesStatisticsCollector()
statistics_collector_mod = ResTypesStatisticsCollector()


for name in FileList:
	filename = name
	pe_instance = PoseEnergies()
	pe_instance.loadFile(filename)
	statistics_collector.add_pose_energies(pe_instance)

for name in FileList_modified:
	filename = name
	pe_instance = PoseEnergies()
	pe_instance.loadFile(filename)
	statistics_collector_mod.add_pose_energies(pe_instance)

#best score terms
#print statistics_collector.restype_av_scores['GLY'].get_best_score('fa_rep')
#print statistics_collector_mod.restype_av_scores['GLY'].get_best_score('fa_rep')

if float(statistics_collector_mod.restype_av_scores['GLY'].get_best_score('fa_rep')[0]) == -3.333333333:
	print 'No error seen, when trying to find the best score!'
else:
	print 'ERROR when trying to find the best score for a specific amino acid and score term!'

