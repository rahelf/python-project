import sys
sys.path.append('../../src')
from ZScoreCalculator import ZScoreCalculator
from constants import *
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
import cPickle





archive_listfile = '../../src/archive_files.txt'
archive_list = []
if archive_listfile !="":
    archhandle = open(archive_listfile, 'r')
    flines = archhandle.readlines()
    archhandle.close()
    for line in flines:
        archive_list.append( line.rstrip('\n'))



statistics_collector_from_archive = ResTypesStatisticsCollector()

#deserialize archives
for archive in archive_list:
    f = file(archive, 'rb')
    statistics_collector_from_archive.add_archived_data( cPickle.load(f) )
    f.close()


#combination of score terms
score_terms_to_be_combined = ['hbond_bb_sc', 'hbond_sc']
for aminoacid in aminoacids:
    statistics_collector_from_archive.restype_av_scores[aminoacid].calculate_sum_of_several_score_terms(score_terms_to_be_combined)


#subtraction of score_terms ----  important for subtraction of fa_dun from total
minuend = 'total'
subtrahend = 'fa_dun'
for aminoacid in aminoacids:
    statistics_collector_from_archive.restype_av_scores[aminoacid].subtract_score_terms(minuend, subtrahend)


error_z = False


pdb_file = '11ba_nohet_1_relax_mod.pdb'
pdb_file_2 = '11ba_nohet_1_relax_mod_for_z_score_calculator.pdb'

reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
z_scores_after_modification = ZScoreCalculator(pdb_file_2, statistics_collector_from_archive)



mutations = 0
for line in reference_z_scores.calculate_differences_in_z_scores(z_scores_after_modification, 'total'):
    print line
    if line[2] == 'mutation':
        mutations += 1
if mutations != 2:
    error_z = True

#print 'total number of mutations:%s' %mutations

if z_scores_after_modification.position_positive_shift != str(57) or z_scores_after_modification.position_negative_shift != str(23):
    error_z = True

if error_z:
    print 'Error in z score calculation!'
else:
    print 'No error in z score calculuation seen'

 

#reference_z_scores.calculate_differences_in_z_scores_combined_plus(z_scores_after_modification, ['hbond_bb_sc', 'hbond_sc'])
reference_z_scores.calculate_differences_in_z_scores_combined_minus(z_scores_after_modification, 'total', 'fa_dun')
 
