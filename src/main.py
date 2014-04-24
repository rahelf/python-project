#!/usr/bin/python
import glob
from ResidueEnergies import ResidueEnergies, PoseEnergies,critical_distance_squared
from ResTypeAverageScores import ResTypeAverageScores
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
from ZScoreCalculator import ZScoreCalculator
from constants import *

import cPickle
import numpy as np
import sys

pdb_listfile = ""
archive_listfile = ""
pdb_file = ''
pdb_file_2 = ''
score_term_z = ''
score_term_minus_z = ''
score_term_plus_z = ''
histogram_location = ''
goodz = ''
badz = ''



#Command Line Arguments
CommandArgs = sys.argv[1:]
for arg in CommandArgs:
    if arg == '-in':
        pdb_listfile = CommandArgs[CommandArgs.index(arg)+1]
        listmode = 1
    elif arg == '-arch':
        archive_listfile = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-pickle':
        pickle_location = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-hist':
        histogram_location = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-z1':
        pdb_file = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-z2':
        pdb_file_2 = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-st':
        score_term_z = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-st+':
        score_term_plus_z = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-st-':
        score_term_minus_z = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-goodz':
        goodz = float(CommandArgs[CommandArgs.index(arg)+1])
    elif arg == '-badz':
        badz = float(CommandArgs[CommandArgs.index(arg)+1])


#data needed (either from archive or calculated from files)
if not (pdb_listfile or archive_listfile):
     sys.exit('Error, please supply name of listfile')



#read given file names and add them to a list (FileList)
if pdb_listfile !="":
    inlist = open(pdb_listfile, 'r')
    liste = inlist.readlines()
    FileList = []
    for item in liste:
        FileList.append(item.rstrip('\n'))
    inlist.close()
else:
    FileList = False

# archive_listfile contains a list of archived files (typically one for each aminoacid).
archive_list = []
if archive_listfile !="":
    archhandle = open(archive_listfile, 'r')
    flines = archhandle.readlines()
    archhandle.close()
    for line in flines:
        archive_list.append( line.rstrip('\n'))


# initialize ResTypesStatisticsCollector
statistics_collector_from_pdb = ResTypesStatisticsCollector()
statistics_collector_from_archive = ResTypesStatisticsCollector()


#initialize PoseEnergies for each file in list
if FileList:
    for filename in FileList:
        pe_instance = PoseEnergies()
        try:
            pe_instance.loadFile(filename)
        except:
            print "Caught exception when trying to read %s" % filename
            continue

        try:
            statistics_collector_from_pdb.add_pose_energies(pe_instance)
        except:
            print "Caught exception when trying to add values from  %s to statistics collector" % filename
            continue


#Serialization
#for aminoacid in aminoacids:
#    statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores(pickle_location+aminoacid+'.txt')


#deserialize archived files
for archive in archive_list:
    f = file(archive, 'rb')
    statistics_collector_from_archive.add_archived_data( cPickle.load(f) )
    #print archive
    f.close()


#Deserialization of one particular file (faster)
#aminoacids = ['GLU']
# = '/home/rahel/uni/wise1314/python/pdbstats/pickle11k/GLU.txt'
#statistics_collector_from_archive.add_archived_data(cPickle.load(file(archive, 'rb')))



#combination of score terms (plus)
score_terms_to_be_combined = ['fa_atr', 'fa_rep']
for aminoacid in aminoacids:
    statistics_collector_from_archive.restype_av_scores[aminoacid].calculate_sum_of_several_score_terms(score_terms_to_be_combined)
score_terms_to_be_combined = ['hbond_bb_sc', 'hbond_sc']
for aminoacid in aminoacids:
    statistics_collector_from_archive.restype_av_scores[aminoacid].calculate_sum_of_several_score_terms(score_terms_to_be_combined)



#subtraction of score_terms ----  important for subtraction of fa_dun from total
minuend = 'total'
subtrahend = 'fa_dun'
for aminoacid in aminoacids:
    statistics_collector_from_archive.restype_av_scores[aminoacid].subtract_score_terms(minuend, subtrahend)


#average and standard deviation
#print statistics_collector_from_archive.get_mean_and_stddev('SER', 'hbond_bb_sc+hbond_sc', 5)
#print statistics_collector_from_archive.get_mean_and_stddev('SER', 'total-fa_dun', 9)


#best score terms
#print '\nfrom archive: Best score is %s. \npdb-identfier of file: %s \nresidue number: %s' %statistics_collector_from_archive.restype_av_scores[aminoacid].get_best_score(score_term)


'''
interesting_score_terms = ['total', 'fa_atr', 'hbond_bb_sc']
for score_term in interesting_score_terms:
    print '\n', score_term
    for score_term in interesting_score_terms:
        print '\n', score_term
        print '%s best score: %s' % (aminoacid, statistics_collector_from_archive.restype_av_scores[aminoacid].get_best_score(score_term))
'''

#frequency of aminoacids and neighbor numbers
#statistics_collector_from_archive.restype_av_scores['GLU'].plot_relative_frequencies_of_numbers_of_neighbors()
#statistics_collector_from_archive.restype_av_scores['GLU'].get_relative_frequency_for_nn(10)
#for aminoacid in aminoacids:
#    statistics_collector_from_archive.restype_av_scores[aminoacid].plot_relative_frequencies_of_numbers_of_neighbors()


aminoacids = ['TYR', 'VAL']
interesting_score_terms = ['fa_dun']

if histogram_location != '' and interesting_score_terms != '':
    for aminoacid in aminoacids:
        if aminoacid in outside_aa:
            main_neighbor_situations = main_neighbor_situations_outside
        elif aminoacid in inside_aa:
            main_neighbor_situations = main_neighbor_situations_inside
        for score_term in interesting_score_terms:
            for neighbor_situation in main_neighbor_situations:
                mean = statistics_collector_from_archive.calculate_averages_and_stddevs_from_subset(aminoacid, score_term, neighbor_situation)[0]
                stddev = statistics_collector_from_archive.calculate_averages_and_stddevs_from_subset(aminoacid, score_term, neighbor_situation)[1]
                number_of_residues = len(statistics_collector_from_archive.restype_av_scores[aminoacid].get_merged_list_for_ncounts(score_term, neighbor_situation))
                statistics_collector_from_archive.restype_av_scores[aminoacid].make_histogram_for_scoreterm_for_ncounts(score_term, neighbor_situation, number_of_residues, histogram_location, mean, stddev)


#Z scores (filenames)
#pdb_file = '/home/rahel/uni/wise1314/python/vergl_designs/3b4x_talcstrlx_0001.pdb'
#pdb_file_2 = '/home/rahel/uni/wise1314/python/vergl_designs/dCM13_talcstrlx_0001.pdb'

score_termsz = ['total', 'fa_atr']

for score_term_z in score_termsz:
    print '\nCalculation of z scores for %s' % score_term_z
    if pdb_file != '' and pdb_file_2 == '':
        print 'z-scores:'
        reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
        zscores = reference_z_scores.calculate_z_scores(score_term_z)
        for key in zscores.keys():
            print '\n%s: %s    %s' %(key, zscores[key][1], zscores[key][0])
    elif pdb_file!= '' and pdb_file_2 != '':
        print 'difference in z scores'
        reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
        z_scores_after_modification = ZScoreCalculator(pdb_file_2, statistics_collector_from_archive)
        delta_z = reference_z_scores.calculate_differences_in_z_scores(z_scores_after_modification, score_term_z)
        for i in range(len(delta_z)):
            print '%s   %s   %s' %(delta_z[i][0], delta_z[i][1], delta_z[i][2])
        print 'position with worst shift: %s (%s) \nposition with best shift: %s (%s)' %(z_scores_after_modification.position_positive_shift, round(z_scores_after_modification.value_pos,4), z_scores_after_modification.position_negative_shift, round(z_scores_after_modification.value_neg,4))
        if goodz != '':
            z_scores_after_modification.get_goodz(goodz)
        if badz != '':
            z_scores_after_modification.get_badz(badz)

score_termsz_plus = ['fa_atr+fa_rep', 'hbond_bb_sc+hbond_sc']
for score_term_plus_z in score_termsz_plus:
    print '\nCalculation of z scores for %s' % score_term_plus_z
    #scoreterms trennen am plus
    score_terms = score_term_plus_z.split("+", 1)
    print score_terms
    for aminoacid in aminoacids:
        statistics_collector_from_archive.restype_av_scores[aminoacid].calculate_sum_of_several_score_terms(score_terms)
    if pdb_file != '' and pdb_file_2 == '':
        print 'z-scores:'
        reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
        zscores = reference_z_scores.calculate_z_scores_combined_plus(score_terms)
        for key in zscores.keys():
            print '\n%s: %s    %s' %(key, zscores[key][1], zscores[key][0])
    elif pdb_file!= '' and pdb_file_2 != '':
        print 'difference in z scores'
        reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
        z_scores_after_modification = ZScoreCalculator(pdb_file_2, statistics_collector_from_archive)
        delta_z = reference_z_scores.calculate_differences_in_z_scores_combined_plus(z_scores_after_modification, score_terms)
        for i in range(len(delta_z)):
            print '%s   %s   %s' %(delta_z[i][0], delta_z[i][1], delta_z[i][2])
        print 'position with worst shift: %s (%s) \nposition with best shift: %s (%s)' %(z_scores_after_modification.position_positive_shift, round(z_scores_after_modification.value_pos,4), z_scores_after_modification.position_negative_shift, round(z_scores_after_modification.value_neg,4))
        if goodz != '':
            z_scores_after_modification.get_goodz(goodz)
        if badz != '':
            z_scores_after_modification.get_badz(badz)

score_termsz_minus = ['total-fa_dun']
for score_term_minus_z in score_termsz_minus:
    print '\nCalculation of z scores for %s' % score_term_minus_z
    #scoreterms trennen am minus
    score_terms = score_term_minus_z.split("-", 1)
    minuend = score_terms[0]
    subtrahend = score_terms[1]
    for aminoacid in aminoacids:
        statistics_collector_from_archive.restype_av_scores[aminoacid].subtract_score_terms(minuend, subtrahend)
    if pdb_file != '' and pdb_file_2 == '':
        print 'z-scores:'
        reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
        zscores = reference_z_scores.calculate_z_scores_combined_minus(minuend, subtrahend)
        for key in zscores.keys():
            print '\n%s: %s    %s' %(key, zscores[key][1], zscores[key][0])
    elif pdb_file!= '' and pdb_file_2 != '':
        print 'difference in z scores'
        reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
        z_scores_after_modification = ZScoreCalculator(pdb_file_2, statistics_collector_from_archive)
        delta_z = reference_z_scores.calculate_differences_in_z_scores_combined_minus(z_scores_after_modification, minuend, subtrahend)
        for i in range(len(delta_z)):
            print '%s   %s   %s' %(delta_z[i][0], delta_z[i][1], delta_z[i][2])
        print 'position with worst shift: %s (%s) \nposition with best shift: %s (%s)' %(z_scores_after_modification.position_positive_shift, round(z_scores_after_modification.value_pos,4), z_scores_after_modification.position_negative_shift, round(z_scores_after_modification.value_neg,4))
        if goodz != '':
            z_scores_after_modification.get_goodz(goodz)
        if badz != '':
            z_scores_after_modification.get_badz(badz)


