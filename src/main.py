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
add_arch = False
pdb_file = ''
pdb_file_2 = ''
interesting_score_terms = ''
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
    elif arg == '-a':
        add_arch = True
    elif arg == '-pickle':
        pickle_location = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-hist':
        histogram_location = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-z1':
        pdb_file = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-z2':
        pdb_file_2 = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-st':
        interesting_score_terms = CommandArgs[CommandArgs.index(arg)+1].split(",")
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
if pickle_location != '':
    for aminoacid in aminoacids:
        statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores(pickle_location+aminoacid+'.txt')


#deserialize archived files
if archive_listfile != '' and pdb_listfile == '':
    if add_arch == True:
        for archive in archive_list:
            f = file(archive, 'rb')
            statistics_collector_from_pdb.add_archived_data( cPickle.load(f) )
            f.close()
            statistics_collector_from_archive = statistics_collector_from_pdb
    elif add_arch ==False:
        for archive in archive_list:
            f = file(archive, 'rb')
            statistics_collector_from_archive.add_archived_data( cPickle.load(f) )
            f.close()



#Deserialization of one particular file (faster)
#aminoacids = ['GLU']
# = '/home/rahel/uni/wise1314/python/pdbstats/pickle11k/GLU.txt'
#statistics_collector_from_archive.add_archived_data(cPickle.load(file(archive, 'rb')))

#average and standard deviation
#print statistics_collector_from_archive.get_mean_and_stddev('SER', 'hbond_bb_sc+hbond_sc', 5)
#print statistics_collector_from_archive.get_mean_and_stddev('SER', 'total-fa_dun', 9)


#best score terms
#print '\nfrom archive: Best score is %s. \npdb-identfier of file: %s \nresidue number: %s' %statistics_collector_from_archive.restype_av_scores[aminoacid].get_best_score(score_term)
#for score_term in interesting_score_terms:
 #   print '\n', score_term
  #  for score_term in interesting_score_terms:
   #     print '\n', score_term
    #    print '%s best score: %s' % (aminoacid, statistics_collector_from_archive.restype_av_scores[aminoacid].get_best_score(score_term))


#frequency of aminoacids and neighbor numbers
#statistics_collector_from_archive.restype_av_scores['GLU'].plot_relative_frequencies_of_numbers_of_neighbors()
#statistics_collector_from_archive.restype_av_scores['GLU'].get_relative_frequency_for_nn(10)
#for aminoacid in aminoacids:
#    statistics_collector_from_archive.restype_av_scores[aminoacid].plot_relative_frequencies_of_numbers_of_neighbors()


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

if not (pdb_file == '' and pdb_file_2 == ''):
    for score_term_z in interesting_score_terms:
        print '\n---------Calculation of z scores for %s-------------------------' % score_term_z
        if '+' in score_term_z or '-' in score_term_z:
            statistics_collector_from_archive.calculate_combined_score_terms(score_term_z)
        if pdb_file != '' and pdb_file_2 == '':
            print 'z-scores:'
            reference_z_scores = ZScoreCalculator(pdb_file, statistics_collector_from_archive)
            zscores = reference_z_scores.calculate_z_scores(score_term_z)
            for i in reference_z_scores.z_sorted:
                print i
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
