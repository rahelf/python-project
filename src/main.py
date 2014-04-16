#!/usr/bin/python
import glob
from ResidueEnergies import ResidueEnergies, PoseEnergies,critical_distance_squared
from ResTypeAverageScores import ResTypeAverageScores
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
from constants import *

import cPickle
import numpy as np
import sys

pdb_listfile = ""
archive_listfile = ""



#Command Line Arguments
CommandArgs = sys.argv[1:]
for arg in CommandArgs:
    if arg == '-pdbfiles':
        pdb_listfile = CommandArgs[CommandArgs.index(arg)+1]
        listmode = 1
    elif arg == '-archived':
        archive_listfile = CommandArgs[CommandArgs.index(arg)+1]

if not pdb_listfile:
     sys.exit('Error, please supply name of listfile')



#read files
inlist = open(pdb_listfile, 'r')
liste = inlist.readlines()
FileList = []
for item in liste:
    FileList.append(item.rstrip('\n'))
inlist.close()

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


#initialize PoseEnergies for eacht file in list
for filename in FileList:
    #filename = '../../pdbdir/'+filename
    pe_instance = PoseEnergies()
    try:
        pe_instance.loadFile(filename)
    except:
        print "Caught exception when trying to read %s" % filename
        continue

    try:
        statistics_collector_from_pdb.add_pose_energies(pe_instance)
    except:
        print "Caught exception when trying to add values from  %s to statistics calculator" % filename
        continue

#combination of score terms
#score_terms_to_be_combined = ['rama', 'fa_atr']
#aminoacid = 'TYR'
#statistics_collector_from_pdb.restype_av_scores[aminoacid].calculate_sum_of_several_score_terms(score_terms_to_be_combined)
#print'statistics for %s and the combined score term %s: ' %(aminoacid, score_terms_to_be_combined), statistics_collector_from_pdb.calculate_averages_and_stdevs(aminoacid, 'rama+fa_atr')

#number_of_entbries = 0
#for nn in number_of_neighbors_list:
    #print nn, ':', len(statistics_collector_from_pdb.restype_av_scores[aminoacid].res_type_all_score_dict['fa_atr'][nn])
    #number_of_entries += len(statistics_collector_from_pdb.restype_av_scores[aminoacid].res_type_all_score_dict['fa_atr'][nn])
#print 'total number of entries for the aminoacid: ', number_of_entries

#print 'statistics for first scoreterm only: ', statistics_collector_from_pdb.calculate_averages_and_stdevs(aminoacid, 'rama')
#print 'statistics for second scoreterm only:', statistics_collector_from_pdb.calculate_averages_and_stdevs(aminoacid, 'fa_atr')



#average and standard deviation
#print statistics_collector_from_pdb.calculate_averages_and_stdevs('SER', 'rama' )



#Serialization
for aminoacid in aminoacids:
    statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores('../pickled-files/'+aminoacid+'.txt')


#deserialize archives
#for archive in archive_list:
#    f = file(archive, 'rb')
#    statistics_collector_from_pdb.add_archived_data( cPickle.load(f) )
#    f.close()


#best score terms
#aminoacid = 'TRP'
#score_term = 'fa_rep'
#print 'Best score is %s. \npdb-identfier of file: %s \nresidue number: %s' %statistics_collector_from_pdb.restype_av_scores[aminoacid].get_best_score(score_term)

#Histograms
#statistics_collector_from_archive.restype_av_scores[aminoacid].make_histogram_for_scoreterm_for_ncounts(score_term, range(0,41))
#statistics_collector_from_pdb.restype_av_scores[aminoacid].make_histogram_for_scoreterm_for_ncounts(score_term, range(0,41))
