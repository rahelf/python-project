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




statistics_collector_from_pdb = ResTypesStatisticsCollector()
statistics_collector_from_archive = ResTypesStatisticsCollector()

for filename in FileList:
    filename = '../../pdbdir/'+filename
    pe_instance = PoseEnergies()
    pe_instance.loadFile(filename)
    statistics_collector_from_pdb.add_pose_energies(pe_instance)

#print 'GLU fa_atr:\n'
#print 'SER with 10-20 neighbors:', statistics_collector_from_pdb.calculate_averages_and_stdevs_from_subset('SER', 'rama', range(10, 21))
#print 'GLU with any number of neighbors:', statistics_collector_from_pdb.calculate_averages_and_stdevs_from_subset('GLU', 'fa_atr', number_of_neighbors_list)
print 'SER with any number of neighbors (other method):', statistics_collector_from_pdb.calculate_averages_and_stdevs('SER', 'rama')



#Serialization
for aminoacid in aminoacids:
    statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores('../pickled-files/'+aminoacid+'.txt')

#deserialize archives
for archive in archive_list:
    f = file(archive, 'rb')
    statistics_collector_from_archive.add_archived_data( cPickle.load(f) )
    f.close()

#Histograms
#statistics_collector_from_archive.restype_av_scores['GLU'].make_histogram_for_scoreterm('fa_rep')
#statistics_collector_from_archive.restype_av_scores['GLU'].make_histogram_for_scoreterm_for_ncounts('fa_atr', range(0,41))
#statistics_collector_from_pdb.restype_av_scores['GLU'].make_histogram_for_scoreterm_for_ncounts('fa_atr', range(0,41))

#statistics_collector_from_pdb.restype_av_scores['GLU'].make_histogram_for_scoreterm_for_ncounts('fa_atr', range(0,12))
#statistics_collector_from_pdb.restype_av_scores['GLU'].make_histogram_for_scoreterm_for_ncounts('fa_atr', range(20,41))

