#!/usr/bin/python
import glob
from ResidueEnergies import ResidueEnergies, PoseEnergies, aminoacids
from ResTypeAverageScores import ResTypeAverageScores
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
import cPickle
import numpy as np
import sys

pdb_listfile = ""
archive_listfile = ""


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



#Serialization
for aminoacid in aminoacids:
    statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores('../pickled-files/'+aminoacid+'.txt')

#Deserialiazation
'''
for aminoacid in aminoacids:
    filename = '../pickled-files/'+aminoacid + '.txt'
    f = file(filename, 'rb')
    loaded_object = cPickle.load(f)
    statistics_collector_from_archive.add_archived_data( loaded_object)
    f.close()
    '''
#deserialize archives
for archive in archive_list:
    f = file(archive, 'rb')
    loaded_object = cPickle.load(f)
    statistics_collector_from_archive.add_archived_data( loaded_object )
    f.close()

#Histograms
statistics_collector_from_archive.restype_av_scores['GLU'].make_histogram_for_scoreterm('fa_rep')
statistics_collector_from_archive.restype_av_scores['GLU'].make_histogram_for_scoreterm('fa_atr')
statistics_collector_from_pdb.restype_av_scores['ASP'].make_histogram_for_scoreterm('rama')
