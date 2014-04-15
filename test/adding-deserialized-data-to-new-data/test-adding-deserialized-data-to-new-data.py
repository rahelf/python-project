#test: adding deserialized data to statistics_collector_from_pdb
import sys
sys.path.append('../../src')
from ResidueEnergies import ResidueEnergies, PoseEnergies,critical_distance_squared
from ResTypeAverageScores import ResTypeAverageScores
from ResTypesStatisticsCollector import ResTypesStatisticsCollector
from constants import *

import cPickle
import numpy as np
import sys

pdb_listfile = 	'shortlist1.txt'

inlist = open(pdb_listfile, 'r')
liste = inlist.readlines()
FileList = []
for item in liste:
    FileList.append(item.rstrip('\n'))
inlist.close()

statistics_collector_from_pdb = ResTypesStatisticsCollector()

for filename in FileList:
    pe_instance = PoseEnergies()
    pe_instance.loadFile(filename)
    statistics_collector_from_pdb.add_pose_energies(pe_instance)

#Serialization
for aminoacid in aminoacids:
    statistics_collector_from_pdb.restype_av_scores[aminoacid].pickle_res_type_average_scores('archived/'+aminoacid+'.txt')

#print 'shortlist1:', statistics_collector_from_pdb.calculate_averages_and_stdevs('SER', 'rama')

####################################################################################


pdb_listfile =  'shortlist2.txt'
archive_listfile = 'test-archived-files.txt'

inlist = open(pdb_listfile, 'r')
liste = inlist.readlines()
FileList2 = []
for item in liste:
    FileList2.append(item.rstrip('\n'))
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

for filename in FileList2:
    pe_instance = PoseEnergies()
    pe_instance.loadFile(filename)
    statistics_collector_from_pdb.add_pose_energies(pe_instance)


#print 'shortlist2:',  statistics_collector_from_pdb.calculate_averages_and_stdevs('SER', 'rama')

#deserialize archives
for archive in archive_list:
    f = file(archive, 'rb')
    #print "trying to unpickle %s" % archive
    #statistics_collector_from_archive.add_archived_data( cPickle.load(f) )
    #print "succesful 1"
    #f.seek(0)
    statistics_collector_from_pdb.add_archived_data( cPickle.load(f) )
    #print "succesful 2"
    f.close()


print 'total from shortlist2 and archived shortlist1', statistics_collector_from_pdb.calculate_averages_and_stdevs('SER', 'rama')


#################################################################################################


pdb_listfile =  'list_total.txt'

inlist = open(pdb_listfile, 'r')
liste = inlist.readlines()
FileList = []
for item in liste:
    FileList.append(item.rstrip('\n'))
inlist.close()

statistics_collector_from_all_pdb = ResTypesStatisticsCollector()

for filename in FileList:
    pe_instance = PoseEnergies()
    pe_instance.loadFile(filename)
    statistics_collector_from_all_pdb.add_pose_energies(pe_instance)

print 'total all pdb', statistics_collector_from_all_pdb.calculate_averages_and_stdevs('SER', 'rama')


error_adding_data_seen = False

all_pdb_result = statistics_collector_from_all_pdb.calculate_averages_and_stdevs('SER', 'rama')
archive_result = statistics_collector_from_pdb.calculate_averages_and_stdevs('SER', 'rama')
if  round( all_pdb_result[0], 5 )!= round( archive_result[0], 5):
    error_adding_data_seen = True

if  round( all_pdb_result[1], 5 )!= round( archive_result[1], 5):
    error_adding_data_seen = True

if error_adding_data_seen:
    print 'Error in adding archived data to calculated data!'
else:
    print 'No error seen in adding archvied data to calculated data!'