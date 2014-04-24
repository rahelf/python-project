import re
import numpy as np
import sys
from ResTypeAverageScores import ResTypeAverageScores
from CalphaCoordinates import CalphaCoordinates, critical_distance_squared
from constants import *


class ResidueEnergies(object):

    def __init__(self, line_string, score_term_list, pdb_identifier ):
        line_string_list = line_string.split()
        for item in line_string_list[1:]:
            item = float(item)
        self.res_type = line_string[0:3]
        self.res_num = re.search('(\d+)', line_string_list[0]).group()
        self.pdb_identifier = pdb_identifier
        self.score_dict = {}
        self.information_dict = {}
        st_counter = 1
        while st_counter < len(score_term_list):
            self.score_dict[score_term_list[st_counter]] = float(line_string_list[st_counter])
            self.information_dict[score_term_list[st_counter]] = (float(line_string_list[st_counter]), self.pdb_identifier, self.res_num)
            #print self.information_dict[score_term_list[st_counter]]
            st_counter += 1

        if len(line_string_list) != len(score_term_list):
            sys.exit('At least one residue has a different amount of score term entries')
        self.number_of_neighbors = 0

    def get_value( self, score_term):
        return float(self.score_dict.get(score_term))

    def get_value_for_combined_score_term_plus(self, list_of_scoreterms):
        return_val = 0.0
        for score_term in list_of_scoreterms:
            return_val = return_val + self.get_value(score_term)
        return return_val

    def get_value_for_combined_score_term_minus(self, minuend, subtrahend):
        return self.get_value(minuend) - self.get_value(subtrahend)

    def get_res_type(self):
       return self.res_type

    def get_nneighbors(self):
        return self.number_of_neighbors


class PoseEnergies(object):

    def __init__(self):
        self.res_e_list = []
        calpha_list = []
        self.pdb_identifier = 0

    def initalize_restype_av_scores( self, score_term_list, pdb_identifier):
        self.restype_av_scores = {}
        for aminoacid in aminoacids:
            self.restype_av_scores[aminoacid] = ResTypeAverageScores(aminoacid, score_term_list, pdb_identifier)
        
    def loadFile(self, filename):
        self.res_e_list = [] #safety
        calpha_list = []
        f = open(filename)
        lines = f.readlines()
        f.close()
        start_read_line = 0
        end_read_line = 0

        for i in range(len(lines)):

            # read atom positions and instantiate CalphaCoordinates
            if lines[i][0:4] == 'ATOM' and lines[i][13:15] == 'CA' and lines[i][17:20] in aminoacids:
                calpha_list.append( CalphaCoordinates (lines[i][17:20], int(lines[i][22:26]), float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54])))

            #find beginning and end of amino acid score terms
            elif re.match('#BEGIN_POSE_ENERGIES_TABLE', lines[i]):
                start_read_line = i+4
                self.pdb_identifier = lines[i].split()[1][0:4]
            elif re.match('#END_POSE_ENERGIES_TABLE', lines[i]):
                end_read_line = i


        #neighbor calculation
        for calpha_1 in range(len(calpha_list)):
            for calpha_2 in range(calpha_1 + 1, len(calpha_list)):
                calpha_list[calpha_1].are_neighbors(calpha_list[calpha_2], critical_distance_squared)

        #initialize ResTypeAverageScores
        current_line_no = start_read_line
        self.score_term_list = lines[start_read_line - 3].split()
        self.initalize_restype_av_scores(self.score_term_list, self.pdb_identifier)
        all_res_counter = 0
        while current_line_no < end_read_line:

            current_line = lines[current_line_no]
            current_line_no += 1
            #do not add non-amino acids to res_e_list
            if not current_line[0:3] in aminoacids:
                #self.res_e_list.pop()r
                #print "skipping non recognized res %s" % current_line[0:3]
                all_res_counter += 1
                continue

            if not (current_line[0:3] == calpha_list[all_res_counter].res_type ):
                sys.exit('ERROR in file %s: residue type of c alpha atoms does not match residue type of ResidueEnergies' %self.pdb_identifier)
            self.res_e_list.append( ResidueEnergies(current_line, self.score_term_list, self.pdb_identifier) )
            last_read_res = len( self.res_e_list) -1
            current_residue_type = self.res_e_list[last_read_res].get_res_type()
            self.res_e_list[last_read_res].number_of_neighbors = calpha_list[all_res_counter].neighbor_counter

            self.restype_av_scores[ current_residue_type ].add_residue_energies( self.res_e_list[ last_read_res ] )
            #print "just added a res of type %s from line %s" %(current_residue_type, current_line_no-1)
            all_res_counter += 1



    def calculate_averages_and_stddevs(self, res_type, score_term):
        cur_mean = self.restype_av_scores[res_type].get_mean_val(score_term)
        cur_stddev = self.restype_av_scores[res_type].get_stddev(score_term)
        return (cur_mean, cur_stddev)

    def get_rtas_for_restype( self, restype):
        return self.restype_av_scores[restype]