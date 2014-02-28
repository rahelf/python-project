import re
import numpy as np
from ResTypeAverageScores import ResTypeAverageScores

aminoacids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']


class ResidueEnergies(object):

    def __init__(self, line_string, score_term_list ):
        line_string_list = line_string.split()
        for item in line_string_list[1:]:
            item = float(item)
        self.res_type = line_string[0:3]
        self.score_dict = {}
        st_counter = 1
        while st_counter < len(score_term_list):
            self.score_dict[score_term_list[st_counter]] = line_string_list[st_counter]
            st_counter += 1

    def get_value( self, score_term):
        return float(self.score_dict.get(score_term))

    def get_res_type(self):
       return self.res_type


class PoseEnergies(object):

    def __init__(self):
        self.res_e_list = []
        self.pdb_identifier = 0


    def initalize_restype_av_scores( self, score_term_list, pdb_identifier):
        self.restype_av_scores = {}
        for aminoacid in aminoacids:
            self.restype_av_scores[aminoacid] = ResTypeAverageScores(aminoacid, score_term_list, pdb_identifier)
            #self.restype_av_scores[aminoacid] = ResTypeAverageScores.from_pdb_file(aminoacid, score_term_list, pdb_identifier)
        #print self.restype_av_scores["ALA"]
        
    def loadFile(self, filename):
        self.res_e_list = [] #safety
        f = open(filename)
        lines = f.readlines()
        f.close()
        start_read_line = 0
        end_read_line = 0
        for i in range(len(lines)):
            if re.match('#BEGIN_POSE_ENERGIES_TABLE', lines[i]):
                start_read_line = i+4
                self.pdb_identifier = lines[i].split()[1][0:4]
            elif re.match('#END_POSE_ENERGIES_TABLE', lines[i]):
                end_read_line = i
        current_line_no = start_read_line
        self.score_term_list = lines[start_read_line - 3].split()
        self.initalize_restype_av_scores(self.score_term_list, self.pdb_identifier)
        while current_line_no < end_read_line:

            current_line = lines[current_line_no]
            self.res_e_list.append( ResidueEnergies(current_line, self.score_term_list) )
            last_read_res = len( self.res_e_list) -1
            current_residue_type = self.res_e_list[last_read_res].get_res_type()
            if current_residue_type == "VRT":
                self.res_e_list.pop()
                continue
            self.restype_av_scores[ current_residue_type ].add_residue_energies( self.res_e_list[ last_read_res ] )
            current_line_no += 1

    def calculate_averages_and_stdevs(self, res_type, score_term):
        cur_mean = self.restype_av_scores[res_type].get_mean_val(score_term)
        cur_stddev = self.restype_av_scores[res_type].get_stddev(score_term)
        return (cur_mean, cur_stddev)

    def get_rtas_for_restype( self, restype):
        return self.restype_av_scores[restype]
