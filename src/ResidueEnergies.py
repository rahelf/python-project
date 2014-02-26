import re
import numpy as np


class ResidueEnergies(object):

    def __init__(self, line_string, score_term_list ):
        line_string_list = line_string.split()
        for item in line_string_list[1:]:
            item = float(item)
        self.res_type = line_string[0:3]
        self.score_dict = dict(zip(score_term_list[1:], line_string_list[1:]))

    def get_value( self, score_term):
        return float(self.score_dict.get(score_term))

    def get_res_type(self):
        return self.res_type



class ResTypeAverageScores(object):

    def __init__( self, residue_type, score_term_list):
        initial_entries = [[]]* (len(score_term_list) -1)
        self.res_type_all_score_dict = dict( zip( score_term_list[1:], initial_entries))
        self.res_type = residue_type
        self.num_entries = 0

    def add_residue_energies(self,  ResidueEnergies_instance):
        if self.res_type != ResidueEnergies_instance.get_res_type():
            print "ERROR: score terms of ResidueEnergies_instance are added to wrong residue_type!"
            sys.exit()

        self.num_entries += 1
        for score_term in self.res_type_all_score_dict.keys():
            self.res_type_all_score_dict[score_term].append( ResidueEnergies_instance.get_value( score_term) )

    def get_mean_val(self, score_term):
        return np.mean( self.res_type_all_score_dict[score_term] )

    def get_stddev(self, score_term):
        return np.std(self.res_type_all_score_dict[score_term])




class PoseEnergies(object):

    def __init__(self):
        self.res_e_list = []

    def initalize_restype_av_scores( self, score_term_list):
        self.restype_av_scores = {}
        aminoacids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
        for aminoacid in aminoacids:
            self.restype_av_scores[aminoacid] = ResTypeAverageScores(aminoacid, score_term_list)
        


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
            elif re.match('#END_POSE_ENERGIES_TABLE', lines[i]):
                end_read_line = i
        current_line_no = start_read_line
        self.score_term_list = lines[start_read_line - 3].split()
        self.initalize_restype_av_scores(self.score_term_list)
        while current_line_no < end_read_line:
            current_line = lines[current_line_no]
            self.res_e_list.append( ResidueEnergies(current_line, self.score_term_list) )
            last_read_res = len( self.res_e_list) -1
            current_residue_type = self.res_e_list[last_read_res].get_res_type()
            self.restype_av_scores[ current_residue_type ].add_residue_energies( self.res_e_list[ last_read_res ] )
            current_line_no += 1
        

    def calculate_averages_and_stdevs(self, res_type, score_term):
        cur_mean = self.restype_av_scores[res_type].get_mean_val(score_term)
        cur_stddev = self.restype_av_scores[res_type].get_stddev(score_term)
        return (cur_mean, cur_stddev)

    #def get_score_term_value_for_residue( self, resnum, score_term): #resunm: position in sequence
        #return self.res_e_list[resnum - 1].get_value(score_term)





pose_energies = PoseEnergies() # creates instance of PoseEnergies
pose_energies.loadFile("../test/testpdb.txt")

#print pose_energies.res_e_list[0].get_value('fa_atr') 
#print pose_energies.get_score_term_value_for_residue(1, "fa_rep")
#print pose_energies.get_number_of_residues('ALA')
#print pose_energies.get_average('ALA', 'fa_atr')
#print pose_energies.get_standard_deviation('ALA', 'fa_atr')



Ala_average_stddev = pose_energies.calculate_averages_and_stdevs("ALA", 'fa_atr')
print Ala_average_stddev[0]
print Ala_average_stddev[1]
print Ala_average_stddev

