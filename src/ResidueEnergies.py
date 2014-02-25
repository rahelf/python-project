import re


###############################################################################
#Residue Energies: jede Zeile der PDB-Datei entspricht einer Instanz

class ResidueEnergies(object):

    def __init__(self, line_string, score_term_list ):
        line_string_list = line_string.split()
        self.score_dict = dict(zip(score_term_list, line_string_list))

    def get_value( self, score_term):
        return self.score_dict.get(score_term)
###############################################################################
#Pose Energies

class PoseEnergies(object):

    def __init__(self):
        self.res_e_list = []

    def loadFile(self, filename):
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
        while current_line_no < end_read_line:
            current_line = lines[current_line_no]
            self.res_e_list.append(ResidueEnergies(current_line, self.score_term_list))
            current_line_no += 1


    def get_score_term_value_for_residue( self, resnum, score_term):
        return self.res_e_list[resnum-1].get_value(score_term)




###############################################################################

pose_energies = PoseEnergies() # instance of PoseEnergies
pose_energies.loadFile("../test/testpdb.txt")

#print pose_energies.res_e_list[0].get_value('fa_atr') 
print pose_energies.get_score_term_value_for_residue(1, "fa_atr")