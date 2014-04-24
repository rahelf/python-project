import numpy as np
from operator import add
import sys
import cPickle
from constants import *
import imp


#import matplotlib.pyplot as matplotlib.pyplot


class ResTypeAverageScores(object):

    def __init__( self, residue_type, score_term_list, pdb_identifier):
        self.res_type_all_score_dict = {}
        self.score_term_list = score_term_list[1:]
        self.best_score_term_dict = {}
        for score_term in self.score_term_list:
            self.best_score_term_dict[score_term] = (10000, "", 0)

        st_counter = 1
        while st_counter < len( score_term_list):
            self.res_type_all_score_dict[score_term_list[st_counter]] = {}
            for n in number_of_neighbors_list:
                (self.res_type_all_score_dict[score_term_list[st_counter]])[number_of_neighbors_list[n]] = []
            st_counter += 1
      
        self.res_type = residue_type
        if( score_term_list == []):
            self.num_entries = 0
        else:
            self.num_entries = 1
        self.pdb_identifier_list = [pdb_identifier]
        #print self.pdb_identifier_list

    @classmethod
    def empty_init( cls, residue_type ):
        return cls( residue_type, [], [])

    def add_other_instance(self, other_instance):

        #print other_instance.best_score_term_dict['rama']
        #1a. safety check whether other instance has same residue type
        if self.res_type != other_instance.res_type:
            sys.exit("ERROR in function add_other_instance: ResTypeAverageScores instance to be added has different residue_type!")

        #in case we got initialized empty we simply copy everything and return
        if self.num_entries == 0: #means we got initialized empty
            self.score_term_list = other_instance.score_term_list
            self.num_entries = other_instance.num_entries
            self.pdb_identifier_list = other_instance.pdb_identifier_list
            self.res_type_all_score_dict = other_instance.res_type_all_score_dict
            self.best_score_term_dict = other_instance.best_score_term_dict
            return
        #1b. safety check whether other instance has same score terms
        if self.score_term_list != other_instance.score_term_list:
            #print 'self.score_term_list:', self.score_term_list, '\nother_instance.score_term_list:', other_instance.score_term_list
            sys.exit("ERROR in function add_other_instance: ResTypeAverageScores_instance to be added has different score_term_list!")
        #2. add num_entries and append lists, also append pdb identifier 
        #= ResTypeAverageScores_instance
        self.num_entries += other_instance.num_entries
        #pdb-identifier should not be added more than once
        for i in range(len(other_instance.pdb_identifier_list)):
            if other_instance.pdb_identifier_list[i] in self.pdb_identifier_list:
                sys.exit('ERROR: pdb identifiers occur multiple times')
        self.pdb_identifier_list.extend(other_instance.pdb_identifier_list)

        for scoreterm in self.score_term_list:
            for n in number_of_neighbors_list:
                self.res_type_all_score_dict[scoreterm][n].extend(other_instance.res_type_all_score_dict[scoreterm][n])

        for score_term in self.score_term_list:
            if self.best_score_term_dict[score_term][0] > other_instance.best_score_term_dict[score_term][0]:
                self.best_score_term_dict[score_term] = other_instance.best_score_term_dict[score_term]



    def add_residue_energies(self,  ResidueEnergies_instance):
        if self.res_type != ResidueEnergies_instance.get_res_type():
            sys.exit("ERROR in function add_residue_energies: score terms of ResidueEnergies_instance are added to wrong residue_type!")
        self.num_entries += 1
        #print ResidueEnergies_instance.res_type, ResidueEnergies_instance.res_num, ResidueEnergies_instance.get_nneighbors()
        #print self.res_type_all_score_dict
        which_neighbor = ResidueEnergies_instance.get_nneighbors()
        if which_neighbor > max_neighbor_count:
            which_neighbor = max_neighbor_count

        for score_term in self.res_type_all_score_dict.keys():
            self.res_type_all_score_dict[score_term][which_neighbor].append( ResidueEnergies_instance.get_value(score_term) )

        for score_term in self.score_term_list:
            if ResidueEnergies_instance.information_dict[score_term][0] < self.best_score_term_dict[score_term][0]:
                self.best_score_term_dict[score_term] = ResidueEnergies_instance.information_dict[score_term]



    def get_merged_list_for_ncounts(self, score_term, nn_list):
        merged_list = []
        for ncount in nn_list:
            merged_list += self.res_type_all_score_dict[score_term][ncount]
        return merged_list


    def get_merged_list_for_all_nn( self, score_term):
        merged_list =[]
        for item in number_of_neighbors_list:
            merged_list += self.res_type_all_score_dict[score_term][item]
        return merged_list

    def get_relative_frequency_for_nn(self, nn):
        nn = [nn]
        relative_frequency = float(len(self.get_merged_list_for_ncounts('fa_atr', nn)))/len(self.get_merged_list_for_all_nn('fa_atr'))
        return relative_frequency


    def get_best_score(self, score_term):
        return self.best_score_term_dict[score_term]


    def get_mean_val(self, score_term):
        return np.mean( self.get_merged_list_for_all_nn( score_term))

    def get_stddev(self, score_term):
        return np.std( self.get_merged_list_for_all_nn( score_term) )

    def get_mean_val_for_ncounts(self, score_term, nn_list): #nnlist should be list of the desired neighbour counts        
        return np.mean(self.get_merged_list_for_ncounts( score_term, nn_list))

    def get_stddev_for_ncounts(self, score_term, nn_list):
        return np.std( self.get_merged_list_for_ncounts( score_term, nn_list))

    def get_scores_for_scoreterm( self, score_term):
        return self.res_type_all_score_dict[score_term]


    def make_histogram_for_scoreterm_for_ncounts(self, score_term, nn_list, number_of_residues, histogram_location, mean, stddev):
        import matplotlib.pyplot as plt
        plt.clf()
        data = self.get_merged_list_for_ncounts(score_term, nn_list)
        minx = int( np.floor(np.min(data)) )
        maxx = int( np.ceil(np.max(data)) )
        plt.title('%s %s: %s to %s neighbors' %(self.res_type, score_term, nn_list[0], nn_list[-1]))
        bins = int( (maxx - minx)/0.25)
        if bins == 0:
            bins = 80
        plt.hist( data, bins, range=[minx, maxx], label=score_term, histtype='stepfilled', normed = True, color = 'g')
        plt.ylim(0, 4)
        plt.xlim(-10,10)
        plt.annotate('number of residues: %s \nmean: %s \nstddev: %s' %(number_of_residues, round(mean, 4), round(stddev, 4)), xy=(0.60, 0.85), xycoords='axes fraction')
        plt.xticks(range(-10,11,2))
        plt.xlabel('score')
        plt.ylabel('relative frequency')
        plt.savefig(histogram_location + '%s_%s_%s_to_%s_neighbors.pdf' %(self.res_type, score_term, nn_list[0], nn_list[-1]))
        #plt.show()

    def plot_relative_frequencies_of_numbers_of_neighbors(self):
        import matplotlib.pyplot as plt
        plt.clf()
        data = []
        for nn in number_of_neighbors_list:
            data.append(self.get_relative_frequency_for_nn(nn))
        plt.title('neighbors of %s' %self.res_type)
        plt.ylim(0, 0.09)
        plt.ylabel('relative frequency')
        plt.xlabel('neighbor number')
        plt.plot(number_of_neighbors_list, data, 'ro')
        plt.savefig('../../neighbor-plots/'+self.res_type+'.pdf' )
        #plt.show()
        
    def pickle_res_type_average_scores(self, filename):
        pickle_file = open(filename, 'w')
        cPickle.dump(self, pickle_file)
        pickle_file.close()

    def calculate_sum_of_several_score_terms(self, score_terms_to_be_combined):
        combined_score_term = combined_score_term_plus(score_terms_to_be_combined)
        #print combined_score_term
        self.score_term_list.append(combined_score_term)
        self.res_type_all_score_dict[combined_score_term] = {}
        for nn in number_of_neighbors_list:
            self.res_type_all_score_dict[combined_score_term][nn] = []
            for i in range(len(self.res_type_all_score_dict[self.score_term_list[1]][nn])):
                sum_of_scores = 0
                for j in range(len(score_terms_to_be_combined)):
                    sum_of_scores += (self.res_type_all_score_dict[score_terms_to_be_combined[j]][nn])[i]
                self.res_type_all_score_dict[combined_score_term][nn].append(sum_of_scores)


    def subtract_score_terms(self, minuend, subtrahend):
        combined_score_term = combined_score_term_minus(minuend, subtrahend)
        #print combined_score_term
        self.score_term_list.append(combined_score_term)
        self.res_type_all_score_dict[combined_score_term] = {}
        for nn in number_of_neighbors_list:
            self.res_type_all_score_dict[combined_score_term][nn] = []
            for i in range(len(self.res_type_all_score_dict[self.score_term_list[1]][nn])):
                diff = self.res_type_all_score_dict[minuend][nn][i] - self.res_type_all_score_dict[subtrahend][nn][i]
                self.res_type_all_score_dict[combined_score_term][nn].append(diff)