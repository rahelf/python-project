from ResidueEnergies import PoseEnergies, number_of_neighbors_list
from ResTypeAverageScores import ResTypeAverageScores
from constants import *

#number_of_neighbors_list = range(0,41)


class ResTypesStatisticsCollector(object):

	def __init__(self):
		self.restype_av_scores = {}
		for aminoacid in aminoacids:
			self.restype_av_scores[aminoacid] = ResTypeAverageScores.empty_init( aminoacid ) 


	def add_pose_energies( self, PoseEnergies_instance):
		for aminoacid in aminoacids:
			self.restype_av_scores[aminoacid].add_other_instance( PoseEnergies_instance.get_rtas_for_restype(aminoacid) )


	def calculate_averages_and_stdevs(self, res_type, score_term):
		cur_mean = self.restype_av_scores[res_type].get_mean_val(score_term)
		cur_stddev = self.restype_av_scores[res_type].get_stddev(score_term)
		return (cur_mean, cur_stddev)
		print cur_mean, cur_stddev

	def calculate_averages_and_stdevs_from_subset(self, res_type, score_term, nn_list):
		cur_mean = self.restype_av_scores[res_type].get_mean_val_for_ncounts( score_term, nn_list)
		cur_stddev = self.restype_av_scores[ res_type].get_stddev_for_ncounts( score_term, nn_list)
		return (cur_mean, cur_stddev)
		

	def add_archived_data(self, archived_res_type_average_scores):
		this_aa = archived_res_type_average_scores.res_type
		self.restype_av_scores[this_aa].add_other_instance(archived_res_type_average_scores)

