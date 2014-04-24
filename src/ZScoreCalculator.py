import numpy as np
from constants import *
from ResidueEnergies import PoseEnergies
from ResTypesStatisticsCollector import ResTypesStatisticsCollector


class ZScoreCalculator(object):

	def __init__(self, pdbfile, statistics_collector):
		self.pe_instance = PoseEnergies()
		self.pe_instance.loadFile(pdbfile)
		self.statistics_collector = statistics_collector
		self.z_score_dict = {}
		self.position_positive_shift = 0
		self.value_pos = 0
		self.position_negative_shift = 0
		self.value_neg = 0


	def calculate_z_scores(self, score_term):
		for re in self.pe_instance.res_e_list:
			score = re.score_dict[score_term]
			neighbor_situation = determine_neighbor_situation(re.res_type, re.number_of_neighbors)
			result = self.statistics_collector.get_mean_and_stddev(re.res_type, score_term, re.number_of_neighbors)
			mean = result[0]
			stddev = result[1]
			zscore = (score - mean)/stddev
			self.z_score_dict[re.res_num] = (zscore, re.res_type)
			#print self.z_score_dict
			#print len(self.z_score_dict.keys())
		#print self.z_score_dict
		return self.z_score_dict

	def calculate_z_scores_combined_minus(self, minuend, subtrahend):
		score_term = combined_score_term_minus(minuend, subtrahend)
		for re in self.pe_instance.res_e_list:
			score = re.score_dict[minuend] - re.score_dict[subtrahend]
			neighbor_situation = determine_neighbor_situation(re.res_type, re.number_of_neighbors)
			result = self.statistics_collector.get_mean_and_stddev(re.res_type, score_term, re.number_of_neighbors)
			mean = result[0]
			stddev = result[1]
			zscore = (score - mean) / stddev
			self.z_score_dict[re.res_num] = (zscore, re.res_type)
		return self.z_score_dict


	def calculate_z_scores_combined_plus(self, score_terms):
		score_term = combined_score_term_plus(score_terms)
		for re in self.pe_instance.res_e_list:
			score = re.score_dict[score_terms[0]] - re.score_dict[score_terms[1]]
			neighbor_situation = determine_neighbor_situation(re.res_type, re.number_of_neighbors)
			result = self.statistics_collector.get_mean_and_stddev(re.res_type, score_term, re.number_of_neighbors)
			mean = result[0]
			stddev = result[1]
			zscore = (score - mean) / stddev
			self.z_score_dict[re.res_num] = (zscore, re.res_type)
		return self.z_score_dict



	def calculate_differences_in_z_scores(self, other_instance, score_term):
		if len(self.pe_instance.res_e_list) != len(other_instance.pe_instance.res_e_list):
			sys.exit('ERROR: The compared sequences do not have the same length!')

		self.calculate_z_scores(score_term)
		other_instance.calculate_z_scores(score_term)
		delta_z_scores = {}

		for res_num in range(len(self.pe_instance.res_e_list)):
			if not str(res_num) in self.z_score_dict.keys():
				#print 'skipped %s' %res_num
				continue
			else: 
				mutation = '-'
				if self.z_score_dict[str(res_num)][1] != other_instance.z_score_dict[str(res_num)][1]:
					mutation = 'mutation'
				delta_z_scores[str(res_num)] =  (other_instance.z_score_dict[str(res_num)][0] - self.z_score_dict[str(res_num)][0]), mutation

		other_instance.position_positive_shift = max(delta_z_scores, key=delta_z_scores.get)
		other_instance.value_pos = delta_z_scores[other_instance.position_positive_shift][0]
		other_instance.position_negative_shift = min(delta_z_scores, key=delta_z_scores.get)
		other_instance.value_neg = delta_z_scores[other_instance.position_negative_shift][0]
		delta_z_scores_sorted = sorted(delta_z_scores, key=delta_z_scores.get, reverse=True)
		#print 'position with biggest positive shift(worst): %s\n (shift: %s)' %(other_instance.position_positive_shift, round(other_instance.value_pos, 3) )
		#print 'position with biggest negative shift (best): %s\n (shift: %s)' %(other_instance.position_negative_shift, round(other_instance.value_neg, 3))
		result = []
		for key in delta_z_scores_sorted:
			result.append([key, round(delta_z_scores[key][0], 4), delta_z_scores[key][1]])
		return result




	def calculate_differences_in_z_scores_combined_minus(self, other_instance, minuend, subtrahend):
		if len(self.pe_instance.res_e_list) != len(other_instance.pe_instance.res_e_list):
			sys.exit('ERROR: The compared sequences do not have the same length!')

		score_term = combined_score_term_minus(minuend, subtrahend)

		self.calculate_z_scores_combined_minus(minuend, subtrahend)
		other_instance.calculate_z_scores_combined_minus(minuend, subtrahend)
		delta_z_scores = {}
		
		for res_num in range(len(self.pe_instance.res_e_list)):
			if not str(res_num) in self.z_score_dict.keys():
				#print 'skipped %s' %res_num
				continue
			else: 
				mutation = '-'
				if self.z_score_dict[str(res_num)][1] != other_instance.z_score_dict[str(res_num)][1]:
					mutation = 'mutation'
				delta_z_scores[str(res_num)] =  (other_instance.z_score_dict[str(res_num)][0] - self.z_score_dict[str(res_num)][0]), mutation

		other_instance.position_positive_shift = max(delta_z_scores, key=delta_z_scores.get)
		other_instance.value_pos = delta_z_scores[other_instance.position_positive_shift][0]
		other_instance.position_negative_shift = min(delta_z_scores, key=delta_z_scores.get)
		other_instance.value_neg = delta_z_scores[other_instance.position_negative_shift][0]
		delta_z_scores_sorted = sorted(delta_z_scores, key=delta_z_scores.get, reverse=True)
		#print 'position with biggest positive shift(worst): %s\n (shift: %s)' %(other_instance.position_positive_shift, round(other_instance.value_pos, 3) )
		#print 'position with biggest negative shift (best): %s\n (shift: %s)' %(other_instance.position_negative_shift, round(other_instance.value_neg, 3))
		result = []
		for key in delta_z_scores_sorted:
			result.append([key, round(delta_z_scores[key][0], 4), delta_z_scores[key][1]])
		return result


	def calculate_differences_in_z_scores_combined_plus(self, other_instance, score_terms):
		if len(self.pe_instance.res_e_list) != len(other_instance.pe_instance.res_e_list):
			sys.exit('ERROR: The compared sequences do not have the same length!')

		score_term = combined_score_term_plus(score_terms)

		self.calculate_z_scores_combined_plus(score_terms)
		other_instance.calculate_z_scores_combined_plus(score_terms)
		delta_z_scores = {}
		
		for res_num in range(len(self.pe_instance.res_e_list)):
			if not str(res_num) in self.z_score_dict.keys():
				#print 'skipped %s' %res_num
				continue
			else: 
				mutation = '-'
				if self.z_score_dict[str(res_num)][1] != other_instance.z_score_dict[str(res_num)][1]:
					mutation = 'mutation'
				delta_z_scores[str(res_num)] =  (other_instance.z_score_dict[str(res_num)][0] - self.z_score_dict[str(res_num)][0]), mutation

		other_instance.position_positive_shift = max(delta_z_scores, key=delta_z_scores.get)
		other_instance.value_pos = delta_z_scores[other_instance.position_positive_shift][0]
		other_instance.position_negative_shift = min(delta_z_scores, key=delta_z_scores.get)
		other_instance.value_neg = delta_z_scores[other_instance.position_negative_shift][0]
		delta_z_scores_sorted = sorted(delta_z_scores, key=delta_z_scores.get, reverse=True)
		#print 'position with biggest positive shift(worst): %s\n (shift: %s)' %(other_instance.position_positive_shift, round(other_instance.value_pos, 3) )
		#print 'position with biggest negative shift (best): %s\n (shift: %s)' %(other_instance.position_negative_shift, round(other_instance.value_neg, 3))
		result = []
		for key in delta_z_scores_sorted:
			result.append([key, round(delta_z_scores[key][0], 4), delta_z_scores[key][1]])
		return result