from ResidueEnergies import PoseEnergies, aminoacids
from ResTypeAverageScores import ResTypeAverageScores

class ResTypesStatisticsCollector(object):

	def __init__(self):
		self.restype_av_scores = {}
		for aminoacid in aminoacids:
			self.restype_av_scores[aminoacid] = ResTypeAverageScores.empty_init( aminoacid ) 


	def add_pose_energies( self, PoseEnergies_instance):
		for aminoacid in aminoacids:
			self.restype_av_scores[aminoacid].add_other_instance( PoseEnergies_instance.get_rtas_for_restype(aminoacid) )
			#print PoseEnergies_instance.get_rtas_for_restype(aminoacid)
			#print self.restype_av_scores[aminoacids]
			#rtas_other_instance = 
			#add_other_instance(rtas_other_instance)


	def calculate_averages_and_stdevs(self, res_type, score_term):
		cur_mean = self.restype_av_scores[res_type].get_mean_val(score_term)
		cur_stddev = self.restype_av_scores[res_type].get_stddev(score_term)
		return (cur_mean, cur_stddev)
		
