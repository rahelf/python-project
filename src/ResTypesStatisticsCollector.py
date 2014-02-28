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
			#rtas_other_instance = 
			#add_other_instance(rtas_other_instance)
		
