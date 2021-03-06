critical_distance_squared = 100
#number_of_neighbors_list = range(0,41)

class CalphaCoordinates(object):

	def __init__(self, res_type, res_number, x, y, z):
		self.res_type = res_type
		self.res_number = res_number
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.neighbor_counter = 0

	def are_neighbors(self, other_Calpha_instance, critical_distance_squared):
		if critical_distance_squared >= (self.x - other_Calpha_instance.x)**2 + (self.y - other_Calpha_instance.y)**2 + (self.z - other_Calpha_instance.z)**2:
			self.neighbor_counter += 1
			other_Calpha_instance.neighbor_counter += 1
			#if self.res_number == 209:
			#	print "Counting res %s%s and %s%s as neighbors" % (self.res_type, self.res_number, other_Calpha_instance.res_type, other_Calpha_instance.res_number)

