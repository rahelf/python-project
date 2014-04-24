max_neighbor_count = 40
number_of_neighbors_list = range(0,max_neighbor_count+1)
aminoacids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

def combined_score_term_plus(score_terms_to_be_combined):
	return "+".join(score_terms_to_be_combined)

def combined_score_term_minus(minuend, subtrahend):
	return '-'.join([minuend, subtrahend])

surface_outside = range(0, 11)
surface_inside = range(0, 14)
buried_outside = range(21, 41)
buried_inside = range(24, 41)
intermediate_outside = range(11, 21)
intermediate_inside = range(14, 24)

outside_aa = ['ASP', 'ARG', 'ASN', 'GLN', 'GLU', 'HIS', 'LYS', 'SER', 'THR', 'GLY']
inside_aa = ['ALA', 'CYS', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL', 'PRO']


def determine_neighbor_situation(res_type, number_of_neighbors):
	neighbor_situation = 0
	if res_type in outside_aa:
		if number_of_neighbors in surface_outside:
			neighbor_situation = surface_outside
		elif number_of_neighbors in intermediate_outside:
			neighbor_situation = intermediate_outside
		elif number_of_neighbors in buried_outside:
			neighbor_situation = buried_outside
	elif res_type in inside_aa:
		if number_of_neighbors in surface_inside:
			neighbor_situation = surface_inside
		elif number_of_neighbors in intermediate_inside:
			neighbor_situation = intermediate_inside
		elif number_of_neighbors in buried_inside:
			neighbor_situation = buried_inside
	return neighbor_situation


any_nn = range(0, 41)
main_neighbor_situations_outside = [surface_outside, intermediate_outside, buried_outside, any_nn]
main_neighbor_situations_inside = [surface_inside, intermediate_inside, buried_inside, any_nn]




#score_terms_to_be_combined = ['hbond_bb_sc', 'hbond_sc']