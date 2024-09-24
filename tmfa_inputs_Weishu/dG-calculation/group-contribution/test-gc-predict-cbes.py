import sys
import csv
from psamm.datasource.native import ModelReader
from psamm.reaction import Reaction
from psamm.reaction import Direction
sys.path.append('/Applications/MarvinSuite/bin')


from group_contribution import group_contribution

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--model", help="the path of model.yaml file")
parser.add_argument("--gc-db",help="the compound list with gc-db mapping results")
parser.add_argument("--smiles", help="the compound list with SMILES information")
parser.add_argument("--pH", type=float, help="the pH used to calculate the delta-G")
parser.add_argument("--IS", type=float, help="the ionic strength used to calculate the delta-G")

args=parser.parse_args()


cid_mapping = {}
for row in csv.reader(open(args.gc_db, mode='rU'), delimiter='\t'):
	cid_mapping[row[0]] = row[1]

# print(cid_mapping)


def compound_name(id):
	return cid_mapping.get(id, id)
	# if id not in [i for i in cid_mapping.keys()]:
	# 	return id
	# else:
	# 	return cid_mapping[id]
#
# example_rxn_formula = ['2pg = 3pg']
# in cases where the compound is not in database, the molstring of the compound (e.g. smiles form) need to be provided

wp2_structure = {}
for row in csv.reader(open(args.smiles, mode='rU'), delimiter='\t'):
	wp2_structure[row[0]] = row[1]
# print(wp2_structure)


mr = ModelReader.reader_from_path(args.model)
nm = mr.create_model()
mm = nm.create_metabolic_model()
rxn_id_list = []

rxn_list = []
rxn_dict = {}
for rxn in sorted(mm.reactions):
	if 'sink' in rxn:
		continue
	elif mm.is_exchange(rxn):
	# 	continue
	# elif rxn == 'OMP_H':
		continue
	elif rxn in ['h2ot']:  # exclude the H2O diffusion (if there's a reaction like R00124, ATP + ADP <=> ADP + ATP, it also need to be exluded)
		continue
	else:
		no_structure = 0
		# print('starting reaction {}'.format(rxn))
		rx = str(mm.get_reaction(rxn))
		rx_original = mm.get_reaction(rxn)
		for cpd, _ in rx_original.compounds:
			cpd = cpd.name
			if cpd not in wp2_structure:
				no_structure = 1

		if no_structure == 0:
			lhs = []
			rhs = []
			for i in rx_original.left:
				cpd = i[0].name
				if cpd != 'C00080':
					lhs.append(i)
			for j in rx_original.right:
				cpd = j[0].name
				if cpd != 'C00080':
					rhs.append(j)
			no_h_rxn = Reaction(Direction.Both, lhs, rhs)
			rx = str(no_h_rxn.translated_compounds(compound_name))
			# rx = str(no_h_rxn)

			rx_trans_string_no_comp = rx.replace('C00080[', '[')
			for comp in nm.compartments:
				rx_trans_string_no_comp = rx_trans_string_no_comp.replace('[%s]' % comp.id, '')
			rx_trans_string_no_comp = rx_trans_string_no_comp.replace('(', '')
			rx_trans_string_no_comp = rx_trans_string_no_comp.replace(')', '')
			rx_trans_string_no_comp = rx_trans_string_no_comp.replace('<=>', '=')
			rx_trans_string_no_comp = rx_trans_string_no_comp.replace('=>', '=')
			rxn_list.append(rx_trans_string_no_comp)
			rxn_dict[rx_trans_string_no_comp] = rxn
		else:
			print('\t'.join([rxn, 'missing at least one compound structure']))


gc = group_contribution()

# test_list = ['cpd_ACP + CHB_15422 + cpd_fa1 = CHB_16027 + cpd_fa1ACP + CHB_18361']
# test_structures = {'cpd_ACP': 'CC(C)(COP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS', 'cpd_fa1': 'CC(C)CCCCCCCCCCC(=O)O', 'cpd_fa1ACP':
#                    'CC(C)CCCCCCCCCCC(=O)OSCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)O'}
#
# new_structures = {'cpd_ACP': 'NCCS', 'cpd_fa1': 'CC(C)CCCCCCCCCCC(=O)O', 'cpd_fa1ACP':
#                    'CC(C)CCCCCCCCCCC(=O)OSCCN'}


dgr_vals = gc.calc_dGr(rxn_list, pH = args.pH, IS = args.IS, T = 298.15, cpd_molstring_dict=wp2_structure)
# dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=new_structures)
# print(dgr_vals)
#
for i, j in enumerate(dgr_vals):
	print('\t'.join([rxn_list[i], rxn_dict[rxn_list[i]], str(j)]))
