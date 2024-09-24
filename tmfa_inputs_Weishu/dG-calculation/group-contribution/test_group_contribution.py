import sys
sys.path.append('/Applications/MarvinSuite/bin')
from group_contribution import group_contribution

gc = group_contribution()

# test case 1: from Keith's script
# test_list = ['cpd_ACP + CHB_15422 + cpd_fa1 = CHB_16027 + cpd_fa1ACP + CHB_18361']
# test_structures = {'cpd_ACP': 'CC(C)(COP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS', 'cpd_fa1': 'CC(C)CCCCCCCCCCC(=O)O', 'cpd_fa1ACP':
#                    'CC(C)CCCCCCCCCCC(=O)OSCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)O'}
#
# new_structures = {'cpd_ACP': 'NCCS', 'cpd_fa1': 'CC(C)CCCCCCCCCCC(=O)O', 'cpd_fa1ACP':
#                    'CC(C)CCCCCCCCCCC(=O)OSCCN'}
#
#
# # dgr_vals = gc.calc_dGr(rxn_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=wp2_structure)
# dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=new_structures)
# print(dgr_vals)


# # test case 2: R00428: C00044 + C00001 => C05922
# test_list2 = ['CHB_15996 + CHB_15377 = C05922']
# # test using smiles from modelseed
# test_structures = {'C05922': 'Nc1nc(N[C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)c(NC=O)c(=O)[nH]1'}
# # test using cannonical smiles converted by above
# # test_structures = {'C05922': 'O=CNc1c(N[C@@H]2O[C@@H]([C@H]([C@H]2O)O)COP(=O)(OP(=O)(OP(=O)([O-])[O-])[O-])[O-])nc([nH]c1=O)N'}
# # test using the cannonical smiles of the uncharged version of C05922 from pubchem
# # test_structures = {'C05922': 'C(C1C(C(C(O1)NC2=C(C(=O)NC(=N2)N)NC=O)O)O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O'}
# # test using the Isomeric smiles of the uncharged version of C05922 from pubchem
# # test_structures = {'C05922': 'C([C@@H]1[C@H]([C@H]([C@@H](O1)NC2=C(C(=O)NC(=N2)N)NC=O)O)O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O'}
#
#
# dgr_vals = gc.calc_dGr(test_list2, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
# print(dgr_vals)


# test case 3: ATPSYN: 4 C00080 + C00002 + C00001 = 5 C00080 + C00008 + C00009 (proton C00080 should be not included)
# test using gc_db data directly
test_list_gcdb = ['CHB_15422 + CHB_15377 = CHB_16761 + CHB_26078']

# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 7.0, IS = 0.0, T = 298.15)
# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 4.0, IS = 0.0, T = 298.15)
# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 9.0, IS = 0.0, T = 298.15)


# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 7.0, IS = 0.25, T = 298.15)
# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 4.0, IS = 0.25, T = 298.15)
# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 9.0, IS = 0.25, T = 298.15)

# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 7.0, IS = 0.036, T = 298.15)
# dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 4.0, IS = 0.036, T = 298.15)
dgr_vals = gc.calc_dGr(test_list_gcdb, pH = 9.0, IS = 0.036, T = 298.15)


# # test using SMILEs to replace some of compound from gc_db mapping (C00002, C00008)
# # test_list = ['C00002 + CHB_15377 = C00008 + CHB_26078']
# # # test
# # test_structures = {'C00002': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N',
# #                     'C00008': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])[O-])O)O)N'}
# test_list = ['C00002 + CHB_15377 = C00008 + CHB_26078']
# # test
# test_structures = {'C00002': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N',
#                     'C00008': 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])[O-])O)O)N'}

#dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)

print(dgr_vals)
