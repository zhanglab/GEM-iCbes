import sys
sys.path.append('/Applications/MarvinSuite/bin')
from group_contribution import group_contribution

gc = group_contribution()


rxn = 'R00575'
test_list = ['2 CHB_15422 + CHB_18050 + C00288 + CHB_15377 = 2 CHB_16761 + CHB_26078 + CHB_29985 + CHB_17672']
test_structures = {'C00288': 'C(=O)(O)[O-]'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)

rxn = 'R00742'
test_list = ['CHB_15422 + CHB_15351 + C00288 = CHB_16761 + CHB_26078 + CHB_15531']
test_structures = {'C00288': 'C(=O)(O)[O-]'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)


rxn = 'R01859'
test_list = ['CHB_15422 + C00100 + C00288 = CHB_16761 + CHB_26078 + CHB_15466']
test_structures = {'C00288': 'C(=O)(O)[O-]',
                    'C00100': 'CCC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)([O-])[O-])O'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)


rxn = 'R05168'
test_list = ['C06376 + CHB_15377 = C06377 + CHB_30089']
test_structures = {'C06376': 'CC(=O)NC1C(C(C(OC1O)COP(=O)([O-])[O-])O)O',
                    'C06377': 'C(C1C(C(C(C(O1)O)[NH3+])O)O)OP(=O)([O-])[O-]'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)

rxn = 'R06590'
test_list = ['C12214 + CHB_52742 = C12215 + MAN_10207']
test_structures = {'C12214': 'C(OP([O-])(=O)[O-])[C@@H]1([C@@H](O)[C@H]([NH3+])[C@@](O)(CO)O1)',
                    'C12215': 'C([C@H](O)[C@H](O)COP([O-])(=O)[O-])=[NH2+]'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)


rxn = 'R07404'
test_list = ['CHB_15422 + CHB_28843 + C00288 = CHB_16761 + CHB_26078 + C15667']
test_structures = {'C00288': 'C(=O)(O)[O-]',
                    'C15667': 'C1=C(N(C=N1)C2C(C(C(O2)COP(=O)([O-])[O-])O)O)NC(=O)[O-]'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)


rxn = 'R08924'
test_list = ['CHB_15422 + C00288 + C00136 = CHB_16761 + CHB_26078 + C18026']
test_structures = {'C00288': 'C(=O)(O)[O-]',
                    'C00136': 'CCCC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)([O-])[O-])O',
                    'C18026': 'CCC(C(=O)[O-])C(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)([O-])[O-])O'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)


rxn = 'hco3'
test_list = ['CHB_16526 + CHB_15377 = C00288']
test_structures = {'C00288': 'C(=O)(O)[O-]'}
dgr_vals = gc.calc_dGr(test_list, pH = 7.0, IS = 0.0, T = 298.15, cpd_molstring_dict=test_structures)
print(rxn, dgr_vals)
