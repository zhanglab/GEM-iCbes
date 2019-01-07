from psamm.datasource.native import NativeModel, ModelReader
import csv
import os
import yaml
from collections import defaultdict, OrderedDict


mr = ModelReader.reader_from_path('/Users/admin/Projects/4_Caldi/GEM-iCbes/kegg-draft-cbes-model/'
                                  'yaml-reaction-databases/kegg-reaction-database/model.yaml')

nm = mr.create_model()

rxn_list = []
with open('/Users/admin/Projects/4_Caldi/GEM-iCbes/kegg-draft-cbes-model/rxns_not_in_filter.tsv', 'r') as f:
    for row in csv.reader(f, delimiter='\t'):
        rxn_list.append(row[0])

# kegg_rxns = set()
# for rxn in nm.reactions:
#     kegg_rxns.add(rxn.id)
# non_filter=set()
# for i in rxn_list:
#     if i not in kegg_rxns:
#         print(i)
#         non_filter.add(i)
# print(len(non_filter))


def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())


def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))


def model_reactions():
    a = 0
    for reaction in sorted(nm.reactions):
        if reaction.id in rxn_list:
            d = OrderedDict()
            a += 1

            d['id'] = encode_utf8(reaction.id)
            enzymes_list = []
            for i in reaction._properties['enzymes']:
                enzymes_list.append(i)

            pathways_list = []
            if reaction._properties['pathways'] is not None:
                for i in reaction._properties['pathways']:
                    pathways_list.append(i)
            if len(pathways_list) == 0:
                pathways_list = None

            orth_list = []
            if reaction._properties.get('orthology') is not None:
                for i in reaction._properties['orthology']:
                    orth_list.append(i)
            if len(orth_list) == 0:
                orth_list = None

            rpairs_list = []
            if reaction._properties['rpairs'] is not None:
                for i in reaction._properties['rpairs']:
                    rpairs_list.append(i)

            if hasattr(reaction, 'name') and reaction.name is not None:
                d['name'] = encode_utf8(reaction.name)
            if reaction._properties.get('names') is not None:
                names_l = []
                for i in reaction._properties['names']:
                    names_l.append(i)
                d['names'] = encode_utf8(names_l)
            if hasattr(reaction, 'equation') and reaction.equation is not None:
                d['equation'] = encode_utf8(str(reaction.equation))
            if reaction._properties.get('enzymes') is not None:
                d['enzymes'] = encode_utf8(enzymes_list)
            if reaction._properties.get('pathways') is not None:
                d['pathways'] = encode_utf8_list(pathways_list)
            if reaction._properties.get('comment') is not None:
                d['comment'] = encode_utf8(str(reaction._properties['comment']))
            if reaction._properties.get('orthology') is not None:
                d['orthology'] = encode_utf8(orth_list)
            if reaction._properties.get('rpairs') is not None:
                d['rpairs'] = encode_utf8(rpairs_list)
            yield d
    print('number of reactions:', a)


def encode_utf8(s):
    if isinstance(s, unicode):
        return s.encode('utf-8')
    return s


def encode_utf8_list(s):
    if isinstance(s, unicode):
        return s[1].encode('utf-8')
    return s


yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                     dict_constructor)
dest = '.'
yaml_args = {'default_flow_style': False,
             'encoding': 'utf-8',
             'allow_unicode': True}

with open(os.path.join(dest, 'cbes-reactions-kegg-generic.yaml'), 'w+') as f:
    yaml.dump(sorted(list(model_reactions())), f, **yaml_args)






