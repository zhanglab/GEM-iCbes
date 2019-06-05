import csv
from psamm.datasource.native import ModelReader, ModelWriter, NativeModel
import argparse

parser = argparse.ArgumentParser(description='Process some inputs.')

parser.add_argument('--model', help='input model file (a model.yaml file)')
parser.add_argument('--output', help='output YAML file, a subset of input')

args = parser.parse_args()

mr = ModelReader.reader_from_path(args.model)
nm = mr.create_model()
# mm = nm.create_metabolic_model()

out_nm = NativeModel()

a, b = 0, 0
rxns = set()
for rxn in nm.reactions:
    if 'enzymes' in rxn.properties:
        # print(rxn.id)
        a += 1
        rxn.properties['ec'] = rxn.properties['enzymes']
        rxn.properties['enzymes'] = None
    else:
        print(rxn.id)

    out_nm.reactions.add_entry(rxn)

print(a)
with open(args.output, 'w') as out:
    out_mw = ModelWriter()
    out_mw.write_reactions(out, sorted(nm.reactions, key=lambda x: x.properties.get('pathway')))


