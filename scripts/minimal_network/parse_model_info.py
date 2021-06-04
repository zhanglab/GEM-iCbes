# this script is used to read the old Cbes model and get gene-ec-rxn-pathway association information, the resulting table will be used to make Cbes master annotation table
import csv
from psamm.datasource import native
from collections import defaultdict
from psamm.expression import boolean
from six import iteritems
import re
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--model', type=str, help='Path to input model: model.yaml')
parser.add_argument('--genename', type=str,
                    help="name of gene property, only allow one name")
parser.add_argument('--output', type=str,
                    help='name of output, the gene-ec-rxn-pathway association table, as well as the equation')
args = parser.parse_args()

model_reader = native.ModelReader.reader_from_path(args.model)
native_model = model_reader.create_model()
mm = native_model.create_metabolic_model()

compounds_name = {}
for compound in native_model.compounds:
    compounds_name[compound.id] = compound.name

gene_ec = defaultdict(list)
gene_rxn = defaultdict(list)
gene_pathway = defaultdict(list)
gene_equation = defaultdict(list)
gene_set = set()
cpd_set = set()
invalid_genes = ['None', 'no gene associated', 'gap', 'spontaneous', 'unknown']
a = 0
for rxn in native_model.reactions:
    if rxn.id in set(mm.reactions):
        a += 1
        for c, _ in rxn.equation.compounds:
            cpd_set.add(c.name)
        if rxn.genes is not None and rxn.genes.lower() not in invalid_genes:
            genes = boolean.Expression(rxn.properties[args.genename])
            for g in genes.variables:
                if 'fig' in str(g):
                    g = re.sub('fig\|31899.10.', '', str(g))
                    g = re.sub('\|SAMN[0-9a-zA-Z_]*]', '', str(g))
                gene_set.add(g)

                # build gene_ec dict
                if 'ec' in rxn.properties:
                    if type(rxn.properties['ec']) == list:
                        for i in list(rxn.properties['ec']):
                            gene_ec[g].append(i)
                    else:
                        gene_ec[g].append(rxn.properties['ec'])

                # build gene_pathway dict
                if 'subsystem' in rxn.properties and \
                        rxn.properties['subsystem'] is not None:
                    if type(rxn.properties['subsystem']) == list:
                        for i in list(rxn.properties['subsystem']):
                            gene_pathway[g].append(i)
                    else:
                        gene_pathway[g].append(rxn.properties['subsystem'])

                # build gene_rxn dict
                gene_rxn[g].append(rxn.id)

                # get translated equation
                rx = rxn.equation
                translated_equation = str(rx.translated_compounds(
                    lambda x: compounds_name.get(x, x)))
                gene_equation[g].append(translated_equation)


print('number of unique genes:', len(gene_set))
print('number of reactions (with exchange):', len(set(mm.reactions)))
print('number of reactions (without exchange):', a)
print('number of compounds (used in the model):', len(cpd_set))


with open(args.output, mode='w') as f:
    for i in gene_set:
        f.write('{}\t{}\t{}\t{}\t{}\n'.format(
            i, '#'.join(gene_ec[i]), '#'.join(gene_rxn[i]),
            '#'.join(gene_pathway[i]), '#'.join(gene_equation[i])))
