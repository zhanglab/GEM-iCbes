import csv
from collections import defaultdict

ec_rxn = defaultdict(list)
with open('/Users/admin/Projects/4_Caldi/GEM-iCbes/kegg-draft-cbes-model/reaction_enzyme.list', 'r') as f:
    for row in csv.reader(f, delimiter='\t'):
        ec_rxn[row[1]].append(row[0])

ko_rxn = defaultdict(list)
with open('/Users/admin/Projects/4_Caldi/GEM-iCbes/kegg-draft-cbes-model/reaction_ko.list', 'r') as f:
    for row in csv.reader(f, delimiter='\t'):
        ko_rxn[row[1]].append(row[0])

gene_ko = defaultdict(list)
with open('/Users/admin/Projects/4_Caldi/GEM-iCbes/kegg-draft-cbes-model/ate/ate_ko.list', 'r') as f:
    for row in csv.reader(f, delimiter='\t'):
        gene_ko[row[0]].append(row[1])

gene_ec = defaultdict(list)
print('{}\t{}\t{}\t{}\t{}'.format('gene_id', 'KO', 'rxn from KO', 'EC', 'rxn from EC'))
with open('/Users/admin/Projects/4_Caldi/GEM-iCbes/kegg-draft-cbes-model/ate/ate_enzyme.list', 'r') as f:
    for row in csv.reader(f, delimiter='\t'):
        gene_ec[row[0]].append(row[1])

gene_set = set()
for i in gene_ec:
    gene_set.add(i)
for j in gene_ko:
    gene_set.add(j)

for g in gene_set:
    ko_list = gene_ko.get(g, '-')
    korxns = []
    for ko in ko_list:
        korxns += ko_rxn.get(ko, '-')
    ec_list = gene_ec.get(g, '-')
    ecrxns = []
    for ec in ec_list:
        ecrxns += ec_rxn.get(ec, '-')
    new_ko = ','.join(ko_list)
    new_korxns = ','.join(korxns)
    new_ec = ','.join(ec_list)
    new_ecrxns = ','.join(ecrxns)


    print('{}\t{}\t{}\t{}\t{}'.format(g, new_ko, new_korxns, new_ec, new_ecrxns))





