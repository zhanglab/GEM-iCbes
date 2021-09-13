import csv
import logging
import argparse
from collections import defaultdict
from six import iteritems

logger = logging.getLogger(__name__)
parser = argparse.ArgumentParser(description='Process some inputs')
parser.add_argument('--input', type=str, help='the combined table of '
                                              'randomsparse result')
parser.add_argument('--repeat', type=str,
                    help='number of randomsparse repeats, for example, the '
                         'inpiut table is summary table of 1000 randomsparse, '
                         'then --repeat is 1000')
parser.add_argument('--out1', type=str,
                    help='summary of gene classification (core, alternative, '
                         'non-essential),   ')
# parser.add_argument('--out2', type=str,
#                     help='numbers of core, alternative, non-essential genes '
#                          'in first x randomsparse simulations, x is a range: '
#                          '[1, maximum of repeats)')
args = parser.parse_args()

invalid_genes = ['gap', 'spontaneous', 'unknown']
gene_class = defaultdict(dict)  # for out1
core, alternative, non_essential = 0, 0, 0
repeat_class = defaultdict(list)
for row in csv.reader(open(args.input, mode='r'), delimiter='\t'):
    if row[0].lower() not in invalid_genes:
        # print('\t'.join([row[0], row[1], row[2]]))
        # create dict for out1
        gene_class[row[0]]['essential'] = row.count('1')
        gene_class[row[0]]['non-essential'] = row.count('0')
        if row.count('1') + row.count('0') != int(args.repeat):
            logger.warning("sum of essential and non-essential cases "
                           "doesn't match the repeat number for gene: "
                           "{}".format(row[0]))
        if row.count('1') == int(args.repeat):
            group = 'Core'
            core += 1
        elif row.count('0') == int(args.repeat):
            group = 'Non-essential'
            non_essential += 1
        elif 0 < row.count('1') < int(args.repeat):
            group = 'Flexible'
            alternative += 1
        else:
            group = 'invalid group'
            logger.warning("invalid group for gene {}".format(row[0]))
        gene_class[row[0]]['group'] = group

        # create dict for out2
        for i in range(1, int(args.repeat)+1):
            index = 'repeat_{}'.format(i)
            repeat_class[index].append(row[i])

logger.info("After {} repeated randomsparse, Core genes: {}; Flexible "
            "genes: {}; Non_essential genes: {}".format(
    args.repeat, core, alternative, non_essential))


# print out1: gene is core, alternative or non-essential
if args.out1 is not None:
    with open(args.out1, mode='w') as f:
        f.write('gene_id\tclassification\tin_minimal_network\t'
                'not_in_minimal_network\n')
        for gene, props in iteritems(gene_class):
            f.write('{}\t{}\t{}\t{}\n'.format(
                gene, gene_class[gene]['group'], gene_class[gene]['essential'],
                gene_class[gene]['non-essential']))

# # print out2: in first x repeats, the number of core, alternative,
# # non-essential genes. x is a range: [1, 1000]
# if args.out2 is not None:
#     with open(args.out1, mode='w') as f:
#         f.write('randomsparse_repeats\tcore\talternative\tnon-essential\n')
#         repeat_class_update = defaultdict(dict)
#         for repeat, props in iteritems(repeat_class):
#             repeat_class_update[repeat]['core'] = 0
#             repeat_class_update[repeat]['alternative'] = 0
#             repeat_class_update[repeat]['non-essential'] = 0
#             for i in range(1, len(props)):
#                 if all(value == '1' for value in props[1:i]):
#                     repeat_class_update[repeat]['core'] += 1
#                 elif all(value == '0' for value in props[1:i]):
#                     repeat_class_update[repeat]['alternative'] += 1
#                 elif all(value in ['0', '1'] for value in props[1:i]):
#                     repeat_class_update[repeat]['alternative'] = 0
