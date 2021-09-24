#!/usr/bin/env python
import argparse
from psamm.datasource import native
from psamm.fluxanalysis import FluxBalanceProblem, FluxBalanceError
from psamm.lpsolver import cplex, lp
import logging
from psamm.expression import boolean
from six import string_types
import math
import csv
from psamm.reaction import Compound
from decimal import Decimal

logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(
    description=('Run robustness-fva by varying single or multiple reaction '
                 'fluxes from min to max'))
parser.add_argument('--model', type=str, help='Path to model.yaml file')
parser.add_argument('--objective', type=str, help='reaction to maximize')
parser.add_argument('--min', type=float,
                    help='min varying reaction flux (default 0)')
parser.add_argument('--max', type=float,
                    help='max varying reaction flux (default maximum value from FBA)')
parser.add_argument('--steps', type=int,
                    help='how many steps to be calculated')
parser.add_argument('--rule', type=str,
                    help='constraint rule to set, "ratio1,rxn1,ratio2,rxn2"',
                    action='append')
parser.add_argument('--gene', nargs='*', action='append',
                    help='IDs of genes to delete')
parser.add_argument('--addrxn', nargs='*', action='append',
                    help='IDs of reactions to added')
parser.add_argument('--exchange-list',
                    help=('the list of exchange constraints, three columns'
                          ': compound ID, lower bound, upper bound'))
parser.add_argument('--rxn-constrain',
                    help=('the list of metabolic reaction constraints, three '
                          'columns: reaction ID, lower bound, upper bound'))
parser.add_argument('--vary', type=str,
                    help=("ID(s) of varied reaction(s), separated by comma. "
                          "e.g. 'rxn1,rxn2,rxn3'"))
args = parser.parse_args()

model_reader = native.ModelReader.reader_from_path(args.model)
native = model_reader.create_model()


def truncate(value, n):
    """setting decimal place range without rounding"""
    if value >= 0:
        return math.floor(value * 10 ** n) / 10 ** n
    else:
        value_temp = -1 * value
        return -(math.floor(value_temp * 10 ** n) / 10 ** n)


# read the exchange list to set exchange constraint
if args.exchange_list is not None:
    input_exchange = []
    with open(args.exchange_list, 'r') as t:
        table = csv.reader(t, dialect='excel-tab')
        for row in table:
            input_exchange.append(row)
            cpd = Compound(row[0], native.extracellular_compartment)
            native.exchange[cpd] = (cpd, 'EX_' + str(cpd),
                                    Decimal(row[1]), Decimal(row[2]))

# read constraint of metabolic reactions
if args.rxn_constrain is not None:
    with open(args.rxn_constrain, 'r') as t:
        table = csv.reader(t, dialect='excel-tab')
        for row in table:
            native.limits[row[0]] = (row[0],
                                     Decimal(row[1]), Decimal(row[2]))

# identify gene-rxn associations
genes = set()
gene_assoc = {}
for reaction in native.reactions:
    assoc = None
    if reaction.genes is None:
        continue
    elif isinstance(reaction.genes, string_types):
        assoc = boolean.Expression(reaction.genes)
    else:
        variables = [boolean.Variable(g) for g in reaction.genes]
        assoc = boolean.Expression(boolean.And(*variables))
    genes.update(v.symbol for v in assoc.variables)
    gene_assoc[reaction.id] = assoc

mm = native.create_metabolic_model()

# remove reactions associated with genes that need to be removed
testing_genes = []
if args.gene is not None:
    for l in args.gene:
        for g in l:
            testing_genes.append(g)
    reactions = set(mm.reactions)
    deleted_reactions = set()

    for reaction in reactions:
        if reaction not in gene_assoc:
            continue
        assoc = gene_assoc[reaction]
        if any(boolean.Variable(gene) in assoc.variables for gene in
               testing_genes):
            new_assoc = assoc.substitute(
                lambda v: v if v.symbol not in testing_genes else False)
            if new_assoc.has_value() and not new_assoc.value:
                logger.info('Deleting reaction {}...'.format(reaction))
                deleted_reactions.add(reaction)
    for rxn in deleted_reactions:
        mm.remove_reaction(rxn)

# add reactions based on command line argument --addrxn
if args.addrxn is not None:
    for l in args.addrxn:
        for r in l:
            mm.add_reaction(r)

# set uo a flux balance problem
solver = cplex.Solver()
p = FluxBalanceProblem(mm, solver)

# add constraints to flux balance problem based on command line argument
if args.rule is not None:
    for rule in args.rule:
        ratio1, rxn1, ratio2, rxn2 = rule.split(',')
        rxn1_var = p.get_flux_var(rxn1)
        rxn2_var = p.get_flux_var(rxn2)
        p.prob.add_linear_constraints(
            rxn1_var >= float(ratio1) / float(ratio2) * rxn2_var)

# determined varied reaction or expression
varying_dict = dict()
if args.vary is not None:
    vary_rxns = args.vary.split(',')
    for rxn in vary_rxns:
        varying_dict[rxn] = 1
varying_expression = p.flux_expr(varying_dict)
logger.info('varied reaction: {}'.format(args.vary))

# determine the minimal and maximum flux of varied reaction/expression
if args.min is None:
    args.min = 0
if args.max is None:
    p.prob.set_objective(varying_expression)
    p._solve()
    args.max = p.prob.result.get_value(varying_expression)
    for r in varying_dict.keys():
        logging.info('{}: {}'.format(r, p.get_flux(r)))
logging.info('min: {}'.format(args.min))
logging.info('max: {}'.format(args.max))
if args.min == args.max:
    raise RuntimeError("Can't vary!")

step_range = [args.min]
num = args.min
step = (args.max - args.min) / float(args.steps)
while num + step <= args.max:
    num = num + step
    step_range.append(num)
if round(step_range[-1], 6) < round(args.max, 6):
    step_range.append(args.max)


logging.info("varied range: {}".format(step_range))

# at each flux point of the varied reaction, print flux range of all reactions
for step in step_range:
    bio_var = p.get_flux_var(args.objective)
    p.prob.set_objective(bio_var)
    edited_step = truncate(step, 9)
    c1, = p.prob.add_linear_constraints(varying_expression == edited_step)
    p._solve()
    edited_biomass = truncate(p.get_flux(args.objective), 9)
    if edited_biomass < 1e-9:
        edited_biomass = 0
    c2, = p.prob.add_linear_constraints(p.get_flux_var(
        args.objective) >= edited_biomass)
    for rxn in list(mm.reactions):
        for i in range(10):
            try:
                lower = p.flux_bound(rxn, -1)
                lower = truncate(lower, 9)
            except FluxBalanceError as e:
                # print(rxn, e)
                lower = '-'
            else:
                # print(rxn, 'solved')
                break
        for i in range(10):
            try:
                upper = p.flux_bound(rxn, 1)
                upper = truncate(upper, 9)
            except FluxBalanceError as e:
                # print(rxn, e)
                upper = '-'
            else:
                # print(rxn, 'solved')
                break
        print('{}\t{}\t{}\t{}'.format(
            rxn, edited_step, lower,
            upper))
    c1.delete()
    c2.delete()
