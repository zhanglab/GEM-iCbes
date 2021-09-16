import argparse
import csv
from decimal import Decimal
from psamm.datasource.native import ModelReader
from psamm.fluxanalysis import FluxBalanceProblem, FluxBalanceError
from psamm.lpsolver.cplex import Solver
from psamm.reaction import Compound
import logging
from psamm.expression import boolean
from six import string_types
import math
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('This script is used to run FVA with setting '
                     'lower and upper bound of target compound(s) in exchange,'
                     ' and optimize the flux of objective'))
    parser.add_argument('-m', '--model', help='path to model')
    parser.add_argument('--objective', help='objective to optimize')
    parser.add_argument('--rule', action='append',
                        help='constraint rule to set, '
                             '"ratio1,rxn1,ratio2,rxn2", means '
                             'flux(rxn1) : flux(rxn2) = '
                             'float(ratio1) : float(ratio2)')
    parser.add_argument('--exchange-list',
                        help=('the list of exchange constraints, three columns'
                              ': compound ID, lower bound, upper bound'))
    parser.add_argument('--gene', nargs='*', action='append',
                        help='IDs of genes to delete')
    parser.add_argument('--addrxn', nargs='*', action='append',
                        help='IDs of reactions to added')
    parser.add_argument('--rxn-constrain', help=(
        'the list of metabolic reaction constraints, '
        'three columns: reaction ID, lower bound, upper bound'))
    parser.add_argument('--ruleEqual', action='append',
                        help='constraint rule to set, for example, '
                             '"rxn1,rxn2,rxn3", means '
                             'flux(rxn1) = flux(rxn2)+flux(rxn3)" ')
    parser.add_argument('--rule-non0',
                        help='constrain reaction flux can NOT be 0. '
                             'Usage: "rxn1,rxn2,rxn3"')
    parser.add_argument('--method', help='choose to export FBA or FVA results',
                        choices=['fva', 'fba'], default='fva')
    parser.add_argument('--decimal', help='determine how many digits will be '
                                          'kept after decimal point', default=4)
    parser.add_argument('--loop-removal', type=str,
                        default=None, help='specify loop removal method')
    args = parser.parse_args()
    model_reader = ModelReader.reader_from_path(args.model)
    native = model_reader.create_model()

    epsilon = 1e-5

    # read the exchange list
    input_exchange = []
    if args.exchange_list is not None:
        with open(args.exchange_list, 'r') as t:
            table = csv.reader(t, dialect='excel-tab')
            for row in table:
                input_exchange.append(row)
                cpd = Compound(row[0], native.extracellular_compartment)
                native.exchange[cpd] = (cpd, 'EX_' + str(cpd),
                                        Decimal(row[1]), Decimal(row[2]))

    if args.rxn_constrain is not None:
        with open(args.rxn_constrain, 'r') as t:
            table = csv.reader(t, dialect='excel-tab')
            for row in table:
                native.limits[row[0]] = (row[0],
                                         Decimal(row[1]), Decimal(row[2]))

    compounds_name = {}
    for compound in native.compounds:
        compounds_name[compound.id] = compound.name

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
        gene_assoc[reaction.id] = assoc

    def truncate(value, n):
        """setting decimal place range without rounding"""
        if value > 0:
            return math.floor(value * 10 ** n) / 10 ** n
        else:
            value_temp = -1 * value
            return -(math.floor(value_temp * 10 ** n) / 10 ** n)

    # below are the constraints of the compound provided by the exchange list
    mm = native.create_metabolic_model()

    # delete rxns associated with deleted genes (args.gene)
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

    # add rxns as command line requires (args.addrxn)
    if args.addrxn is not None:
        for l in args.addrxn:
            for r in l:
                mm.add_reaction(r)

    # FVA
    p = FluxBalanceProblem(mm, Solver())
    if args.rule is not None:
        # varying_dict = dict()
        for rule in args.rule:
            ratio1, rxn1, ratio2, rxn2 = rule.split(',')
            rxn1_var = p.get_flux_var(rxn1)
            rxn2_var = p.get_flux_var(rxn2)
            # varying_dict[rxn1] = 1
            # varying_dict[rxn2] = 1
            p.prob.add_linear_constraints(
                rxn1_var == float(ratio1) / float(ratio2) * rxn2_var)
        # # constrains for multiple reactions
        # varying_expression = p.flux_expr(varying_dict)
    if args.ruleEqual is not None:
        for rule in args.ruleEqual:
            rxn_list = rule.split(',')
            rxn1_var = p.get_flux_var(rxn_list[0])
            rxnTotal_var = 0
            for i in range(1, len(rxn_list)):
                rxnTotal_var += p.get_flux_var(rxn_list[i])
            p.prob.add_linear_constraints(rxn1_var == rxnTotal_var)

    if args.rule_non0 is not None:
        for rxn in args.rule_non0.split(','):
            rxn_var = p.get_flux_var(rxn)
            p.prob.add_linear_constraints(rxn_var >= 1e-6)

    if args.objective is not None:
        obj = args.objective.split(',')
        if len(obj) == 1:
            obj_rxn = args.objective
            bio_var = p.get_flux_var(obj_rxn)
            bio_flux = truncate(p.flux_bound(obj_rxn, 1), 6)
            p.prob.add_linear_constraints(bio_var >= bio_flux)
        else:
            varying_dict = dict()
            for rxn in obj:
                varying_dict[rxn] = 1
            varying_expression = p.flux_expr(varying_dict)
            p.prob.set_objective(varying_expression)
            p._solve()
            args.max = p.prob.result.get_value(varying_expression)
            objective_rxns = '+'.join(obj)
            print('maximum flux of objective {}: {}'.format(objective_rxns, args.max))
            c1, = p.prob.add_linear_constraints(
                varying_expression >= truncate(args.max, 6))

    else:
        obj = [native.biomass_reaction]
        obj_rxn = native.biomass_reaction
        bio_var = p.get_flux_var(obj_rxn)
        bio_flux = truncate(p.flux_bound(obj_rxn, 1), 6)
        p.prob.add_linear_constraints(bio_var >= bio_flux)

    # Maximize reaction flux
    if args.method == 'fba' and args.loop_removal == 'l1min':
        if len(obj) == 1:
            try:
                p.maximize(obj_rxn)
            except FluxBalanceError as e:
                logger.error(e)
            fluxes = {r: p.get_flux(r) for r in mm.reactions}
            p.minimize_l1()
        else:
            logger.info('objective contains more than one reaction, cannot apply l1min in current script')
    p._solve()

    print('exchange input:', input_exchange)
    count = 0
    for rxn in mm.reactions:
        rx = mm.get_reaction(rxn)
        translated_equation = str(rx.translated_compounds(
            lambda x: compounds_name.get(x, x)))
        # print(rxn, p.flux_bound(rxn, -1), p.flux_bound(rxn, 1), sep='\t')
        if args.method == 'fva':
            flux_range = '[{}, {}]'.format(round(p.flux_bound(rxn, -1), int(args.decimal)),
                                           round(p.flux_bound(rxn, 1), int(args.decimal)))
            print('{}\t{}\t{}\t{}\t{}'.format(
                rxn, p.flux_bound(rxn, -1), p.flux_bound(rxn, 1), flux_range,
                translated_equation))
        elif args.method == 'fba':
            if args.loop_removal is None:
                print('{}\t{}\t{}'.format(
                    rxn, p.get_flux(rxn), translated_equation))
            elif args.loop_removal == 'l1min':
                flux = p.get_flux(rxn)
                if abs(flux - fluxes[rxn]) > epsilon:
                    count += 1
                print('{}\t{}\t{}'.format(
                    rxn, p.get_flux(rxn), translated_equation))
    logger.info('Minimized reactions: {}'.format(count))
