import argparse
import csv
from decimal import Decimal
from psamm.datasource.native import ModelReader
from psamm.fluxanalysis import FluxBalanceProblem, FluxBalanceError
from psamm.lpsolver.cplex import Solver
from psamm.reaction import Compound
import logging
from psamm.expression import boolean
from six import string_types, iteritems
import math
from psamm import moma
from psamm.lpsolver import lp

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('This script is used to run FBA or Moma with setting '
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
                        help='constraint rule to set, fpr example, '
                             '"rxn1,rxn2,rxn3", means '
                             'flux(rxn1) = flux(rxn2)+flux(rxn3)" ')
    parser.add_argument(
        '--method', metavar='method',
        choices=['fba', 'lin_moma', 'lin_moma2', 'moma', 'moma2'],
        default='fba', help='Select which method to use. '
                            '(fba, lin_moma, lin_moma2, moma, moma2)')
    args = parser.parse_args()
    model_reader = ModelReader.reader_from_path(args.model)
    native = model_reader.create_model()

    # read the exchange list
    input_exchange = []
    with open(args.exchange_list, 'r') as t:
        table = csv.reader(t, dialect='excel-tab')
        for row in table:
            input_exchange.append(row)
            cpd = Compound(row[0], native.extracellular_compartment)
            native.exchange[cpd] = (cpd, 'EX_' + str(cpd),
                                    Decimal(row[1]), Decimal(row[2]))
    print('exchange input:', input_exchange)

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
    deleted_reactions = set()
    testing_genes = []
    if args.gene is not None:
        for l in args.gene:
            for g in l:
                testing_genes.append(g)
        reactions = set(mm.reactions)

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

        # for rxn in deleted_reactions:
        #     mm.remove_reaction(rxn)

    # add rxns as command line requires (args.addrxn)
    if args.addrxn is not None:
        for l in args.addrxn:
            for r in l:
                mm.add_reaction(r)

    def add_constrain(p, rules, ruleEqual):
        if rules is not None:
            for rule in rules:
                print('test rule:', rule)
                ratio1, rxn1, ratio2, rxn2 = rule.split(',')
                rxn1_var = p.get_flux_var(rxn1)
                rxn2_var = p.get_flux_var(rxn2)
                p.prob.add_linear_constraints(
                    rxn1_var == float(ratio1) / float(ratio2) * rxn2_var)
        if ruleEqual is not None:
            for rule in ruleEqual:
                rxn_list = rule.split(',')
                rxn1_var = p.get_flux_var(rxn_list[0])
                rxnTotal_var = 0
                for i in range(1, len(rxn_list)):
                    rxnTotal_var += p.get_flux_var(rxn_list[i])
                p.prob.add_linear_constraints(rxn1_var == rxnTotal_var)

        # if objective is not None:
        #     obj = args.objective.split(',')
        #     if len(obj) == 1:
        #         obj_rxn = args.objective
        #         # bio_var = p.get_flux_var(obj_rxn)
        #         # bio_flux = truncate(p.flux_bound(obj_rxn, 1), 9)
        #         # p.prob.add_linear_constraints(bio_var >= bio_flux)
        #     else:
        #         varying_dict = dict()
        #         for rxn in obj:
        #             varying_dict[rxn] = 1
        #         varying_expression = p.flux_expr(varying_dict)
        #         obj_rxn = varying_expression
        #
        #         # p._solve()
        #         # args.max = p.prob.result.get_value(varying_expression)
        #         # print(args.max)
        #         # c1, = p.prob.add_linear_constraints(
        #         #     varying_expression >= truncate(args.max, 6))
        #         # p._solve()
        #
        # else:
        #     obj_rxn = native.biomass_reaction

        return p

    # def flux_expr_moma(p, reaction):
    #     """Get LP expression representing the reaction flux."""
    #     if isinstance(reaction, dict):
    #         return p._v.expr(iteritems(reaction))
    #     return p._v(reaction)
    #
    # def get_minimal_fba_flux_objDict(self, objective, wt_obj_flux):
    #     """Find the FBA solution that minimizes all the flux values.
    #
    #     Maximize the objective flux then minimize all other fluxes
    #     while keeping the objective flux at the maximum.
    #
    #     Args:
    #         objective: The objective reaction that is maximized.
    #
    #     Returns:
    #         A dictionary of all the reactions and their minimized fluxes.
    #     """
    #     # Define constraints
    #     vs_wt = self._v_wt.set(self._model.reactions)
    #     zs = self._z.set(self._model.reactions)
    #
    #     # wt_obj_flux = self.get_fba_obj_flux(objective)
    #
    #     self._solve()
    #     with self.constraints() as constr:
    #         constr.add(
    #             zs >= vs_wt, vs_wt >= -zs,
    #             self.prob.result.get_value(objective) >= wt_obj_flux)
    #         self._prob.set_objective(self._z.sum(self._model.reactions))
    #         result = self._solve(lp.ObjectiveSense.Minimize)
    #
    #     fba_fluxes = {}
    #     for key in self._model.reactions:
    #         fba_fluxes[key] = result.get_value(self._v_wt[key])
    #     return fba_fluxes

    def get_opt_flux(rxn_set, p, obj_rxn):
        if isinstance(obj_rxn, str):
            p.maximize(obj_reaction)
            opt_flux = p.get_flux(obj_reaction)
        else:
            p.prob.set_objective(obj_reaction)
            p._solve()
            opt_flux = p.prob.result.get_value(obj_reaction)
        return opt_flux

    # optimization
    if args.method == 'fba':
        logger.info('Solving using FBA...')
        print('Solving using FBA...')
        prob = FluxBalanceProblem(mm, Solver())

        prob = add_constrain(prob, args.rule, args.ruleEqual)

        if args.objective is not None:
            obj = args.objective.split(',')
            if len(obj) == 1:
                obj_reaction = args.objective
            else:
                varying_dict = dict()
                for rxn in obj:
                    varying_dict[rxn] = 1
                varying_expression = prob.flux_expr(varying_dict)
                obj_reaction = varying_expression
        else:
            obj_reaction = native.biomass_reaction

        try:
            wild = get_opt_flux(set(mm.reactions), prob, obj_reaction)
        except FluxBalanceError as e:
            logger.error(e)

        if len(deleted_reactions) > 0:
            for reaction in deleted_reactions:
                flux_var = prob.get_flux_var(reaction)
                prob.prob.add_linear_constraints(flux_var == 0)

        opt_obj_flux = get_opt_flux(set(mm.reactions), prob, obj_reaction)

    elif args.method in ['lin_moma', 'lin_moma2', 'moma', 'moma2']:
        prob = moma.MOMAProblem(mm, Solver())
        prob = add_constrain(prob, args.rule, args.ruleEqual)

        if args.objective is not None:
            obj = args.objective.split(',')
            if len(obj) == 1:
                obj_reaction = args.objective
            else:
                varying_dict = dict()
                for rxn in obj:
                    varying_dict[rxn] = 1
                obj_reaction = varying_dict
        else:
            obj_reaction = native.biomass_reaction

        # if isinstance(obj_reaction, dict):
        #     # p._v.expr(iteritems(obj_reaction))
        #     # varying_expression = p._v(obj_reaction)
        #     varying_expression = flux_expr_moma(p, obj_reaction)
        #     obj_reaction = varying_expression
        #     p.prob.set_objective(varying_expression)
        #     p._solve()
        #     # wt_fluxes = p.get_minimal_fba_flux(varying_expression)
        #     # wild = wt_fluxes[varying_expression]
        #     obj_flux = p.prob.result.get_value(varying_expression)
        #     wt_fluxes = get_minimal_fba_flux_objDict(p, varying_expression, obj_flux)
        #
        #     # print('wild:', varying_expression, wild)
        # else:

        wt_fluxes = prob.get_minimal_fba_flux(obj_reaction)
        wild = wt_fluxes[obj_reaction]

        for reaction in deleted_reactions:
            flux_var = prob.get_flux_var(reaction)
            prob.prob.add_linear_constraints(flux_var == 0)

        try:
            if args.method == 'moma':
                logger.info('Solving using MOMA...')
                print('Solving using MOMA...')
                prob.moma(wt_fluxes)
            elif args.method == 'lin_moma':
                logger.info('Solving using linear MOMA...')
                print('Solving using linear MOMA...')
                prob.lin_moma(wt_fluxes)
            elif args.method == 'moma2':
                logger.info('Solving using combined-model MOMA...')
                print('Solving using combined-model MOMA...')
                prob.moma2(obj_reaction, wild)
            elif args.method == 'lin_moma2':
                logger.info('Solving using combined-model linear MOMA...')
                print('Solving using combined-model linear MOMA...')
                prob.lin_moma2(obj_reaction, wild)
        except moma.MOMAError:
            logger.error('Error computing the MOMA result.')

        opt_obj_flux = prob.get_flux(obj_reaction)

    else:
        print('invalid method: {}'.format(args.method))

    for rxn in mm.reactions:
        if rxn in prob._v:
            rx = mm.get_reaction(rxn)
            translated_equation = str(rx.translated_compounds(
                lambda x: compounds_name.get(x, x)))
            print('{}\t{}\t{}'.format(rxn, prob.get_flux(rxn), translated_equation))

    print('optimized objective function "{}" flux: {}'.format(obj_reaction, opt_obj_flux))
