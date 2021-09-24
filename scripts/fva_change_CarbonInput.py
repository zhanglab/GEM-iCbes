import argparse
import csv
from decimal import Decimal
import math

from psamm.datasource.native import ModelReader
from psamm.fluxanalysis import FluxBalanceProblem
from psamm.lpsolver.cplex import Solver
from psamm.reaction import Compound, Reaction, Direction
import logging
from psamm.expression import boolean
from six import string_types

logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO)

def truncate(value, n):
    """setting decimal place range without rounding"""
    return math.floor(value * 10 ** n) / 10 ** n


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('This script is used to run FVA with setting '
                     'lower and upper bound of target compound(s) in exchange,'
                     ' and optimize the flux of objective')
    )
    parser.add_argument('-m', '--model', help='path to model')
    parser.add_argument('--objective', help='objective to optimize')
    parser.add_argument('--exchange-list',
                        help=('the list of carbon input, three columns: '
                              'compound ID (no compartment), lower bound, '
                              'upper bound'))
    args = parser.parse_args()

    model_reader = ModelReader.reader_from_path(args.model)
    native = model_reader.create_model()

    cpd_dict = {}
    for compound in native.compounds:
        cpd_dict[compound.id] = compound

    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
        'carbon_source', 'compound_name', 'formula', 'objective', 'lower_bound_of_objective',
        'upper_bound_of_objective', 'lower_bound_of_carbon',
        'upper_bound_of_carbon', 'carbon_used_up?'))

    # read the exchange list
    for row in csv.reader(open(args.exchange_list, 'r'), delimiter='\t'):
        cpd = Compound(row[0], 'e')
        # native.exchange[cpd] = (cpd, 'EX_' + str(cpd) + '[e]',
        #                         Decimal(row[1]), Decimal(row[2]))
        mm = native.create_metabolic_model()
        reaction_id = 'EX_{}[e]'.format(str(cpd))
        mm.database.set_reaction(
            reaction_id, Reaction(Direction.Both, {cpd: -1}))
        mm.add_reaction(reaction_id)
        mm.limits[reaction_id].lower = Decimal(row[1])
        mm.limits[reaction_id].upper = Decimal(row[2])

        # below are the constraints of the compound provided by the exchange list

        # FVA
        p = FluxBalanceProblem(mm, Solver())

        bio_var = p.get_flux_var(args.objective)
        # print('test:', bio_var)
        bio_flux = p.flux_bound(args.objective, 1)
        edited_bioflux = truncate(bio_flux, 6)
        # print('test:', bio_flux, edited_bioflux)
        p.prob.add_linear_constraints(
            bio_var == edited_bioflux)

        formula = 'NA'
        if row[0] in cpd_dict:
            formula = cpd_dict[row[0]].formula
            name = cpd_dict[row[0]].name
        else:
            name = row[3]

        lower_carbon = p.flux_bound(reaction_id, -1)
        upper_carbon = p.flux_bound(reaction_id, 1)
        if round(lower_carbon, 2) == round(upper_carbon, 2) and lower_carbon == float(row[1]):
            used_up = 'Yes'
        else:
            used_up = 'No'

        # print(objective)
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            row[0], name, formula, args.objective, p.flux_bound(args.objective, -1),
            p.flux_bound(args.objective, 1), lower_carbon, upper_carbon,
            used_up))
