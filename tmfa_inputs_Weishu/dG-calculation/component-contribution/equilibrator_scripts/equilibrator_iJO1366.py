# The MIT License (MIT)
#
# Copyright (c) 2013 Weizmann Institute of Science
# Copyright (c) 2018 Institute for Molecular Systems Biology,
# ETH Zurich
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import argparse
import csv
import logging
import sys

from numpy import nan, sqrt

from equilibrator_api import (
    ComponentContribution,
    Q_,
    Reaction,
    default_physiological_p_h,
    default_physiological_ionic_strength
)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Calculate potentials for a number of reactions.')
    parser.add_argument(
        'outfile', type=argparse.FileType('w'),
        help='path to output file')
    parser.add_argument(
        '--i',
        type=float,
        help='ionic strength in M',
        default=default_physiological_ionic_strength.m_as("M")
    )
    parser.add_argument(
        '--ph',
        type=float,
        help='pH level',
        default=default_physiological_p_h.m_as("")
    )
    logging.getLogger().setLevel(logging.WARNING)

    args = parser.parse_args()

    sys.stderr.write('pH = %.2f\n' % args.ph)
    sys.stderr.write('I = %.3g M\n' % args.i)

    cc = ComponentContribution()
    cc.p_h = Q_(args.ph)
    cc.ionic_strength = Q_(args.i, "M")

    ids = []
    reactions = []
    with open('iJO1366_reactions.csv', 'r') as f:
        for row in csv.DictReader(f):
            ids.append(row['bigg.reaction'])
            try:
                rxn = cc.parse_reaction_formula(row['formula'])
                cc.standard_dg_prime(rxn)
            except Exception as e:
                print('warning: cannot parse reaction %s because of %s' %
                      (row['bigg.reaction'], str(e)))
                rxn = Reaction({})
            reactions.append(rxn)

    dG0_prime, U = cc.standard_dg_prime_multi(reactions)

    writer = csv.writer(args.outfile)
    header = ['reaction', 'pH', 'ionic strength [M]', 'dG\'0 [kJ/mol]',
              'uncertainty [kJ/mol]', 'ln(Reversibility Index)', 'comment']
    writer.writerow(header)
    for s, r, dg0, u in zip(ids, reactions,
                            dG0_prime.flat, U.diagonal().flat):
        row = [s, args.ph, args.i]
        if r.is_empty():
            row += [nan, nan, nan, 'reaction is empty']
        elif r.is_balanced():
            ln_RI = cc.ln_reversibility_index(r)
            row += [
                '%.2f' % dg0.m_as("kJ/mol"),
                '%.2f' % sqrt(u).m_as("kJ/mol"),
                '%.2f' % ln_RI.nominal_value,
                ''
            ]
        else:
            row += [nan, nan, nan, 'reaction is not chemically balanced']
        writer.writerow(row)

    args.outfile.flush()
