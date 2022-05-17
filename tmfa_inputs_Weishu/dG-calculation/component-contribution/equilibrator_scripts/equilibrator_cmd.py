"""Command-line script for estimating dG of single reactions."""
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


# A stand-alone version of Component Contribution that can calculate the
# Delta-Gr'0 of any reaction (with KEGG notation, i.e. whose reactants are
# already cached in our database), at a given pH and I.

import argparse
import logging
import sys

from equilibrator_api import Q_, ComponentContribution, ureg


def MakeParser():
    """Parser for equilibrator_cmd."""
    parser = argparse.ArgumentParser(
        description=('Estimate the Gibbs energy of a reaction. For example,'
                     'the following calculates dGr0 for ATP hydrolysis '
                     'at pH 6: equlibrator_cmd.py --ph 6 "kegg:C00002 + '
                     'kegg:C00001 = kegg:C00008 + kegg:C00009"'))
    parser.add_argument('--ph', type=str,
                        help='pH level (default = 7.0)',
                        default="7.0")
    parser.add_argument('--pmg', type=str,
                        help='pMg level (default = 3.0)',
                        default="3.0")
    parser.add_argument('--i', type=str,
                        help='ionic strength (in molar, default 0.25 M)',
                        default="0.25M")
    parser.add_argument('--t', type=str,
                        help='temperature (in kalvin, default 298.15 K)',
                        default="298.15K")
    parser.add_argument('reaction', type=str,
                        help='chemical formula using accessions')
    return parser


def main(args):
    """Run main script, calculates the reaction Gibbs energy change."""
    logging.getLogger().setLevel(logging.WARNING)
    ureg.default_format = ".2f"

    p_h = Q_(args.ph)
    assert p_h.check(None)

    p_mg = Q_(args.pmg)
    assert p_mg.check(None)

    ionic_strength = Q_(args.i)
    if ionic_strength.check(None):
        ionic_strength *= Q_("M")
    else:
        assert ionic_strength.check("[concentration]")

    temperature = Q_(args.t)
    if temperature.check(None):
        temperature *= Q_("K")
    else:
        assert temperature.check("[temperature]")

    sys.stderr.write(f"pH = {p_h}\n")
    sys.stderr.write(f"pMg = {p_mg}\n")
    sys.stderr.write(f"I = {ionic_strength:.2g}\n")
    sys.stderr.write(f"T = {temperature}\n")
    sys.stderr.flush()

    equilibrator = ComponentContribution()
    equilibrator.p_h = p_h
    equilibrator.p_mg = p_mg
    equilibrator.ionic_strength = ionic_strength
    equilibrator.temperature = temperature

    # parse the reaction
    try:
        reaction = equilibrator.parse_reaction_formula(args.reaction)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(-1)

    n_e = reaction.check_half_reaction_balancing()
    if n_e is None:
        sys.stderr.write(f"This reaction is not chemically balanced:\n")
        atom_bag = reaction._get_reaction_atom_bag()
        for atom, count in atom_bag.items():
            side = "left" if count < 0 else "right"
            sys.stderr.write(
                f"* {abs(count):g} extra '{atom}'s on the {side}-hand side\n")
        sys.exit(-1)
    elif n_e == 0:
        standard_dg_prime = equilibrator.standard_dg_prime(reaction)
        sys.stdout.write(f"ΔG'° = {standard_dg_prime}\n")

        ln_RI = equilibrator.ln_reversibility_index(reaction)
        sys.stdout.write(f"ln(Reversibility Index) = {ln_RI}\n")
    else:  # treat as a half-reaction
        logging.warning("This reaction isn't balanced, but can still be "
                        "treated as a half-reaction")
        standard_e_prime = equilibrator.standard_e_prime(reaction)
        sys.stdout.write(f"E'° = {standard_e_prime}\n")

    sys.stdout.flush()
    sys.stderr.write(r"* the range represents the 95% confidence interval"
                     " due to Component Contribution estimation uncertainty\n")


if __name__ == "__main__":
    parser = MakeParser()
    args = parser.parse_args()
    main(args)
