"""Command-line script for estimating dGs in batch."""
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
import logging
import sys

import numpy as np
import pandas as pd

from equilibrator_api import Q_, ComponentContribution, ureg


def MakeParser():
    parser = argparse.ArgumentParser(
        description='Calculate potentials for a number of reactions.')
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='path to file containing reactions (TXT format)')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='path to result file (CSV format)')
    parser.add_argument('--plaintext', action='store_true',
                        help="indicate that reactions are given in plain text"
                             " (not accessions)")
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

    logging.getLogger().setLevel(logging.WARNING)

    return parser

def main(args):
    """Run main script, calculates the reaction Gibbs energy change."""
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

    infile_lines = list(filter(None, map(str.strip, args.infile.readlines())))

    if args.plaintext:
        raise NotImplementedError("ReactionMatcher needs to be reimplemented")
    else:
        reactions = list(map(equilibrator.parse_reaction_formula, infile_lines))

    result_df = pd.DataFrame(data=list(zip(infile_lines, reactions)),
                             columns=["reaction (original)",
                                      "reaction (parsed)"])

    result_df['pH'] = p_h
    result_df['pMg'] = p_mg
    result_df['I'] = ionic_strength
    result_df['T'] = temperature
    standard_dg_prime, dg_sigma = equilibrator.standard_dg_prime_multi(reactions)
    result_df["ΔG'°"] = list(standard_dg_prime.flat)
    result_df["σ(ΔG'°)"] = list(np.sqrt(dg_sigma.diagonal()).flat)
    result_df["ln(Reversibility index)"] = [
        equilibrator.ln_reversibility_index(r).m_as("") for r in reactions
    ]
    result_df["comment"] = ""

    non_balanced = result_df["reaction (parsed)"].apply(
        lambda r: not r.is_balanced(raise_exception=False)
    )
    result_df.loc[non_balanced, "ln(Reversibility index)"] = np.nan
    result_df.loc[non_balanced, "comment"] = ("reaction is not chemically "
                                              "balanced")
    result_df.to_csv(args.outfile)
    args.outfile.flush()


if __name__ == '__main__':
    parser = MakeParser()
    args = parser.parse_args()
    main(args)
