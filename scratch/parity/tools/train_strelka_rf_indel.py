#!/usr/bin/env python
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT OTHER WARRANTY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""Train the StrelkaRFIndel model (somatic tree), py2/py3 shared driver.

``evs_learn.py`` cannot select this model in the somatic tree:
``evs/__init__`` imports only ``somatic_rf`` (pre-existing at HEAD, identical
under py2 and py3: "argument -m/--model: invalid choice: 'strelka.rf.indel'").
This driver runs the same training flow as ``evs_learn.py`` (same data
handling, same parameter file) with an explicit import of the model module,
so both tool chains exercise the converted ``strelka_rf_indel`` code.

Usage:
    train_strelka_rf_indel.py --libdir <somatic/lib> --parameters FILE -o OUT INPUT_CSV
"""

import argparse
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="labeled indel feature CSV (with NT/INDELTYPE)")
    parser.add_argument("--libdir", required=True,
                        help="the somatic lib directory (holds the evs package)")
    parser.add_argument("-p", "--parameters",
                        help="training parameter JSON file")
    parser.add_argument("-o", "--output", required=True, help="output pickle")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        raise Exception("Can't find input features CSV: '%s'" % args.input)
    sys.path.insert(0, os.path.abspath(args.libdir))

    import evs  # noqa: E402
    import evs.features  # noqa: E402
    import evs.strelka_rf_indel  # noqa: E402,F401  (registers strelka.rf.indel)
    import pandas  # noqa: E402

    fset = evs.features.FeatureSet.make("somatic.indel")
    features = fset.trainingfeatures()

    pars = {}
    if args.parameters:
        pars = json.load(open(args.parameters))
        print("Using custom learning parameters: %s" % str(pars))

    # same data handling as evs_learn.py getDataSet (no subsampling/balancing)
    df = pandas.read_csv(args.input)
    df = df[df["tag"] != "FN"]
    df["weight"] = 1
    tpdata = df[df["tag"] == "TP"]
    fpdata = df[df["tag"] == "FP"]
    print("TP: %d FP: %d" % (tpdata.shape[0], fpdata.shape[0]))

    model = evs.EVSModel.createNew("strelka.rf.indel")
    model.train(tpdata, fpdata, features, **pars)
    model.save(args.output)
    print("wrote", args.output)


if __name__ == "__main__":
    main()
