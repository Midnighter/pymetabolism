#!/usr/bin/env python
# -*- coding: utf-8 -*-


import operator
import warnings


def run(model):
    deleted_edges = list()
    ignored_edges = list()
    last_rxn = None
    old_model = None
    while True:
        model._model._model.reset() # caching seems to be a problem atm
        status = model.fba()

        if status == 3:
            print "infeasible model"
            if not last_rxn:
                print "could not revover biomass"
                break
#            model.modify_reaction_bounds({last_rxn: (0, 999999)})
            model = old_model
            ignored_edges.append(last_rxn)
            last_rxn = None
            continue
        elif status > 2:
            raise Exception("Gurobi status = %d" % (status))
        else:
            value = model.get_objective_value()
            if(value <= 0):
                print "no biomass"
                if not last_rxn:
                    print "could not revover biomass"
                    break
#                model.modify_reaction_bounds({last_rxn: (0, 999999)})
                model = old_model
                ignored_edges.append(last_rxn)
                last_rxn = None
                continue

            if last_rxn:
                deleted_edges.append(last_rxn)

# maybe ignore transports, drains, and exchange reactions?
# ignore or knock-out all reactions with 0 flux at once? ignoring right now
            candidates = [c for c in model.get_flux_distribution()
                    if not (c[0] in deleted_edges or c[0] in ignored_edges) and
                    c[1] > 0.0]

            if len(candidates) == 0:
                warnings.warn("no more candidates left")
                break

            candidates.sort(key=operator.itemgetter(1))

            if len(candidates) > 0:
                print "biomass", value, "| ignored", len(ignored_edges),\
                        "| deleted", len(deleted_edges), "| candidates",\
                        len(candidates), "| current", candidates[0]
            else:
                print "biomass", value, "| ignored", len(ignored_edges),\
                        "| deleted", len(deleted_edges), "| candidates",\
                        len(candidates), "| current --"

            old_model = model.copy()
            last_rxn = candidates[0][0]
            model.knockout_reaction(last_rxn)


if __name__ == "__main__":
    from pymetabolism.builders import MetabolicModelBuilder
    builder = MetabolicModelBuilder()
#    (model, known_fluxes) = builder.generate_model("pymetabolism/tests/data/Ec_core_flux1.xml")
    (model, known_fluxes) = builder.generate_model("pymetabolism/tests/data/Ec_iAF1260_flux1.xml")
    print list(model.get_objective_reaction())[0]
    outer = [cmpd for cmpd in model.get_compounds() if cmpd.endswith("_b")]
    model.add_compound_drain(outer, bounds=(0, 1000))
    model.add_compound_transport(outer, bounds=(0, 1000))
    model._model._model.setParam("OutputFlag", 0)
    run(model)
