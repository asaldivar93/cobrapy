from __future__ import absolute_import

from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero
from copy import deepcopy


def build_aggregate_model(model):
    ag_model = deepcopy(model)
    prob = ag_model.problem
    # empty all constraints
    for constraint in ag_model.constraints:
        ag_model.remove_cons_vars(constraint)

    abundance_var = prob.Variable(
        'abundance_{}'.format(ag_model.name),
        lb=0,
    )
    growth_rate = prob.Variable(
        'growth_rate',
        lb=0,
    )
    ag_model.add_cons_vars([abundance_var, growth_rate])

    rxn = ag_model.reactions.get_by_id(
        'BIOMASS_{}'.format(ag_model.name)
    )
    growth_cons = prob.Constraint(
        1.0 * rxn.forward_variable * abundance_var - growth_rate * abundance_var,
        name='growth_{}'.format(ag_model.name),
        lb=0,
        ub=0
    )
    ag_model.add_cons_vars(growth_cons)
    ag_model.solver.update()

    cons = Zero
    all_constraints = []
    for met in ag_model.metabolites:
        for rxn in met.reactions:
            coeff = float(rxn.metabolites[met])
            cons += coeff * rxn.forward_variable * abundance_var
            - coeff * rxn.reverse_variable * abundance_var

        # Create constraint
        met_cons = prob.Constraint(
            cons,
            name=met.id,
            lb=0,
            ub=0
        )
        all_constraints.append(met_cons)

    ag_model.add_cons_vars(all_constraints)
    ag_model.solver.update()

    return ag_model
