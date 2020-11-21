from __future__ import absolute_import

from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero


def build_aggregate_model(model):
    prob = model.problem
    # empty all constraints
    for constraint in model.constraints:
        model.remove_cons_vars(constraint)

    abundance_var = prob.Variable(
        'abundance_{}'.format(model.name),
        lb=0,
    )
    growth_rate = prob.Variable(
        'growth_rate',
        lb=0,
    )
    model.add_cons_vars([abundance_var, growth_rate])

    rxn = model.reactions.get_by_id(
        'BIOMASS_{}'.format(model.name)
    )
    growth_cons = prob.Constraint(
        1.0 * rxn.forward_variable * abundance_var - growth_rate * abundance_var,
        name='growth_{}'.format(model.name),
        lb=0,
        ub=0
    )
    model.add_cons_vars(growth_cons)
    model.solver.update()

    for met in model.metabolites:
        # Create constraint
        met_cons = prob.Constraint(
            Zero,
            name=met.id,
            lb=0,
            ub=0
        )

        model.add_cons_vars(met_cons)
        rxn_coeffs = []
        # Get stoichiometrich coefficients
        # for all reactions of met
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            rxn_coeffs.append(
                [rxn.forward_variable, coeff * abundance_var]
            )
            rxn_coeffs.append(
                [rxn.reverse_variable, -coeff * abundance_var]
            )

        # Add constraint to model
        rxn_coeffs = dict(rxn_coeffs)
        model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)
