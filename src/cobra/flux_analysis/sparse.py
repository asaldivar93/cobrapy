from __future__ import absolute_import

from cobra.flux_analysis.helpers import normalize_cutoff
from optlang.symbolics import Zero


# %% codecell
def find_leaks(model):
    """Following [ref] finds all the posible leaks of a model
    by maximizing the following program:

    max_y,v ||y||_0
    s.t. S*v - y = 0
         0 <= y for met in metabolites
         lb <= v <= ub

    Since ||y||_0 is an NP-Hard problem, this is re casted
    as a LP using L1-norm as an approximation
    Parameters
    ----------
    model: cobra.core.Model
        The cobra model to find leaks from
    Returns
    -------
    leaks: list
        The list of meatabolites considered leaks.
    References
    -------
    [ref]
    """
    w_model = model.copy()
    prob = w_model.problem
    zero_cutoff = normalize_cutoff(model, None)

    # keep only internal reactions
    for rxn in w_model.reactions:
        if rxn.boundary:
            w_model.remove_reactions(rxn)

    # empty all constraints
    for const in w_model.constraints:
        w_model.remove_cons_vars(const)

    obj_vars = []
    for met in w_model.metabolites:
        # Create auxilliary variable y
        met_var = prob.Variable("aux_{}".format(met.id),
                                lb=0
                                )
        # Create constraint
        const = prob.Constraint(Zero,
                                name=met.id,
                                lb=0,
                                ub=0
                                )
        w_model.add_cons_vars([met_var, const])

        rxn_coeffs = []
        # Get stoichiometrich coefficients
        # for all reactions of met
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            rxn_coeffs.append([rxn.forward_variable, coeff])
            rxn_coeffs.append([rxn.reverse_variable, -coeff])
        # Add auxiliary var to constraint
        rxn_coeffs.append([met_var, -1])
        # Add constraint to model
        rxn_coeffs = dict(rxn_coeffs)
        w_model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)

        obj_vars.extend([met_var])

    # Set objective to max(sum(met_vars))
    w_model.objective = prob.Objective(Zero, direction='max')
    w_model.objective.set_linear_coefficients({o: 1 for o in obj_vars})

    w_model.optimize()
    leaks = []
    # if y_i > 0 met_i is a leak in the model
    for var in w_model.variables:
        if 'aux' in var.name:
            if var.primal > zero_cutoff:
                leaks.append([var.name.replace('aux_', ''), var.primal])

    return leaks


def find_leak_mode(model, leaks=[], cutoff_mult=1):
    """ The minimal set of reactions needed for production
    of leak metabolites (leak_mode) cand be found by min of
    the following program:

    min_y,v ||v||_0
    s.t. S*v - y = 0
         y >= 1 for met in leaks
         y >= 0 for met not in leaks
         lb <= v <= ub

    To approximate ||v||_0, we add an indicator variable z
    for every reaction and minimized L_2 norm of vector z in
    the following QP:

    min_z,y,v (sum_r(z_r**2))**(1/2)
    s.t. S*v - y = 0
         v - z = 0
         z unbounded
         y >= 1 for met in leaks
         y >= 0 for met not in leaks
         lb <= v <= ub

    Parameters
    ----------
    model: cobra.core.Model
        The cobra model to find leaks from
    leaks: list
        A list of leaking metabolites and fluxes
    Returns
    -------
    leak_modes: dict
        A dictionary of active reactions for every
        leak metabolite
    References
    -------
    [ref]
    """
    zero_cutoff = normalize_cutoff(model, None)
    w_model = model.copy()
    prob = w_model.problem
    objective = Zero
    rxn_vars_and_cons = []

    # Keep only internal reactions
    for rxn in w_model.reactions:
        if rxn.boundary:
            w_model.remove_reactions(rxn)
    # Empty all constraints
    for const in w_model.constraints:
        w_model.remove_cons_vars(const)

    # Add y variables and associated constraints
    for met in w_model.metabolites:
        met_var = prob.Variable("aux_{}".format(met.id),
                                lb=0
                                )

        const = prob.Constraint(Zero,
                                name=met.id,
                                lb=0,
                                ub=0
                                )
        w_model.add_cons_vars([met_var, const])

        rxn_coeffs = []
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            rxn_coeffs.append([rxn.forward_variable, coeff])
            rxn_coeffs.append([rxn.reverse_variable, -coeff])

        rxn_coeffs.append([met_var, -1])
        rxn_coeffs = dict(rxn_coeffs)
        w_model.constraints.get(met.id).set_linear_coefficients(rxn_coeffs)

    # Add z variables and associated constraints
    for rxn in w_model.reactions:
        rxn_var = prob.Variable("rxn_{}".format(rxn.id))

        const = prob.Constraint(rxn.forward_variable + rxn.reverse_variable - rxn_var,
                                name="rxn_{}".format(rxn.id),
                                lb=0,
                                ub=0
                                )
        rxn_vars_and_cons.extend([rxn_var, const])
        objective += rxn_var**2

    w_model.add_cons_vars(rxn_vars_and_cons)

    # Add objective. since solvers only take linear or quadratic objectives
    # the objective function is left as sum_r(z_r**2)
    w_model.objective = prob.Objective(objective, direction='min')

    # find active reactions for every leak individually
    leak_modes = {}
    for leak, flux in leaks:
        rxns_in_mode = []

        met_var = w_model.variables.get("aux_{}".format(leak))
        met_var.lb = 1
        sol = w_model.optimize()
        met_var.lb = 0

        rxns_in_mode = [[rxn.id, sol.fluxes[rxn.id]] for rxn in w_model.reactions
                        if abs(sol.fluxes[rxn.id]) >= cutoff_mult * zero_cutoff]
        leak_modes[leak] = rxns_in_mode
        print(leak, sol.status, len(rxns_in_mode), met_var.primal)

    return leak_modes
