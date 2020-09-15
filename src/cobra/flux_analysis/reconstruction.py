from cobra.flux_analysis.helpers import normalize_cutoff
from cobra.core import Configuration
from multiprocessing import Pool
import time
CONFIGURATION = Configuration()


def _init_worker(model, carbon_source, type):
    global _model
    global _carbon_source
    global _type
    _model = model
    _carbon_source = carbon_source
    _type = type


def _test_met_parallel(metabolite_id):
    global _model
    global _carbon_source
    global _type

    if _type == 'can_produce':
        lb = 0
        ub = 1000
    elif _type == 'can_consume':
        lb = -1000
        ub = 0

    zero_cutoff = normalize_cutoff(_model, None)
    blocked_met = []
    produced_from_nothing = []
    available_met = []

    with _model as model:
        rxn_id = 'SK_' + metabolite_id
        if rxn_id in model.reactions:
            met_to_test = model.reactions.get_by_id(rxn_id)
        else:
            met_to_test = model.add_boundary(
                model.metabolites.get_by_id(metabolite_id),
                type='sink',
                lb=lb,
                ub=ub
            )
        model.objective = met_to_test
        solution = model.optimize()

        if abs(solution.objective_value) < zero_cutoff:
            blocked_met = [
                metabolite_id
            ]
        else:
            if abs(solution.fluxes[_carbon_source]) <= 0:
                produced_from_nothing = [
                    metabolite_id,
                    solution.objective_value
                ]
            else:
                available_met = [
                    metabolite_id,
                    solution.objective_value,
                    solution.fluxes[_carbon_source]
                ]

    return blocked_met, available_met, produced_from_nothing


def _test_met(model, metabolite_id, carbon_source, type):
    if type == 'can_produce':
        lb = 0
        ub = 1000
    elif type == 'can_consume':
        lb = -1000
        ub = 0

    zero_cutoff = normalize_cutoff(model, None)
    blocked_met = []
    produced_from_nothing = []
    available_met = []

    with model as model:
        rxn_id = 'SK_' + metabolite_id
        if rxn_id in model.reactions:
            met_to_test = model.reactions.get_by_id(rxn_id)
        else:
            met_to_test = model.add_boundary(
                model.metabolites.get_by_id(metabolite_id),
                type='sink',
                lb=lb,
                ub=ub
            )
        model.objective = met_to_test
        solution = model.optimize()

        if abs(solution.objective_value) < zero_cutoff:
            blocked_met = [
                metabolite_id
            ]
        else:
            if abs(solution.fluxes[carbon_source]) <= 0:
                produced_from_nothing = [
                    metabolite_id,
                    solution.objective_value
                ]
            else:
                available_met = [
                    metabolite_id,
                    solution.objective_value,
                    solution.fluxes[carbon_source]
                ]

    return blocked_met, available_met, produced_from_nothing, solution


def find_blocked_mets(model, demands=[], carbon_source='', type='can_produce', parallel=True, n_workers=None):
    if n_workers is None:
        n_workers = CONFIGURATION.processes

    if demands:
        demands_to_test = demands
    else:
        demands_to_test = [
            met.id for met in model.metabolites
        ]

    carbon_source = 'EX_' + carbon_source
    blocked_mets = []
    available_mets = []
    produced_from_nothing = []

    # print('Start testing mets')
    # start_time = time.time()

    if parallel and len(demands_to_test) > 1:
        chunk_size = len(demands_to_test) // n_workers
        with Pool(
            n_workers, initializer=_init_worker, initargs=(model, carbon_source, type)
        ) as pool:
            for bcmet, avmet, pnmet in pool.imap(
                _test_met_parallel, demands_to_test, chunksize=chunk_size
            ):
                if bcmet:
                    blocked_mets.extend(bcmet)
                if avmet:
                    available_mets.append(avmet)
                if pnmet:
                    produced_from_nothing.append(pnmet)

        return blocked_mets, available_mets, produced_from_nothing
    else:
        bcmet, avmet, pnmet, sol = _test_met(
            model, demands_to_test[0], carbon_source, type
        )
        if bcmet:
            blocked_mets.extend(bcmet)
        if avmet:
            available_mets.append(avmet)
        if pnmet:
            produced_from_nothing.append(pnmet)

        return blocked_mets, available_mets, produced_from_nothing, sol
        # stop_time = time.time()
        # print(f'Testing {len(demands_to_test)} metabolites took {stop_time - start_time} secs')


def find_mets_to_connect(model, blocked_mets, carbon_source='', n_workers=None):
    mets_to_connect = []
    start_time = time.time()
    with model as model:
        for met in blocked_mets:
            rxn_id = 'SK_' + met
            if rxn_id not in model.reactions:
                model.add_boundary(
                    model.metabolites.get_by_id(met), type='sink'
                )
            blocked_mets.remove(met)
            demands = blocked_mets
            bmnew, amnew, pnnew = find_blocked_mets(
                model,
                demands=demands,
                carbon_source=carbon_source,
                n_workers=n_workers
            )
            mets_to_connect.append([len(blocked_mets) - len(bmnew) + 1, met])
            for m in blocked_mets:
                if m not in bmnew:
                    blocked_mets.remove(m)
    stop_time = time.time()
    mets_to_connect.sort(reverse=True)
    print(f'Searching mets took {(stop_time - start_time) / 60} min')

    return mets_to_connect


def find_dead_ends(model, carbon_source='', n_workers=None):
    one_mets = []
    for met in model.metabolites:
        if len(met.reactions) == 1 and '_ex' not in met.id:
            one_mets.extend([met.id])

    start_time = time.time()
    print('Searching can_produce metabolites')
    can_produce = []
    bm, av, pn = find_blocked_mets(
        model,
        demands=one_mets,
        carbon_source=carbon_source,
        n_workers=n_workers,
        type='can_produce'
    )
    for met in av:
        can_produce.extend([met[0]])
    for met in pn:
        can_produce.extend([met[0]])

    print('Searching can_consume metabolites')
    can_consume = []
    bm, av, pn = find_blocked_mets(
        model,
        demands=one_mets,
        carbon_source=carbon_source,
        n_workers=n_workers,
        type='can_consume'
    )
    for met in av:
        can_consume.extend([met[0]])
    for met in pn:
        can_consume.extend([met[0]])

    dead_ends = []
    for met in one_mets:
        if met not in can_produce and met not in can_consume:
            dead_ends.extend([met])
    stop_time = time.time()
    print(f'Process took {stop_time - start_time} secs')
    return dead_ends
