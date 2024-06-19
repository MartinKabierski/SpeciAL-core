from pm4py.algo.simulation.playout.petri_net import algorithm as simulator
from tqdm import tqdm


def simulate_model(net, im, fm, repetitions=200, traces=5000):
    logs = []
    for i in tqdm(range(repetitions), desc='Simulating model'):
        log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: traces})
        #ensure no zero-length traces are in final log
        for j in range(len(log)):
            tr = log[j]
            if len(tr)==0:
                length = 0
                while length == 0:
                    tr = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 1})[0]
                    length = len(tr)
                log[j] = tr
        logs.append(log)
    return logs
