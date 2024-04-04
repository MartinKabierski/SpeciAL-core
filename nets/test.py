import os
import pm4py
import pm4py.stats

from pm4py.algo.simulation.playout.petri_net import algorithm as simulator

if __name__ == "__main__":
    net, im, fm = pm4py.read_pnml("net1.pnml")
    pm4py.view_petri_net(net, im, fm)
    
    simulated_log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 1})
    variants = pm4py.get_variants(simulated_log, activity_key='concept:name', case_id_key='case:concept:name', timestamp_key='time:timestamp')
    print(variants.keys())
    
    simulated_log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 10})
    variants = pm4py.get_variants(simulated_log, activity_key='concept:name', case_id_key='case:concept:name', timestamp_key='time:timestamp')
    print(variants.keys())
    
    simulated_log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 20})
    variants = pm4py.get_variants(simulated_log, activity_key='concept:name', case_id_key='case:concept:name', timestamp_key='time:timestamp')
    print(variants.keys())
    
    simulated_log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 50})
    variants = pm4py.get_variants(simulated_log, activity_key='concept:name', case_id_key='case:concept:name', timestamp_key='time:timestamp')
    print(variants.keys())
    
    
    simulated_log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 100})
    variants = pm4py.get_variants(simulated_log, activity_key='concept:name', case_id_key='case:concept:name', timestamp_key='time:timestamp')
    print(variants.keys())
    
    
    simulated_log = simulator.apply(net, im, variant=simulator.Variants.BASIC_PLAYOUT, parameters={simulator.Variants.BASIC_PLAYOUT.value.Parameters.NO_TRACES: 1000})
    variants = pm4py.get_variants(simulated_log, activity_key='concept:name', case_id_key='case:concept:name', timestamp_key='time:timestamp')
    print(variants.keys())
    #print(simulated_log[0])
