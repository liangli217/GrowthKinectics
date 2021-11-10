import logging
from collections import defaultdict


def run_tool(args):
    logging.debug('Arguments {}'.format(args))
    import data.Patient as Patient
    from .GrowthKineticsEngine import GrowthKineticsEngine

    patient_data = Patient.Patient(indiv_name=args.indiv_id)
    mcmc_trace_cell_abundance, num_itertaions = load_mcmc_trace_abundances(args.abundance_mcmc_trace)
    times_sample, times, wbc_dfd = load_wbc_file(args.wbc_time_file)
    gk_engine = GrowthKineticsEngine(patient_data, wbc_dfd, times, times_sample)
    gk_engine.estimate_growth_rate(mcmc_trace_cell_abundance, n_iter=min(num_itertaions, args.n_iter))

    # Output and visualization
    import output.PhylogicOutput
    phylogicoutput = output.PhylogicOutput.PhylogicOutput()
    phylogicoutput.write_growth_rate_tsv(gk_engine.growth_rates, args.indiv_id)
    phylogicoutput.plot_growth_rates(gk_engine.growth_rates, args.indiv_id)


def load_mcmc_trace_abundances(in_file):
    iterations = set()
    cell_abundance_mcmc_trace = defaultdict(lambda: defaultdict(list))
    with open(in_file, 'r') as reader:
        for line in reader:
            values = line.strip('\n').split('\t')
            if line.startswith('Patient_ID'):
                header = {k: v for v, k in enumerate(values)}
            else:
                sample_id = values[header['Sample_ID']]
                cluster_id = int(values[header['Cluster_ID']])
                abundance = int(values[header['Abundance']])
                iterations.add(values[header['Iteration']])
                cell_abundance_mcmc_trace[sample_id][cluster_id].append(abundance)
    return cell_abundance_mcmc_trace, len(iterations)


### Parse wbc file
def load_wbc_file(file):
    times_sample =[]
    times = []
    wbc =[]
    sample_wbc = []

    with open(file) as wbc_file:
        next(wbc_file)
        for line in wbc_file:
            values = line.strip('\r\n').split('\t')
            wbc.append(float(values[-1]))
            times.append(int(values[1]))
            if line.startswith('JB'):
                times_sample.append(int(values[1]))
                sample_wbc.append(values[0])

    return times_sample, times, wbc