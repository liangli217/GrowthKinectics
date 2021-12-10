import numpy as np
import bisect
from collections import defaultdict
from scipy.stats import linregress
from scipy import interpolate


class GrowthKineticsEngine:

    def __init__(self, patient, wbc, times, times_sample, sample_wbc):
        self._patient = patient
        self._wbc = wbc
        self._growth_rates = defaultdict(list)
        self.times = times
        self.times_sample = times_sample
        ## convert days to years
        self.times_in_year = [x / 365 for x in self.times]
        self.times_sample_in_year = [x / 365 for x in self.times_sample]

        self.sample_wbc = sample_wbc

    @property
    def growth_rates(self):
        return self._growth_rates

    @property
    def wbc(self):
        return self._wbc

    def estimate_growth_rate(self, mcmc_trace_cell_abundance, n_iter=100, conv=1e-4):
        '''
        
        '''
        # Number of samples for this Patient
        # sample_list = list(mcmc_trace_cell_abundance.keys())
        sample_list = self.sample_wbc
        time_points = len(sample_list)
        # n_clusters = len(list(mcmc_trace_cell_abundance[sample_list[0]]))
        #cluster_rates = defaultdict(list)
        # If times of samples are not provided, treat as an integers
        if not self.times:
            times = np.array(range(time_points)) + 1
        for n in range(n_iter):
            adj_wbc = self._wbc * (1 + np.array([(np.random.random() - 0.5) / 100. for x in range(len(self._wbc))]))
            for cluster_id in list(mcmc_trace_cell_abundance[sample_list[0]]):
                cluster_abundances = []
                ## iterate through the samples in the wbc file to make sure the order is correct
                for sample_name in sample_list:
                    sample_abundances = mcmc_trace_cell_abundance[sample_name]
                    cluster_abundances.append(float(sample_abundances[cluster_id][n] ))

                # for sample_name, sample_abundances in mcmc_trace_cell_abundance.items():
                #     if sample_name in sample_list:
                #         print(sample_abundances[cluster_id][n])
                #         cluster_abundances.append(sample_abundances[cluster_id][n] + conv)

                ## Interpolation of cluster_abundances
                scaled_cluster_abundances = [i /100+ conv for i in cluster_abundances]

                interpolate_func = interpolate.interp1d(self.times_sample_in_year, scaled_cluster_abundances)
                cluster_abundances_interpolate = interpolate_func(self.times_in_year)
                ## added the log transformationq
                cluster_slope = linregress(self.times_in_year, np.log(cluster_abundances_interpolate * adj_wbc)).slope
                adj_slope = (np.exp(cluster_slope) - 1)
                if cluster_id ==6:
                    print adj_slope
                ## adjust the slope calculation

                self._growth_rates[cluster_id].append(adj_slope)

    def line_fit(self, x, c_idx, fb_x_vals, len_pre_tp, adj_dens):
        """ """
        slope, intercept = x
        y_domain = [np.log(self.grid * self._wbc[tp_idx] + 1e-40) for tp_idx in range(len_pre_tp)]
        y_weights = [adj_dens[tp_idx][c_idx] for tp_idx in range(len_pre_tp)]
        line_y_vals = slope * fb_x_vals + intercept

        selected_weight = [
            min(
                sum(y_weights[tp_idx][:min(bisect.bisect(y_domain[tp_idx], line_y_vals[tp_idx]), 100) + 1]),
                sum(y_weights[tp_idx][min(bisect.bisect(y_domain[tp_idx], line_y_vals[tp_idx]), 100):])
            )
            for tp_idx in range(len_pre_tp)]

        return selected_weight

    def line_fit_err(self, x, c_idx, wbc, fb_x_vals, len_pre_tp, adj_dens):
        slope, intercept = x
        grid = np.arange(101)
        y_domain = [np.log(grid * wbc[tp_idx] + 1e-40) for tp_idx in range(len_pre_tp)]
        y_weights = [adj_dens[tp_idx][c_idx] for tp_idx in range(len_pre_tp)]

        line_y_vals = slope * fb_x_vals + intercept

        selected_weight = []
        for tp_idx in range(len_pre_tp):
            sum0 = sum(y_weights[tp_idx][min(bisect.bisect(y_domain[tp_idx], line_y_vals[tp_idx]), 100):])
            sum1 = sum(y_weights[tp_idx][:min(bisect.bisect(y_domain[tp_idx], line_y_vals[tp_idx]), 100) + 1])
            selected_weight.append(min(sum0, sum1))

        return -sum(selected_weight)

    def line_fit_pval(self, x, c_idx, wbc, fb_x_vals, len_pre_tp, adj_dens):
        slope, intercept = x
        grid = np.arange(101)
        y_domain = [np.log(grid * wbc[tp_idx] + 1e-40) for tp_idx in range(len_pre_tp)]
        y_weights = [adj_dens[tp_idx][c_idx] for tp_idx in range(len_pre_tp)]

        line_y_vals = slope * fb_x_vals + intercept

        selected_weight = [
            min(
                sum(y_weights[tp_idx][:min(bisect.bisect(y_domain[tp_idx], line_y_vals[tp_idx]), 100) + 1]),
                sum(y_weights[tp_idx][min(bisect.bisect(y_domain[tp_idx], line_y_vals[tp_idx]), 100):])
            )
            for tp_idx in range(len_pre_tp)]

        return min(selected_weight)
