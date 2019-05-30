from scipy import stats
from random import shuffle
import collections
import numpy as np
import itertools
import operator
import logging

# add as command line parameter
np.random.seed()

logging.basicConfig(filename='cell_population_engine.log',
                    filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=getattr(logging, "INFO"))


class CellPopulationEngine:

    def __init__(self, patient):
        self._patient = patient
        if patient.ClusteringResults:
            self._clustering_results = patient.ClusteringResults
        else:
            logging.error('Clustering results are not found, need to run Clustering module prior to Cell Population')
        if patient.TopTree:
            self._top_tree = patient.TopTree
        else:
            logging.error('Build Tree results are not found, need to run BuildTree prior to running Cell Population')

    def get_random_cluster(self):
        return np.random.choice(self._clusters)

    @staticmethod
    def sample_ccf(xk, pk):
        """

        :param xk:
        :param pk:
        :return:
        """
        if sum(pk) == 0:
            return 0
        else:
            custm = stats.rv_discrete(name='custm', values=(xk, pk))
            return custm.rvs(size=1)[0]

    @staticmethod
    def normalized(distribution):
        """

        :param distribution:
        :return:
        """
        total_sum = float(sum(distribution))
        if total_sum > 0.0:
            return [d / total_sum for d in distribution]
        else:
            return distribution

    def compute_node_constrained_distribution(self, node, cluster_ccf, iter_ccf, hist):
        parent = node.parent
        if parent:
            if parent.identifier in iter_ccf:
                parent_ccf = iter_ccf[parent.identifier]
                siblings_total_ccf = sum([iter_ccf[sibling] for sibling in node.siblings if sibling in iter_ccf])
                leftover_ccf = int(parent_ccf - siblings_total_ccf)
                logging.debug('Node {} has parent {} with ccf {}'.format(node.identifier,
                                                                         node.parent.identifier,
                                                                         parent_ccf))
                logging.debug('Node {} has siblings {} with total ccf {}'.format(node.identifier,
                                                                                 node.siblings,
                                                                                 siblings_total_ccf))
                logging.debug('Node {} has leftover ccf {}'.format(node.identifier,
                                                                   leftover_ccf))
                constrained_ccf = list(cluster_ccf[:leftover_ccf + 1]) + [0.0] * (100 - leftover_ccf)
                return self.normalized(constrained_ccf)
            else:
                logging.error('Parent {} ccf was not assigned'.format(parent.identifier))
                return None
        else:
            return cluster_ccf

    def sample_node_ccf(self, node, cluster_ccf, iter_ccf, hist):
        constrained_cluster_distribution = self.compute_node_constrained_distribution(node, cluster_ccf, iter_ccf, hist)
        logging.debug(
            'Constrained distribution for node {} is \n{}'.format(node.identifier, constrained_cluster_distribution))
        if constrained_cluster_distribution is not None:
            return self.sample_ccf(hist, constrained_cluster_distribution)
        else:
            logging.warn('Constrained ccf for node {} is None'.format(node.identifier))
            return 0.0

    @staticmethod
    def _get_most_frequent_configuration(constrained_ccfs_configs):
        sorted_config = []
        for d in constrained_ccfs_configs:
            sorted_config.append(tuple(sorted(d.items(), key=operator.itemgetter(1))))
        counts = collections.Counter(sorted_config)
        most_frequent = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)[0]
        return sorted(most_frequent[0], key=operator.itemgetter(0)), most_frequent[1]

    def clusters_constrained_ccf(self, sample_clusters_ccf, tree_levels, n_iter=100):
        '''

        Args:
            sample_clusters_ccf:
            tree_levels:
            n_iter:
        Returns:
        '''
        hist = range(101)
        all_configurations = []
        for i in range(n_iter):
            logging.debug('Iteration {}'.format(i))
            iter_ccf = dict.fromkeys(itertools.chain(self._top_tree.nodes.keys()), 0.0)
            # Traverse tree from root to it's leaves
            for level in tree_levels:
                level_nodes = tree_levels[level]
                shuffle(level_nodes)
                # For each node in the level
                for node_id in level_nodes:
                    node = self._top_tree.nodes[node_id]
                    node_ccf = self.sample_node_ccf(node, sample_clusters_ccf[node_id], iter_ccf, hist)
                    logging.debug('Node {} has constrained ccf {}'.format(node_id, node_ccf))
                    iter_ccf[node_id] = node_ccf
                all_configurations.append(iter_ccf)
        most_frequent_config, most_frequent_count = self._get_most_frequent_configuration(all_configurations)
        logging.debug('Most frequent constrained ccf configuration with count {} \n{} '.format(most_frequent_count,
                                                                                               most_frequent_config))
        return most_frequent_config

    @staticmethod
    def _get_sample_clusters_densities(sample_id, clusters_ccf):
        sample_cluster_densities = {}
        for c in clusters_ccf:
            sample_cluster_densities[c] = clusters_ccf[c].densities[sample_id]
        return sample_cluster_densities

    def samples_average_constrained_ccf(self, n_iter=10000):
        """
        For each sample, iterates over clusters and computes average constrained ccf for that cluster
        :param n_iter:
        :return:
        """
        tree_levels = self._top_tree.get_tree_levels()
        logging.debug('Loaded top tree with edges {}'.format(self._top_tree.edges))
        clusters_ccf = self._clustering_results.clusters
        sample_names = self._clustering_results.samples
        logging.debug('Samples for this patient {}'.format(sample_names))
        # iterate over samples to get average cell abundance in each sample
        sample_constrained_ccf = {sample_name: [] for sample_name in sample_names}
        for sample_id in sample_names:
            sample_clusters_ccf = self._get_sample_clusters_densities(sample_id, clusters_ccf)
            sample_constrained_ccf[sample_id] = self.clusters_constrained_ccf(sample_clusters_ccf, tree_levels,
                                                                              n_iter)
        return sample_constrained_ccf

    def get_cell_abundance(self, constrained_ccf):
        """

        :param constrained_ccf:
        :return:
        """
        cell_abundances = {}
        for sample_id, clusters_constrained_ccf in constrained_ccf.items():
            sample_cell_abundance = {key: 0 for key, value in clusters_constrained_ccf}
            for node_id, node_ccf in clusters_constrained_ccf:
                node = self._top_tree.nodes[node_id]
                parent = node.parent
                if parent:
                    logging.debug('Node {} has parent {} with abundance {}'.format(node_id, parent.identifier,
                                                                                   sample_cell_abundance[
                                                                                       parent.identifier]))
                    sample_cell_abundance[parent.identifier] -= node_ccf
                else:
                    logging.debug('Node {} has no parent'.format(node_id))
                sample_cell_abundance[node_id] += node_ccf
            # check that cancer cell population in the sample sums up to 100%
            assert sum([a for cl, a in sample_cell_abundance.items()]) <= 100.0
            cell_abundances[sample_id] = sample_cell_abundance
            logging.debug('Cell abundance for sample {} \n {}'.format(sample_id, sample_cell_abundance))
        return cell_abundances