import pandas as pd
import numpy as np
import sklearn.cluster
import networkx as nx
import scipy.stats as sts
import os
import os.path
import tarfile
import requests

__author__ = 'Rasmus Magnusson, Eduard Ghemes'
__COPYRIGHT__ = 'Copyright (C) 2023 Rasmus Magnusson, Eduard Ghemes'
__contact__ = 'rasma774@gmail.com'

PATH = os.path.dirname(os.path.abspath(__file__)).replace('utilities', '')

np.random.seed(0) # Set a seed for reproducibility

def load_nw(
    nw: str, 
    cutoff: float = 0
    ):
    """
    Load network data from a specified source, applying a cutoff filter to the scores.

    Parameters:
    nw : The name of the network to load. Currently supports 'STRINGdb'.
    cutoff [Default = 0]: The minimum score for edges to be included. 

    Returns:
    pandas.DataFrame: Contains the filtered network data.
    """
    if nw == 'STRINGdb':
        tmppath = PATH + 'data/STRINGdb/clean_string.csv'
        # Check if the file exists, if not, extract it from a compressed archive
        if not os.path.isfile(tmppath):
            file = tarfile.open(PATH + 'data/STRINGdb/clean_string.tar.gz')
            file.extractall(PATH + 'data/STRINGdb/', filter=lambda tarinfo, _: tarinfo)
        df = pd.read_csv(tmppath)
        df = df.set_index('p1') # Set the index to the first protein identifier
        df = df[df.score >= cutoff] # Filter rows based on the score cutoff
    return df

def _find_alternative_ids(missing_genes):
    """Attempt to find alternative IDs for missing genes using the STRINGdb API."""
    alternative_identifiers = {}
    no_identifier = []

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    
    for gene in missing_genes:
        params = {
            "identifiers": gene,
            "species": 9606, # Human
            "limit": 1,
            "echo_query": 1,
            "caller_identity": "Best bioinformatician"
        }
        request_url = "/".join([string_api_url, output_format, method])
        results = requests.post(request_url, data=params)
        if results.text.strip():
            for line in results.text.strip().split("\n"):
                l = line.split("\t")
                input_identifier, string_identifier, common_name = l[0], l[2], l[5]
                alternative_identifiers[input_identifier] = (string_identifier, common_name)
        else:
            no_identifier.append(gene)
    
    return alternative_identifiers, no_identifier

def map2nw(
    start_nodes: np.array, 
    nw: pd.DataFrame, 
    augment_factor: int = 2, 
    norm: bool = True):
    """
    Map a list of genes to a network, optionally augmenting the list based on network scores.

    Parameters:
    start_nodes: The initial list of nodes to map.
    nw: The network data as a pandas DataFrame.
    augment_factor [Default = 2]: The factor by which to augment the list of nodes.
    norm [Default = True]: Whether to normalize scores by background values.

    Returns:
    np.array: An array of the augmented list of nodes.
    dict: A report containing information about missing genes, alternative IDs, and added genes.
    """
    # Ensure start_nodes is an np.array
    if isinstance(start_nodes, list):
        start_nodes = np.array(start_nodes)
    elif not isinstance(start_nodes, np.ndarray):
        raise ValueError('start_nodes must be np.array or list')

    # Identify nodes present in the network
    in_nw = np.isin(start_nodes, np.unique(nw.index))
    missing_genes = start_nodes[~in_nw]

    # Attempt to find alternative IDs for missing genes
    alternative_identifiers, no_identifier = _find_alternative_ids(missing_genes)

    # Print summary information about the mapping process
    print(f"{len(missing_genes)} out of {len(start_nodes)} genes not found.")
    if alternative_identifiers:
        print(f"{len(alternative_identifiers)} out of {len(start_nodes)} input IDs were mapped to STRING IDs and alternative IDs.")
    if no_identifier:
        print(f"{len(no_identifier)} out of {len(start_nodes)} have not been found with STRING IDs or alternative IDs.")

    # Filter and augment the list of starting genes based on the network scores
    start_nodes_filtered = start_nodes[in_nw]
    score_sum = nw.loc[start_nodes_filtered].groupby(nw.columns[0]).sum()

    # Normalize scores if requested. Supports only STRINGdb at the moment
    if norm:
        norm_vals = pd.read_csv(PATH + 'data/network_background/STRINGdb_background.csv', index_col=0)
        score_sum = score_sum / norm_vals.reindex(score_sum.index)

    # Sort and filter the scores, then select genes to augment the list
    score_sum = score_sum.sort_values('score', ascending=False)
    score_sum = score_sum[~score_sum.index.isin(start_nodes)]
    # Adjust the number of genes to be added based on the actual number of found genes
    # Does not include the genes that needed alternative IDs or have no known ID
    gene_exp = score_sum.index[:int(augment_factor * len(start_nodes_filtered))]

    print(f"{len(gene_exp)} genes added.")

    # Compile a report of the mapping process
    report = {
        "missing_genes": missing_genes,
        "alternative_IDs": alternative_identifiers,
        "no_IDs": no_identifier,
        "added_genes": gene_exp
    }

    return np.array(gene_exp), report

def cluster_nw(
    genes: np.array, 
    nw: pd.DataFrame, 
    n_clusters: int = None
    ):
    """
    Cluster genes within a network using KMeans clustering.

    Parameters:
    genes: The list of genes to cluster.
    nw: The network data as a pandas DataFrame.
    n_clusters [Default = None]: The number of clusters to form. If None, 'auto' is used.

    Returns:
    np.array: The labels of the clusters.
    index: The index of the genes that were clustered.
    """
    # Filter genes to those present in the network
    genes = genes[np.isin(genes, np.unique(nw.index))]
    subnw = nw.loc[genes]
    subnw = subnw[subnw.iloc[:, 0].isin(genes)]
    subnw = subnw.pivot(columns=nw.columns[0])
    subnw = subnw.fillna(0)

    # Perform KMeans clustering
    if n_clusters is None:
        kmeans = sklearn.cluster.KMeans(n_init='auto')
    else:
        kmeans = sklearn.cluster.KMeans(n_clusters=n_clusters)

    kmeans.fit(subnw.values)
    return kmeans.labels_, subnw.index

def filter_walk(walked_genes, cutoff=2, test=None):
    """NOT USED"""
    # possibly add other tests in future releases
    gene_name, count = np.unique(walked_genes, return_counts=True)
    gene_name = gene_name[count > cutoff]
    return gene_name

def _build_graph(nw):
    """Build a networkx graph from network data. """
    
    print('Building graph in networkx')
    graph = nx.Graph(list(zip(nw.index, nw.iloc[:, 0])))
    nx.connected_components(graph)

    # Extract all connected components
    components = [graph.subgraph(comp).copy() for comp in nx.connected_components(graph)]
    
    # Find the largest component
    biggest_component = max(components, key=len)
    
    # Print the fraction of genes lost by selecting the largest component
    print('Lost ', 1 - (len(biggest_component.nodes)/len(graph.nodes)),
          'of genes by extracting largest connected component of the network')
    
    return biggest_component, components

def _get_dists(
    network: nx.Graph, 
    gene_list1: list, 
    gene_list2: list, 
    n_rand: int = 100):
    """
    Calculate distances between random pairs of genes from two lists in a network.

    This function randomly selects pairs of genes from two lists and calculates
    the shortest path length between each pair in the given network.
    """
    rand_lens = []

    # Generate random positions to select genes from each list
    randpos_1 = np.random.randint(0, len(gene_list1), (n_rand, 1))
    randpos_2 = np.random.randint(0, len(gene_list2), (n_rand, 1))

    # Select the genes at the random positions
    genes1 = np.array(gene_list1)[randpos_1[:, 0]]
    genes2 = np.array(gene_list2)[randpos_2[:, 0]]

    # Iterate over the selected genes to calculate distances
    for i in range(len(randpos_1)):
        # Skip calculation if the selected genes are the same
        if genes1[i] == genes2[i]:
            continue

        # Calculate and append the shortest path length between the genes
        rand_lens.append(nx.shortest_path_length(network,
                                                 source=genes1[i],
                                                 target=genes2[i])
                        )
    return np.array(rand_lens)

def _cluster_dist_until_convergence(
    network: nx.Graph, 
    genes1: list, 
    genes2: list):
    """
    Repeatedly calculate the average shortest pathway length between two lists of genes
    by random sampling until the results converge.

    Convergence is defined as the point where the ratio of the mean of the current distances
    to the mean of the new distances (including the newly calculated ones) is close to 1.
    """

    # Initial calculation of distances
    lens = _get_dists(network, genes1, genes2)
    q = np.inf # Initialize q to infinity for the while loop condition
 
    # Loop until the ratio of means converges to 1
    while not (0.99 < q < 1.01):
        print(q)
        lens_new = _get_dists(network, genes1, genes1) # New distances
        new_lens = np.concatenate((lens, lens_new)) # Old + new distances
        q = lens.mean()/new_lens.mean() # Ratio of means
        lens = new_lens
    return new_lens

def compare_distance(
    cluster_res: dict, 
    nw: pd.DataFrame, 
    max_size: int = 300):
    """
    Compare the average distance between genes within clusters to the background network.

    This function calculates the mean shortest path length between genes within each cluster
    and compares it to the mean shortest path length in the entire network. It uses a Monte Carlo
    method to estimate distances when the cluster size is large (>300).

    Parameters:
    cluster_res: A dictionary where keys are cluster identifiers and values are lists of genes in each cluster.
    nw: The network data as a pandas DataFrame.
    max_size [Default = 300]: The maximum size of a cluster for direct calculation. Clusters larger than this will use Monte Carlo estimation.

    Returns:
    tuple: Contains two np.array: mean_fold and p_mannwhitney. 
        mean_fold is the fold change of the mean distance within clusters compared to the background network. 
        p_mannwhitney is the p-value from a Mann-Whitney U test comparing distances within clusters to the background.
    """
    biggest_component, components = _build_graph(nw.copy())

    # Similar to the previous function, but now we compare the mean distance within clusters
    lens = _get_dists(biggest_component, biggest_component.nodes, biggest_component.nodes)
    q = np.inf
    
    while not (0.99 < q < 1.01):
        lens_new = _get_dists(biggest_component, biggest_component.nodes, biggest_component.nodes)
        new_lens = np.concatenate((lens, lens_new))
        q = lens.mean()/new_lens.mean()
        lens = new_lens
    mean_background = new_lens.mean()

    # Initialize matrices to store results
    keys = list(cluster_res.keys())
    res_mu = np.zeros((len(keys), len(keys)))
    res_mu[:] = np.nan # Initialize with NaNs
    p_mannwhitney = res_mu.copy() # Copy the structure for p-values
    
    # Compare each pair of clusters
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            genes1 = cluster_res[keys[i]][1]
            genes2 = cluster_res[keys[j]][1]

            # Use Monte Carlo estimation for large clusters
            if max(len(genes1), len(genes2)) > max_size:
                lens_tmp = _cluster_dist_until_convergence(biggest_component,
                                                           genes1,
                                                           genes2)
            # Direct calculation for smaller clusters
            else:
                lens_tmp = []
                for g1 in genes1:
                    for g2 in genes2:
                        lens_tmp.append(
                            nx.shortest_path_length(biggest_component,
                                                    source=g1,
                                                    target=g2))

            # Calculate mean distance and p-value
            res_mu[i, j] = np.mean(lens_tmp)
            res_mu[j, i] = res_mu[i, j] # Symmetric matrix
            p_mannwhitney[i, j] =  sts.mannwhitneyu(new_lens, lens_tmp)[1]
            p_mannwhitney[j, i] = p_mannwhitney[i, j]
    
    # Calculate fold change of mean distances
    mean_fold = res_mu/mean_background
    
    return mean_fold, p_mannwhitney
