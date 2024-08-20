"""
Functions for building a residue interaction graph and finding paths
between specific residues
"""

import networkx as nx
import numpy as np

from copy import deepcopy
from warnings import warn

from Bio.PDB import PDBParser

############################
# Building distance matrix #
############################

def calc_intra_dist_matrix(chain):
    """
    Given one protein chain from a PDB file, calculate a matrix containing all
    pairwise Euclidean distances between residues within the chain.

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object

    Returns
    -------
    dist_mtx: array-like, matrix containing pairwise Euclidean distances between
              residues, of dimensions (chain, chain)
    """

    dist_mtx = np.zeros((len(chain), len(chain)))
    for i, residue_a in enumerate(chain):
        for j, residue_b in enumerate(chain):
            if i != j: # Skip measuring distances to self
                dist_mtx[i, j] = calc_residue_distance(residue_a, residue_b)
    return dist_mtx

def calc_residue_distance(residue_one, residue_two):
    """
    Given two residues, calculate the Euclidean distance between them.
    The distance is measured between the beta carbons (or, in the case of
    glycine, with respect to the alpha carbon).

    Arguments
    ---------
    residue_one: Bio.PDB.Residue.Residue object
    residue_two: Bio.PDB.Residue.Residue object

    Returns
    -------
    euclid_dist:  float, Euclidean distance between the two residues
    """
    is_one_glycine = (residue_one.resname == 'GLY')
    is_two_glycine = (residue_two.resname == 'GLY')

    try:
        if is_one_glycine and is_two_glycine:
            diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
        elif is_one_glycine and not is_two_glycine:
            diff_vector = residue_one["CA"].coord - residue_two["CB"].coord
        elif not is_one_glycine and is_two_glycine:
            diff_vector = residue_one["CB"].coord - residue_two["CA"].coord
        else:
            diff_vector = residue_one["CB"].coord - residue_two["CB"].coord
    except KeyError: # Missing atom
        euclid_dist = np.nan
        return euclid_dist

    euclid_dist = np.sqrt(np.sum(diff_vector * diff_vector))
    return euclid_dist

############################
# Graph building functions #
############################

def calc_intra_dist_matrix(chain):
    """
    Given one protein chain from a PDB file, calculate a matrix containing all
    pairwise Euclidean distances between residues within the chain.

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object

    Returns
    -------
    dist_mtx: array-like, matrix containing pairwise Euclidean distances between
              residues, of dimensions (chain, chain)
    """

    dist_mtx = np.zeros((len(chain), len(chain)))
    for i, residue_a in enumerate(chain):
        for j, residue_b in enumerate(chain):
            if i != j: # Skip measuring distances to self
                dist_mtx[i, j] = calc_residue_distance(residue_a, residue_b)
    return dist_mtx



def distance_matrix_to_graph(distance_matrix, residue_idxs, threshold=8):
    """
    Uses a protein structure distance matrix to create a residue interaction network,
    where two residues are connected by an edge if the distance between them is below
    a specified threshold.
    
    Arguments
    ---------
    distance_matrix: array, square matrix representing pairwise distances 
                     between residues
    residue_idxs:    list, residue indexes corresponding to the protein structure
    threshold:       float, the distance threshold below which an edge
                     connects two residues
    
    Returns
    -------
    G: networkx.Graph, residue interaction network
    """
    # Initialize graph
    G = nx.Graph()
    num_residues = len(residue_idxs)
    G.add_nodes_from(residue_idxs)
    
    # Iterate over upper triangle of the distance matrix
    for i in range(num_residues):
        for j in range(i + 1, num_residues):
            if distance_matrix[i, j] < threshold:
                G.add_edge(residue_idxs[i], residue_idxs[j])    
    return G


def trim_graph_to_manifold(G, manifold_residue_idxs):
    """
    Remove nodes from the graph that 
    
    Arguments
    ---------
    G: networkx Graph, residue interaction network
    residue_idxs: list, residue indexes corresponding to the protein structure
    manifold_residue_idxs: list, residue indexes within the high strain manifold
    
    Returns
    -------
    G: networkx.Graph, residue interaction network    
    """
    valid_residues = set(manifold_residue_idxs)
    nodes_to_remove = [node for node in G.nodes if node not in valid_residues]
    trimmed_G = deepcopy(G)
    # Automatically removes adjacent edges too
    trimmed_G.remove_nodes_from(nodes_to_remove)
    return trimmed_G

def get_number_of_connected_components(G):
    no_connected_components = nx.number_connected_components(G)
    print(f"{no_connected_components} connected components in the graph")
    return no_connected_components

##################
# Path functions #
##################

def find_path_in_graph(G, source, destination):
    """
    Finds if there is a path in the graph G between the source and destination nodes.
    
    Arguments
    ---------
    G: networkx.Graph, the graph to search
    source: node, the source node to start the path search
    destination: node, the destination node to end the path search
    
    Returns
    -------
    path_exists: bool, True if a path exists between source and destination,
                 False otherwise.
    """
    # Check if both source and destination are in the graph
    if source not in G:
        warn(f'Source node {source} not found in the graph')
        return False
    if destination not in G:
        warn(f'Destination node {destination} not found in the graph')
        return False
    
    # Check if there is a path between source and destination
    if nx.has_path(G, source, destination):
        return True
    else:
        return False


if __name__ == "__main__":

    parser = PDBParser(QUIET=True)
    pdb_path = "../data/processed/pdb_pairs/chains_by_protein/Q16539/6hwv_A.pdb"
    structure = parser.get_structure("6hwv_A",pdb_path)
    chain = structure[0]["A"]

    distance_matrix = calc_intra_dist_matrix(chain)
    # Get residue indexes directly from the structure
    residue_indexes = [res._id[1] for res in chain]
    G = distance_matrix_to_graph(distance_matrix, residue_indexes)

    # Random bunch of residues to stand in as the manifold
    valid_residues = list(range(240,260))
    trimmed_G = trim_graph_to_manifold(G, valid_residues)

    path_1 = find_path_in_graph(trimmed_G, 240, 300)
    path_2 = find_path_in_graph(trimmed_G, 1, 250)
    path_3 = find_path_in_graph(trimmed_G, 250, 500)

    print('Et voila!')