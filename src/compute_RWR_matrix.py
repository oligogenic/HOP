"""Compute the transition matrix for the RWR algorithms from the BOCK graph"""
import os
import pickle
from graph_tool.all import remove_parallel_edges, remove_self_loops, load_graph, adjacency
import configparser

config = configparser.ConfigParser()
config.read('config.ini')
PATH_TO_DATA = config['BASE']['PATH_TO_DATA']
PATH_TO_BOCK = config['BASE']['PATH_TO_BOCK']

def load_and_filter_graph(graph_file):
    """Load the BOCK graph from the graphml format and removes the OligogenicInteraction and Disease nodes"""
    g = load_graph(graph_file)
    filtered_vertices = g.new_vertex_property("bool")
    filtered_edges = g.new_edge_property("bool")
    # Filter nodes
    unique_vertex_labels = []
    for v in g.iter_vertices():
        unique_vertex_labels.append(g.vp.labels[v])
        filtered_vertices[v] = (g.vp.labels[v] not in [':Disease', ':OligogenicCombination'] and g.get_out_degrees([v])[0] > 0)

    # Filter edges
    unique_edge_labels = []
    for e in g.edges():
        unique_edge_labels.append(g.ep.label[e])
        filtered_edges[e] = g.ep.label[e] not in ['involves', 'causes']
    g.set_edge_filter(filtered_edges)
    g.set_vertex_filter(filtered_vertices)
    # Remove parallel edges and self loops
    remove_parallel_edges(g)
    remove_self_loops(g)
    return g

def get_adjacency_matrix(graph, property_map):
    adj_matrix = adjacency(graph, weight=None, vindex=property_map)
    return adj_matrix


def get_instances_to_index(g):
    """
    Get mapping between id of the instance in the graph and the row number in the adjacency matrix
    :param g: The graph
    :return: the graph,
    """
    instances_to_index = {}
    new_node_indices = g.new_vertex_property('int')
    count = 0
    for v in g.iter_vertices():
        if g.vp.id[v] != '':
            v_id = g.vp.id[v]
        elif g.vp.name[v] != '':
            v_id = g.vp.name[v]
        new_node_indices[v] = count
        instances_to_index[v_id] = count
        count += 1
    return g, new_node_indices, instances_to_index


def remove_duplicates(matrix):
    """Converts all non zero entries of a matrix to ones.
    Allows to remove the additions when there are several links between two nodes.w"""
    indices = matrix.nonzero()
    for i,j in zip(indices[0], indices[1]):
            matrix[i,j] = 1
    return matrix


def compute_adjacency_matrix(g):
    """Compute adjacency matrix and get the mapping between the instances in the graph and the index in the matrix"""
    g, new_node_indices, instances_to_index = get_instances_to_index(g)
    adj_matrix = get_adjacency_matrix(g, new_node_indices)
    adj_matrix = remove_duplicates(adj_matrix)
    return adj_matrix, instances_to_index


def get_all_gene_names(g):
    """Get unique gene names"""
    all_genes ={}
    ensg_to_genename = {}
    for v in g.iter_vertices():
        if g.vp.labels[v] == ':Gene' and g.vp.id[v] != '':
            v_id = g.vp.id[v]
            ensg_to_genename[v_id] = g.vp.name[v]
            all_genes[v_id] = 1
    return all_genes, ensg_to_genename


def compute_rw_matrix():
    """Create the transition matrix that will be used by the RWR algorithm and a dictionary mapping each
    node in the graph to the index of the row in the matrix it corresponds to"""
    g = load_and_filter_graph(PATH_TO_BOCK)
    all_genes, ensg_to_geneName = get_all_gene_names(g)
    #print(len(all_genes.keys()))
    os.makedirs(PATH_TO_DATA + 'rwr_data')
    with open(PATH_TO_DATA + 'rwr_data/all_genes.p', 'wb') as file:
        pickle.dump(all_genes, file)
    with open(PATH_TO_DATA + 'rwr_data/ensg_to_geneName.p', 'wb') as file:
        pickle.dump(ensg_to_geneName, file)
    matrix_whole_kg, instances_to_index = compute_adjacency_matrix(g)

    with open(PATH_TO_DATA + 'rwr_data/final_wholeKG_matrix.p', 'wb') as file:
        pickle.dump(matrix_whole_kg, file)
    with open(PATH_TO_DATA + 'rwr_data/final_instances_to_index.p', 'wb') as file:
        pickle.dump(instances_to_index, file)
    return None


