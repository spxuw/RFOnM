import numpy as np
import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor
import rfonm.utils as utils
from functools import partial
import networkx as nx

def postprocess_theta(theta, adj, gene_name_v1, gene_name_v2, upper_size, thres):
    sin_values = np.sin(theta)
    cos_values = np.cos(theta)

    mean_zz = []
    max_thre = []

    for upper in upper_size:
        mean_z = []
        for t in thres:
            # Count how many genes exceed the threshold in either sin or cos
            counts = np.sum((sin_values > t) | (cos_values > t), axis=1)
            Hc = np.argmin(np.abs(counts - upper))  # Find H closest to target size

            disease_gene = np.where((sin_values[Hc, :] > t) | (cos_values[Hc, :] > t))[0]

            g = nx.from_numpy_array(adj)
            subg = g.subgraph(disease_gene)

            components = list(nx.connected_components(subg))
            if not components:
                mean_z.append(0)
                continue

            largest_component = max(components, key=len)
            node_ids = np.array(list(largest_component))

            # Calculate mean over v1 and v2 for selected genes
            val1 = np.sum(gene_name_v1[node_ids]) / len(node_ids)
            val2 = np.sum(gene_name_v2[node_ids]) / len(node_ids)
            mean_z.append(max(val1, val2))

        mean_zz.append(np.mean(mean_z))
        max_thre.append(thres[np.argmax(mean_z)])

    optimal_size = upper_size[np.argmax(mean_zz)]
    best_thre = max_thre[np.argmax(mean_zz)]
    counts = np.sum((sin_values > best_thre) | (cos_values > best_thre), axis=1)
    Hc = np.argmin(np.abs(counts - optimal_size))

    disease_gene = np.where((sin_values[Hc, :] > best_thre) | (cos_values[Hc, :] > best_thre))[0]
    g = nx.from_numpy_array(adj)
    subg = g.subgraph(disease_gene)

    components = list(nx.connected_components(subg))
    if not components:
        final_genes = np.array([], dtype=int)
    else:
        largest_component = max(components, key=len)
        final_genes = np.array(list(largest_component))

    print("=== SUMMARY ===")
    print(f"target upper_size array: {upper_size}")
    print(f"theta array: {thres}")
    print(f"Selected optimal size={optimal_size}, best_thre={best_thre}")
    print(f"Counts at best_threshold across H: min={counts.min()}, max={counts.max()}")
    print(f"Selected Hc index={Hc} with counts[Hc]={counts[Hc]}")
    print(f"Genes passing threshold before graph={disease_gene.size}")
    print(f"Largest connected component size (final)={final_genes.size}")
    
    return final_genes

def run(max_iterations, n_workers, adj_path, feature_path, top_size, bottom_size, result_path):
    discreated_action = np.linspace(-np.pi,np.pi,num=51)

    adj = np.loadtxt(adj_path, delimiter=',')
    num_nodes = adj.shape
    degrees = np.sum(adj, axis=1)
    row_indices, col_indices = np.nonzero(adj)
    edge_weights = np.mean(degrees) / np.sqrt(degrees[row_indices] * degrees[col_indices])

    external_field = np.loadtxt(feature_path, delimiter=',')
    external_field = np.where(external_field == -np.inf, np.min(external_field[external_field != -np.inf]), external_field)

    HH = np.arange(-10, 10.5, 0.5)

    sa_task_partial = partial(
        sa_task,
        discreated_action=discreated_action,
        num_nodes=num_nodes[0],
        external_field=external_field,
        row_indices=row_indices,
        col_indices=col_indices,
        edge_weights=edge_weights,
        max_iterations=max_iterations
    )

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        state_final = np.array(list(executor.map(sa_task_partial, HH)))

    theta = state_final  # from SA output
    upper_size = list(range(int(bottom_size), int(top_size) + 1, 100))
    thres = np.linspace(0.85, 0.95, 11)
    
    # Output
    final_gene_indices = postprocess_theta(
        theta, adj, external_field[:, 0], external_field[:, 1], upper_size, thres
    )
    
    os.makedirs(result_path, exist_ok=True)
    output_file = os.path.join(result_path, f'final_module.csv')
    np.savetxt(output_file, final_gene_indices, fmt='%d', delimiter=',')

    output_file_state = os.path.join(result_path, f'gene_state.csv')
    np.savetxt(output_file_state, state_final, fmt='%f', delimiter=',')


def sa_task(H_value, discreated_action, num_nodes, external_field, row_indices, col_indices, edge_weights, max_iterations):
    import rfonm.utils as utils
    initial_state = np.random.choice(discreated_action, size=num_nodes)
    best_configuration, _ = utils.simulated_annealing(
        initial_state,
        external_field + H_value,
        row_indices,
        col_indices,
        edge_weights,
        max_iterations,
        discreated_action
    )
    print(f"--- {H_value} external field ---")
    return best_configuration.reshape(-1)

def main():
    parser = argparse.ArgumentParser(description='RFOnM runner')
    parser.add_argument('--max_iterations', type=int, default=300000, help='Maximum number of iterations')
    parser.add_argument('--n_workers', type=int, default=4, help='Maximum number of workers')
    parser.add_argument('--adj_path', required=True, help='Path to adjacency matrix CSV')
    parser.add_argument('--feature_path', required=True, help='Path to node feature matrix CSV')
    parser.add_argument('--top_size', required=True, help='Top size to determine module')
    parser.add_argument('--bottom_size', required=True, help='Bottom size to determine module')
    parser.add_argument('--result_path', required=True, help='Directory to save output CSV')
    args = parser.parse_args()
    run(
        max_iterations = args.max_iterations,
        n_workers = args.n_workers,
        adj_path=args.adj_path,
        feature_path=args.feature_path,
        top_size = args.top_size,
        bottom_size = args.bottom_size,
        result_path=args.result_path
    )
