import numpy as np
import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor
import RFOnM.utils as utils
from functools import partial

def run(max_iterations, n_workers, adj_path, feature_path, result_path):
    discreated_action = np.linspace(0, 2 * np.pi, 10, endpoint=False)

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

    os.makedirs(result_path, exist_ok=True)
    output_file = os.path.join(result_path, f'gene_state.csv')
    np.savetxt(output_file, state_final, delimiter=',')


def sa_task(H_value, discreated_action, num_nodes, external_field, row_indices, col_indices, edge_weights, max_iterations):
    import RFOnM.utils as utils
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
    parser.add_argument('--max_iterations', type=int, default=300001, help='Maximum number of iterations')
    parser.add_argument('--n_workers', type=int, default=4, help='Maximum number of workers')
    parser.add_argument('--adj_path', required=True, help='Path to adjacency matrix CSV')
    parser.add_argument('--feature_path', required=True, help='Path to node feature matrix CSV')
    parser.add_argument('--result_path', required=True, help='Directory to save output CSV')
    args = parser.parse_args()
    run(
        max_iterations = args.max_iterations,
        n_workers = args.n_workers,
        adj_path=args.adj_path,
        feature_path=args.feature_path,
        result_path=args.result_path
    )
