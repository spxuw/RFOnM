import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.sparse import coo_matrix
import pdb
import networkx as nx
import argparse
import os.path
import time
import random
from matplotlib.colors import ListedColormap
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor



parser = argparse.ArgumentParser(description='RFOnM')
parser.add_argument('--dataset', default=None, help='dataset')

args = parser.parse_args()
dataset = (args.dataset)

discreated_action = np.linspace(-np.pi,np.pi,num=51)


def simulated_annealing(initial_state, external_field, max_iterations=2000001, start_temp=10.0, end_temp=0.01, alpha=0.99):
    current_state = np.array(initial_state)
    current_energy = compute_energy(current_state, external_field)
    best_state = np.array(current_state)
    best_energy = current_energy
    energies = []
    temperature = start_temp
    for iteration in range(max_iterations):
        node_index = np.random.randint(0, len(current_state))
        proposed_state = np.random.choice(discreated_action)
        new_state = np.array(current_state)
        new_state[node_index] = proposed_state
        new_energy = compute_energy(new_state, external_field)
        delta_energy = new_energy - current_energy
        if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy / temperature):
            current_state, current_energy = new_state, new_energy
            if current_energy < best_energy:
                best_state, best_energy = np.array(current_state), current_energy

        if iteration % 100 ==0:
            energies.append(best_energy)
        temperature = max(end_temp, temperature * alpha)
    return best_state, energies

def compute_energy(state, external_field):
    J = 1  # Assuming J = 1 for simplicity
    energy_interaction = -J * np.sum(np.cos(state[row_indices] - state[col_indices]))
    state_vectors = np.array([np.cos(state), np.sin(state)]).T
    energy_external = -np.sum(np.einsum('ij,ij->i', external_field, state_vectors))
    return energy_interaction + energy_external


# Main code
adj = np.loadtxt('../data/GEO/graph_'+str(dataset)+'.csv',delimiter=',')
num_nodes = adj.shape

row_indices, col_indices = np.nonzero(adj)

external_field = np.loadtxt('../data/GEO/features_'+str(dataset)+'.csv',delimiter=',')  # Load external field values for each node from input file
external_field = np.where(external_field == -np.inf, np.min(external_field[external_field!=-np.inf]), external_field)

external_field[:,0] = external_field[:,0] * np.sum(adj, axis=1)
external_field[:,1] = external_field[:,1] * np.sum(adj, axis=1)

HH = np.arange(-10, 10.5, 0.5)
state_final = np.zeros((num_nodes[0],len(HH)))

def sa_task(H_value):
    start_time = time.time()
    initial_state = np.random.choice(discreated_action, size=num_nodes[0],)
    best_configuration, sa_energies = simulated_annealing(initial_state, external_field + H_value)
    print("--- %s external field ---" % H_value)

    # Plot energies
    plot_energies(sa_energies)

    return best_configuration.reshape(-1)

# Use multiprocessing to parallelize the for loop
with ProcessPoolExecutor(max_workers=5) as executor:
    state_final = np.array(list(executor.map(sa_task, HH)))

np.savetxt('../results/GEO_v4/State_'+str(dataset)+'.csv', state_final, delimiter=',')


