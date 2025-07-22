import numpy as np

def simulated_annealing(initial_state, external_field, row_indices, col_indices, edge_weights, max_iterations, discreated_action, start_temp=10.0, end_temp=0.01, alpha=0.99):
    current_state = np.array(initial_state)
    current_energy = compute_energy(current_state, external_field, row_indices, col_indices, edge_weights)
    best_state = np.array(current_state)
    best_energy = current_energy
    energies = []
    temperature = start_temp

    for iteration in range(max_iterations):
        node_index = np.random.randint(0, len(current_state))
        proposed_state = np.random.choice(discreated_action)
        new_state = np.array(current_state)
        new_state[node_index] = proposed_state
        new_energy = compute_energy(new_state, external_field, row_indices, col_indices, edge_weights)
        delta_energy = new_energy - current_energy
        if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy / temperature):
            current_state, current_energy = new_state, new_energy
            if current_energy < best_energy:
                best_state, best_energy = np.array(current_state), current_energy
        if iteration % 100 == 0:
            energies.append(best_energy)
        temperature = max(end_temp, temperature * alpha)

    return best_state, energies

def compute_energy(state, external_field, row_indices, col_indices, edge_weights):
    interaction_terms = np.cos(state[row_indices] - state[col_indices])
    energy_interaction = -np.sum(edge_weights * interaction_terms)
    state_vectors = np.array([np.cos(state), np.sin(state)]).T
    energy_external = -np.sum(np.einsum('ij,ij->i', external_field, state_vectors))
    return energy_interaction + energy_external
