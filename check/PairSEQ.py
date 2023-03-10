
import numpy as np

def PairSEQ_Simulation(target_sequence_occupancy, n_permutations = 10000, n_wells = 96):
    """
    Simulate null P value distributions
    
    Parameters:
        target_sequence_occupancy (list): A list of target sequences with different well occupancies
        
        n_permutations (integer): Number of permutations to run
        
        n_wells (integer): Number of wells in the experiment

    Returns:
        minimum_p_value (dict): A dictionary with minimum P value for each number of query wells j
    """
    # Initialize an empty dictionary to store minimum P value for each number of query wells j
    minimum_p_value = {}
    for permutation in range(n_permutations):
        # Initialize an empty dictionary to store P value for each number of query wells j in this permutation
        p_value = {}
        for i in range(1, n_wells + 1):
            # Get the number of target sequences that occupy i wells
            t = target_sequence_occupancy[i]
            # Sample Ti random numbers in [0,1] and find the largest number gi
            g = np.max(np.random.uniform(0, 1, size=t))
            for j in range(1, n_wells + 1):
                # Use gi to determine the number of shared wells Nij in a cumulative distribution function for sequences that occupy i and j wells
                nij = np.sum(np.random.uniform(0, g, size=target_sequence_occupancy[j]))
                # Compute a P value dj for occupancy level j
                dj = nij / target_sequence_occupancy[j]
                # If dj is smaller than the smallest P value seen so far at level j, store it
                if j not in p_value or dj < p_value[j]:
                    p_value[j] = dj
        for j in p_value:
            # If dj is smaller than the smallest P value seen so far for number of query wells j, store it
            if j not in minimum_p_value or p_value[j] < minimum_p_value[j]:
                minimum_p_value[j] = p_value[j]
    # Return the minimum P value for each number of query wells j
    return minimum_p_value
