import os
import sys
current_directory = os.getcwd()
sys.path.insert(0, current_directory)
 
from no_demo_simu.simulation import simulation
from multiprocessing import Pool
import numpy as np
from no_demo_config import Rt_sce,simu_times,initial_number,import_sets,num_processors
import random

def run_simulation(params):
    index,sce, import_set = params
    
    # Run the simulation
    result_matrix = simulation(initial_number,sce, import_set)
    
    # Define the filename format
    random_number = random.randint(0, 10000)
    filename = f"../test_results/simulation_{index}_{sce}_{initial_number}_{import_set}_{random_number}.npy"
        
    # Save the result matrix to .npy file
    np.save(filename, result_matrix)

    return filename

if __name__ == "__main__":
    # Generate parameter combinations
    print("start no demo simulation.")
    params = [(index,sce, import_set) 
              for index in range(1, simu_times+1)
              for sce in Rt_sce
              for import_set in import_sets]
    
    with Pool(num_processors) as pool:
        # Run simulations in parallel with the generated parameters
        results = pool.map(run_simulation, params)

    print("Simulation results saved.")
