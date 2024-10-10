import os
import sys
current_directory = os.getcwd()
sys.path.insert(0, current_directory)
import argparse
from simulation import simulation
from multiprocessing import Pool
import numpy as np
from config import *
import random
import json


def run_simulation(param,config):
    index,sce, imu, import_set = param
    
    # Run the simulation
    result_path = config['result_path']
    result_matrix, inf_preg_dic = simulation(initial_number,sce, imu,import_set,config,ifVac=False)
    
    # Define the filename format
    random_number = random.randint(0, 10000)
    npyfile = f"{result_path}/simulation_{index}_{sce}_{initial_number}_{import_set}_{imu}_{random_number}.npy"
    np.save(npyfile, result_matrix)

    jsonfile = f"{result_path}/simulation_{index}_{sce}_{initial_number}_{import_set}_{imu}_{random_number}.json"

    with open(jsonfile, 'w') as file:
        json.dump(inf_preg_dic, file)


    return 

if __name__ == "__main__":
    # make the modification explicit 
    # python main.py -v=test; single argument does not need 
    parser = argparse.ArgumentParser(description="Run the program with the specified configuration.")
    parser.add_argument('-v', default='real world', choices=['test', 'real world'], help="Specify the configuration to use.\
                         Options are 'test' or 'real world'. Default is 'real world'.")
    args = parser.parse_args()

    if args.v == 'test':
        config = test_config
    else:
        config = actual_config

    # build result folder
    if not os.path.exists(config['result_path']):
        os.makedirs(config['result_path'])
        print(f"Folder created at {config['result_path']}")
    else:
        print(f"Folder already exists at {config['result_path']}")
    
    # unpack config print(config)
    simu_times = config['simu_times']
    num_processors = config['num_processors']

    params = [(index,sce,imu, import_set) 
              for index in range(1, simu_times+1)
              for sce in Rt_sce
              for imu in herd_imus
              for import_set in import_sets]
    
    with Pool(num_processors) as pool:
        # Run simulations in parallel with the generated parameters
        results = pool.starmap(run_simulation, [(param, config) for param in params])

    print("Simulation results saved.")
