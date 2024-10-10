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
    #print(param)

    index,sce, imu, import_set,eff,cov,vac_sce = param
    
    # Run the simulation
    result_path = config['result_path']
    result_matrix, inf_preg_dic = simulation(initial_number,sce, imu,import_set,config,eff, cov, vac_sce)
    
    # Define the filename format
    random_number = random.randint(0, 10000)
    npyfile = f"{result_path}/simulation_{index}_{sce}_{initial_number}_{import_set}_{imu}_{eff}_{cov}_{vac_sce}_{random_number}.npy"
    np.save(npyfile, result_matrix)

    jsonfile = f"{result_path}/simulation_{index}_{sce}_{initial_number}_{import_set}_{imu}_{eff}_{cov}_{vac_sce}_{random_number}.json"

    with open(jsonfile, 'w') as file:
        json.dump(inf_preg_dic, file)


    return 

if __name__ == "__main__":
    # make the modification explicit 
    # python main.py -v=test; single argument does not need 
    parser = argparse.ArgumentParser(description="Run the program with the specified configuration.")
    parser.add_argument('-v', default='real world', choices=['test', 'real world'], help="Specify the configuration to use.\
                         Options are 'test' or 'real world'. Default is 'real world'.")
    
    parser.add_argument('-f', default='no vac', choices=['no vac', 'vac', 'test vac'], help="Specify the configuration to use.\
                         Options are 'no vac', 'vac', 'test vac'. Default is 'no vac'.")
    args = parser.parse_args()

    if args.v == 'test':
        config = test_config
    else:
        config = actual_config

    #vaccination related arguments
    if args.f == 'vac':
        vac_config = actual_vac_config
        config = actual_config
    elif args.f == 'test vac':
        vac_config = test_vac_config
        config = test_config


    # build result folder
    if not os.path.exists(config['result_path']):
        os.makedirs(config['result_path'])
        print(f"Folder created at {config['result_path']}")
    else:
        print(f"Folder already exists at {config['result_path']}")
    
    # unpack config print(config)
    simu_times = config['simu_times']
    num_processors = config['num_processors']

    # if  vaccination, specify... 

    if  args.f == 'no vac':
            params = [(index,sce,imu, import_set) 
              for index in range(1, simu_times+1)
              for sce in Rt_sce
              for imu in herd_imus
              for import_set in import_sets]
        
    else:
        print( 'do the vaccination simulations')
        Rt_sce = ["sustained"]
        import_sets = [(40, 50)]
        
        # vaccination related parameter
        vac_effs = vac_config['VAC_EFFICACY_LIST']
        vac_covs = vac_config['VAC_COVERAGE_LIST']
        vac_sces = vac_config['vac_sces']
        params = [(index,sce,imu, import_set, eff, cov,vac_sce) 
              for index in range(1, simu_times+1)
              for sce in Rt_sce
              for imu in herd_imus
              for import_set in import_sets
              for eff in vac_effs
              for cov in vac_covs
              for vac_sce in  vac_sces
              ]



    
    with Pool(num_processors) as pool:
        # Run simulations in parallel with the generated parameters
        results = pool.starmap(run_simulation, [(param, config) for param in params])

    print("Vaccination Simulation results saved.")
