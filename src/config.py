import json
common_config = {

    # Add other shared parameters here
}

# Specific configuration for the test scenario
test_specific_config = {
    'num_year': 10,
    'simu_times': 1,
    'num_processors': 16,
    'num_time_steps': 3650,  # num_year is 10
    'sampling_fraction': 6052713 / 4159788,
    'initial_popu': 6052713,
    'demo_path': '../population_data',
    # Add other real-world specific parameters here
    'result_path': '../test/results_617'
}

# Specific configuration for the real-world scenario
actual_specific_config = {
    'num_year': 10,
    'simu_times': 100,
    'num_processors': 64,
    'num_time_steps': 365*10,  # num_year is 10
    'sampling_fraction': 6052713 / 4159788,
    'initial_popu': 6052713,
    'demo_path': '../population_data',
    # Add other real-world specific parameters here
    'result_path': '../results'
}

# 2024 6052709/  2023 5917653

# Additional shared parameters
year_0 = 2024
year_end = 2033

zika_params = {
    "mu": 1 / 80 / 365,
    "sigma": 1 / 5.9,
    "gamma": 1 / 5.5,
}

mu = 1 / 80 / 365
initial_number = 0
herd_imus = [0, 0.02, 0.05]
import_sets = [(0, 3), (40, 50)] # 
Rt_sce = ["sustained",'dying'] # 

p = zika_params    
f =  (p["gamma"] + p["mu"]) * (p["sigma"] + p["mu"]) / p["sigma"]





# vaccination 
# VAC_EFFICACY_LIST = [0.5,0.8,1.0]
# VAC_COVERAGE_LIST = [0.2,0.5,0.8]

vac_sce = {
    1: (18,100,1),
    2: (20,100,1),
    3: (25,100,1),
    4: (18,5,5),# 42
    5: (20,5,4),#39
    6: (25,5,3),# 39
    7:(27,2,4), # 34
    8:("special pre-pregnancy vaccination ") # 1 year 
}
# 1,2,3,4,5,6 actually no difference 

# Merge common configuration with specific configurations
test_config = {**common_config, **test_specific_config}
actual_config = {**common_config, **actual_specific_config}

actual_vac_config = {
    'VAC_EFFICACY_LIST' : [0.5,0.8,1.0],
    'VAC_COVERAGE_LIST' : [0.2,0.5,0.8],
    'vac_sces' : [1,2,3,4,5,6,7,8]
}


test_vac_config ={
    'VAC_EFFICACY_LIST' : [0.5,0.8],
    'VAC_COVERAGE_LIST' : [0.2,0.5],
    'vac_sces' : [1,4,7,8]#[1,5,7]
}

def load_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)
    

female_age_dic_path = '../population_data/female_age_dic.json'
add_v_all_sce_path = '../population_data/vac_data/add_v_all_sce.json'
delete_v_all_sce_path = '../population_data/vac_data/delete_v_all_sce.json'

female_age_dic =  (load_json(female_age_dic_path))
add_v_all_sce =  (load_json(add_v_all_sce_path))
delete_v_all_sce = ( load_json(delete_v_all_sce_path))


