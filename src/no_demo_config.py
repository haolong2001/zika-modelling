
# year_0 = 2024 
# year_end = 2038

# modify the code here 
num_year = 1
simu_times = 10
num_processors = 10
num_time_steps = 50 * num_year

example_params = {
    "mu": 1 / 80 /365,
    # "omega": 1 / (365*2),
    "sigma": 1 / 5.9,
    "gamma": 1 / 5.5,
}
mu = 1/80/365
initial_number = 10
import_sets = [(0,3)]

Rt_sce = ["sustained"]
### get_beta(p): "dying"
p = example_params
f = (p["gamma"] + p["mu"] ) * (p["sigma"] + p["mu"]) / p["sigma"]



