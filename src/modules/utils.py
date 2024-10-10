
import numpy as np
import random
import sys
sys.path.insert(0, '/Users/haolong/Documents/GitHub/zika_new/zika_project(reformatted)/src')
from no_demo_config import *

def generate_R0(day,sce):
    sigma2 = 0.23
    # Function with sine wave and additional noise N(0, 0.2^2)

    if( sce == "dying"):
        Rt = np.random.uniform(0.7, 0.9)
    if ( sce == "sustained"):
        # R0 = max(0.216,np.random.normal(1.,0.4))
        # R0 = min(R0, 1.78)
        # mu = 1
        seasonal_variation = 1 + 0.064 * np.sin(2 * np.pi * day / 365 + 0.028)
        noise = np.random.normal(0, sigma2)
        Rt = seasonal_variation + noise


    return Rt

def get_beta(p):
    global f
    # if p["sigma"] > 0 or p["mu"] > 0:
    #     f = (p["gamma"] + p["mu"] ) * (p["sigma"] + p["mu"]) / p["sigma"]
    # else:
    #     f = p["gamma"]
    beta = p["R0"] * f
    return beta

def do_exposed(S,E,I,R,date,persons_dict,sce,example_params,IfVac):
    
    num_S = len(S)
    num_I = len(I)
    num_E = len(E)
    num_R = len(R)
    num_N = len(persons_dict)
    
    # get beta
    example_params["R0"] = generate_R0(date,sce)
    beta = get_beta(example_params)
    exposed_set =set()
    exposed_today_num = np.random.binomial(num_S, beta * num_I/ num_N)
    # print(f"Day {date}: Exposed today - {exposed_today_num}")
    # example_params["beta"]
    indices_to_exposed = set(random.choices(tuple(S),k = exposed_today_num))
    # logging.info(f"Day {day}: Exposed today - {exposed_today_num}")
    for index in indices_to_exposed:
        person = persons_dict[index]
        person.contact_with_infected(date)
        # == -1 means that the person is resistent to the infection 
        if(person.expected_exposed_time!= -1):
            exposed_set.add(index)

    return exposed_set

def initial_infected(initial_number,persons_dict,N):

    randomly_chosen_indices =  set(random.choices(tuple(N),k = initial_number))
    for random_index in randomly_chosen_indices:
        persons_dict[random_index].infection_direct(0)
    S = N - randomly_chosen_indices
    E = set()
    I = randomly_chosen_indices
    R = set()

    return  S, E, I, R
