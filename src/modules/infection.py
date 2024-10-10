

from collections import defaultdict 
import numpy as np
import random
import os
import sys
import random
current_directory = os.path.dirname(os.path.abspath(__file__))  
# Get the parent directory
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)
# Insert the parent directory into sys.path
from modules.Person import Person
from config import *



def generate_new_Rt(day, sce):
    


def generate_Rt(day,sce):
    
  
    # Function with sine wave and additional noise N(0, 0.2^2)

    if( sce == "dying"):
        Rt = 0.8 + 0.1 * np.sin(2 * np.pi * day / 365 - 0.18)
   



    if ( sce == "sustained"):
        # R0 = max(0.216,np.random.normal(1.,0.4))
        # R0 = min(R0, 1.78)
        # mu = 1
        # sigma2 = 0.19
        # seasonal_variation =  0.071 * np.sin(2 * np.pi * day / 365 -0.033)

        # seasonal_variation =  0.064 * np.sin(2 * np.pi * day / 365 + 0.028)
        # noise = np.random.normal(0, 0.23)
        # Rt = 1 + seasonal_variation + noise

        Rt = max(0.2,Rt)
        Rt = min(Rt, 2.0)
        sigma2 = 0.264
        seasonal_variation = 0.072 * np.sin(2 * np.pi * day / 365 - 0.18)
        noise = np.random.normal(0, sigma2)
        Rt = 1 + seasonal_variation + noise

        Rt = max(0.2,Rt) # make sure it's larger than 0
        Rt = min(Rt,1.8)
        


   


    return Rt



def do_exposed(S,I, num_N, date,infected_dic,sce,IfVac):
    '''
    sce: scenarios;
    params: store parameters  
    '''
    global f 
    global Rt_list
    num_S = len(S)
    num_I = len(I)  
    # get beta
    Rt = generate_Rt(date,sce)
    #
    beta =  f * Rt
    exposed_set =set()
    exposed_today_num = np.random.binomial(num_S, beta * num_I/ num_N)
    indices_to_exposed = set(random.choices(tuple(S),k = exposed_today_num))

    for index in indices_to_exposed:
        infected_dic[index] = Person()
        infected_dic[index].contact_with_infected(date)
        # == -1 means that the person is resistent to the infection 
        if (infected_dic[index].expected_exposed_time!= -1):
            exposed_set.add(index)

    return exposed_set




def initial_infected(infected_dic,initial_number,N):

    randomly_chosen_indices =  set(random.choices(tuple(N),k = 0))
    for random_index in randomly_chosen_indices:
        # initialize the person and added to infected_dic
        infected_dic[random_index] = Person()
        infected_dic[random_index].infection_direct(0)
        
    S = N - randomly_chosen_indices
    E = set()
    I = randomly_chosen_indices
    R = set()

    return  S, E, I, R


def get_inf_preg(infected_dic, prg_dates_dic):
        '''
    This function identifies females who were infected during their pregnancy period.
    
    Parameters:
    infected_dic (dict): Dictionary where keys are female IDs and values are objects with attributes 'exposed_date', 
                         'expected_exposed_time', and 'expected_infected_time'.
    prg_dates_dic (dict): Dictionary where keys are female IDs and values are lists of birth dates.
    
    Returns:
    list: Dictionary of every year's new born 
        '''
   
        relevant_females = set(infected_dic.keys()).intersection(prg_dates_dic.keys())
        # Step 1: Construct arrays
        relevant_indices = np.array(sorted(relevant_females))  # Sorting to ensure consistent order

        infected_preg_per_year = {y:0 for y in range(2024,2035)}
        if len(relevant_indices) == 0:
            return infected_preg_per_year
        
        exposed_dates, expected_exposed_times, expected_infected_times = zip(*[(infected_dic[index].exposed_date,
                                                                            infected_dic[index].expected_exposed_time,
                                                                            infected_dic[index].expected_infected_time)
                                                                            for index in relevant_indices])

        # Convert lists to numpy arrays
        exposed_dates = np.array(exposed_dates)
        expected_exposed_times = np.array(expected_exposed_times)
        expected_infected_times = np.array(expected_infected_times)
        # Constructing birth_dates and flattened_indices to keep track of which birth_date corresponds to which female
        birth_dates = []
        flattened_indices = []
        pos_ls = [] # positional operation
        for pos, index in enumerate(relevant_indices):
            dates = prg_dates_dic[index]
            birth_dates.extend(dates)
            flattened_indices.extend([index] * len(dates))
            pos_ls.extend([pos] * len(dates))
             # mother index 
        birth_dates = np.array(birth_dates)
        flattened_indices = np.array(flattened_indices)

        # Step 2: Calculate end_day, conception
        end_days = exposed_dates + expected_exposed_times + expected_infected_times - 1
        conceptions = birth_dates - 279

        # Map end_days and exposed_dates to the flattened indices
        end_days_mapped = end_days[pos_ls]
        exposed_dates_mapped = exposed_dates[pos_ls]

        # Step 3: Calculate infected column, use 
        infected = np.maximum(-1, np.minimum(end_days_mapped, birth_dates) - np.maximum(exposed_dates_mapped, conceptions))

        # Step 4: Filter rows where infected != -1
        valid_rows = infected != -1
        valid_birth_dates = birth_dates[valid_rows]

        # Step 5: Calculate birth years and count frequencies
        birth_years = (valid_birth_dates // 365) + 2024
        
        #  store results 
        

        for bir_y in birth_years:
            infected_preg_per_year[bir_y] += 1

        return infected_preg_per_year











    # for index in relevant_females:
    #     female = infected_dic[index]
    #     exposed_date = female.exposed_date
    #     expected_exposed_time = female.expected_exposed_time
    #     expected_infected_time = female.expected_infected_time

    #     carry_days = (exposed_date, exposed_date + expected_exposed_time + expected_infected_time - 1)
        
    #     birth_days = prg_dates_dic[index]
        
    #     for birth_date in birth_days:
    #         preg_days = (birth_date - 279, birth_date)
    #         if get_overlap(carry_days, preg_days) != -1:
    #             infants_per_year[int(birth_date)//365] += 1
    #             break
  
    # return infants_per_year

# mother - index- date 
# relevant_females
# mother - date construct  x*2 array from prg_dates_dic dic , female_index:[date1,date2], preg_ary
# add mother - exposed_date - expected_exposed_time-expected_infected_time array from infected_dic
# new column : end_day  exposed_date + expected_exposed_time + expected_infected_time - 1
# new column: conception =  date - 279
# new column   infected column: max(-1, min(  end_day , date) - max( exposed_date , date - 279));  # rewrite using numpy to get faster
# get the rows where infected column !=0
# use collection to count int(birth_date)//365 frequencies; should be a dictionary




if __name__ == "__main__":


    # Example usage, with no problem  
    infected_dic = {i: Person() for i in range(6)}
    for index in infected_dic:
        infected_dic[index].contact_with_infected(300)

    prg_dates_dic = {
        1: [100 ],  # List of birth dates for female ID 1
        2: [300], 
         3:[300,305] # List of birth dates for female ID 2
        # Add more entries as needed
    }
    for index in infected_dic:
        print(  infected_dic[index].__dict__ )

    infants_per_year = get_inf_preg(infected_dic, prg_dates_dic)
    print(infants_per_year)
  

            