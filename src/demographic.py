'''
initialize Singapore demographic
save in the json file;
avoid memory leaks 
'''
from collections import defaultdict
import json
import numpy as np
import pandas as pd 
import argparse
from modules.birth_function import GetPregnancyDates, getGender
from modules.death_function import getDeathAge,batchDeathRate0, batchDeathRate_general
from config import * 
import time 
import itertools
import time 


def save_dict_to_json(dic, file_path):
    with open(file_path, 'w') as f:
        json.dump(dic, f)

def load_dict_from_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)
    
def revert_dic(dic):
    '''
    convert the default dic of list 
    from string to np.int32 format
    '''
    return {np.int32(key): [np.int32(value) for value in values] for key, values in dic.items()}


def convert_dic(dic):
    '''
    convert the int32 to int for storing purpose 
    '''
    return {int(key): [int(value) for value in values] for key, values in dic.items()}


def get_2023_popu(sampling_fraction):
    """"
    Returns:    
    np array: consists of the demographic information/population structure of Singapore 2023 
    male / female of age group 0 to 85 and 85+
    """

    male_2023 = pd.read_csv("../original_data/Male_2023.csv", index_col=0) 
    female_2023 = pd.read_csv("../original_data/Female_2023.csv", index_col=0)
    # only get  0 - 85 years old and 85+ , target population 5917653
    male_vec = np.array(male_2023.iloc[list(range(1,87)) + [95],0]) * sampling_fraction
    female_vec = np.array(female_2023.iloc[list(range(1,87)) + [95],0]) * sampling_fraction
    
    male_vec_rounded = np.round(male_vec).astype(int)
    female_vec_rounded = np.round(female_vec).astype(int)
    # print( sum(np.round(male_vec)) + sum(np.round(female_vec)))
    return male_vec_rounded, female_vec_rounded




def initialize_demo(male_vec, female_vec,config):
    '''
    Input:
        male_vec (np.array): Array of male population counts by age
        female_vec (np.array): Array of female population counts by age
    Output:
        np.array: Population data as of 2024.01.01 approximated using 2023.6 data
    The output array has the following columns:
        0 - index
        1 - yearOfAdd
        2 - AgeOfAdd (age + 1)
        3 - age
        4 - gender (1 for male, 2 for female)
        5 - deathAge
    '''
    cur_index = 0
    year_of_add = 2024
    total_population = np.sum(male_vec) + np.sum(female_vec)
    
    # Initialize the population array with the correct shape and data type
    # only count the death rate within these 10 years
    popu_array = np.zeros((total_population, 6), dtype=np.int32)
    popu_array[:,0] = np.arange(total_population)
    popu_array[:,1] = year_of_add

    # generate age of age
    ages = np.arange(1, 88) # generate 1 to 87
    male_ages = np.repeat(ages, male_vec)
    female_ages = np.repeat(ages, female_vec)
    flattened_ages = np.concatenate((male_ages, female_ages))
    popu_array[:,2] = flattened_ages

    # age , same as age of add
    popu_array[:,3] = flattened_ages

    # gender
    male_genders = np.repeat(1, np.sum(male_vec))
    female_genders = np.repeat(2, np.sum(female_vec))
    gender_vector = np.concatenate((male_genders, female_genders))
    popu_array[:,4] = gender_vector

    # death age
    popu_array[:,5] = batchDeathRate_general(flattened_ages)

    female_start_index = np.sum(male_vec)
    female_end_index = total_population
    cur_index = female_end_index

    return popu_array,female_start_index,female_end_index,cur_index

def get_birth_dic(popu_array, female_start_index, female_end_index, cur_index, config):
    '''
    input: 

    2024 population array; indices of all the female, the last index

    output:
    
    infant_array: stores the attributes of a baby
    prg_dates_dic: mother index: dates of giving out a baby
    mother_infant_dic: mother index: child indices
    prg_dates_dic_by_dates:    date: indices of infants born on this day

        0 - index
        1 - yearOfAdd
        2 - AgeOfAdd (age + 1)
        3 - age
        4 - gender (1 for male, 2 for female)
        5 - deathAge

    '''

    FR_matrix = np.load("../original_data/FRBY.npy") 
    
    mother_infant_dic = defaultdict(list) # not useful

    # Filter the array for relevant female indices and age <= 49
    relevant_females = popu_array[female_start_index:female_end_index]
    childbearing_females = relevant_females[relevant_females[:, 2] <= 49]

    year_of_add = year_0

    female_indices = childbearing_females[:,0]
    ages_of_add = childbearing_females[:,2] # the 2nd col from 0 
    ages_of_death = childbearing_females[:,5]


    # construct birth matrix according to females 
    num_year = 10 
    output_matrix = np.zeros((len(female_indices), num_year+1))
    for i, age in enumerate(ages_of_add):
            output_matrix[i, :] = FR_matrix[age, :]

    
    random_matrix = np.random.rand( len(female_indices), num_year+1)
    birth_matrix = (random_matrix <output_matrix ).astype(np.int32)

    endpoints = np.minimum(ages_of_death - ages_of_add, 11)
    col_indices = np.arange(num_year +1 )
    mask = col_indices >= endpoints[:, None]
    birth_matrix[mask] = 0

    # input is the birth matrix ; female_indices[i] will be the mother index
    #  stand for fertility condition from 2024 - 2034 
    # output 
    # infant_array, prg_dates_dic, mother_infant_dic,prg_dates_dic_by_dates
    # mother_infant_dic = {} # put it to empty 
    # prg_dates_dic_by_dates   # dates: index on that date 
    # prg_dates_dic # mother_index to dates of birth

    prg_dates_dic = defaultdict(list)

    # construct a list of list 
    birth_years = []
    birth_dates = []
    for row in birth_matrix:
        years = np.nonzero(row)[0] # will return array([])
        dates = years * 365 + np.random.randint(0,365)
        birth_years.append(years + 2024)
        birth_dates.append(dates)
    
    # empty_count = sum(dates.size == 0 for dates in birth_dates)
    # print(f"Number of empty arrays in birth_dates: {empty_count}")
    prg_dates_dic = {index: date for index, date in zip(female_indices, birth_dates) if date.size != 0}


    # construct infants array 
    # flatten the list of list 
    flattened_years = list(itertools.chain.from_iterable(birth_years))
    flattened_dates = list(itertools.chain.from_iterable(birth_dates))
    infants_begin = cur_index
    num_infants =  len(flattened_years)
    infants_end = cur_index + num_infants
    infant_array = np.zeros((num_infants, 6), dtype=np.int32)
    infant_array[:,0] = range(infants_begin,infants_end)
    infant_array[:,1] = flattened_years
    random_vector = np.random.rand(num_infants)
    infant_array[:,4]  = np.where(random_vector >= 0.485, 1, 2)
    infant_array[:,5] = batchDeathRate0( num_infants)

    # construct birth dates 
    prg_dates_dic_by_dates = defaultdict(list)
    for key, value in zip(flattened_dates, infant_array[:,0]):
        prg_dates_dic_by_dates[key].append(value)

    # add test here 
    #print(prg_dates_dic_by_dates.values)
    # np.save('../population_data/2024_infants_data.npy',infant_array )
 
    # save_dict_to_json( convert_dic(prg_dates_dic_by_dates), f'../population_data/test_infants_dates_dic.json')
    # print("Saved pregnancy dates dictionary to test_infants_dates_dic.json.")


    return infant_array, prg_dates_dic, mother_infant_dic,prg_dates_dic_by_dates





def get_delete_dic(popu_array):
    '''
    input: popu_array; output delete_dic,delete_year_dic
        0 - index
        1 - yearOfAdd
        2 - AgeOfAdd (age + 1)
        3 - age
        4 - gender (1 for male, 2 for female)
        5 - deathAge
    '''
    delete_dic = defaultdict(list)
    delete_year_dic = defaultdict(list)
    
    # Filter the array to exclude rows where deathAge is 999
    filtered_popu_array = popu_array[popu_array[:, 5] != 999]

    # Extract the columns of interest
    indices = filtered_popu_array[:, 0]
    deathAges = filtered_popu_array[:, 5]
    agesOfAdd = filtered_popu_array[:, 2]
    yearsOfAdd = filtered_popu_array[:, 1]

    # Calculate death_years and death_dates
    death_years = yearsOfAdd - agesOfAdd + deathAges
    death_dates = (death_years - year_0) * 365 + np.random.randint(0, 365, size=death_years.shape)

    delete_year_dic = defaultdict(list)

    # Populate the defaultdicts using vectorized operations
    for  (index,death_date) in (zip(indices,death_dates)):

        delete_year_dic[death_date].append(index)
            
    return  delete_year_dic






def merge_dicts(delete_year_dic, delete_infant):
    for key, value in delete_infant.items():
        # Check if the key exists in delete_year_dic
        if key in delete_year_dic:
            # If the key exists, extend the list with the values from delete_infant
            delete_year_dic[key].extend(value)
        else:
            # If the key does not exist, add it to delete_year_dic
            delete_year_dic[key] = value


def read_demo(config):
    '''
     output the variables we need 
     {num of population 
     delete_year_dic
     prg_dates_dic,prg_dates_dic_by_dates}
    '''
    # load the storing population and death  
    demo_path = config['demo_path']
    delete_year_dic_path = f'{demo_path}/delete_year_dic.json'
    delete_year_dic = revert_dic( load_dict_from_json(delete_year_dic_path) )

    with open(f'{demo_path}/2024_population_data.npy', 'rb') as file:
        population_array = np.load(file)


    # read the parameters 
    indices_array = np.load(f"{demo_path}/indices.npy")
    female_start_index = indices_array[0]
    female_end_index = indices_array[1]
    cur_index = indices_array[2]

    # construct the new birth related file 
    infant_array, prg_dates_dic, mother_infant_dic,prg_dates_dic_by_dates = get_birth_dic(
            population_array, female_start_index, female_end_index, cur_index,config
            )
    
    # combine the two population matrix 
    population_data = np.concatenate((population_array, infant_array), axis=0)

    # construct the delete dictionary,and merge it into delete_year_dic
    delete_infant = get_delete_dic(infant_array)
    merge_dicts(delete_year_dic, delete_infant) # mutable, no need for return 

    return population_data,delete_year_dic,prg_dates_dic,prg_dates_dic_by_dates, female_end_index
    






def test(config):
    start_time = time.time()
    num_time_steps = config['num_time_steps'] 
    demo_path = config['demo_path']   
    popu_array,delete_year_dic,prg_dates_dic,prg_dates_dic_by_dates,female_end_index = read_demo(config)

    end_time = time.time()
    print(f" time:{end_time -start_time}")
    np.save(f'{demo_path}/test_popu_array.npy', popu_array)

    save_dict_to_json( convert_dic(prg_dates_dic), f'{demo_path}/prg_dates_dic.json')
    print("Saved pregnancy dates dictionary to prg_dates_dic.json.")

    save_dict_to_json( convert_dic(prg_dates_dic_by_dates), f'{demo_path}/prg_dates_dic_by_dates.json')
    print("Saved pregnancy dates dictionary to prg_dates_dic_by_dates.json.")

    print(female_end_index)


def main(config):
        demo_path = config['demo_path']
        sampling_fraction = config['sampling_fraction']
        
        # initialize the base population  
        male_vec, female_vec = get_2023_popu(sampling_fraction)
        popu_array,female_start_index,female_end_index,cur_index  = initialize_demo(male_vec, female_vec,config)
        delete_year_dic = get_delete_dic(popu_array)
        indices_array = np.array([female_start_index, female_end_index, cur_index])
        # save 2024 population
        np.save(f'{demo_path}/2024_population_data.npy', popu_array)
        # Save the parameters for female construction to a .npy file
        np.save(f"{demo_path}/indices.npy", indices_array)
        # save the delete 
        save_dict_to_json( convert_dic(delete_year_dic), f'{demo_path}/delete_year_dic.json')


        print("Saved")




if __name__ == "__main__":
        # the parser parts
        parser = argparse.ArgumentParser(description="Run the program with the specified configuration.")
        parser.add_argument('-v', default='real world', choices=['test', 'real world'], help="Specify the configuration to use.\
                            Options are 'test' or 'real world'. Default is 'real world'.")
        args = parser.parse_args()
        if args.v == 'test':
            config = test_config
        else:
            config = actual_config

        test(config)
        #main(config)
        


        ###################################################################################
        # start_time = time.time()
        # infant_array, prg_dates_dic, mother_infant_dic,prg_dates_dic_by_dates = get_birth_dic(
        # popu_array, female_start_index, female_end_index, cur_index
        # )

        # population_data = np.concatenate((popu_array, infant_array), axis=0)


        # np.save(f'{result_path}/population_data.npy', population_data)

        # Save the dictionaries to JSON files
        # save_dict_to_json( convert_dic(delete_dic), f'{result_path}/delete_year_dic.json')
        # print("Saved delete dictionary to delete_year_dic.json.")

        

        # save_dict_to_json( convert_dic(prg_dates_dic), f'{result_path}/prg_dates_dic.json')
        # print("Saved pregnancy dates dictionary to prg_dates_dic.json.")

        # save_dict_to_json( convert_dic(mother_infant_dic), f'{result_path}/mother_infant_dic.json')
        # print("Saved mother-infant dictionary to mother_infant_dic.json.")

        # save_dict_to_json( convert_dic(prg_dates_dic_by_dates), f'{result_path}/prg_dates_dic_by_dates.json')
        # print("Saved pregnancy dates by dates dictionary to prg_dates_dic_by_dates.json.")

    # check female beginning index 
    
    # expected_end = np.sum( female_vec ) + cur_index
    # print(f'expected ending point{expected_end}')
    # female_start_index = cur_index 
    
    # for age in range(0,86):
    #     # age 86 stands for 85+ 
    #     i = age
    #     num_per_age = female_vec[i]
        
    #     for _ in range(num_per_age):
    #         popu_array[cur_index,] = ( cur_index, year_of_add,age+1, age+1,2, getDeathAge(yearOfadd,age) )
    
    #         cur_index = cur_index + 1
    
    # # add 85+ to the population dictionary
    # for _ in range(female_vec[86]):
    #         age = 86
    #         popu_array[cur_index,] = ( cur_index, year_of_add,age,age 2, getDeathAge(yearOfadd,age) )
    #         # for age of 86, simply add them to 85+ group

    #         # get fertility and a dictionary 
    #         pass
    #         # also get year: fertlity dictionary 
    
    #         cur_index = cur_index + 1
    # female_end_index = cur_index
    # print(f"end of new birth index {cur_index}")
    # return cur_index, popu_array, prg_dates_dic, prg_dates_dic_by_dates, death_by_year_dic 
   

