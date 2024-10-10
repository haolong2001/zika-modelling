import pandas as pd
import numpy as np
import os



average_death_rates = np.load("../original_data/average_death_rates.npy")

# write a batch processing

def getDeathRates(ages):
    '''
    just for age 0- 10
    '''
    global average_death_rates
    cats = ages // 5
    return average_death_rates[cats]



def batchDeathRate_general(ages):
    '''
    input ages, suppose np.array([15,20,20,20...21])
    output deathages

    generate an len(ages) * 10 matrix ;
            the new matrix each row  is from age to age + 9
            
            get a matrix where each row  cats = min(ages // 5, len(average_death_rates)-1) for each row
            and then each row will be average_death_rates[cats]
            use numpy to make it faster

    boundary problems: 


    '''
    global average_death_rates
    num = len(ages)
    len_death = len(average_death_rates)

    # Generate the len(ages) * 10 matrix
    age_matrix = np.arange(10) + ages[:, None]

    # Compute the category for each age
    # Step 2: Divide age_matrix by 5 to categorize the ages
    category_matrix = age_matrix // 5

    # Step 3: Bound the values by len(average_death_rates) - 1
    category_matrix = np.minimum(category_matrix, len_death - 1)

 
    # Step 4: Map these bounded categories to the average death rates
    death_age_matrix = average_death_rates[category_matrix]

    # get the death age
    random_matrix = np.random.rand(num,10)

    comparison_matrix = (random_matrix < death_age_matrix).astype(int)


    row_sums = np.sum(comparison_matrix, axis=1)
    result = np.full_like(row_sums, fill_value=999)
    non_zero_sum_indices = np.where(row_sums != 0)[0]
    first_nonzero_indices = np.argmax(comparison_matrix[non_zero_sum_indices] != 0, axis=1)
    
    result[non_zero_sum_indices] = first_nonzero_indices + ages[non_zero_sum_indices]

    return result







def batchDeathRate0(num):
    """
    Compares each column of the random_matrix with a single row obtained from getDeathRates function.
    
    Returns:
    - warning: 
    - If the row sum == 0, returns 999.
    - Otherwise, returns a numpy array containing the first location of the 1 in each row.
    """
    death_rate = getDeathRates(np.array([range(10)]))

    # construct death matrix 
    random_matrix = np.random.rand(num,10)
    comparison_matrix = (random_matrix < death_rate).astype(int)

    # initialization
    
    # Check if the sum of each row is zero, means no death 
    row_sums = np.sum(comparison_matrix, axis=1)
    result = np.full_like(row_sums, fill_value=999)
    non_zero_sum_indices = np.where(row_sums != 0)[0]
    first_nonzero_indices = np.argmax(comparison_matrix[non_zero_sum_indices] != 0, axis=1)
    
    result[non_zero_sum_indices] = first_nonzero_indices
    
    return result



def getDeathAge(yearOfadd,age_start):
    """
    Calculate the age at which a person will experience death based on 
    the year of addition and the starting age.

    Parameters:
    - yearOfadd (int): The year when the person's information is recorded.
    - age_start (int): The initial age of the person when their information is added.

    Returns:
    - int: The age at which the person is expected to experience death. 
    If the person does not experience death within the simulated time frame (up to the year 2050), 
    it returns 999.

    Note:
    - The function simulates the aging process from the starting age until the year 2050.
    - The death age is determined probabilistically based on age-specific death rates 
    obtained from the `getAgeSpecificDeathRate` function.
    - If a random value generated is less than the calculated death probability, 
    the function returns the age at which death occurs.
    - If no death occurs within the simulated period, the function returns 999.
    - ...

    """
    

    age_2050 = age_start + year_end - yearOfadd
    
    for age in range(age_start,age_2050+1):
        death_prob = getAgeSpecificDeathRate(age)
        r = np.random.rand()
        if r < death_prob:
            return age   
        
    return 999 # no death


def getAgeSpecificDeathRate(age):
    """
    Input: age

    Method: The death rates in the dataframe are provided per 5 age intervals, 
    so it uses age//5 to determine the corresponding interval.

    Special Case:
    The last entry in the average_death_rates dataframe 
    (i.e., average_death_rates.iloc[len(average_death_rates)-1]) represents the age group "85+".
    If the input age is larger than the upper limit of the last age group (85+), 
    the method assumes the person's age is greater than 85, and it uses the same death rate for the 
    "85+" group.
    """

    global average_death_rates
    cat = age // 5
    return average_death_rates[min(cat, len(average_death_rates)-1)] # 85+ are all the same 





import sys
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
parent_dir = os.path.dirname(script_dir) 
sys.path.insert(0, parent_dir)
from config import year_0,year_end

def saveFile():
    # Read the CSV file, use the path relative to the current file


    data_path = os.path.join(parent_dir, "../original_data/Death_table.csv")

    #print(data_path)
    df = pd.read_csv(data_path,index_col=0)
    df.replace('-', 0.1, inplace=True)
    df = df.apply(pd.to_numeric, errors='coerce')

    selected_ages = df.iloc[list(range(1,18)) + [22]] 

    # Calculate average death rates for each age group across years 2015-2019
    average_death_rates = selected_ages.mean(axis=1) /1000.
    np.save(os.path.join("../original_data/average_death_rates.npy"), average_death_rates)

if __name__ == "__main__":
    # result = getDeathAge(2020, 70)
    # print(result)
    #age = np.random.randint(10, 15, size=(20, 0))
    result = batchDeathRate0(2000)
    print(result)

    # 87 - 97 


