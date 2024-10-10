import os
import sys
import random
current_directory = os.path.dirname(os.path.abspath(__file__))  # Assuming this code is in a script
# Get the parent directory
parent_directory = os.path.dirname(current_directory)

# Insert the parent directory into sys.path
sys.path.insert(0, parent_directory)
from config import year_0,year_end
import numpy as np
import pandas as pd

ASFR = pd.read_csv("../original_data/aggregated_ASFR.csv") 
# size: 28 *35

# def getFertilityRate(year, age):
#     global ASFR
#     # fertility begins with year 2023
#     # with age 15
#     year_idx = year - 2023
#     age_idx = age - 15 
#     # print(year_idx, age_idx)
#     fer_prob = ASFR.iloc[year_idx, age_idx] 
#     return fer_prob
def getGender():
    r = random.random()
    threshold =  0.485 # 1/(1 + 1.059)
    if(r >= threshold):
        return 1
    else:
        return 2
    
def GetPregnancyDates(yearOfadd, AgeOfAdd,age_of_death ):
        '''
        get pregnancy within the range

        start ( year_0, age of 15)
        end (year_end , age of 49, or age of death)

        output: list of pregnancy dates 
        '''
        pregnancy_start = max(year_0, yearOfadd - AgeOfAdd + 15)

        death_year = yearOfadd  - AgeOfAdd + age_of_death 
        pregnancy_end = min(year_end +1, yearOfadd  - AgeOfAdd + 49, death_year)

        if (pregnancy_start <= pregnancy_end ):
            age_of_start = pregnancy_start - yearOfadd + AgeOfAdd

            preg_list = getPregnancyList(pregnancy_start,
                            pregnancy_end,
                            age_of_start
                            )
        else:
            preg_list = []

        return preg_list
        


def getPregnancyList(pregnancy_start,
                          pregnancy_end,
                          age_of_start):
    """
    pregnancy_start  year of 15 or 2024
    pregnancy_end    year of 49 or 2050
    age_of_start     age of pregnancy_start year 

    we only consider the pregnancy_start and pregnancy_end in the range of 2024 to year_end 
    special case:
    1. pregnancy_start > pregnancy_end 
    """

    years = np.arange(pregnancy_start, pregnancy_end + 1)
    ages = np.arange(age_of_start, age_of_start + len(years))

    # Vectorized retrieval of fertility rates
    col_index = ages - 15
    y_index = years -2023

    fertility_rates = [ASFR.iloc[year, col] for year, col in zip(y_index, col_index)]
    # Determine pregnancies based on fertility rates
    binary_preg = np.random.random(len(fertility_rates)) <= fertility_rates
    # corresponding to the years 
    # Create periods_list
    periods_list = []

    for year, is_pregnant in zip(years, binary_preg):
        if is_pregnant:
            period = (year - year_0) * 365 + np.random.randint(0, 365) 
            periods_list.append(  period )

    return periods_list
        


if __name__ == "__main__":
        # Get the current directory of the script
   
    year_of_add = 2024
    age_of_add = 25


    pregnancy_dates = GetPregnancyDates(year_of_add, age_of_add,55)
    print(pregnancy_dates)
    # test death at 28 
    pregnancy_dates = GetPregnancyDates(year_of_add, age_of_add,28)
    print(pregnancy_dates)
     



# def generatePeriod(year):
#     # random
#     end_of_conception = (year - year_0) * 365 + np.random.randint(0, 365) + 1
#     start_of_conception = end_of_conception - 280 + 1
#     return start_of_conception, end_of_conception


# def getChildbearing_Years_and_Periods(pregnancy_start,
#                           pregnancy_end,
#                           age_of_start):
#     """
#     pregnancy_start  year of 15 or 2024
#     pregnancy_end    year of 49 or 2050
#     age_of_start     age of pregnancy_start year 

#     we only consider the pregnancy_start and pregnancy_end in the range of 2024 to 2050 
   
    
#     special case:
#     1. pregnancy_start > pregnancy_end 
#     """

#     if pregnancy_start > pregnancy_end:
#             return [], []
    

    
#     age = age_of_start
#     list_of_fertility_rate = []
#     list_of_years = []

#     for  year in range(pregnancy_start,pregnancy_end+1):
#             # if age < 15:
                
#             #     age += 1
#             #     continue
            
#             # if age > 49:
#             #     break
            
#             list_of_fertility_rate.append(getFertilityRate(year,age))
#             list_of_years.append(year)

#             age += 1
    
#     # print(list_of_fertility_rate,list_of_years)
#     binary_preg = [1 if np.random.random() < rate else 0 for rate in list_of_fertility_rate]

#     # corresponding to the years 
#     # Create periods_list
#     periods_list = []
#     pregnancy_years_list = []

#     if sum(binary_preg) == 0:
#         return [], []
#     else:
#         for i, binary_value in enumerate(binary_preg):
#             if binary_value == 1:
#                 year = list_of_years[i]
#                 period = generatePeriod(year)
#                 periods_list.append(period)
#                 pregnancy_years_list.append(year)
            
#         # print(periods_list)    

#         return pregnancy_years_list, periods_list