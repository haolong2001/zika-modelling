import numpy as np
import json
from collections import defaultdict
from config import vac_sce

def convert_dic(dic):
    '''
    convert the int32 to int for storing purpose 
    '''
    return {int(key): [int(value) for value in values] for key, values in dic.items()}


def get_age_indices():
    '''
    Dictionary mapping age to a tuple of (begin index, end index) for females.
    '''
    population_data = np.load('../population_data/2024_population_data.npy')
    
    female_age_dic = defaultdict(list)
    female_popu = population_data[population_data[:, 4] == 2]
    for age in range(1, 50):
        rows = female_popu[female_popu[:, 2] == age]
        begin_idx = rows[0, 0]
        # Since it's sorted, we only need to find the first and last row, column 0

        end_idx = rows[-1, 0]
        female_age_dic[age] = [begin_idx, end_idx]
    
    # Save female_age_dic to population_data/female_age_dic.json 
    with open('../population_data/female_age_dic.json', 'w') as f:
        json.dump(convert_dic(female_age_dic), f)
    
    return female_age_dic

def V_construct(sce):
    '''
    Constructs vaccination schedule dictionaries for adding and deleting ages.
    Keys are years, values are lists of ages.
    Does not include newborns.
    '''
    start_age, vac_year, expected_vac_freq = vac_sce[sce]
    end_age = start_age + vac_year * expected_vac_freq
    ages = list(range(1, np.minimum( 49,end_age )))
    
    start_year = np.maximum(0, start_age - np.array(ages)) + 2024
    
    # Use 2040 as a cut-off point 
    add_V = defaultdict(list)
    for i in range(len(start_year)):
        if start_year[i] < 2034:
            add_V[start_year[i]].append(ages[i])
    
    total_freq = expected_vac_freq - np.maximum(0, (np.array(ages) - start_age) // vac_year)
    total_period = total_freq * vac_year
    end_years = start_year + total_period
    
    delete_V = defaultdict(list)
    for i in range(len(end_years)):
        if end_years[i] < 2034:
            delete_V[end_years[i]].append(ages[i])
    
    return  convert_dic(add_V), convert_dic(delete_V)

def main():
    

    add_v_all_sce = {}
    delete_v_all_sce = {}
    for sce in range(1, 8):
        add_V, delete_V = V_construct(sce)
        add_v_all_sce[sce] = add_V
        if sce > 3:
            delete_v_all_sce[sce] = delete_V
    
    saving_path = '../population_data/vac_data'
    
    # Save the two dictionaries to JSON files
    with open(f'{saving_path}/add_v_all_sce.json', 'w') as f:
        json.dump(add_v_all_sce, f)
    with open(f'{saving_path}/delete_v_all_sce.json', 'w') as f:
        json.dump(delete_v_all_sce, f)
    




if __name__ == "__main__":
    pass

    #get_age_indices()
























# import numpy as np
# from config import vac_sce


# def get_age_indices():
#     '''
#     dictionary age: (begin index,end index )
#     ''' 
#     population_data = np.load('../population_data/2024_population_data.npy')
    
#     female_age_dic = {}
#     female_popu = population_data[ population_data[,4] == 2]
#     for age in range(18:50):
#         rows = female_popu[ female_popu[,2] ==age ]
#         begin_idx =rows[0,0]
#     # since it's sorted, we only need to find first row and last row, column 1 
#         end_idx = rows[rows.shape[0]-1,0]
#         female_age_dic[age] = (begin_idx,end_idx) 
    
#     # save female_age_dic to population_data/female_age_dic.json 
#     return 
# def V_construct(sce):
#     '''
#     add_V dictionary
#     the keys are years, the values are list of addof ages
#     currently does not touch with newborns 
#     '''
#     start_age, vac_year,expected_vac_freq = vac_sce[sce]
#     end_age = start_age + vac_year * expected_vac_freq
#     ages = range(1,end_age)
#     #  
#     # if ages above start_age, then it will be 2024 
#     # if ages below start_age(e.g. 25), then we get the years when it reaches 25
#     start_year = - np.min(0,start_age - ages) + 2024
    
#     # use 2040 as a cut point 
#     add_V = defaultdic( [])
#     for i in range( len(start_year)):
#         if (start_year < 2034):
#             add_V[ start_year[i] ].append( ages[i] )
    
#     # if ages > start_age , and this value is larger than 5, then means that the duration
#     # years are smaller 
#     total_freq = vac_freq - np.max(0, ages -start_age )  // vac_year
#     total_period = total_freq * vac_year 
#     end_years = start_year + total_period 

#     delete_V = defaultdic([])
#     for i in range( len(start_year)):
#             if (start_year < 2034):
#                 delete_V[ end_years[i] ].append( ages[i] )
    
#     return add_V, delete_V




# def main():
#      add_v_all_sce = {}
#      delete_v_all_sce = {}
#      for sce in range(1,8):
#           add_V, delete_V = V_construct(sce)
     
#         add_v_all_sce[sce] = add_V
#         delete_v_all_sce[sce] = delete_V
    
#         saving_path = ../population_data/vac_data
#         save the two dictionary to json files 
          
          
          



# if main == main :
#      main()
#      get_age_indices()