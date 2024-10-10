


import numpy as np

from config import add_v_all_sce,delete_v_all_sce

def update_vaccination(day, Vac, perc, S, V, female_age_dic):

    if day % 365 == 0 and Vac != 0 and Vac != 8:
        vac_set = set()
        cur_year = str(day // 365 + 2024)
        str_Vac = str(Vac)
        
        # Get current year vaccination
        fem_list = []
        if cur_year in add_v_all_sce[str_Vac].keys():
            for age in add_v_all_sce[str_Vac][cur_year]:
                str_age = str(age)
                index_begin, index_end = female_age_dic[str_age]
                fem_list.extend(range(index_begin, index_end + 1))

            # Convert set to a numpy array for sampling
            vac_array = np.array(fem_list)
            
            # Calculate the number of samples to take
            sample_size = int(len(vac_array) * perc)
            
            # Randomly sample 20% of the elements
            vac_set = set(np.random.choice(vac_array, size=sample_size, replace=False))
            
            go_v = S.intersection(vac_set)
            V.update(go_v)
            S -= go_v

        # Only do exit with sce 4, 5, 6, 7
        if Vac > 3:
            out_fem_list = []
            if cur_year in delete_v_all_sce[str_Vac].keys():
                for age in delete_v_all_sce[str_Vac][cur_year]:
                    str_age = str(age)
                    index_begin, index_end = female_age_dic[str_age]
                    out_fem_list.extend(range(index_begin, index_end + 1))

                # Remove the set from S
                out_vac_set = V.intersection(set(out_fem_list))
                V -= out_vac_set
                
