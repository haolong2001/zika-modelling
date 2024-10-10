

import json
import numpy as np
from config import *
from modules.infection import *
import os
from demographic import get_delete_dic,get_birth_dic,read_demo
from modules.Person import Person 
from modules.vac_function import pre_vac_schedule


    
def simulation(initial_number,sce,imu, import_set,config,eff, cov, Vac):
    '''
    the principle of this version is to create dynamic program
    only when infected do we initialize the infection Person 
    and then 
    '''
    num_time_steps = config['num_time_steps']   
    total_popu = config['initial_popu']   
    popu_array,delete_year_dic,prg_dates_dic,prg_dates_dic_by_dates,_ = read_demo(config)



    # get all the dates of delete_years, and preg_dates
    death_list = sorted(delete_year_dic.keys())
    birth_list = sorted(prg_dates_dic_by_dates.keys())

    # track the location of birth dates/ death dates
    i = 0 
    j = 0 

    # initialize resistent percentage
    perc = eff * cov

    # special set up for vaccination 8 
    if ( Vac == 8):
        vac_in,vac_out = pre_vac_schedule(prg_dates_dic,perc)

        vac_in_list = sorted(vac_in.keys())
        vac_out_list = sorted(vac_out.keys())

        i_vac = 0
        j_vac = 0

        # print(vac_in_list)
        # print( len(vac_in_list))
        # print( len(vac_out_list))



    # only use dictionary to store the person
    # total_popu is the length of all the populaton in 2024
    infected_dic = {}
    # num of participated cases 
    num_of_parti = int(total_popu * (1-imu))
    #N = set( np.arange(num_of_parti) )
    N_parti = np.random.choice(np.arange(total_popu), num_of_parti, replace=False)
    N = set(N_parti)

    S, E, I, R = initial_infected(infected_dic,0,N) # no infection at beginning  
    V = set()

    # initialize storing results 
    total_infected = 0 
    matrix_shape = (num_time_steps, 7)
    result_matrix = np.zeros(matrix_shape, dtype=int)

    

    for day in range(num_time_steps):


        
        # 1. Infection
        num_S = len(S)
        num_I = len(I)
        num_E = len(E)
        num_R = len(R)
        num_N_with_imu = total_popu

        # put in the V
        # first year is special, and then campaign 
        # three doses at first 
        # update  v, ignore those who are in E,I,R state 



        # do daily update vaccination
        if ( Vac == 8):
            if  i_vac < len(vac_in_list):
                if day == vac_in_list[i_vac]   :
                    vac_set = set( (vac_in[vac_in_list[i_vac]]) )
                    i_vac += 1
                    # update S,V; ignore those go into E,I,R...
                    go_v = S.intersection(vac_set)    
                    V.update(go_v)
                    S-=(go_v)

            if j_vac < len(vac_out_list):
                if day == vac_out_list[j_vac] :
                    exp_out_set = set( vac_out[vac_out_list[j_vac]] )
                    j_vac += 1
                    # update S,V
                    out_vac_set = V.intersection(exp_out_set)
                    V-=(out_vac_set)
                    S.update(out_vac_set)
                    


        # do yearly update vaccination

        if ( day % 365 ==0 and Vac != 0 and Vac != 8):
            
            vac_set = ()
            cur_year = str(day // 365 + 2024) 
            str_Vac = str(Vac)
            # get current year vaccination 
            fem_list = []
            if cur_year in add_v_all_sce[str_Vac].keys():
                for age in add_v_all_sce[str_Vac][cur_year]: # [27,28,29]
                    
                    str_age = str(age)
                    index_begin, index_end = female_age_dic[str_age] # both included 
                    # randomw choose 
                    fem_list.extend( range( index_begin, index_end+1) )

                # Convert set to a numpy array for sampling
                vac_array = np.array(fem_list)
                
                
                # Calculate the number of samples to take
                sample_size = int(len(vac_array) * perc)
                
                # Randomly sample 20% of the elements
                vac_set = set( np.random.choice(vac_array, size=sample_size, replace=False))
        
                go_v = S.intersection(vac_set)    
                
                V.update(go_v)
                S-=(go_v)


            # only do exit with sce 4,5,6,7
            if Vac > 3 :
                out_fem_list = []
                if cur_year in delete_v_all_sce[str_Vac].keys():
                    for age in delete_v_all_sce[str_Vac][cur_year]: # [27,28,29]
                        
                        str_age = str(age)
                        index_begin, index_end = female_age_dic[str_age] # both included 
                        # randomw choose 
                        out_fem_list.extend( range( index_begin, index_end+1) )

                    # remove the set to S
                    # intersection of females with V, to know who are truly protected 
                    # O(1) but may cost a bit time 
                    out_vac_set = V.intersection(set(out_fem_list))
                    V-=(out_vac_set)
                    S.update(out_vac_set)
            
            # special pre pregnancy scenario


        

        exposed_set =do_exposed(S,I, num_N_with_imu, day, infected_dic,sce,False)
        
        # Move those indices from S to E,update the state
        # end of day, update E I R
        infected_set = set()
        for index in E:
            person = infected_dic[index]
            person.update_by_day()
            if (person.days >= person.expected_exposed_time) :
                person.make_infection()
                infected_set.add(index)
        new_infected_num = len(infected_set)

        # recovery set
        recovery_set = set()
        for index in I:
            person = infected_dic[index]
            person.update_by_day()
            if (person.days >= person.expected_infected_time ):
                person.make_recovery()
                recovery_set.add(index)
        
        S -= exposed_set
        E.update(exposed_set)
        
        E-= infected_set
        I.update(infected_set)
                
        I-= recovery_set
        R.update(recovery_set)
        #new_infected_num = 0 #delete later 

        num_S = len(S)
        num_I = len(I)
        num_E = len(E)
        num_R = len(R)

        # update population, here the new born does not have immunity 
        # set operation will ignore the immunity parts  
        # limit the day, because it's just 3650, so no problem at all   
        delete_set =set()
        if day ==  death_list[i]:
            delete_set =set(  delete_year_dic[day] )
            i+=1
        
        add_set = set()
        if day == birth_list[j]:
            add_set = set( prg_dates_dic_by_dates[day] )
            j+=1
            
        S -= delete_set                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        S.update(add_set)
        N.update(add_set)
        E -= delete_set
        I -= delete_set
        R -= delete_set
        N -= delete_set
        V -= delete_set

        #monthly imports, simply assume they belong to the same day
        if day % 30 == 0:
                num_add = np.random.randint( import_set[0],import_set[1]+1)
                import_to_infected  = set(random.choices(tuple(S),k = num_add))
                for random_index in import_to_infected:
                        infected_dic[random_index] = Person()
                        infected_dic[random_index].contact_with_infected(day)
                        infected_dic[random_index].make_infection()
                S.difference_update(import_to_infected)
                I.update(import_to_infected)
                # count the new_infected_num 
                new_infected_num+= num_add
        
        #store the result in matrix

        total_infected += new_infected_num 
        num_I +=new_infected_num

        result_matrix[day, :] = [day,num_S, num_E, num_I, num_R, new_infected_num, total_infected]
    # start infected pregnancy calculation
    infants_per_year =  get_inf_preg(infected_dic, prg_dates_dic)
     
    return result_matrix,infants_per_year # json

# do a test  
if __name__ == "__main__":
    simulation( 1,'sustained', 0.05, (40, 50), 0.5, 0.2, 8)
       






# # Example usage:
# delete_year_dic = {
#     '2019': ['data1', 'data2'],
#     '2020': ['data3', 'data4']
# }

# delete_infant = {
#     '2019': ['data5', 'data6'],
#     '2021': ['data7', 'data8']
# }

# merge_dicts(delete_year_dic, delete_infant)

# print(delete_year_dic)
