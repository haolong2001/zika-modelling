
#current_directory = os.getcwd()
# Add the parent directory to sys.path if not already present
import os
import sys
# the cwd should be src
current_directory = os.getcwd()
sys.path.insert(0, current_directory)

from modules.Person import Person 
from no_demo_config import *
import random
import numpy as np
import os
import sys
from modules.infection import *

def simulation(initial_number,sce, import_set):
    
    # initialization  
    infected_dic = {}
    num_persons = 5917648
    N = set( range(num_persons))
    S, E, I, R = initial_infected(infected_dic,initial_number,N)
    
    # store the results
    total_infected = 0 
    matrix_shape = (num_time_steps, 7)
    result_matrix = np.zeros(matrix_shape, dtype=int)

    #print(num_time_steps)
    cwd_ID = num_persons
    for day in range(num_time_steps):
        
        # 1. Infection
        num_S = len(S)
        num_I = len(I)
        num_E = len(E)
        num_R = len(R)
        num_N = len(N)

        exposed_set =do_exposed(S,I, num_N, day, infected_dic,sce,zika_params,False)
        # Move those indices from S to E,update the state
        # end of day, update E I R
        infected_set = set()
        for index in E:
            person = infected_dic[index]
            person.update_by_day()
            if (person.days >= person.expected_exposed_time ) :
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

        # back to S, no going back; consistent with others
        # back_to_S = set()
        # for index in R:
        #     person = persons_dict[index]
        #     person.update_by_day()
        #     # 0 to 365*2 -1 is the two year
        #     if (person.days >= 100*2):
        #         person.back_to_S()
        #         back_to_S.add(index)

        # R-= back_to_S
        # S.update(back_to_S)
        num_S = len(S)
        num_I = len(I)
        num_E = len(E)
        num_R = len(R)

        # delete people 
        delete_set =set()
        death_people = np.random.binomial(num_N, mu)
        #indices_to_delete = set(random.choices(tuple(S),k = exposed_today_num))

        # add new people to replace the dead people 
        indices_to_delete = set(random.choices(tuple(N),k = death_people))
        indices_to_add = set( range( cwd_ID, cwd_ID + death_people ))
        cwd_ID = cwd_ID + death_people

        # enter birth/death set 
        add_set = indices_to_add
        delete_set = indices_to_delete

        S -= delete_set                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        S.update(add_set)
        E -= delete_set
        I -= delete_set
        R -= delete_set

        N.update(add_set)
        N-= delete_set

        # monthly updates
        if day % 30 == 0:
                num_add = np.random.randint( import_set[0],import_set[1]+1)
                #num_add = np.random.binomial(num_N, new/5000000)
                import_to_infected  = set(random.choices(tuple(S),k = num_add))
                for random_index in import_to_infected:
                        infected_dic[random_index] = Person()
                        infected_dic[random_index].contact_with_infected(day)
                        infected_dic[random_index].make_infection()
                S.difference_update(import_to_infected)
                I.update(import_to_infected)
        
        # store the result in matrix
        total_infected += new_infected_num
        result_matrix[day, :] = [day,num_S, num_E, num_I, num_R, new_infected_num, total_infected]
    
    return result_matrix


if __name__ == "__main__":
    pass   