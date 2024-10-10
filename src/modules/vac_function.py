import numpy as np
import json
from collections import defaultdict


def pre_vac_schedule(prg_dates_dic, perc):
    '''
    '''
    # construct a date:index format array 


    preg_list = []
    for index, dates in prg_dates_dic.items():
        for date in dates:
            preg_list.append([index, date])

    preg_array = np.array(preg_list)

    percentage = perc
    num_rows = len(preg_array)
    num_selected = int(num_rows * percentage)
    indices = np.random.choice(np.arange(num_rows), num_selected, replace=False)


    vac_preg_array = preg_array[indices]

    vac_indices = vac_preg_array[:,0]
    vaccination_dates = np.maximum(0, vac_preg_array[:, 1] - 280)
    lose_protect_dates = np.maximum(0, vaccination_dates + 365)

    vac_in = defaultdict(list)
    vac_out = defaultdict(list)

    for i in range( len(vac_preg_array)):
        index, vaccination_date,lose_protect_date= vac_indices[i],vaccination_dates[i],lose_protect_dates[i]

        vac_in[vaccination_date].append(index)
        vac_out[lose_protect_date].append(index)

    return vac_in,vac_out