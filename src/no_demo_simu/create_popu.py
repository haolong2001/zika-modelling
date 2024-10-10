
import sys
import os
current_directory = os.getcwd()
sys.path.insert(0, current_directory)
from modules.Person import Person
import pickle

num_persons = 6052709
persons_dict = {}
for i in range(num_persons):
    age = -1
    gender = -1 
    yearOfadd = 2024

    person_instance = Person(index=i, 
                             age=age, 
                             yearOfadd = yearOfadd,                           
                             gender=gender)
    persons_dict[i] = person_instance

with open('persons.pkl', 'wb') as file:
    pickle.dump(persons_dict, file) 