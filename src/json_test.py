import json

# Your dictionary
my_dict = {1: [1, 2]}

# Specify the file name
file_name = 'my_dict.json'

# Open the file in write mode and use json.dump() to write the dictionary to the file
with open(file_name, 'w') as file:
    json.dump(my_dict, file)

