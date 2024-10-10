
import psutil
import pickle
import sys
import os
import multiprocessing
# current_directory = os.getcwd()
# sys.path.insert(0, current_directory)

# # Print system memory
# def print_system_memory():
#     memory = psutil.virtual_memory()
#     print(f"Total memory: {memory.total / (1024 ** 3):.2f} GB")
#     print(f"Available memory: {memory.available / (1024 ** 3):.2f} GB")
#     print(f"Used memory: {memory.used / (1024 ** 3):.2f} GB")
#     print(f"Percentage used: {memory.percent}%")

# # Print memory usage before loading the file
# print("Memory usage before loading the file:")
# print_system_memory()

# # Load the dictionary
# with open("population_data\persons_no_demo.pkl", 'rb') as file:
#     persons_dict = pickle.load(file)

# # Print memory usage after loading the file
# print("\nMemory usage after loading the file:")
# print_system_memory()

# # Get the size of the dictionary
# dict_size = sys.getsizeof(persons_dict)
# print(f"\nMemory size of persons_dict: {dict_size / (1024 ** 2):.2f} MB")

def read_file(file_path, result_queue):
    # Read the file
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    
    # Put the loaded data into the result queue
    result_queue.put(data)

if __name__ == "__main__":
    # Define the file path
    file_path = 'population_data/persons_no_demo.pkl'
    
    # Get initial memory usage
    initial_memory = psutil.virtual_memory().used
    
    # Number of processes to spawn
    num_processes = 5
    
    # Create a multiprocessing Queue to collect results
    result_queue = multiprocessing.Queue()
    
    # Create and start processes
    processes = []
    for _ in range(num_processes):
        p = multiprocessing.Process(target=read_file, args=(file_path, result_queue))
        p.start()
        processes.append(p)
    
    # Wait for all processes to finish
    for p in processes:
        p.join()
    
    # Get final memory usage
    final_memory = psutil.virtual_memory().used
    
    # Calculate memory usage change
    memory_change = final_memory - initial_memory
    
    # Get the size of the loaded data
    loaded_data_size = os.path.getsize(file_path) / (1024 * 1024)  # Convert bytes to MB
    
    # Get the average memory usage change per process
    avg_memory_change_per_process = memory_change / num_processes
    
    print("Memory usage before loading the file:")
    print(f"Initial memory: {initial_memory / (1024 * 1024):.2f} MB")
    print("\nMemory usage after loading the file in parallel:")
    print(f"Final memory: {final_memory / (1024 * 1024):.2f} MB")
    print(f"Memory change: {memory_change / (1024 * 1024):.2f} MB")
    print(f"Size of the loaded data: {loaded_data_size:.2f} MB")
    print(f"Average memory change per process: {avg_memory_change_per_process / (1024 * 1024):.2f} MB")
