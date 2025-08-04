#utils.py
import logging
import time
import os
from functools import wraps
import psutil

def track_usage(label=None):
    cpu_usage = psutil.cpu_percent(interval=None)
    memory_usage = psutil.virtual_memory().percent
    logging.info(f"Process: {label}, CPU: {cpu_usage}%, Memory: {memory_usage}%")

def setup_logging(filename):
    logging.basicConfig(
        filename=filename,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filemode='a'
    )

def profile_iteration(func):
    """Decorator to profile CPU, memory, and time usage for each iteration."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        pid = os.getpid()
        process = psutil.Process(pid)
        start_time = time.time()
        start_memory = process.memory_info().rss / (1024 * 1024)  # in MB
        start_cpu = process.cpu_percent(interval=None)
        result = func(*args, **kwargs)  # Run the function
        end_time = time.time()
        end_memory = process.memory_info().rss / (1024 * 1024)  # in MB
        end_cpu = process.cpu_percent(interval=None)
        logging.info(f"Iteration {args[0]}: Time={end_time - start_time:.2f}s, "
              f"Memory Used={end_memory - start_memory:.2f}MB, CPU={end_cpu - start_cpu}%")
        return result
    return wrapper

def get_chromosome(insert_position, chromosome_dir):
    '''
    Checks within which chromosome the insertion landed.
    '''
    for chromosome, coordinates in chromosome_dir.items():
        start, end = coordinates
        if start <= insert_position <= end:
            return chromosome
    return None

def check_barcoding(n_barcodes):
    return n_barcodes >= 1

def barcode_genome(chromosome_dir, barcode):
    '''
    Adds prefix to chromosome dict
    '''
    chromosome_dir = {f'Barcode_{barcode}_{k}': v for k, v in chromosome_dir.items()}
    return chromosome_dir

def roi_barcoding(roi_dict, n_barcodes):
    '''
    Add _i suffix to each key in the dictionary for each index i in the specified range.
    '''
    new_dict = {}
    print(roi_dict)
    for key, value in roi_dict.items():
        for i in range(n_barcodes):
            new_key = f"{key}_{i}"
            new_dict[new_key] = value
    return new_dict