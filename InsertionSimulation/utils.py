#utils.py

import time #wrapper
import psutil #wrapper
import os #wrapper

def elapsed_since(start):
	return time.strftime("%H:%M:%S", time.gmtime(time.time() - start))

def get_process_memory():
	process = psutil.Process(os.getpid())
	return process.memory_info().rss

def profile(func):
	'''
	Wrapper to monitor the resources used by key steps.
	'''
	def wrapper(*args, **kwargs):
		mem_before = get_process_memory()
		start = time.time()
		result = func(*args, **kwargs)
		elapsed_time = elapsed_since(start)
		mem_after = get_process_memory()
		print("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
			func.__name__,
			mem_before, mem_after, mem_after - mem_before,
			elapsed_time))
		return result
	return wrapper

def get_chromosome(insert_position, chromosome_dir):
	'''
	Checks within which chromosome the random insertion landed.
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
	for key, value in roi_dict.items():
		for i in range(n_barcodes):
			new_key = f"{key}_{i}"
			new_dict[new_key] = value
	return new_dict