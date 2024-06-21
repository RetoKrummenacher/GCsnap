import os
from typing import Callable, Union

from threading import Thread, Lock

from multiprocessing import Pool as ProcessPool
from multiprocessing.pool import ThreadPool

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed

import dask
from dask.distributed import Client
from dask import config


# Parallel wrappers for any function:
# ------------------------------------------------------
def sequential_wrapper(cores: int, parallel_args: list, func: Callable) -> list:
    
    result_list = []
    
    for arg in parallel_args:
        result_list.append(func(arg))   
        
    return result_list


def threading_wrapper(n_threads: int, parallel_args: list, func: Callable) -> list:
    def worker(index):
        result = func(parallel_args[index])
        with lock:
            results_list[index] = result

    results_list = [None] * len(parallel_args)
    threads = []
    lock = Lock()

    for i in range(len(parallel_args)):
        t = Thread(target=worker, args=(i,))
        threads.append(t)
        t.start()

    for t in threads:
        t.join()

    return results_list


def processpool_wrapper(n_processes: int, parallel_args: list, func: Callable, return_results: bool = True) -> Union[list,None]:       
    pool = ProcessPool(processes = n_processes)
    
    # imap and map take just one argument, hence unpacking within function
    # imap is the lazy version of map,
    # _unordered is without controlling the result order. (map_async)
    # multiple arguments use startmap(),
    # but starmap() and map() may be inefficient with large lists:
        # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.map
    # great explanation about:
        # https://stackoverflow.com/questions/26520781/multiprocessing-pool-whats-the-difference-between-map-async-and-imap
        
    results = pool.map_async(func, parallel_args)
    
    pool.close() # close pool for further processes
    pool.join() # wait until all have finished.
    
    result_list = results.get()
    if return_results:
        # Convert to a list of results from the MapResult object
        return result_list


def threadpool_wrapper(n_threads: int, parallel_args: list, func: Callable) -> list:
    pool = ThreadPool(n_threads)

    # Use map_async to apply the function asynchronously
    results = pool.map_async(func, parallel_args)

    pool.close() # close pool for further processes
    pool.join() # wait until all have finished.

    # Convert to a list of results from the MapResult object
    result_list = results.get()
    return result_list


def futures_thread_wrapper(n_threads: int, parallel_args: list, func: Callable) -> list:
    # executor.submit() to submit each task to the executor. 
    # This returns a Future object for each task. You then use as_completed() 
    # to get an iterator that yields futures as they complete. 
    # Create a ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [executor.submit(func, arg) for arg in parallel_args]
        # list comprehension with future.result() to get the results of the futures. 
        # This means that your code will start processing the results as soon as they become available, in the order they finish.
        result_list = [future.result() for future in as_completed(futures)]
        # but it will wait till all are done to go on

    return result_list

def futures_process_wrapper(n_processes: int, parallel_args: list, func: Callable) -> list: 
    # Same as with futures and threads
    # Create a ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [executor.submit(func, arg) for arg in parallel_args]
        result_list = [future.result() for future in as_completed(futures)]

    return result_list


def daskthread_wrapper(n_threads: int, parallel_args: list, func: Callable) -> list:    
    # n_workers: number of processes (Defaults to 1)
    # threads_per_worker: threads in each process (Defaults to None)
        # i.e. it uses all available cores.        

    # list of delayed objects to compute
    delayed_results = [dask.delayed(func)(*arg) for arg in parallel_args]

    # dask.compute, Dask will automatically wait for all tasks to finish before returning the results   
    # even though a delayed object is usesd, the computation starts right away when using copmpute() 
    with Client(threads_per_worker=n_threads,n_workers=1) as client:
        futures = client.compute(delayed_results)  # Start computation in the background
        result_list = client.gather(futures)  # Block until all results are ready
    
    return result_list


def daskprocess_wrapper(n_processes: int, parallel_args: list, func: Callable) -> list:   
    # n_workers: number of processes (Defaults to 1)
    # threads_per_worker: threads in each process (Defaults to None)
        # i.e. it uses all available cores.        

    # list of delayed objects to compute
    delayed_results = [dask.delayed(func)(*arg) for arg in parallel_args]

    # dask.compute, Dask will automatically wait for all tasks to finish before returning the results   
    # even though a delayed object is usesd, the computation starts right away when using copmpute() 
    with Client(n_workers=n_processes, threads_per_worker=1) as client:
        futures = client.compute(delayed_results)  # Start computation in the background
        result_list = client.gather(futures)  # Block until all results are ready
    
    return result_list
# ------------------------------------------------------



# Downloading
# -----------
# pip install wget
import wget
import requests
# pip install asyncio
import asyncio # single threaded coroutines
# pip install aiohttp
import aiohttp  # asynchronous asyncio

# asyncio
async def download_asyncio(url: str, file_path: str, file_name: str) -> None:
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as resp:
            with open(os.path.join(file_path, file_name), 'wb') as f:
                f.write(await resp.read())    

async def wrapper_download_asyncio(args_list: list[tuple[str,str,str]]) -> None:
    tasks = [download_asyncio(args[0], args[1], args[2]) for args in args_list]
    await asyncio.gather(*tasks)
    
# wget
def download_wget(url: str, file_path: str, file_name: str) -> None:
    wget.download(url, out=os.path.join(file_path, file_name))
        
def wrapper_download(args_list: list[tuple[str,str,str]]) -> None:
    for args in args_list:
        download_wget(args[0], args[1], args[2])    
        
# request
def download_request(url: str, file_path: str, file_name: str) -> None:
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(os.path.join(file_path, file_name), 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

def wrapper_download_request(args_list: list[tuple[str,str,str]]) -> None:
    for args in args_list:
        download_request(args[0], args[1], args[2])
        


# Styling & logging
# -----------------
import logging # already sets the loggin process

logging.basicConfig(
    filename='gcsnap.log',
    filemode='w',
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('gcsnap')


