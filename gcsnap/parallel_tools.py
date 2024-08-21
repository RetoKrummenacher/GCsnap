
from typing import Callable

from dask.distributed import Client

from concurrent.futures import as_completed
from mpi4py.futures import MPIPoolExecutor

from gcsnap.slurm_control import DaskSlurmCluster
from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

class ParallelTools:

    """
    Methods and attributes to handle parallel processing. The class is initialized with a Configuration object.

    Attributes:
        config (Configuration): The Configuration object with the parsed arguments.
        n_nodes (int): The number of nodes to use for parallel processing.
        n_cpu (int): The number of CPUs per node to use for parallel processing.
        workers (int): The number of workers to use for parallel processing with mpoi4py.
        tool (str): The parallel processing tool to use.
        dask_scheduler (str): The Dask scheduler to use.
        memory_per_node (int): The memory per node to use for parallel processing.
        cluster (DaskSlurmCluster): The Dask cluster object.
        console (RichConsole): The RichConsole object.
    """

    def __init__(self, config: Configuration):
        """
        Initialize the ParallelTools object.

        Args:
            config (Configuration): The Configuration object with the parsed arguments.
        """        
        self.config = config
        self.n_nodes = config.arguments['n_nodes']['value']
        self.n_cpu = config.arguments['n_cpu_per_node']['value']
        self.tool = config.arguments['parallel-tool']['value']
        self.dask_scheduler = config.arguments['dask_scheduler']['value']
        self.memory_per_node = config.arguments['memory_per_node']['value']
        self.cluster = None

        # set up the cluster in case Dask is requested
        if self.tool == 'dask' and self.dask_scheduler is None:
            self.cluster = DaskSlurmCluster(self.n_nodes, self.memory_per_node, self.n_cpu, self.n_cpu)
        elif self.tool == 'dask':
            self.cluster = self.dask_scheduler
        elif self.tool == 'mpi':
            # for mpi, we actually use - 1, as one is running the main thread
            self.workers = (self.nodes * self.n_cpu) - 1
            
        self.console = RichConsole()

        # Assign this instance to the class-level variable
        # Now self is available in the static method as on class level
        # usually instance is not on class level
        ParallelTools._instance = self

    @staticmethod
    def process_wrapper(parallel_args: list[tuple], func: Callable) -> list:
        """
        A static method that calls the process_wrapper method of the stored instance.
        This is called from the modules using parallel processing.

        Args:
            parallel_args (list[tuple]):A list of tuples, where each tuple contains the arguments for the function.
            func (Callable): The function to apply to the arguments.

        Returns:
            list: A list of results from the function applied to the arguments in the order they are provided.            
        """
        if ParallelTools._instance is None:
            raise RuntimeError("ParallelTools instance has not been initialized.")
        
        return ParallelTools._instance.process_wrapper(parallel_args, func)        

    def process_wrapper(self, parallel_args: list[tuple], func: Callable) -> list:
        """
        Process wrapper to select the desired parallel tool.

        Args:
            parallel_args (list[tuple]):A list of tuples, where each tuple contains the arguments for the function.
            func (Callable): The function to apply to the arguments.

        Returns:
            list: A list of results from the function applied to the arguments in the order they are provided.
        """
        if self.tool == 'mpi':
            result_list = self.mpiprocess_wrapper(parallel_args, func)
        elif self.tool == 'dask':
            result_list = self.daskprocess_wrapper(parallel_args, func)
        else:
            self.console.print_error('Specified parallel tool {} not supported'.format(tool))

        return result_list

    def daskprocess_wrapper(self, parallel_args: list[tuple], func: Callable) -> list: 
        """
        Apply a function to a list of arguments using Dask with processes. The arguments are passed as tuples
        and are unpacked within the function. Dask is executed asynchronusly, but with the order of the results guaranteed.

        Args:
            parallel_args (list[tuple]): A list of tuples, where each tuple contains the arguments for the function.
            func (Callable): The function to apply to the arguments.

        Returns:
            list: A list of results from the function applied to the arguments in the order they are provided.
        """      
        # For a Dask Local Cluster on one machine:
        # n_workers: number of processes (Defaults to 1)
        # threads_per_worker: threads in each process (Defaults to None)
            # i.e. it uses all available cores.    
            #  with Client(n_workers=n_processes, threads_per_worker=1) as client:
                # futures = client.compute(delayed_results)  # Start computation in the background
                # result_list = client.gather(futures)  # Block until all results are ready    

        if isinstance(self.cluster, DaskSlurmCluster):
            # get the cluster object from the cluster instance
            client = Client(self.cluster.get_cluster())
        else:
            client = Client(scheduler_file = self.cluster)

        # # Check the number of workers
        # print(client.scheduler_info())
        # # this will not work: https://dask.discourse.group/t/how-to-retrieve-the-requested-number-of-cores/798/6
        # print('Cores {}'.format(sum(w['cores'] for w in client.scheduler_info()['workers'].values())))
        # print('Threads {}'.format(sum(w['nthreads'] for w in client.scheduler_info()['workers'].values())))
        
        # list of delayed objects to compute
        futures = [client.submit(func, arg) for arg in parallel_args]
        #futures = client.map(func, parallel_args)

        # Collect results as they complete: https://docs.dask.org/en/stable/futures.html
        # result_list = [future.result() for future in as_completed(futures)]
        # client gather should be more efficient as done concurrently
        # but as we collect them as completed, this is not sure and the documentation is not clear
        result_list = client.gather(futures)

        client.close()
        
        return result_list

    def mpiprocess_wrapper(self, parallel_args: list[tuple], func: Callable) -> list:
        """
        Apply a function to a list of arguments using ProcessPoolExecutor. The arguments are passed as tuples
        and are unpacked within the function. As completed is used to get the results in the order they finish.

        Args:
            parallel_args (list[tuple]): A list of tuples, where each tuple contains the arguments for the function.
            func (Callable): The function to apply to the arguments.

        Returns:
            list: A list of results from the function applied to the arguments in the order they finish.
        """    
        # Same as with ProcessPoolExecutor from cuncurrent.futures
        # https://mpi4py.readthedocs.io/en/stable/mpi4py.futures.html#parallel-tasks
        # with MPIPoolExecutor(max_workers = workers) as executor:
        with MPIPoolExecutor(max_workers = self.workers) as executor:        
            # get numbers of worker
            #print('Number of workers given: {}'.format(executor.num_workers))
            #print('Number of workers asked: {}'.format(workers))
            futures = [executor.submit(func, arg) for arg in parallel_args]
            result_list = [future.result() for future in as_completed(futures)]

        return result_list
    
    def close_cluster(self):
        """
        Close the Dask cluster.
        """
        self.cluster.close_cluster()