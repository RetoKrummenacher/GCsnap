import os

from dask.distributed import LocalCluster
from dask_jobqueue import SLURMCluster

from gcsnap.rich_console import RichConsole 

class DaskSlurmCluster():


    def __init__(self, nodes: int, memory_per_node: int, processes_per_node: int, 
                 tool: str = 'dask', identifier: str = None, result_path: str = None):
        """
        Initialize the DaskSlurmCluster object.

        Args:
            nodes (int): Number of nodes to use.
            memory_per_node (int): Available memory per node in GB.
            processes_per_node (int): Processes per node.
            tool (str, optional): The parallel processing tool to use. Defaults to 'dask'.
            identifier (str, optional): Identifier for the SLURM output. Defaults to None.
            result_path (str, optional): Where to store the SLURM job scripts and the output. Defaults to None.
        """        
        self.nodes = nodes
        self.memory_per_node = '{}GB'.format(memory_per_node)
        self.processes_per_node = processes_per_node
        self.tool = tool
        self.identifier = identifier
        self.result_path = result_path

        self.console = RichConsole()

        # set output identifier for Dask workers (each produces its own out file)
        if identifier is None:
            self.identifier = 'dask_worker'
        if result_path is None:
            self.result_path = os.getcwd()
        
        if not os.path.exists(os.path.join(self.result_path, 'dask_workers')):
            os.makedirs(os.path.join(self.result_path, 'dask_workers'))

        self.worker_out = os.path.join(self.result_path, 'dask_workers', self.identifier)
        
        # when using only one node, it would be possible to use LocalCluster
        if self.tool == 'dask_local':
            self.start_local_cluster()
        elif self.tool == 'dask':
            # read the configuration file
            self.cluster_kwargs = self.read_slurm_env_vars()
            
            # start and return the cluster object
            self.start_cluster()
        else:
            self.console.print_error('The tool {} is not supported. Select from "dask, mpi, dask local"'.format(self.tool))
            exit(1)

    def get_cluster(self) -> SLURMCluster:
        """
        Get the SLURM cluster object

        Returns:
            SLURMCluster: The SLURM cluster object.
        """        
        return self.cluster

    def read_slurm_env_vars(self) -> dict:
        """
        SLURM job script configuration. This is submitted when using cluster.scale() in start_cluster().

        Returns:
            dict: The SLURM job script configuration.
        """        
        # SLURM configuration that was produced with this class with nodes=3, cores_per_node=4, processes_per_node=2
        #SBATCH -J dask-worker
        #SBATCH -n 1
        #SBATCH --cpus-per-task=4
        #SBATCH --mem=30G
        #SBATCH -t 00:30:00
        #SBATCH --hint=nomultithread
        #SBATCH --exclusive
        #SBATCH --partition=xeon
        #ml Python
        #/opt/apps/easybuild/software/Python/3.10.8-GCCcore-12.2.0/bin/python3 -m distributed.cli.dask_worker tcp://10.34.59.1:42774 
        #--name dummy-name --nthreads 2 --memory-limit 14.90GiB --nworkers 2 --nanny --death-timeout 60
        
        # summary:
            # cluster.scale(jobs=3): request 3 times the slurm_cluster_kwargs settings from SLURM, we end up with 3 dask worker nodes
                # this in contrast to cluster.scale(n=3), which asks for 3 dask workers <> number of worker nodes
                # https://dask.discourse.group/t/dask-jobqueue-slurmcluster-multi-threaded-workloads-and-the-effect-of-setting-cores/2310
            # 'cores': Is the number of threads (dask --nthreads) 
            # 'processes': How many processes/workers (dask --nworkers) are on thop of those --nthreads. Here its 2, which is also the n_workers argument from LocalCluster
                # As we have 4 CPU cores, we have 2 threads (threads_per_worker argument from LocalCluster): 2 process * 2 threads = 4 CPU cores
            # 'job_cpu': This is actually the number of CPU Cores assigend. It sets --cpus-per-taks, 
                # despite --nthreads and --nworkers together might lead to a different value for --cpus-per-taks: we suggest no using it
            # The total number of actual nodes that are used is arbitrary and determined from Dask      
            
            # Using Threads in gerneral is an issues, as Dask does not release the GIL:
                # https://dask.discourse.group/t/how-does-dask-avoid-gil-in-multi-threading/1846

        slurm_cluster_kwargs = {
            'cores': self.processes_per_node, # threads per job (if cores < processes --nthreads=1)
            'processes': self.processes_per_node, # cut up job in this many processes: 
                # Default process ~= sqrt(cores) 
                # so that the number of processes and the number of threads per process is roughly the same.
                # here we use as many processes as cores assigned, so each process running 1 thread as those ar bound by GIL
            'memory': self.memory_per_node, # memory per node, 0GB doesn't work despite SLURM would react to it, but dask fails
            # https://discourse.pangeo.io/t/when-using-dask-slurmcluster-can-i-avoid-passing-the-memory-argument/2362/7
            'walltime': '00:30:00',
            'job_extra_directives': ['--hint=nomultithread', 
                                     '--exclusive', 
                                     '--partition=xeon',
                                     f'--output={self.worker_out}_%j.out',
                                     ],
            'job_script_prologue' : [
                'ml Python'
                'source ~/miniconda3/etc/profile.d/conda.sh',
                'source activate GCsnap'  # Activate your conda environment
            ]            
        }

        return slurm_cluster_kwargs
    
    def start_cluster(self) -> SLURMCluster:
        """
        Start the SLURM cluster.

        Returns:
            SLURMCluster: The SLURM cluster object.
        """        
        cluster = SLURMCluster(**self.cluster_kwargs)
        # scale the cluster to as many nodes as desired. #self.nodes jobs are submitted forming a worker pool
        cluster.scale(jobs=self.nodes)
        # export the job script that Dask produces
        with open(os.path.join(self.result_path, 'dask_workers', 'dask_job_script.job'), 'w') as f:
            f.write(cluster.job_script())

        self.cluster = cluster
    
    def close_cluster(self) -> None:
        """
        Close the SLURM cluster.
        """        
        # Close the cluser (client is closed in with statement)
        self.cluster.close()    
        
    def start_local_cluster(self):
        """
        Start a Dask Local Cluster to be used on one node.
        """        
        cluster = LocalCluster(n_workers=self.processes_per_node,  # number of workers
                               threads_per_worker=1 # threads per worker 
                               )
        self.cluster = cluster