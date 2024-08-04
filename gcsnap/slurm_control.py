import os

from dask.distributed import Client
from dask_jobqueue import SLURMCluster

class SLURMcontrol():

    def __init__(self, total_cores: int):
        self.total_cores = total_cores
        self.cluster_kwargs = self.read_slurm_env_vars()

    def read_slurm_env_vars(self) -> dict:
        # SLURM configuration in run.job file
        # example, 2 nodes of the sei: 128 cores
        #!/bin/bash
        #SBATCH --exclusive
        #SBATCH --time=0:30:00
        #SBATCH --qos=30min      # Quality of service
        #SBATCH --ntasks-per-node=128
        #SBATCH --cpus-per-task=1        
        #SBATCH --hint=nomultithread        
        #SBATCH --partition=sei 
        #SBATCH --nodelist=sei[50-51]        # not needed, SLURM handles this

        # defined in script calling sbatch run.job for the experiments
        #--job-name=your_job_name  
        #--output=your_job_name.out
        #--nodes=2

        # number of CPU cores is set in config.yaml (or in experiments with CLI --n-cpu)
        # SLURM handles the rest

        # Read SLURM environment variables
        partition = os.getenv('SLURM_JOB_PARTITION')
        job_name = os.getenv('SLURM_JOB_NAME')
        nodes = int(os.getenv('SLURM_JOB_NUM_NODES'))
        ntasks_per_node = int(os.getenv('SLURM_NTASKS_PER_NODE'))
        walltime = os.getenv('SLURM_TIMELIMIT')

        slurm_cluster_kwargs = {
            'partition': partition,
            'account': job_name,
            'cores': 1,
            'memory': '4GB', # memory per worker, a node has 512, meaning 4GB per CPU core (128)
            'walltime': walltime,
            'job_extra': ['--hint=nomultithread', '--exclusive', '--qos=30min']
        }

        self.total_cores = nodes * ntasks_per_node

        return slurm_cluster_kwargs
    
    def start_cluster(self) -> SLURMCluster:
        cluster = SLURMCluster(**self.cluster_kwargs)
        cluster.scale(self.total_cores)

        return cluster