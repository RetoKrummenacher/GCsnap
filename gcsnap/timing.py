import time
import pandas as pd

class Timer:
    """ 
    Methods and attributes to time the execution of a code block.

    Attributes:
        name (str): The name of the timer.
        start_time (float): The start time of the timer.
        end_time (float): The end time of the timer.
        elapsed_time (float): The elapsed time of the timer.
    """

    def __init__(self, name: str):
        """
        Initialize the Timer object.

        Args:
            name (str): The name of the timer.
        """        
        self.name = name
        self.start_time = None
        self.end_time = None
        self.elapsed_time = None

    def start(self) -> None:
        """
        Start the timer by setting the start time.
        """        
        self.start_time = time.time()
        self.end_time = None  # Reset end time when restarted

    def stop(self) -> None:
        """
        Stop the timer by setting the end time and calculating the elapsed time.
        """        
        if self.start_time is not None:
            self.end_time = time.time()
            self.elapsed_time = self.end_time - self.start_time

    def get_elapsed_time(self) -> float:
        """
        Get the elapsed time of the timer.

        Returns:
            float: The elapsed time of the timer in seconds.
        """        
        if self.elapsed_time is not None:
            return self.elapsed_time
        elif self.start_time is not None:
            return time.time() - self.start_time
        else:
            return 0

class Timing:
    """
    Methods and attributes to time the execution of multiple code blocks.

    Attributes:
        timers (dict): The dictionary with the timers.
        counter (int): The counter to keep track of the timers.
    """

    def __init__(self):
        """
        Initialize the Timing object.
        """        
        self.timers = {}
        self.counter = 0

    def timer(self, name: str) -> Timer:
        """
        Initialize a new Timer object, start it and add it to the timers dictionary.

        Args:
            name (str): Name of the new Timer object.

        Returns:
            Timer: The new Timer object.
        """        
        self.counter += 1
        timer = Timer(name)
        self.timers[self.counter] = timer
        timer.start()
        return timer
    
    def get_timing_data(self) -> list:
        """
        Get the information about name and elapsed time of all timers in the timers dictionary.

        Returns:
            list: The timing data of all timers.
        """            
        data = []
        for counter, timer in self.timers.items():
            data.append({
                "Step": timer.name,
                # "Start Time": timer.start_time,
                # "End Time": timer.end_time,
                "Time (sec) ": timer.get_elapsed_time()
            })
        return data

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert the timing data to a pandas DataFrame.

        Returns:
            pd.DataFrame: The timing data as a pandas DataFrame.
        """        
        
        data = self.get_timing_data()
        return pd.DataFrame(data)

    def to_csv(self, file_name: str) -> None:
        """
        Save the timing data to a CSV file.

        Args:
            file_name (str): Name of the CSV file.
        """        
        df = self.to_dataframe()
        df.to_csv(file_name, index=False)