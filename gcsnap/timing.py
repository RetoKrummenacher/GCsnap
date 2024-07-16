import time
import pandas as pd

class Timer:
    def __init__(self, name):
        self.name = name
        self.start_time = None
        self.end_time = None
        self.elapsed_time = None

    def start(self):
        self.start_time = time.time()
        self.end_time = None  # Reset end time when restarted

    def stop(self):
        if self.start_time is not None:
            self.end_time = time.time()
            self.elapsed_time = self.end_time - self.start_time

    def get_elapsed_time(self):
        if self.elapsed_time is not None:
            return self.elapsed_time
        elif self.start_time is not None:
            return time.time() - self.start_time
        else:
            return 0

class Timing:
    def __init__(self):
        self.timers = {}
        self.counter = 0

    def timer(self, name):
        self.counter += 1
        timer = Timer(name)
        self.timers[self.counter] = timer
        timer.start()
        return timer
    
    def get_timing_data(self):
        data = []
        for counter, timer in self.timers.items():
            data.append({
                "Step": timer.name,
                # "Start Time": timer.start_time,
                # "End Time": timer.end_time,
                "Time (sec) ": timer.get_elapsed_time()
            })
        return data

    def to_dataframe(self):
        data = self.get_timing_data()
        return pd.DataFrame(data)

    def to_csv(self, file_name):
        df = self.to_dataframe()
        df.to_csv(file_name, index=False)