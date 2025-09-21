import csv
import os
import time
from pathlib import Path

class LiveResultsWriter:
    def __init__(self, output_path, buffer_size=100):
        self.output_path = output_path
        self.buffer_size = buffer_size
        self.buffer = []
        self.total_written = 0
        
        # Create output file with headers if it doesn't exist
        if not os.path.exists(output_path):
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sequence_id', 'predicted_activity', 'sequence_length', 'timestamp'])
        
        # Create progress file
        self.progress_file = output_path.replace('.csv', '_progress.txt')
        
    def add_result(self, sequence_id, activity, seq_length):
        result = {
            'sequence_id': sequence_id,
            'predicted_activity': activity,
            'sequence_length': seq_length,
            'timestamp': time.time()
        }
        self.buffer.append(result)
        
        if len(self.buffer) >= self.buffer_size:
            self.flush()
    
    def flush(self):
        if not self.buffer:
            return
            
        # Append to CSV
        with open(self.output_path, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['sequence_id', 'predicted_activity', 'sequence_length', 'timestamp'])
            writer.writerows(self.buffer)
        
        self.total_written += len(self.buffer)
        self.buffer = []
        
        # Update progress
        self.update_progress()
    
    def update_progress(self):
        with open(self.progress_file, 'w') as f:
            f.write(f"Processed: {self.total_written}\n")
            f.write(f"Last update: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            if hasattr(self, 'start_time') and self.total_written > 0:
                elapsed = time.time() - self.start_time
                rate = self.total_written / elapsed
                f.write(f"Rate: {rate:.2f} sequences/second\n")
                f.write(f"Elapsed: {elapsed/3600:.2f} hours\n")
    
    def start_timing(self):
        self.start_time = time.time()
    
    def finish(self):
        self.flush()
        print(f"Completed! Total processed: {self.total_written}")