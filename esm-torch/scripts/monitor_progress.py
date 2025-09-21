# Create: scripts/monitor_progress.py
import time
import os
import pandas as pd

def monitor_progress(output_file, progress_file, target_total=1000000):
    """Monitor live prediction progress"""
    
    while True:
        try:
            # Read current progress
            if os.path.exists(progress_file):
                with open(progress_file, 'r') as f:
                    content = f.read()
                print(f"\n=== Progress Update ===")
                print(content)
            
            # Read current results
            if os.path.exists(output_file):
                df = pd.read_csv(output_file)
                current_count = len(df)
                
                if current_count > 0:
                    print(f"Results so far: {current_count}/{target_total} ({100*current_count/target_total:.1f}%)")
                    print(f"Latest predictions:")
                    print(df.tail(3)[['sequence_id', 'predicted_activity']].to_string(index=False))
                    
                    # Basic stats
                    activities = df['predicted_activity']
                    valid_activities = activities[activities != -999]  # Exclude error markers
                    if len(valid_activities) > 0:
                        print(f"Activity stats: mean={valid_activities.mean():.2f}, "
                              f"range=[{valid_activities.min():.2f}, {valid_activities.max():.2f}]")
            
            time.sleep(30)  # Update every 30 seconds
            
        except KeyboardInterrupt:
            print("Monitoring stopped")
            break
        except Exception as e:
            print(f"Monitor error: {e}")
            time.sleep(30)

if __name__ == "__main__":
    import sys
    output_file = sys.argv[1] if len(sys.argv) > 1 else "predictions.csv"
    progress_file = output_file.replace('.csv', '_progress.txt')
    monitor_progress(output_file, progress_file)