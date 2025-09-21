import matplotlib.pyplot as plt
import numpy as np

from src.data.data_loader import load_esm_embeddings_with_activity

# Load activities
embeddings, activities, sequence_ids = load_esm_embeddings_with_activity()

# Filter out proteins with 0 activity
non_zero_mask = activities > 0
filtered_activities = activities[non_zero_mask]
filtered_embeddings = embeddings[non_zero_mask]

# Handle sequence_ids based on its type
if isinstance(sequence_ids, np.ndarray):
    filtered_sequence_ids = sequence_ids[non_zero_mask]
else:
    # If sequence_ids is a list, convert mask to indices
    non_zero_indices = np.where(non_zero_mask)[0]
    filtered_sequence_ids = [sequence_ids[i] for i in non_zero_indices]

# Log transform the non-zero activities
log_activities = np.log(filtered_activities)

print(f"Original dataset size: {len(activities)}")
print(f"After removing zero activities: {len(filtered_activities)}")
print(f"Removed {len(activities) - len(filtered_activities)} proteins with zero activity")

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(filtered_activities, bins=30)
plt.title('Non-zero Activity Distribution')
plt.xlabel('Activity')
plt.ylabel('Frequency')

plt.subplot(1, 2, 2)
plt.hist(log_activities, bins=30)
plt.title('Log-transformed Activity (Non-zero)')
plt.xlabel('Log(Activity)')
plt.ylabel('Frequency')

plt.tight_layout()
plt.savefig('filtered_activity_distribution.png')
# plt.show()  # Comment out or remove this line