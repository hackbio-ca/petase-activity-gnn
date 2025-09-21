import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer
import umap
import matplotlib.pyplot as plt


activity_file = r"C:\Users\divya\OneDrive\Desktop\metadata_project\activity.csv"
activity_df = pd.read_csv(activity_file)

activity_df.columns = [col.strip().lower() for col in activity_df.columns]


col_map = {}
if 'srr_id' in activity_df.columns:
    col_map['srr_id'] = 'SRR_id'
elif 'srrid' in activity_df.columns:
    col_map['srrid'] = 'SRR_id'
if 'gene' in activity_df.columns:
    col_map['gene'] = 'gene'
if 'readout_value' in activity_df.columns:
    col_map['readout_value'] = 'readout_value'

activity_df.rename(columns=col_map, inplace=True)

# Verify columns
expected_cols = ['SRR_id', 'gene', 'readout_value']
for col in expected_cols:
    if col not in activity_df.columns:
        raise ValueError(f"activity.csv must have column '{col}'")

gene_matrix = activity_df.pivot(index='SRR_id', columns='gene', values='readout_value').fillna(0)

avg_activity = activity_df.groupby('SRR_id')['readout_value'].mean()
avg_activity = avg_activity.loc[gene_matrix.index]

model = SentenceTransformer('all-MiniLM-L6-v2')
gene_strings = gene_matrix.astype(str).agg(' '.join, axis=1)
embeddings = model.encode(gene_strings.tolist(), show_progress_bar=True)

reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
umap_embeds = reducer.fit_transform(embeddings)

plt.figure(figsize=(10, 8))
sc = plt.scatter(umap_embeds[:,0], umap_embeds[:,1], c=avg_activity.values, cmap='viridis', s=50)
plt.colorbar(sc, label="Average Readout Value")
plt.title("UMAP of SRRs Colored by Average Activity")
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.tight_layout()
plt.show()

