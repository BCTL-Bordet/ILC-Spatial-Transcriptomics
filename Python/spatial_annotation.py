import numpy as np
import pandas as pd
import json
import cv2
import matplotlib.pyplot as plt
from PIL import Image

# Allow processing of large images
Image.MAX_IMAGE_PIXELS = None

# Define sample name
sample = "ST41"

# Define file paths. All inpute files are present on Zenodo repository ("spacerange" outputs, link to Zenodo in README file)
# Define file names (relative paths for portability)
image_file = "tissue_hires_image.png"
mask_file = f"{sample}.png"
positions_file = "tissue_positions_list.csv"
scalefactors_file = "scalefactors_json.json"
out_file = "tissue_positions_list_annotation.csv"
out_adjusted_file = f"percentages_tissue_{sample}_percentages.csv"

# Define color mapping (RGB -> tissue class)
color_map = {
    "11201": "Tumor",
    "000": "Necrosis",
    "00128": "Fat_tissue",
    "2332330": "High_TILs_stroma",
    "255153128": "Cellular stroma",
    "232209187": "Acellular stroma",
    "22000": "Vessels",
    "110360": "Artefact",
    "153128230": "Canal_galactophore",
    "12812826": "Nodule_lymphoid",
    "204255204": "In_situ",
    "77128128": "Nerve",
    "19665127": "Lymphocyte",
    "64229246": "Hole",
    "655276": "Microcalcification",
    "255255255": "Out",
    "2550255": "Apocrine_metaplasia"
}

# Function to compute bounding box for a spot
def get_bb_from_spot(spot, diameter):
    radius = diameter // 2
    return (spot[0] - radius, spot[1] - radius), (spot[0] + radius, spot[1] + radius)

# Load images and metadata
image = np.array(Image.open(image_path))
mask = np.array(Image.open(mask_path).convert("RGB"))
df_positions = pd.read_csv(positions_path, header=None)
with open(scalefactors_path, 'r') as f:
    scalefactors = json.load(f)

# Extract scale factor and spot diameter
scale_factor = scalefactors['tissue_hires_scalef']
diam = scalefactors['spot_diameter_fullres']

# Resize mask to match tissue image dimensions
new_shape = (image.shape[1], image.shape[0])
mask_resized = cv2.resize(mask, new_shape, interpolation=cv2.INTER_NEAREST)

# Overlay mask on tissue image for visualization
plt.figure(figsize=(10, 10))
plt.imshow(image)
plt.imshow(mask_resized, alpha=0.5)
plt.show()

# Compute per-spot tissue composition
annotation_array = []
for i in range(len(df_positions)):
    spot = (round(df_positions.iloc[i, 4] * scale_factor), round(df_positions.iloc[i, 5] * scale_factor))
    top_left, bottom_right = get_bb_from_spot(spot, round(diam * scale_factor))
    tile_mask = mask_resized[top_left[1]:bottom_right[1], top_left[0]:bottom_right[0]]
    
    # Get unique colors in the mask tile with their frequencies
    colors, counts = np.unique(tile_mask.reshape(-1, 3), axis=0, return_counts=True)
    counts_fraction = counts / sum(counts)
    
    # Initialize tissue fractions dictionary
    color_counts = {key: 0 for key in color_map.keys()}
    
    # Assign fractions to corresponding tissue class
    for j in range(len(colors)):
        color_code = ''.join(str(x) for x in colors[j])
        if color_code in color_counts:
            color_counts[color_code] = counts_fraction[j]
    
    annotation_array.append(color_counts)

# Create annotated DataFrame
df_annotations = pd.DataFrame(annotation_array).rename(columns=color_map)
df_combined = pd.concat([df_positions, df_annotations], axis=1)
df_combined.to_csv(out_path, header=None, index=False)

# Compute total tissue composition
colors, counts = np.unique(mask_resized.reshape(-1, 3), axis=0, return_counts=True)
counts_fraction = counts / sum(counts)

color_counts_total = {key: 0 for key in color_map.keys()}
for j in range(len(colors)):
    color_code = ''.join(str(x) for x in colors[j])
    if color_code in color_counts_total:
        color_counts_total[color_code] = counts_fraction[j]

df_total = pd.DataFrame([color_counts_total]).rename(columns=color_map)

# Remove artefacts and normalize, and write the results
subtract = df_total.loc[0, ['Artefact', 'Hole', 'Out']].sum()
df_clean = df_total.drop(columns=['Artefact', 'Hole', 'Out']) / (1 - subtract)
df_clean = df_clean.transpose()
df_clean.columns = [sample]
df_clean.to_csv(out_adjusted_path, index=True)
