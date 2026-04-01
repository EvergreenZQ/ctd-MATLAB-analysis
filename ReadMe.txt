Project Overview
This project contains two MATLAB scripts for processing and analyzing STORM microscopy data:
Crop Script
Extracts regions of interest (ROIs) from .bin files using manually annotated masks.
Analysis Script
Performs quantitative analysis on cropped regions, including mask size, colocalization, distance, and sub-dot detection.

1. Crop Script
Purpose: The crop script extracts specific regions (e.g., dendritic regions) from STORM .bin files using manually annotated PNG masks, and generates new .bin files for each segmented region.

Workflow
Step 1: Manual Annotation (IMPORTANT)
Open images in PowerPoint or any drawing tool
Select regions of interest
Paint selected regions in pure blue (RGB: 0, 0, 255)

Save images as:

*_crop.png
Place them in the same folder as the corresponding .bin files

⚠️ Strict requirement:

Each .bin file must correspond to exactly one *_crop.png
Files must be in matching order

Step 2: Run Crop Script

The script will:

Load .bin localization data
Read the corresponding _crop.png
Convert blue regions into binary masks
Label connected regions
Assign molecules to each region
Rotate each region to align orientation
Save each region as a new .bin file

2. Analysis Script
Purpose: The analysis script processes cropped .bin files to quantify:

Mask size (Channel 1 & 2)
Pearson Correlation Coefficient (PCC)
Channel ratio
Distance between clusters
Overlap area
Number of sub-dots

Workflow
Step 1: Load Data
Extract localization coordinates for:
Channel 1
Channel 2
Step 2: Preprocessing
Remove noise using DBSCAN clustering
Keep only clustered points
Step 3: Rendering
Convert localization data into images using Gaussian rendering
Normalize and resize channels
Step 4: Mask Generation
Background subtraction
Otsu thresholding
Remove small objects
Extract largest connected component
Step 5: Quantification

The script computes:

Mask size (Ch1 & Ch2)
Centroid distance between channels
Overlap area
PCC (correlation between channels)
Step 6: Manual Mask Selection (IMPORTANT)

A GUI will appear:

User manually selects valid region using polygon tool (roipoly)
Final mask is refined based on user input

Step 7: Sub-dot Detection
Difference of Gaussian filtering
Thresholding
Connected component analysis

Outputs:

Number of sub-dots in each channel
Step 8: Visualization Output

A 4-panel image is generated:

Original rendered image
Largest connected regions
Final mask
Sub-dot detection