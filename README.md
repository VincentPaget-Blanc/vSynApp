# vSynApp

## Introduction

The Synaptosome Imaging Analysis App (vSynApp) is designed to facilitate the analysis of synaptosomal proteins imaged using immunofluorescence.

This app provides a suite of tools to analyze immunofluorescence microscopy images of synaptosomes. The primary goal is to analyse fluorescence on a per-structure basis, enabling researchers to assess signal thresholds, potential antibody interactions, co-localizations, and the presence of bleedthrough between channels. Bleedthrough refers to the phenomenon where the signal from one fluorescence channel is detected in another channel, which can lead to false-positive interpretations of co-localization.

The app incorporates a variety of statistical and computational methods to analyze the data. These methods help distinguish between specific staining and background noise or non-specific staining, ensuring accurate interpretation of the microscopy data.

### Citation

If you publish work using this software please use the following citation: "Paget-Blanc, V. (2025). vSynApp (0.094051). Zenodo. https://doi.org/10.5281/zenodo.17533244"

### Purpose of the App

The app is tailored to address specific research questions related to synaptosome imaging, including:
- Intensity thresholding
- Identifying and correcting for bleedthrough between immunofluorescence channels
- Differentiating specific signals from background noise or non-specific staining

By providing a structured and automated approach to these analyses, the app aims to enhance the reproducibility and reliability of synaptosome imaging studies.

### Data to use

Prepare data analysed using IJ-Toolset_SynaptosomesMacro by pooling Result.csv files containing x and y coordinates from the same immunofluorescence.

For a full history of changes between versions, see [CHANGELOG.md](CHANGELOG.md).

---

## Analysis Methods and Corresponding Plots

### 1. Pearson Correlation Analysis

**Purpose**:
Quantifies the degree of linear relationship between the intensities of two selected channels. This helps determine if changes in one protein's intensity correlate linearly with changes in another protein's intensity on a structure-by-structure basis.

**How it Works**:
1. Calculate the Pearson correlation coefficient (r) to measure the strength and direction of the linear relationship between the intensities of the two channels.
2. Determine the significance of the correlation using a p-value.

**Parameters**:
- **Channel 1 and 2 Intensities**: Background-corrected intensity values.
- **Correlation Coefficient (r)**: Indicates the strength and direction of the linear relationship, ranging from -1 to 1.
- **P-value**: Indicates the significance of the correlation, with a typical threshold of 0.05.

**Plot Description**:
Scatter plot where each point represents a single synaptosome, positioned by the intensity values of the two channels. Axes use a symmetric log scale to handle the wide dynamic range typical of immunofluorescence data.

**Example Image**:
![Scatter Plot Example](path/to/scatter_plot_example.png)

**Interpretation**:
- A high positive r value suggests a direct interaction between the two proteins or bleedthrough of one channel into the other.
- A significant p-value implies that the observed correlation is unlikely due to chance, which is crucial for validating potential interactions or bleedthrough effects.

---

### 2. Sub-cluster Correlation Analysis

**Purpose**:
Divides the data into clusters using K-means clustering and performs Pearson correlation analysis within each cluster. This can reveal localized patterns of protein interaction or proximity that might not be evident at the population level.

**How it Works**:
1. Use K-means clustering (k=3) to divide the data into clusters based on intensity values.
2. Perform Pearson correlation analysis within each cluster to identify localized correlation patterns.

**Parameters**:
- **Number of Clusters (k)**: Fixed at 3.
- **Intensity Values**: Background-corrected intensities from the selected channels.
- **Correlation Coefficient and P-value**: Computed per cluster to identify localized patterns of interaction or bleedthrough.

**Plot Description**:
Three scatter subplots, one per cluster, each annotated with the Pearson r and p-value for that subpopulation.

**Example Image**:
![Sub-cluster Scatter Plot Example](path/to/subcluster_scatter_example.png)

**Interpretation**:
- Different clusters may reveal distinct interaction patterns or levels of bleedthrough.
- This method can uncover nuanced relationships that overall analysis might overlook, providing insights into specific synaptosome subpopulations.

---

### 3. Background Comparison

**Purpose**:
Compares the background intensities between the two channels to assess if there is a significant difference in background signal levels. This is essential for identifying non-specific staining and assessing whether stainings are comparable.

**How it Works**:
Perform a two-sample t-test and a Mann-Whitney U test to compare background intensities between channels.

**Parameters**:
- **Background Intensities**: Background intensity values for both channels.
- **T-test Statistic and P-value**: Measures the difference in means and its significance.
- **U Statistic and P-value**: Non-parametric measure of distribution differences.

**Plot Description**:
A violin plot overlaid with individual data points (strip plot, subsampled to 500 points per channel for readability) showing the full distribution of background intensity values for each channel. This is more informative than a boxplot for the right-skewed distributions typical of immunofluorescence data. Statistical annotations (t-test and Mann-Whitney p-values) are displayed on the plot.

**Example Image**:
![Violin/Strip Plot Example](path/to/background_comparison_example.png)

**Interpretation**:
- A significant difference in background intensities may indicate differential non-specific staining or bleedthrough, which would need to be accounted for in subsequent analyses.

---

### 4. Vector Analysis (Rayleigh p-value)

**Purpose**:
Analyzes vectors between centroids of the two channels to assess the uniformity of vector directions. This can help determine if there is a systematic spatial relationship between the two proteins on synaptosomes.

**How it Works**:
1. Calculate vectors between centroids and convert them into angles and magnitudes.
2. Use the Rayleigh test to determine if directions are uniformly distributed or if there is a preferred direction.
3. Vectors are incrementally added in order of magnitude and the Rayleigh p-value is tracked cumulatively, with Bonferroni correction applied at each step.

The Rayleigh test p-value is computed using the exact chi-squared form `chi2.sf(2 * n * R_bar², df=2)`, which is statistically standard and numerically stable across all sample sizes.

**Parameters**:
- **Centroid Positions**: X and Y coordinates of centroids for both channels.
- **Alpha Value**: Significance level for the Rayleigh test.
- **Sorting Order**: Order in which vectors are added by magnitude (ascending or descending).
- **MRL (Mean Resultant Length)**: Measures the concentration of vector directions. The MRL threshold is fully user-configurable and applied consistently across all analysis pathways.

**Plot Description**:
Four-panel figure: polar histograms for all vectors, vectors below the magnitude threshold, and vectors above the magnitude threshold; plus a Rayleigh p-value progression plot (log scale) tracking significance as vectors are added in order of magnitude, with the threshold point annotated.

When a threshold produces zero vectors in a panel, the plot clearly labels "No vectors below/above threshold" rather than displaying a silent empty panel.

**Example Image**:
![Vector Plot Example](path/to/vector_plot_example.png)

**Interpretation**:
- A significant p-value combined with an MRL above the user-defined threshold suggests a preferred direction in the spatial relationship between proteins, which might indicate a specific interaction or consistent bleedthrough effect.

---

### 5. Vector Analysis (MRL Threshold)

**Purpose**:
Focuses on using the Mean Resultant Length (MRL) as a threshold to determine significant vector directions. This can help identify subsets of synaptosomes with significant spatial relationships between proteins.

**How it Works**:
1. Sort vectors by magnitude (ascending or descending, user-controlled).
2. Iteratively calculate the MRL for growing subsets of vectors.
3. The magnitude at which MRL first drops below the user-defined threshold is used as the cutoff, separating below-threshold and above-threshold populations.

**Parameters**:
- **Centroid Positions**: X and Y coordinates of centroids.
- **Sorting Order**: Determines if vectors are sorted by ascending or descending magnitude.
- **MRL Threshold**: User-defined threshold to determine significance.
- **Alpha Value**: Significance level for statistical tests.

**Plot Description**:
Four-panel figure: polar histograms for all vectors, below-threshold vectors, and above-threshold vectors; plus an MRL progression plot tracking how MRL evolves as vectors are added in order of magnitude, with the threshold crossing annotated.

**Example Image**:
![MRL Threshold Vector Plot Example](path/to/mrl_vector_plot_example.png)

**Interpretation**:
- Vectors below the determined threshold may exhibit significant directionality, suggesting specific protein interactions or consistent bleedthrough patterns in certain subsets of synaptosomes.

---

### 6. t-SNE and kNN Analysis

**Purpose**:
t-SNE (t-distributed Stochastic Neighbor Embedding) combined with K-Means clustering identifies natural subpopulations within the data. This method helps visualize high-dimensional intensity data in a lower-dimensional space and group synaptosomes with similar properties.

**How it Works**:
1. Use t-SNE (perplexity=30) to reduce the dimensionality of the data while preserving local structures.
2. Apply K-Means clustering to identify natural groupings within the reduced-dimensional space.
3. Project cluster labels back onto the original intensity scatter plot for interpretation.

t-SNE results are cached after the first computation. Subsequent calls with the same channels and `k` reuse the embedding rather than recomputing it.

**Parameters**:
- **Intensity Values**: Background-corrected intensities from the selected channels.
- **Number of Clusters (k)**: The number of clusters to identify (minimum 2).
- **Perplexity**: t-SNE parameter balancing local and global structure (default 30).

**Plot Description**:
- **t-SNE Clustering Plot**: Data points in the reduced-dimensional space, color-coded by cluster.
- **Clusters in Original Data Plot**: Original intensity scatter plot (symlog scale) with points color-coded by t-SNE-derived cluster labels.
- **Individual Clusters Plot**: Per-cluster scatter subplots with consistent axes for direct visual comparison across clusters.

**Example Image**:
![t-SNE Clustering Plot Example](path/to/tsne_clustering_example.png)

**Interpretation**:
- Clusters identified by K-Means can reveal distinct groups of synaptosomes with similar intensity profiles.
- The spatial distribution in the t-SNE plot reflects similarity between data points, with closer points being more alike.

---

### 7. t-SNE and kNN Thresholding

**Purpose**:
Uses the clusters identified by K-Means (applied to t-SNE embeddings) to set intensity thresholds. This method helps distinguish between different states or conditions of synaptosomes based on their intensity values.

**How it Works**:
1. Perform t-SNE and K-Means analysis to identify clusters (or reuse cached results).
2. Use the minimum intensity of the user-selected cluster as the threshold for each channel independently.
3. Apply these thresholds to segment the data and pass double-positive particles into a polar vector plot.

**Parameters**:
- **Cluster for X Threshold**: Cluster whose minimum intensity sets the Channel 1 threshold.
- **Cluster for Y Threshold**: Cluster whose minimum intensity sets the Channel 2 threshold.

**Plot Description**:
- Per-channel intensity histograms with threshold lines.
- Symlog-scaled scatter plot with both threshold lines overlaid.
- Polar vector histogram for synaptosomes above both channel thresholds simultaneously.

**Example Image**:
![t-SNE and kNN Thresholding Plot](path/to/tsne_thresholding_example.png)

**Interpretation**:
The determined thresholds segment particles based on both intensity and spatial relationships between channels, providing insights into protein interactions and potential bleedthrough effects.

---

### 8. Triangle Thresholding

**Purpose**:
Determines intensity thresholds for each channel based on the shape of the intensity histogram. This helps segment images to distinguish between specific and non-specific staining.

**How it Works**:
1. Construct a histogram of intensity values for each channel.
2. Identify the histogram peak and find the point where the perpendicular distance from the histogram to the line connecting the peak to the tail is maximal. This point is the threshold.

**Parameters**:
- **Intensity Values**: Background-corrected intensities for each channel.
- **Threshold Value**: Intensity threshold for each channel.

**Plot Description**:
Per-channel intensity histograms with threshold lines, a symlog-scaled scatter plot with threshold lines, and a polar vector plot for double-positive particles.

**Example Image**:
![Triangle Thresholding Histogram Example](path/to/triangle_threshold_histogram_example.png)

**Interpretation**:
- The threshold separates significant signal from background noise or bleedthrough, which is crucial for accurate quantification of protein proximity.

---

### 9. K-means Thresholding

**Purpose**:
Uses K-means clustering to determine optimal intensity thresholds based on natural groupings of intensity values. This is more flexible than histogram-based thresholding and can reveal multiple intensity clusters corresponding to different states of protein recruitment or varying levels of bleedthrough.

**How it Works**:
1. Combine background-corrected intensity values with raw background intensities for each channel.
2. Use the elbow method to determine the optimal number of clusters.
3. Perform K-means clustering and set thresholds at the midpoints between adjacent cluster centers.

**Parameters**:
- **Intensity and Background Intensity Values**: For each channel.
- **Optimal Number of Clusters**: Determined automatically using the elbow method.
- **Threshold Values**: Midpoints between adjacent cluster centers.

**Plot Description**:
Per-channel intensity histograms with threshold lines, a symlog-scaled scatter plot, and a polar vector plot for double-positive particles.

**Example Image**:
![K-means Thresholding Cluster Plot Example](path/to/kmeans_cluster_plot_example.png)

**Interpretation**:
- Different intensity clusters may represent distinct states of protein recruitment or different levels of bleedthrough, providing a more nuanced understanding of the experimental data.

---

### 10. Threshold Vector Plots

**Purpose**:
All thresholding methods (Triangle, K-means, GMM, t-SNE/kNN) generate a polar vector plot for double-positive particles — those above threshold in both channels simultaneously. This indicates whether synaptosomes passing the dual threshold exhibit a preferred directionality, potentially revealing a genuine spatial apposition between the two markers rather than random co-localization.

**Plot Description**:
Polar histogram of vector angles for double-positive particles. If directionality is significant (Rayleigh p < 0.05 and MRL above the user-defined threshold), the mean angle is annotated with a red line.

**Example Image**:
![Vector Plot for double-positive particles](path/to/threshold_vector_plot_example.png)

---

### 11. Vector Magnitude Intensity Threshold

**Purpose**:
Determines intensity thresholds for each channel based on the directional coherence of inter-centroid vectors. By integrating spatial and intensity information, this method identifies the intensity cutoff that best separates synaptosomes with coherent spatial apposition from those without.

**How it Works**:
1. Calculate vectors between centroids and determine their magnitudes and angles.
2. Scan candidate intensity thresholds across the data range. At each candidate, evaluate the Rayleigh p-value and MRL of the above-threshold population.
3. Select the threshold where the above-threshold group first achieves significance (p < alpha, MRL >= MRL threshold) and maximize a score combining p-value and MRL.
4. Apply Bonferroni correction for the number of thresholds tested. If no significant threshold is found, fall back to the mean, median, or minimum of the intensity distribution.

**Parameters**:
- **Intensity Values**: Background-corrected intensities for each channel.
- **Centroid Positions**: X and Y coordinates for both channels.
- **Alpha Value (for magnitude)**: Significance level for this analysis (configurable separately from the global alpha).
- **Threshold Method**: Fallback strategy if no significant threshold is found — `mean`, `median`, or `min`.
- **MRL Threshold (for magnitude)**: MRL threshold for this analysis (configurable separately from the global MRL threshold).

**Plot Description**:
Four-panel figure: per-channel intensity histograms with threshold lines and adjusted p-values annotated; vector magnitude distribution with magnitude threshold; symlog-scaled scatter plot with both channel thresholds overlaid.

**Example Image**:
![Vector Magnitude Intensity Threshold Plot Example](path/to/vector_magnitude_threshold_plot_example.png)

**Interpretation**:
- The determined thresholds segment particles based on both intensity and the spatial coherence of their inter-channel vectors, providing a principled separation of specifically apposed pairs from background co-incidence.

---

## Example Usage in the Context of Synaptosomes

1. **Load Data**
   - **Option A — Single file**: Click **Load Data** to open a CSV file directly via a file browser.
   - **Option B — Multi-file search**: Tick **Enable Multi-File Search**, browse to your experiment root folder, type a filename substring (e.g. `Bassoon`) into the "Filename contains" field, and click **Search**. The app will recursively scan all subfolders and list every matching CSV by filename in the dropdown. Select the file you want and click **Load Selected File**.
   - The CSV should contain columns for background-corrected intensity, background intensity, and centroid X/Y positions for both channels, corresponding to validated synaptosome structures of interest.

2. **Select Channels and Columns**
   - Select the appropriate channels from the dropdowns; columns should auto-fill with the corresponding data columns. Manual selection is available if auto-detection does not match your column naming convention.

3. **Choose Analysis Options**
   - Select the desired analysis methods to address specific research questions.

4. **Run Analysis**
   - Execute the selected methods. Numerical results are collected internally and can be exported to CSV.

5. **Generate or Export Plots**
   - Click **Generate Plots** to open plots in interactive windows, or **Export Plots** to save them directly to a folder in PNG, SVG, or PDF format.

6. **Export Results**
   - Click **Export Results to CSV** to save all numerical outputs (thresholds, p-values, MRL values, etc.) to a single CSV file for further analysis or reporting.

---

# Installation and Running Guide for miniConda and Spyder

This guide will help you install miniConda, set up a conda environment, install the necessary packages, and run the script using Spyder on Linux, Windows, and macOS.

---

## 1. Install miniConda

### Windows
1. Download the miniConda installer for Windows from the [official Conda website](https://docs.conda.io/en/latest/miniconda.html).
2. Run the installer and follow the on-screen instructions to complete the installation.
3. Open the Anaconda Prompt from the Start menu.

### macOS
1. Download the miniConda installer for macOS from the [official Conda website](https://docs.conda.io/en/latest/miniconda.html).
2. Open the Terminal and navigate to the directory where the installer was downloaded.
3. Run the installer with the following command:

```bash
bash Miniconda3-latest-MacOSX-x86_64.sh
```

4. Follow the on-screen instructions to complete the installation.
5. Restart the Terminal to start using conda.

### Linux
1. Download the miniConda installer for Linux from the [official Conda website](https://docs.conda.io/en/latest/miniconda.html).
2. Open the Terminal and navigate to the directory where the installer was downloaded.
3. Run the installer with the following command:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

4. Follow the on-screen instructions to complete the installation.
5. Restart the Terminal to start using conda.

---

## 2. Create a Conda Environment

Open a terminal (Linux/macOS) or Anaconda Prompt (Windows) and run the following commands to create a new conda environment with Python 3.9:

```bash
conda create --name immuno_env python=3.9
conda activate immuno_env
```

---

## 3. Install Necessary Packages

Activate the conda environment and install the required packages:

```bash
conda install -c conda-forge pandas numpy matplotlib seaborn scipy scikit-learn tk
```

---

## 4. Install Spyder

Install Spyder in the conda environment to provide an IDE for running the script:

```bash
conda install -c conda-forge spyder
```

---

## 5. Run the Script

1. Open Spyder by running the following command in the terminal (Linux/macOS) or Anaconda Prompt (Windows):

```bash
spyder
```

2. In Spyder, open the script file by going to `File` > `Open` and selecting `SynaptoThresholdingInteraction_Vector_V1.py`.
3. Run the script by pressing the green "Run" button or by pressing `F5`.

---

## Example Script to Test Installation

Here is a simple script to ensure everything is installed correctly. Save this as `test_script.py` and run it in Spyder:

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
import tkinter as tk

print("All packages imported successfully!")
```

If the script runs without errors, your environment is set up correctly.
