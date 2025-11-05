# vSynApp

## Introduction

The Synaptosome Imaging Analysis App (vSynApp) is designed to facilitate the analysis of synatposomal proteins imaged using immunofluorescence. 

This app provides a suite of tools to analyze immunofluorescence microscopy images of synaptosomes. The primary goal is to analyse fluorescence on a per-structure basis, enabling researchers to assess signal thresholds, potential antibodies interactions, co-localizations, and the presence of bleedthrough between channels. Bleedthrough refers to the phenomenon where the signal from one fluorescence channel is detected in another channel, which can lead to false-positive interpretations of co-localization.

The app incorporates a variety of statistical and computational methods to analyze the data. These methods help distinguish between specific staining and background noise or non-specific staining, ensuring accurate interpretation of the microscopy data.

### Citation

If you use this software please use the following citation  “Paget-Blanc, V. (2025). vSynApp (0.094051). Zenodo. https://doi.org/10.5281/zenodo.17533244 “

### Purpose of the App

The app is tailored to address specific research questions related to synaptosome imaging, including:
- Intensity thresholding
- Identifying and correcting for bleedthrough between immunofluorescence channels.
- Differentiating specific signals from background noise or non-specific staining.

By providing a structured and automated approach to these analyses, the app aims to enhance the reproducibility and reliability of synaptosome imaging studies.

### Data to use

Prepare data analysed using  IJ-Toolset_SynaptosomesMacro by pooling Result.csv files containing x and y coordinates from the same immunofluorescence.
   
### Analysis Methods and Corresponding Plots

Possible analysis methods and their corresponding plots.

#### 1. Pearson Correlation Analysis

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
These plots display the relationship between the intensities of two channels. Each point on the scatter plot represents a single synaptosome, with its position determined by the intensity values in the two channels.

**Example Image**:
![Scatter Plot Example](path/to/scatter_plot_example.png)

**Interpretation**:
- A high positive r value suggests a direct interaction between the two proteins or bleedthrough of one channel into the other.
  
- A significant p-value implies that the observed correlation is unlikely due to chance, which is crucial for validating potential interactions or bleedthrough effects.

#### 2. Sub-cluster Correlation Analysis

**Purpose**:
Divides the data into clusters using K-means clustering and performs Pearson correlation analysis within each cluster. This can reveal localized patterns of protein interaction or proximity that might not be evident at the population level.

**How it Works**:
1. Use K-means clustering to divide the data into clusters based on intensity values.
2. Perform Pearson correlation analysis within each cluster to identify localized correlation patterns.

**Parameters**:
- **Number of Clusters (k)**: Typically set to 3, but adjustable based on data characteristics.
- **Intensity Values**: Background-corrected intensities from the selected channels.
- **Correlation Coefficient and P-value**: For each cluster, to identify localized patterns of interaction or bleedthrough.

**Plot Description**:
Scatter plots are used here as well.

**Example Image**:
![Sub-cluster Scatter Plot Example](path/to/subcluster_scatter_example.png)

**Interpretation**:
- Different clusters may reveal distinct interaction patterns or levels of bleedthrough.
- This method can uncover nuanced relationships that overall analysis might overlook, providing insights into specific synaptosome subpopulations.

#### 3. Background Comparison

**Purpose**:
Compares the background intensities between the two channels to assess if there is a significant difference in background signal levels. This is essential for identifying non-specific staining, and assess whether stainings are comparable.

**How it Works**:
Perform a two-sample t-test and a Mann-Whitney U test to compare background intensities between channels.

**Parameters**:
- **Background Intensities**: Background intensity values for both channels.
- **T-test Statistic and P-value**: Measures the difference in means and its significance.
- **U Statistic and P-value**: Non-parametric measure of distribution differences.

**Plot Description**:
Histograms showing the distribution of background intensity values for each channel. The x-axis represents the intensity values, while the y-axis represents the frequency of occurrences of each intensity value.

**Example Image**:
![Histogram Example](path/to/histogram_example.png)

**Interpretation**:
- A significant difference in background intensities may indicate differential non-specific staining or bleedthrough, which would need to be accounted for in subsequent analyses.

#### 4. Vector Analysis (Rayleigh p-value)

**Purpose**:
Analyzes vectors between centroids of the two channels to assess the uniformity of vector directions. This can help determine if there is a systematic spatial relationship between the two proteins on synaptosomes.

**How it Works**:
1. Calculate vectors between centroids and convert them into angles and magnitudes.
2. Use the Rayleigh test to determine if directions are uniformly distributed or if there is a preferred direction.

**Parameters**:
- **Centroid Positions**: X and Y coordinates of centroids for both channels.
- **Alpha Value**: Significance level for the Rayleigh test.
- **MRL (Mean Resultant Length)**: Measures the concentration of vector directions.
- **Rayleigh p-value**: Indicates whether vector directions are uniformly distributed.

**Plot Description**:
Vector plots that display the vectors between the centroids of the two channels. Each vector represents the spatial relationship between the two proteins within a single synaptosome.

**Example Image**:
![Vector Plot Example](path/to/vector_plot_example.png)

**Interpretation**:
- A significant p-value suggests a preferred direction in the spatial relationship between proteins, which might indicate a specific interaction or consistent bleedthrough effect.

#### 5. Vector Analysis (MRL Threshold)

**Purpose**:
Focuses on using the Mean Resultant Length (MRL) as a threshold to determine significant vector directions. This can help identify subsets of synaptosomes with significant spatial relationships between proteins.

**How it Works**:
1. Sort vectors by magnitude and iteratively calculate the MRL for subsets of vectors to determine a threshold for significant directionality.

**Parameters**:
- **Centroid Positions**: X and Y coordinates of centroids.
- **Sorting Order**: Determines if vectors are sorted by ascending or descending magnitude.
- **MRL Threshold**: User-defined threshold to determine significance.
- **Alpha Value**: Significance level for statistical tests.

**Plot Description**:
Vector plots, potentially color-coded or filtered to show vectors below the determined MRL threshold, indicating significant directionality.

**Example Image**:
![MRL Threshold Vector Plot Example](path/to/mrl_vector_plot_example.png)

**Interpretation**:
- Vectors below the determined threshold may exhibit significant directionality, suggesting specific protein interactions or consistent bleedthrough patterns in certain subsets of synaptosomes.

#### 6. t-SNE and kNN Analysis

**Purpose**:

t-SNE (t-distributed Stochastic Neighbor Embedding) and k-NN (k-Nearest Neighbors) analysis is designed to identify natural clusters within the data. This method helps in visualizing high-dimensional data in a lower-dimensional space and identifying groups of synaptosomes with similar properties.

**How it Works**:
1. Use t-SNE to reduce the dimensionality of the data while preserving local structures.
2. Apply k-NN clustering to identify natural groupings within the reduced-dimensional data.
3. Determine the optimal number of clusters and analyze the characteristics of each cluster.

**Parameters**:
- **Intensity Values**: Background-corrected intensities from the selected channels.
- **Number of Clusters (k)**: The number of clusters to be identified by the k-NN algorithm.
- **Perplexity**: A parameter for t-SNE that balances local and global aspects of the data.


**Plot Description**:

- **t-SNE Clustering Plot**: This plot shows the data points in the reduced-dimensional space obtained using t-SNE. Each point represents a synaptosome, and the points are color-coded based on the cluster they belong to, as determined by the k-NN algorithm. The plot helps visualize the clusters and understand the relationships between data points in the reduced-dimensional space.
- **Clusters in Original Data Plot**: This plot shows the original data points, with each point representing a synaptosome. The points are color-coded based on the cluster they belong to, as determined by the k-NN algorithm. This plot helps visualize how the clusters identified in the t-SNE space correspond to the original data.
- **Individual Clusters Plot**: This plot shows the data points for each cluster separately. Each subplot represents a different cluster, and the points are plotted using their original intensity values. This plot helps understand the characteristics of each cluster in the original data space.

**Example Image**:

![t-SNE Clustering Plot Example](path/to/mrl_vector_plot_example.png)
![Clusters in Original Data Plot](path/to/mrl_vector_plot_example.png)
![Individual Clusters Plot](path/to/mrl_vector_plot_example.png)

**Interpretation**:
- Clusters identified by k-NN can reveal distinct groups of synaptosomes with similar properties.
- The spatial distribution of points in the t-SNE plot can indicate the similarity between data points, with closer points being more similar.

#### 6. t-SNE and kNN Thresholding

**Purpose**:

t-SNE and kNN thresholding uses the clusters identified by the k-NN algorithm to set intensity thresholds. This method helps in distinguishing between different states or conditions of synaptosomes based on their intensity values.

**How it Works**:
1. Perform t-SNE and k-NN analysis to identify clusters.
2. Use the identified clusters to set intensity thresholds for each cluster.
3. Apply these thresholds to segment the data and analyze the characteristics of each segment.

**Parameters**:
- **Cluster Number**: Cluster selected for thresholding.
- **Threshold Values**: Intensity thresholds determined for each cluster.


**Plot Description**:

- **t-SNE and kNN Thresholding histogram Plot**: Histogram plots for representating the threshold for each channel

- **t-SNE and kNN Thresholding Scatter Plot**: Scatter plots showing the thresholding applied to the clusters, with intensity thresholds indicated.
- **Vector Polar Plot**: Vector polar plots for data points above both channels thresholds.

**Example Image**:

![t-SNE and kNN Thresholding histogram Plot](path/to/mrl_vector_plot_example.png)
![t-SNE and kNN Thresholding Scatter Plot](path/to/mrl_vector_plot_example.png)
![Vector Polar Plot](path/to/mrl_vector_plot_example.png)

**Interpretation**:

The determined thresholds can help segment images based on both intensity and spatial relationships between channels, providing insights into protein interactions and potential bleedthrough effects.
#### 7. Triangle Thresholding

**Purpose**:
Determines intensity thresholds for each channel based on the histogram of intensity values. This helps segment images to distinguish between specific and non-specific staining or to identify bleedthrough between channels.

**How it Works**:
1. Construct a histogram of intensity values for each channel.
2. Identify the peak of the histogram and determine the threshold where the distance from the histogram peak to the line connecting the peak to the highest intensity is maximal.

**Parameters**:
- **Intensity Values**: Background-corrected intensities for each channel.
- **Threshold Value**: Intensity threshold for each channel.

**Plot Description**:
Histograms showing the distribution of intensity values for each channel. The threshold is often identified as a peak or valley in the histogram.

**Example Image**:
![Triangle Thresholding Histogram Example](path/to/triangle_threshold_histogram_example.png)

**Interpretation**:
- The threshold value can be used to binarize images, separating significant signal from background noise or bleedthrough, which is crucial for accurate quantification of protein proximity.

#### 8. K-means Thresholding

**Purpose**:
Uses K-means clustering to determine optimal intensity thresholds based on natural groupings of intensity values. This is more flexible than histogram-based thresholding and can reveal multiple intensity clusters, which might correspond to different states of protein recruitment or varying levels of bleedthrough.

**How it Works**:
1. Combine intensity values with background intensities for clustering.
2. Use the elbow method to determine the optimal number of clusters.
3. Perform K-means clustering to find natural groupings and determine thresholds at midpoints between cluster centers.

**Parameters**:
- **Intensity and Background Intensity Values**: For each channel.
- **Optimal Number of Clusters**: Determined using the elbow method.
- **Threshold Values**: Calculated based on midpoints between cluster centers.

**Plot Description**:
Histograms or cluster plots showing the natural groupings of intensity values, with thresholds indicated between cluster centers.

**Example Image**:
![K-means Thresholding Cluster Plot Example](path/to/kmeans_cluster_plot_example.png)


**Interpretation**:
- Different intensity clusters may represent distinct states of protein recruitment or different levels of bleedthrough, providing a more nuanced understanding of the experimental data.

#### 9. Threshold Vector plots
**Purpose**:
Both thresholding methods will also generate vector plots for double-positive particles, which can indicate whether synaptosomes that fall above both thresholds exhibit a preferred directionality, potentially indicating bleedthrough.

**Plot Description**:
Vector plots, indicating significant directionality

**Example Image**:
![Vector Plot for double-positive particles](path/to/kmeans_cluster_plot_example.png)
#### 10. Vector Magnitude Intensity Threshold

**Purpose**:
Determines intensity thresholds for channels based on vector magnitudes and angles, focusing on separating vectors with significant directionality. This integrates spatial and intensity information to identify meaningful thresholds, which can help distinguish specific from non-specific staining or detect consistent bleedthrough patterns.

**How it Works**:
1. Calculate vectors between centroids, determine their magnitudes and angles, and use statistical tests to determine intensity thresholds based on vector directionality and magnitude.

**Parameters**:
- **Intensity Values**: Background-corrected intensities for each channel.
- **Centroid Positions**: X and Y coordinates for both channels.
- **Alpha Value**: Significance level for statistical tests.
- **Threshold Method**: Method for determining the threshold ('mean', 'median', or 'min').
- **MRL Threshold**: User-defined threshold for MRL to determine significance.

**Plot Description**:
Vector plots or scatter plots that integrate spatial and intensity information, showing significant thresholds that help distinguish specific from non-specific staining.

**Example Image**:
![Vector Magnitude Intensity Threshold Plot Example](path/to/vector_magnitude_threshold_plot_example.png)

**Interpretation**:
- The determined thresholds can help segment images based on both intensity and spatial relationships between channels, providing insights into protein interactions and potential bleedthrough effects.

### Example Usage in the Context of Synaptosomes

1. **Load Data**
   - Load a CSV file containing columns for intensity, background intensity, and centroid positions for both channels. Ensure the data corresponds to validated synaptosome structures of interest.

2. **Select Channels and Columns**
   - Select the appropriate channels, columns should be auto-filled with the corresponding background-corrected intensity, background intensity, and centroid positions for further analysis if not you can manually do so.

3. **Choose Analysis Options**
   - Select the desired analysis methods to address specific questions.

4. **Run Analysis**
   - Execute the selected methods and visualize the results using plots such as scatter plots for correlation analysis, histograms for thresholding methods, and vector plots for vector analysis.
5. **Generate or export plots**
   - Generate plots for visualization or export them to keep a track of the quality of your immunofluorescence.

6. **Export Results**
   - Save plots and analysis results in desired formats (e.g., PNG, CSV) for further analysis, reporting, or additional quality control steps.


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

   (Code block will be provided later)

2. In Spyder, open your script file by going to `File` > `Open` and selecting your script.
3. Run the script by pressing the green "Run" button or by pressing `F5`.

---

## Example Script to Test Installation

Here's a simple script to ensure everything is installed correctly. Save this as `test_script.py` and run it in Spyder:

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
