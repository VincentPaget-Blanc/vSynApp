#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 20:10:57 2025

@author: vincentpb
"""

import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.stats import pearsonr, ttest_ind, mannwhitneyu, circmean, circstd, chi2
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.manifold import TSNE
import seaborn as sns

class ThresholdingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Immunofluorescence Thresholding App")
        self.root.geometry("800x600")
        self.root.minsize(600, 400)
        self.data = None
        self.data_file_path = None
        self.output_dir = "output"
        self.columns = []
        self.analysis_results = []
        self.export_format = tk.StringVar(value="png")
        self.channel1_var = tk.StringVar()
        self.channel2_var = tk.StringVar()
        self.channel_names = []
        self.threshold_method = tk.StringVar(value='mean')
        self.window_size = tk.IntVar(value=10)
        self.mrl_threshold = tk.DoubleVar(value=0.3)
        self.alpha_value = tk.DoubleVar(value=0.05)
        self.vector_magnitude_alpha = tk.DoubleVar(value=0.05)
        self.vector_magnitude_mrl_threshold = tk.DoubleVar(value=0.3)
        self.sorting_order = tk.StringVar(value='ascending')
        self.gmm_thresh_var = tk.BooleanVar()
        self.tsne_knn_var = tk.BooleanVar()  # New variable for t-SNE and kNN analysis
        self.tsne_knn_thresh_var = tk.BooleanVar()
        self.k_var = tk.IntVar(value=3)
        self.cluster_x_var = tk.StringVar()
        self.cluster_y_var = tk.StringVar()
        self.plots = {}
        self.X_tsne = None
        self.labels_tsne = None
        self.k = None
        self.create_widgets()

    def create_widgets(self):
        # Create a main frame to hold the canvas and scrollbar
        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(fill=tk.BOTH, expand=True)
    
        # Create a canvas object and a vertical scrollbar for scrolling it
        self.canvas = tk.Canvas(self.main_frame)
        self.scrollbar = tk.Scrollbar(self.main_frame, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = tk.Frame(self.canvas)
    
        # Create a window item inside the canvas that contains the scrollable frame
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(
                scrollregion=self.canvas.bbox("all")
            )
        )
    
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
    
        # Pack the canvas and scrollbar
        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")
    
        # Load Data Button
        self.load_button = tk.Button(self.scrollable_frame, text="Load Data", command=self.load_data)
        self.load_button.pack(pady=10)
    
        # Channel Selection Frame
        self.channel_frame = tk.LabelFrame(self.scrollable_frame, text="Channel Selection", padx=10, pady=10)
        self.channel_frame.pack(pady=10)
        tk.Label(self.channel_frame, text="Select Channel 1:").grid(row=0, column=0, padx=5, pady=5, sticky='w')
        self.channel1_dropdown = ttk.Combobox(self.channel_frame, textvariable=self.channel1_var, state="disabled")
        self.channel1_dropdown.grid(row=0, column=1, padx=5, pady=5)
        self.channel1_dropdown.bind("<<ComboboxSelected>>", self.update_column_dropdowns)
        tk.Label(self.channel_frame, text="Select Channel 2:").grid(row=1, column=0, padx=5, pady=5, sticky='w')
        self.channel2_dropdown = ttk.Combobox(self.channel_frame, textvariable=self.channel2_var, state="disabled")
        self.channel2_dropdown.grid(row=1, column=1, padx=5, pady=5)
        self.channel2_dropdown.bind("<<ComboboxSelected>>", self.update_column_dropdowns)
    
        # Column Selection Frame
        self.column_frame = tk.LabelFrame(self.scrollable_frame, text="Column Selection", padx=10, pady=10)
        self.column_frame.pack(pady=10)
        self.column_labels = [
            "Background Corrected Intensity:",
            "Background Intensity:",
            "Centroid X:",
            "Centroid Y:"
        ]
        self.column_vars = []
        self.column_dropdowns = []
        for i, label in enumerate(self.column_labels):
            tk.Label(self.column_frame, text=f"{label} Channel 1:").grid(row=i, column=0, padx=5, pady=5, sticky='w')
            var = tk.StringVar()
            dropdown = ttk.Combobox(self.column_frame, textvariable=var, state="disabled")
            dropdown.grid(row=i, column=1, padx=5, pady=5)
            self.column_vars.append(var)
            self.column_dropdowns.append(dropdown)
        for i, label in enumerate(self.column_labels):
            tk.Label(self.column_frame, text=f"{label} Channel 2:").grid(row=i, column=2, padx=5, pady=5, sticky='w')
            var = tk.StringVar()
            dropdown = ttk.Combobox(self.column_frame, textvariable=var, state="disabled")
            dropdown.grid(row=i, column=3, padx=5, pady=5)
            self.column_vars.append(var)
            self.column_dropdowns.append(dropdown)
    
        # Analysis Options Frame
        self.analysis_frame = tk.LabelFrame(self.scrollable_frame, text="Analysis Options", padx=10, pady=10)
        self.analysis_frame.pack(pady=10)
    
        # Correlation options
        self.corr_var = tk.BooleanVar()
        self.corr_check = tk.Checkbutton(self.analysis_frame, text="Pearson Correlation Analysis", variable=self.corr_var)
        self.corr_check.pack(anchor='w')
    
        self.subcluster_corr_var = tk.BooleanVar()
        self.subcluster_corr_check = tk.Checkbutton(self.analysis_frame, text="Sub-cluster Correlation Analysis", variable=self.subcluster_corr_var)
        self.subcluster_corr_check.pack(anchor='w')
    
        self.bg_comp_var = tk.BooleanVar()
        self.bg_comp_check = tk.Checkbutton(self.analysis_frame, text="Background Comparison", variable=self.bg_comp_var)
        self.bg_comp_check.pack(anchor='w')
    
        # Vector Analysis (Rayleigh p-value) section
        self.vector_rayleigh_frame = tk.LabelFrame(self.analysis_frame, text="Vector Analysis (Rayleigh p-value)", padx=5, pady=5)
        self.vector_rayleigh_frame.pack(pady=5, fill='x')
    
        self.vector_analysis_var = tk.BooleanVar()
        self.vector_analysis_check = tk.Checkbutton(
            self.vector_rayleigh_frame,
            text="Enable Vector Analysis (Rayleigh p-value)",
            variable=self.vector_analysis_var
        )
        self.vector_analysis_check.grid(row=0, column=0, columnspan=3, sticky='w', padx=5, pady=2)
    
        tk.Label(self.vector_rayleigh_frame, text="Sorting Order:").grid(row=2, column=0, sticky='w', padx=5, pady=2)
        sorting_dropdown_rayleigh = ttk.Combobox(
            self.vector_rayleigh_frame,
            textvariable=self.sorting_order,
            values=['ascending', 'descending'],
            state="readonly",
            width=10
        )
        sorting_dropdown_rayleigh.grid(row=2, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.vector_rayleigh_frame, text="(Order for magnitude sorting)").grid(row=2, column=2, sticky='w', padx=5, pady=2)
        tk.Label(self.vector_rayleigh_frame, text="Alpha Value:").grid(row=1, column=0, sticky='w', padx=5, pady=2)
        self.alpha_entry = tk.Entry(self.vector_rayleigh_frame, textvariable=self.alpha_value, width=5)
        self.alpha_entry.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.vector_rayleigh_frame, text="(Significance level for Rayleigh test, typical: 0.05)").grid(row=1, column=2, sticky='w', padx=5, pady=2)
    
        # Vector Analysis (MRL) section
        self.vector_mrl_frame = tk.LabelFrame(self.analysis_frame, text="Vector Analysis (MRL)", padx=5, pady=5)
        self.vector_mrl_frame.pack(pady=5, fill='x')
    
        self.vector_analysis_mrl_var = tk.BooleanVar()
        self.vector_analysis_mrl_check = tk.Checkbutton(
            self.vector_mrl_frame,
            text="Enable Vector Analysis (MRL Threshold)",
            variable=self.vector_analysis_mrl_var
        )
        self.vector_analysis_mrl_check.grid(row=0, column=0, columnspan=3, sticky='w', padx=5, pady=2)
        tk.Label(self.vector_mrl_frame, text="Sorting Order:").grid(row=2, column=0, sticky='w', padx=5, pady=2)
        sorting_dropdown_mrl = ttk.Combobox(
            self.vector_mrl_frame,
            textvariable=self.sorting_order,  # same variable as above
            values=['ascending', 'descending'],
            state="readonly",
            width=10
        )
        sorting_dropdown_mrl.grid(row=2, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.vector_mrl_frame, text="(Order for magnitude sorting)").grid(row=2, column=2, sticky='w', padx=5, pady=2)
        tk.Label(self.vector_mrl_frame, text="MRL Threshold:").grid(row=1, column=0, sticky='w', padx=5, pady=2)
        self.mrl_threshold_entry = tk.Entry(self.vector_mrl_frame, textvariable=self.mrl_threshold, width=5)
        self.mrl_threshold_entry.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.vector_mrl_frame, text="(Typical range: 0.3-0.9, default 0.7)").grid(row=1, column=2, sticky='w', padx=5, pady=2)
    
        # t-SNE and k-NN Analysis Frame (now inside Analysis Options)
        self.tsne_knn_frame = tk.LabelFrame(self.analysis_frame, text="t-SNE and k-NN Analysis", padx=10, pady=5)
        self.tsne_knn_frame.pack(pady=5, fill='x')
        
        self.tsne_knn_analysis_check = tk.Checkbutton(
            self.tsne_knn_frame,
            text="Enable t-SNE and k-NN Analysis",
            variable=self.tsne_knn_var,
            command=self.toggle_tsne_knn_analysis_widgets
        )
        self.tsne_knn_analysis_check.grid(row=0, column=0, columnspan=3, sticky='w', padx=5, pady=2)
        
        # Number of clusters for k-NN
        tk.Label(self.tsne_knn_frame, text="Number of Clusters (k):").grid(row=1, column=0, padx=5, pady=2, sticky='w')
        self.k_entry = tk.Entry(self.tsne_knn_frame, textvariable=self.k_var, width=5)
        self.k_entry.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.tsne_knn_frame, text="(Typical range: 2-10)").grid(row=1, column=2, sticky='w', padx=5, pady=2)

        # Thresholding Frame
        self.thresholding_frame = tk.LabelFrame(self.scrollable_frame, text="Thresholding", padx=10, pady=10)
        self.thresholding_frame.pack(pady=10)


        # t-SNE and kNN Thresholding Options sub-frame
        self.tsne_knn_thresholding_frame = tk.LabelFrame(self.thresholding_frame, text="t-SNE and k-NN Thresholding Options", padx=5, pady=5)
        self.tsne_knn_thresholding_frame.pack(pady=5, fill='x')
        
        self.tsne_knn_thresh_check = tk.Checkbutton(
            self.tsne_knn_thresholding_frame,
            text="t-SNE and k-NN Thresholding",
            variable=self.tsne_knn_thresh_var,
            command=self.toggle_tsne_knn_thresholding_widgets
        )
        self.tsne_knn_thresh_check.grid(row=0, column=0, columnspan=3, sticky='w', padx=5, pady=2)
        
        # Dropdown for selecting clusters for thresholding
        tk.Label(self.tsne_knn_thresholding_frame, text="Cluster for X Threshold:").grid(row=1, column=0, padx=5, pady=2, sticky='w')
        self.cluster_x_dropdown = ttk.Combobox(self.tsne_knn_thresholding_frame, textvariable=self.cluster_x_var, state="disabled")
        self.cluster_x_dropdown.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        
        tk.Label(self.tsne_knn_thresholding_frame, text="Cluster for Y Threshold:").grid(row=2, column=0, padx=5, pady=2, sticky='w')
        self.cluster_y_dropdown = ttk.Combobox(self.tsne_knn_thresholding_frame, textvariable=self.cluster_y_var, state="disabled")
        self.cluster_y_dropdown.grid(row=2, column=1, padx=5, pady=2, sticky='w')
    
        # Triangle thresholding checkbox
        self.triangle_thresh_var = tk.BooleanVar()
        self.triangle_thresh_check = tk.Checkbutton(self.thresholding_frame, text="Triangle Thresholding", variable=self.triangle_thresh_var)
        self.triangle_thresh_check.pack(anchor='w')
    
        # K-means thresholding checkbox
        self.kmeans_thresh_var = tk.BooleanVar()
        self.kmeans_thresh_check = tk.Checkbutton(self.thresholding_frame, text="K-means Thresholding", variable=self.kmeans_thresh_var)
        self.kmeans_thresh_check.pack(anchor='w')
    
        # GMM Thresholding
        self.gmm_thresh_check = tk.Checkbutton(self.thresholding_frame, text="GMM Thresholding", variable=self.gmm_thresh_var)
        self.gmm_thresh_check.pack(anchor='w')
    
        # Vector threshold options frame
        self.vector_threshold_frame = tk.LabelFrame(self.thresholding_frame, text="Vector Threshold Options", padx=5, pady=5)
        self.vector_threshold_frame.pack(pady=5, fill='x')
    
        # Configure grid columns to expand properly
        self.vector_threshold_frame.grid_columnconfigure(2, weight=1)
    
        tk.Label(self.vector_threshold_frame, text="Threshold Method:").grid(row=0, column=0, sticky='w', padx=5, pady=2)
        self.threshold_method_dropdown = ttk.Combobox(
            self.vector_threshold_frame,
            textvariable=self.threshold_method,
            values=['min', 'mean', 'median'],
            state="readonly",
            width=10
        )
        self.threshold_method_dropdown.grid(row=0, column=1, padx=5, pady=2, sticky='w')
    
        # Vector Magnitude Intensity Threshold checkbox
        self.vector_magnitude_intensity_var = tk.BooleanVar()
        self.vector_magnitude_intensity_check = tk.Checkbutton(
            self.vector_threshold_frame,
            text="Vector Magnitude Intensity Threshold",
            variable=self.vector_magnitude_intensity_var
        )
        self.vector_magnitude_intensity_check.grid(row=1, column=0, columnspan=3, sticky='w', padx=5, pady=2)
    
        # Add new controls for Vector Magnitude Intensity Threshold parameters
        tk.Label(self.vector_threshold_frame, text="Alpha Value (for magnitude):").grid(row=2, column=0, sticky='w', padx=5, pady=2)
        self.vector_magnitude_alpha_entry = tk.Entry(self.vector_threshold_frame, textvariable=self.vector_magnitude_alpha, width=5)
        self.vector_magnitude_alpha_entry.grid(row=2, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.vector_threshold_frame, text="(Significance level for magnitude analysis)").grid(row=2, column=2, sticky='w', padx=5, pady=2)
    
        tk.Label(self.vector_threshold_frame, text="MRL Threshold (for magnitude):").grid(row=3, column=0, sticky='w', padx=5, pady=2)
        self.vector_magnitude_mrl_threshold_entry = tk.Entry(self.vector_threshold_frame, textvariable=self.vector_magnitude_mrl_threshold, width=5)
        self.vector_magnitude_mrl_threshold_entry.grid(row=3, column=1, padx=5, pady=2, sticky='w')
        tk.Label(self.vector_threshold_frame, text="(Typical range: 0.3-0.9)").grid(row=3, column=2, sticky='w', padx=5, pady=2)
    
        # Button Frame
        self.button_frame = tk.Frame(self.scrollable_frame)
        self.button_frame.pack(pady=10)
        self.run_button = tk.Button(self.button_frame, text="Run Analysis", command=self.run_analysis)
        self.run_button.grid(row=0, column=0, padx=5, pady=5)
        self.plot_button = tk.Button(self.button_frame, text="Generate Plots", command=self.generate_plots)
        self.plot_button.grid(row=0, column=1, padx=5, pady=5)
        self.export_plot_button = tk.Button(self.button_frame, text="Export Plots", command=self.export_plots)
        self.export_plot_button.grid(row=0, column=2, padx=5, pady=5)
        self.export_button = tk.Button(self.button_frame, text="Export Results to CSV", command=self.export_results_to_csv)
        self.export_button.grid(row=0, column=3, padx=5, pady=5)
    
        # Export Format Frame
        self.export_format_frame = tk.LabelFrame(self.scrollable_frame, text="Export Format", padx=10, pady=10)
        self.export_format_frame.pack(pady=5)
        self.export_format_dropdown = ttk.Combobox(self.export_format_frame, textvariable=self.export_format, values=["png", "svg", "pdf"], state="readonly")
        self.export_format_dropdown.pack()
        
    def perform_tsne_knn_analysis(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
        X = np.column_stack((ch1_intensity, ch2_intensity))
        # Apply t-SNE
        tsne = TSNE(n_components=2, random_state=0, perplexity=30)
        X_tsne = tsne.fit_transform(X)
        # Apply k-NN
        k = self.k_var.get()
        kmeans_tsne = KMeans(n_clusters=k, random_state=0)
        labels_tsne = kmeans_tsne.fit_predict(X_tsne)
        # Store the results
        self.X_tsne = X_tsne
        self.labels_tsne = labels_tsne
        self.k = k
        # Update dropdowns for cluster selection
        self.cluster_x_dropdown['values'] = [f'Cluster {i}' for i in range(k)]
        self.cluster_x_dropdown.set(f'Cluster 0')
        self.cluster_y_dropdown['values'] = [f'Cluster {i}' for i in range(k)]
        self.cluster_y_dropdown.set(f'Cluster 2')
        return X_tsne, labels_tsne, k
    
    def _create_tsne_cluster_plot(self, X_tsne, labels_tsne, k):
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(X_tsne[:, 0], X_tsne[:, 1], c=labels_tsne, cmap='viridis', alpha=0.5)
    
        # Create a legend with a single dot for each cluster
        legend_elements = []
        for i in range(k):
            color = plt.cm.viridis(i / (k - 1))  # Use the same colormap as the scatter plot
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=f'Cluster {i}'))
    
        ax.legend(handles=legend_elements, title='Clusters')
        ax.set_title('t-SNE Clustering')
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        return fig
    
    def toggle_tsne_knn_analysis_widgets(self):
        if self.tsne_knn_var.get():
            if self.data is not None and hasattr(self, 'bkgd_corr_ch1') and hasattr(self, 'bkgd_corr_ch2'):
                self.perform_tsne_knn_analysis(self.bkgd_corr_ch1, self.bkgd_corr_ch2)
            # Remove the warning message when ticking the checkbox
            # The warning will only appear when running the analysis without data

    def toggle_tsne_knn_thresholding_widgets(self):
        state = "readonly" if self.tsne_knn_thresh_var.get() else "disabled"
        self.cluster_x_dropdown.config(state=state)
        self.cluster_y_dropdown.config(state=state)

    def plot_tsne_clusters(self, X_tsne, labels_tsne, k):
        fig = self._create_tsne_cluster_plot(X_tsne, labels_tsne, k)
        self.display_plot_in_new_window(fig, "t-SNE Clustering")
        plt.close(fig)
    
    def export_tsne_clusters_plot(self, directory, X_tsne, labels_tsne, k, bkgd_corr_ch1, bkgd_corr_ch2):
        fig = self._create_tsne_cluster_plot(X_tsne, labels_tsne, k)
        filename = os.path.join(directory, f"t-SNE_Clusters_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
    
    def _create_original_cluster_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k):
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(ch1_intensity, ch2_intensity, c=labels_tsne, cmap='viridis', alpha=0.5)
    
        # Create a legend with a single dot for each cluster
        legend_elements = []
        for i in range(k):
            color = plt.cm.viridis(i / (k - 1))  # Use the same colormap as the scatter plot
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=f'Cluster {i}'))
    
        ax.legend(handles=legend_elements, title='Clusters')
        ax.set_xscale('symlog')
        ax.set_yscale('symlog')
        ax.set_title('Clusters in Original Data')
        ax.set_xlabel(bkgd_corr_ch1)
        ax.set_ylabel(bkgd_corr_ch2)
        return fig
    
    def plot_original_clusters(self, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k):
        fig = self._create_original_cluster_plot(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
        self.display_plot_in_new_window(fig, "Clusters in Original Data")
        plt.close(fig)
    
    def export_original_clusters_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k):
        fig = self._create_original_cluster_plot(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
        filename = os.path.join(directory, f"Original_Clusters_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
    
    def toggle_tsne_knn_analysis_widgets(self):
        if self.tsne_knn_var.get():
            if self.data is not None and hasattr(self, 'bkgd_corr_ch1') and hasattr(self, 'bkgd_corr_ch2'):
                self.perform_tsne_knn_analysis(self.bkgd_corr_ch1, self.bkgd_corr_ch2)
            # Remove the warning message when ticking the checkbox
            # The warning will only appear when running the analysis without data

    def _create_individual_clusters_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k):
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
    
        # Determine the global min and max for x and y axes
        x_min, x_max = np.min(ch1_intensity), np.max(ch1_intensity)
        y_min, y_max = np.min(ch2_intensity), np.max(ch2_intensity)
    
        # Add a small buffer to the min and max values
        x_buffer = (x_max - x_min) * 0.05
        y_buffer = (y_max - y_min) * 0.05
    
        fig, axes = plt.subplots(1, k, figsize=(5*k, 5))
    
        # Use the same colormap as in the "Clusters in Original Data" plot
        cmap = plt.cm.viridis
    
        for i in range(k):
            ax = axes[i]
            cluster_mask = labels_tsne == i
            color = cmap(i / (k - 1))  # Use the same color as in the "Clusters in Original Data" plot
            ax.scatter(
                ch1_intensity[cluster_mask],
                ch2_intensity[cluster_mask],
                color=color,
                alpha=0.5
            )
            ax.set_xscale('symlog')
            ax.set_yscale('symlog')
    
            # Set the same x and y limits for all subplots
            ax.set_xlim(x_min - x_buffer, x_max + x_buffer)
            ax.set_ylim(y_min - y_buffer, y_max + y_buffer)
    
            ax.set_title(f'Cluster {i}')
            ax.set_xlabel(bkgd_corr_ch1)
            ax.set_ylabel(bkgd_corr_ch2)
    
        plt.tight_layout()
        return fig

    def plot_individual_clusters(self, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k):
        fig = self._create_individual_clusters_plot(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
        self.display_plot_in_new_window(fig, "Individual Clusters")
        plt.close(fig)
    
    def export_individual_clusters_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k):
        fig = self._create_individual_clusters_plot(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
        filename = os.path.join(directory, f"Individual_Clusters_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
        
    def extract_thresholds(self, ch1_intensity, ch2_intensity, labels_tsne):
        cluster_x = self.cluster_x_var.get()
        cluster_y = self.cluster_y_var.get()
        if not cluster_x or not cluster_y:
            return None, None
        cluster_x = int(cluster_x.split()[-1])
        cluster_y = int(cluster_y.split()[-1])
    
        x_threshold = np.min(ch1_intensity[labels_tsne == cluster_x])
        y_threshold = np.min(ch2_intensity[labels_tsne == cluster_y])
    
        return x_threshold, y_threshold

    def _get_tsne_knn_thresholding_data(self, bkgd_corr_ch1, bkgd_corr_ch2):
        # Use stored results if available
        if self.X_tsne is None or self.labels_tsne is None or self.k is None:
            messagebox.showwarning("Warning", "Please run t-SNE and kNN analysis first.")
            return None, None, None, None, None, None, None
    
        # Extract thresholds based on selected clusters
        x_threshold, y_threshold = self.extract_thresholds(
            self.data[bkgd_corr_ch1].values,
            self.data[bkgd_corr_ch2].values,
            self.labels_tsne
        )
    
        # Generate plots using these thresholds
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(
            self.data[bkgd_corr_ch1].values,
            self.data[bkgd_corr_ch2].values,
            x_threshold,
            y_threshold
        )
    
        # Get column variables
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
    
        return fig, x_threshold, y_threshold, cx_ch1, cy_ch1, cx_ch2, cy_ch2
    
    def generate_tsne_knn_thresholding_plots(self, bkgd_corr_ch1, bkgd_corr_ch2):
        fig, x_threshold, y_threshold, cx_ch1, cy_ch1, cx_ch2, cy_ch2 = self._get_tsne_knn_thresholding_data(bkgd_corr_ch1, bkgd_corr_ch2)
        self.display_plot_in_new_window(fig, "t-SNE and kNN Thresholding")
        if x_threshold is not None and y_threshold is not None:
            self.generate_single_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, x_threshold, y_threshold)
    
    def export_tsne_knn_thresholding_plots(self, directory, bkgd_corr_ch1, bkgd_corr_ch2):
        fig, x_threshold, y_threshold, cx_ch1, cy_ch1, cx_ch2, cy_ch2 = self._get_tsne_knn_thresholding_data(bkgd_corr_ch1, bkgd_corr_ch2)
        filename = os.path.join(directory, f"tSNE_kNN_Thresholding_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
        if x_threshold is not None and y_threshold is not None:
            self.export_single_polar_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, x_threshold, y_threshold, "tSNE_kNN_Thresholding")

    def add_plot_text(self, ax, text, position=(0.05, 0.05), fontsize=10, bbox_style=None):
        if bbox_style is None:
            bbox_style = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.7)
        ax.text(position[0], position[1], text,
                transform=ax.transAxes,
                verticalalignment='top',
                fontsize=fontsize,
                bbox=bbox_style)

    def load_data(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            try:
                self.data = pd.read_csv(file_path)
                self.data_file_path = file_path
                self.columns = list(self.data.columns)
                self.channel_names = self.detect_channels(self.columns)
                self.channel1_dropdown['values'] = self.channel_names
                self.channel1_dropdown.config(state="readonly")
                self.channel2_dropdown['values'] = self.channel_names
                self.channel2_dropdown.config(state="readonly")
                self.update_column_dropdowns()
                messagebox.showinfo("Success", "Data loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load data: {e}")

    def detect_channels(self, columns):
        common_parts = ["Bkgd_Intensity_", "BkgdCorrectedIntensity_", "Centroid_X_", "Centroid_Y_"]
        channel_names = set()
        for col in columns:
            for part in common_parts:
                if part in col:
                    channel_name = col.replace(part, "").rstrip("0123456789_")
                    channel_names.add(channel_name)
        return sorted(channel_names)

    def update_column_dropdowns(self, event=None):
        for var in self.column_vars:
            var.set("")
        channel1 = self.channel1_var.get()
        channel2 = self.channel2_var.get()
        column_patterns = {
            "Background Corrected Intensity:": "Bkgd_Corr_Intensity_{}",
            "Background Intensity:": "Bkgd_Intensity_{}",
            "Centroid X:": "X_{}_microns",
            "Centroid Y:": "Y_{}_microns"
        }

        def find_matching_column(pattern, channel):
            if not channel:
                return None
            search_pattern = pattern.format(channel)
            matches = [col for col in self.columns if search_pattern.lower() in col.lower()]
            return max(matches, key=len) if matches else None

        num_column_types = len(self.column_labels)
        if channel1:
            for i, label in enumerate(self.column_labels):
                pattern = column_patterns.get(label, "")
                if pattern:
                    matched_column = find_matching_column(pattern, channel1)
                    if matched_column:
                        self.column_vars[i].set(matched_column)
                        self.column_dropdowns[i].set(matched_column)
        if channel2:
            for i, label in enumerate(self.column_labels):
                pattern = column_patterns.get(label, "")
                if pattern:
                    matched_column = find_matching_column(pattern, channel2)
                    if matched_column:
                        index = i + num_column_types
                        self.column_vars[index].set(matched_column)
                        self.column_dropdowns[index].set(matched_column)

        for dropdown in self.column_dropdowns:
            dropdown['values'] = self.columns if self.columns else []
            dropdown.config(state="readonly" if self.columns else "disabled")

    def run_analysis(self):
        # Get the current values from the dropdowns
        bkgd_corr_ch1 = self.column_vars[0].get()
        bkgd_corr_ch2 = self.column_vars[4].get()
        bg_int_ch1 = self.column_vars[1].get()
        bg_int_ch2 = self.column_vars[5].get()
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
    
        if None in (bkgd_corr_ch1, bkgd_corr_ch2, bg_int_ch1, bg_int_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
            tk.messagebox.showerror("Error", "Please select all required columns")
            return
    
        self.analysis_results = []
    
        # Perform correlation analysis if selected
        if self.corr_var.get():
            self.collect_correlation_analysis(bkgd_corr_ch1, bkgd_corr_ch2)
    
        # Perform subcluster correlation if selected
        if self.subcluster_corr_var.get():
            self.collect_subcluster_correlation_analysis(bkgd_corr_ch1, bkgd_corr_ch2)
    
        # Perform background comparison if selected
        if self.bg_comp_var.get():
            self.collect_background_comparison(bg_int_ch1, bg_int_ch2)
    
        # Perform vector analysis if selected
        if self.vector_analysis_var.get():
            self.collect_vector_analysis(cx_ch1, cy_ch1, cx_ch2, cy_ch2)
    
        # Perform vector analysis (MRL) if selected
        if self.vector_analysis_mrl_var.get():
            self.collect_vector_analysis_mrl(cx_ch1, cy_ch1, cx_ch2, cy_ch2)
    
        # Perform triangle thresholding if selected
        if self.triangle_thresh_var.get():
            self.collect_triangle_thresholding(bkgd_corr_ch1, bkgd_corr_ch2)
    
        # Perform k-means thresholding if selected
        if self.kmeans_thresh_var.get():
            self.collect_kmeans_thresholding(bkgd_corr_ch1, bkgd_corr_ch2, bg_int_ch1, bg_int_ch2)
    
        # Perform GMM thresholding if selected
        if self.gmm_thresh_var.get():
            self.collect_gmm_thresholding(bkgd_corr_ch1, bkgd_corr_ch2, bg_int_ch1, bg_int_ch2)
    
        # Perform vector magnitude intensity threshold if selected
        if self.vector_magnitude_intensity_var.get():
            self.collect_vector_magnitude_intensity_threshold(
                bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2,
                alpha=self.vector_magnitude_alpha.get(),  # Use the dedicated alpha value
                threshold_method=self.threshold_method.get(),
                mrl_threshold=self.vector_magnitude_mrl_threshold.get()  # Use the dedicated MRL threshold
            )
    
        # Perform t-SNE and k-NN analysis if selected
        if self.tsne_knn_var.get():
            X_tsne, labels_tsne, k = self.perform_tsne_knn_analysis(bkgd_corr_ch1, bkgd_corr_ch2)
            self.plot_tsne_clusters(X_tsne, labels_tsne, k)
            self.plot_original_clusters(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
            self.plot_individual_clusters(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
    
        # Perform t-SNE and k-NN thresholding if selected
        if self.tsne_knn_thresh_var.get():
            if self.X_tsne is None or self.labels_tsne is None or self.k is None:
                messagebox.showwarning("Warning", "Please run t-SNE and kNN analysis first.")
            else:
                fig, x_threshold, y_threshold, cx_ch1, cy_ch1, cx_ch2, cy_ch2 = self._get_tsne_knn_thresholding_data(bkgd_corr_ch1, bkgd_corr_ch2)
                if fig is not None:
                    self.display_plot_in_new_window(fig, "t-SNE and kNN Thresholding")
                    if x_threshold is not None and y_threshold is not None:
                        self.generate_single_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, x_threshold, y_threshold)

    def generate_plots(self):
        if self.data is None:
            messagebox.showerror("Error", "No data loaded!")
            return
        channel1 = self.channel1_var.get()
        channel2 = self.channel2_var.get()
        if not channel1 or not channel2:
            messagebox.showerror("Error", "Please select both channels!")
            return
        bkgd_corr_ch1 = self.column_vars[0].get()
        bkgd_ch1 = self.column_vars[1].get()
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        bkgd_corr_ch2 = self.column_vars[4].get()
        bkgd_ch2 = self.column_vars[5].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        if not all([bkgd_corr_ch1, bkgd_ch1, cx_ch1, cy_ch1, bkgd_corr_ch2, bkgd_ch2, cx_ch2, cy_ch2]):
            messagebox.showerror("Error", "Please select all columns!")
            return
        if self.corr_var.get():
            self.correlation_analysis_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        if self.subcluster_corr_var.get():
            self.subcluster_correlation_analysis_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        if self.bg_comp_var.get():
            self.background_comparison_plot(bkgd_ch1, bkgd_ch2)
        if self.triangle_thresh_var.get():
            self.triangle_thresholding_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        if self.kmeans_thresh_var.get():
            self.kmeans_thresholding_plot(bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        if self.gmm_thresh_var.get():
            self.gmm_thresholding_plot(bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        if self.vector_analysis_var.get():
            self.vector_analysis_plot(cx_ch1, cy_ch1, cx_ch2, cy_ch2)
        if self.vector_analysis_mrl_var.get():
            self.vector_analysis_mrl_plot(cx_ch1, cy_ch1, cx_ch2, cy_ch2)
            
        if self.tsne_knn_var.get():
            X_tsne, labels_tsne, k = self.perform_tsne_knn_analysis(bkgd_corr_ch1, bkgd_corr_ch2)
            self.plot_tsne_clusters(X_tsne, labels_tsne, k)
            self.plot_original_clusters(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
            self.plot_individual_clusters(bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
    
        if self.tsne_knn_thresh_var.get():
            if self.X_tsne is None or self.labels_tsne is None or self.k is None:
                messagebox.showwarning("Warning", "Please run t-SNE and kNN analysis first.")
            else:
                fig, x_threshold, y_threshold, cx_ch1, cy_ch1, cx_ch2, cy_ch2 = self._get_tsne_knn_thresholding_data(bkgd_corr_ch1, bkgd_corr_ch2)
                if fig is not None:
                    self.display_plot_in_new_window(fig, "t-SNE and kNN Thresholding")
                    if x_threshold is not None and y_threshold is not None:
                        self.generate_single_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, x_threshold, y_threshold)
                
        if self.vector_magnitude_intensity_var.get():
            fig = self.vector_magnitude_intensity_threshold_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2)
            self.display_plot_in_new_window(fig, "Vector Magnitude Intensity Threshold Analysis")
            plt.close(fig)  # Close the figure to prevent memory leaks
    
    def export_plots(self):
        if self.data is None:
            messagebox.showerror("Error", "No data loaded!")
            return
        channel1 = self.channel1_var.get()
        channel2 = self.channel2_var.get()
        if not channel1 or not channel2:
            messagebox.showerror("Error", "Please select both channels!")
            return
        bkgd_corr_ch1 = self.column_vars[0].get()
        bkgd_ch1 = self.column_vars[1].get()
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        bkgd_corr_ch2 = self.column_vars[4].get()
        bkgd_ch2 = self.column_vars[5].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        if not all([bkgd_corr_ch1, bkgd_ch1, cx_ch1, cy_ch1, bkgd_corr_ch2, bkgd_ch2, cx_ch2, cy_ch2]):
            messagebox.showerror("Error", "Please select all columns!")
            return
        directory = filedialog.askdirectory()
        if not directory:
            return
        if self.corr_var.get():
            self.export_correlation_analysis_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2)
        if self.subcluster_corr_var.get():
            self.export_subcluster_correlation_analysis_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2)
        if self.bg_comp_var.get():
            self.export_background_comparison_plot(directory, bkgd_ch1, bkgd_ch2)
        if self.triangle_thresh_var.get():
            self.export_triangle_thresholding_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2)
        if self.kmeans_thresh_var.get():
            self.export_kmeans_thresholding_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        if self.gmm_thresh_var.get():
            self.export_gmm_thresholding_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        if self.vector_analysis_var.get():
            self.export_vector_analysis_plot(directory, cx_ch1, cy_ch1, cx_ch2, cy_ch2)
        if self.vector_analysis_mrl_var.get():
            self.export_vector_analysis_mrl_plot(directory, cx_ch1, cy_ch1, cx_ch2, cy_ch2)
            
        if self.tsne_knn_var.get():
            X_tsne, labels_tsne, k = self.perform_tsne_knn_analysis(bkgd_corr_ch1, bkgd_corr_ch2)
            self.export_tsne_clusters_plot(directory, X_tsne, labels_tsne, k, bkgd_corr_ch1, bkgd_corr_ch2)
            self.export_original_clusters_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)
            self.export_individual_clusters_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, labels_tsne, k)

    
        if self.tsne_knn_thresh_var.get():
            if self.X_tsne is None or self.labels_tsne is None or self.k is None:
                messagebox.showwarning("Warning", "Please run t-SNE and kNN analysis first.")
            else:
                fig, x_threshold, y_threshold, cx_ch1, cy_ch1, cx_ch2, cy_ch2 = self._get_tsne_knn_thresholding_data(bkgd_corr_ch1, bkgd_corr_ch2)
                if fig is not None:
                    filename = os.path.join(directory, f"tSNE_kNN_Thresholding_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
                    fig.savefig(filename, format=self.export_format.get())
                    plt.close(fig)
                    if x_threshold is not None and y_threshold is not None:
                        self.export_single_polar_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, x_threshold, y_threshold, "tSNE_kNN_Thresholding")
        
        if self.vector_magnitude_intensity_var.get():
            self.export_vector_magnitude_intensity_threshold_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2)
        messagebox.showinfo("Success", "Plots exported successfully!")


    def export_results_to_csv(self):
        if not self.analysis_results:
            messagebox.showerror("Error", "No analysis results to export!")
            return
        results_df = pd.DataFrame(self.analysis_results)
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile=f"{os.path.splitext(os.path.basename(self.data_file_path))[0]}_Analysis.csv"
        )
        if file_path:
            results_df.to_csv(file_path, index=False)
            messagebox.showinfo("Success", f"Results saved successfully to {file_path}")

    def display_plot_in_new_window(self, fig, title):
        new_window = tk.Toplevel(self.root)
        new_window.title(title)
        canvas = FigureCanvasTkAgg(fig, master=new_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def triangle_threshold(self, intensities):
        intensities = intensities.flatten()
        hist, bin_edges = np.histogram(intensities, bins=256)
        peak_index = np.argmax(hist)
        end_index = len(hist) - 1
        line = np.linspace(hist[peak_index], hist[end_index], end_index - peak_index)
        for i in range(peak_index, end_index):
            if hist[i] < line[i - peak_index]:
                threshold = bin_edges[i]
                break
        return threshold
    
    def to_superscript(self, n):
        superscripts = {
            "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
            "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"
        }
        return ''.join(superscripts.get(c, c) for c in str(int(n)))

    def log_tick_formatter(self, val, pos):
        if val == 0:
            return "0"
        abs_val = abs(val)
        if abs_val == 0:
            return "0"
        exponent = int(np.floor(np.log10(abs_val)))
        coefficient = val / (10**exponent)

        if abs(coefficient) == 1:
            sign = "-" if coefficient < 0 else ""
            return f"{sign}10{self.to_superscript(exponent)}"
        else:
            if coefficient == int(coefficient):
                coeff_str = f"{int(coefficient)}"
            else:
                coeff_str = f"{coefficient:.1f}".rstrip('0').rstrip('.')
            return f"{coeff_str}×10{self.to_superscript(exponent)}"

    def _setup_symlog_scale(self, ax, data1, data2=None):
        if data2 is not None:
            data_combined = np.concatenate([data1.values if hasattr(data1, 'values') else data1,
                                           data2.values if hasattr(data2, 'values') else data2])
        else:
            data_combined = data1.values if hasattr(data1, 'values') else data1
    
        # Filter out any NaN or infinity values
        data_combined = data_combined[~np.isnan(data_combined)]
        data_combined = data_combined[~np.isinf(data_combined)]
    
        if len(data_combined) == 0:
            linthresh = 1
        else:
            non_zero_data = data_combined[data_combined != 0]
            if len(non_zero_data) > 0:
                min_abs_val = np.min(np.abs(non_zero_data))
                # Set linthresh to the smallest power of 10 less than or equal to min_abs_val
                linthresh = 10 ** int(np.floor(np.log10(min_abs_val))) if min_abs_val != 0 else 1
            else:
                linthresh = 1
    
        ax.set_xscale('symlog', linthresh=linthresh)
        if data2 is not None:
            ax.set_yscale('symlog', linthresh=linthresh)
    
        formatter = FuncFormatter(self.log_tick_formatter)
        ax.xaxis.set_major_formatter(formatter)
        if data2 is not None:
            ax.yaxis.set_major_formatter(formatter)
    
        # Determine the range of exponents needed
        if len(data_combined) > 0:
            non_zero_data = data_combined[data_combined != 0]
    
            if len(non_zero_data) > 0:
                min_abs_val = np.min(np.abs(non_zero_data))
                max_abs_val = np.max(np.abs(non_zero_data))
    
                min_exponent = int(np.floor(np.log10(min_abs_val))) if min_abs_val != 0 else 0
                max_exponent = int(np.ceil(np.log10(max_abs_val))) if max_abs_val != 0 else 0
    
                exponents = list(range(min_exponent, max_exponent + 1))
            else:
                exponents = [0]
    
            ticks = [0]  # Always include zero
            for exp in exponents:
                if exp != 0:
                    ticks.append(10**exp)
                    ticks.append(-10**exp)
    
            # Remove duplicates and sort
            ticks = sorted(list(set(ticks)))
    
            # Set the ticks on both axes if data2 is not None, otherwise just on x-axis
            ax.set_xticks(ticks)
            if data2 is not None:
                ax.set_yticks(ticks)
    
        return linthresh


    
    def collect_correlation_analysis(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1]
        ch2_intensity = self.data[bkgd_corr_ch2]
        corr_coeff, p_value = pearsonr(ch1_intensity, ch2_intensity)
        p_value_rounded = f"{p_value:.5f}" if p_value >= 0.0001 else "p < 0.0001"
        status = "Correlated" if (p_value < 0.05 and abs(corr_coeff) > 0.5) else "Not Correlated"
        self.analysis_results.append({
            "Analysis Type": "Pearson Correlation",
            "Correlation Coefficient": corr_coeff,
            "P-value": p_value_rounded,
            "Status": status
        })
    
    def _create_correlation_plot(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1]
        ch2_intensity = self.data[bkgd_corr_ch2]
        corr_coeff, p_value = pearsonr(ch1_intensity, ch2_intensity)
        p_value_rounded = f"{p_value:.5f}" if p_value >= 0.0001 else "p < 0.0001"
        status = "Correlated" if (p_value < 0.05 and abs(corr_coeff) > 0.5) else "Not Correlated"
        fig, ax = plt.subplots(figsize=(6, 6))
        self._setup_symlog_scale(ax, ch1_intensity, ch2_intensity)
        ax.scatter(ch1_intensity, ch2_intensity, alpha=0.3)
        ax.set_xlabel(bkgd_corr_ch1)
        ax.set_ylabel(bkgd_corr_ch2)
        ax.set_title('Pearson Correlation Analysis')
    
        text_content = (
            f'Correlation Coefficient: {corr_coeff:.2f}\n'
            f'p-value: {p_value_rounded}\n'
            f'Status: {status}\n'
            f'Criteria: p < 0.05 and |r| > 0.5'
        )
    
        # Place text box in the bottom-left corner
        ax.text(
            0.02, 0.02, text_content,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )
    
        plt.tight_layout()
        return fig
    
    def correlation_analysis_plot(self, bkgd_corr_ch1, bkgd_corr_ch2):
        fig = self._create_correlation_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        self.display_plot_in_new_window(fig, "Pearson Correlation Analysis")
    
    def export_correlation_analysis_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2):
        fig = self._create_correlation_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        filename = os.path.join(directory, f"Pearson_Correlation_Analysis_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)

    
    def collect_vector_analysis_mrl(self, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        # Get the current MRL threshold from the UI
        mrl_threshold = self.mrl_threshold.get()
    
        # Calculate vectors
        centroids_ch1_x = self.data[cx_ch1]
        centroids_ch1_y = self.data[cy_ch1]
        centroids_ch2_x = self.data[cx_ch2]
        centroids_ch2_y = self.data[cy_ch2]
    
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_angles = np.arctan2(vectors_y, vectors_x)
        vector_magnitudes = np.sqrt(vectors_x**2 + vectors_y**2)
    
        # Combine magnitudes and angles, and sort by magnitude (ascending)
        combined = list(zip(vector_magnitudes, vector_angles))
        sort_order = self.sorting_order.get()
        combined_sorted = sorted(combined, key=lambda x: x[0], reverse=(sort_order == 'descending'))
    
        # Find threshold based on MRL < threshold value (using UI value)
        threshold = 0
        threshold_index = 0
        mrl_values = []  # To store MRL values for potential plotting
    
        for i in range(1, len(combined_sorted) + 1):
            current_subset = combined_sorted[:i]
            current_angles = np.array([x[1] for x in current_subset])
            mrl = self.calculate_mrl(current_angles)
            mrl_values.append(mrl)
            if mrl < mrl_threshold:  # Use the UI-specified MRL threshold
                threshold = combined_sorted[i-1][0]  # magnitude of the last added vector
                threshold_index = i
                break
    
        # If we never find MRL < threshold, set threshold to max magnitude
        if threshold == 0:
            if len(combined_sorted) > 0:
                threshold = combined_sorted[-1][0]
            else:
                threshold = 0
    
        # Separate vectors based on threshold
        below_threshold_angles = []
        above_threshold_angles = []
    
        for mag, angle in combined_sorted:
            if mag <= threshold:
                below_threshold_angles.append(angle)
            else:
                above_threshold_angles.append(angle)
    
        # Convert to numpy arrays for easier calculations
        below_threshold_angles = np.array(below_threshold_angles)
        above_threshold_angles = np.array(above_threshold_angles)
    
        # Calculate statistics for all vectors
        mean_angle = circmean(vector_angles)
        std_angle = circstd(vector_angles)
        mrl_all = self.calculate_mrl(vector_angles)
    
        # Calculate Rayleigh p-value for all vectors
        rayleigh_p_value_all = self.rayleigh_test_p_value(vector_angles)
    
        # Calculate statistics for below threshold
        below_stats = {}
        if len(below_threshold_angles) > 0:
            below_stats["mean_angle_deg"] = np.degrees(circmean(below_threshold_angles))
            below_stats["std_angle_deg"] = np.degrees(circstd(below_threshold_angles))
            below_stats["mrl"] = self.calculate_mrl(below_threshold_angles)
            below_stats["rayleigh_p_value"] = self.rayleigh_test_p_value(below_threshold_angles)
        else:
            below_stats["mean_angle_deg"] = None
            below_stats["std_angle_deg"] = None
            below_stats["mrl"] = None
            below_stats["rayleigh_p_value"] = None
    
        # Calculate statistics for above threshold
        above_stats = {}
        if len(above_threshold_angles) > 0:
            above_stats["mean_angle_deg"] = np.degrees(circmean(above_threshold_angles))
            above_stats["std_angle_deg"] = np.degrees(circstd(above_threshold_angles))
            above_stats["mrl"] = self.calculate_mrl(above_threshold_angles)
            above_stats["rayleigh_p_value"] = self.rayleigh_test_p_value(above_threshold_angles)
        else:
            above_stats["mean_angle_deg"] = None
            above_stats["std_angle_deg"] = None
            above_stats["mrl"] = None
            above_stats["rayleigh_p_value"] = None
    
        # Determine significance based on MRL threshold from UI
        status_all = "Significant" if mrl_all > mrl_threshold else "Not Significant"
        status_below = "Significant" if below_stats["mrl"] is not None and below_stats["mrl"] > mrl_threshold else "Not Significant"
        status_above = "Significant" if above_stats["mrl"] is not None and above_stats["mrl"] > mrl_threshold else "Not Significant"
    
        # Add results to analysis_results
        self.analysis_results.append({
            "Analysis Type": "Vector Analysis (MRL Threshold)",
            "Mean Angle (degrees)": np.degrees(mean_angle),
            "Standard Deviation of Angles (degrees)": np.degrees(std_angle),
            "Rayleigh Test p-value (all)": rayleigh_p_value_all,
            "Mean Resultant Length (MRL) (all)": mrl_all,
            "Directionality Status (all)": status_all,
            "MRL Threshold Used": mrl_threshold,
            "Threshold for Magnitudes": threshold,
            "Below Threshold Count": len(below_threshold_angles),
            "Above Threshold Count": len(above_threshold_angles),
            "Below Threshold Mean Angle (degrees)": below_stats["mean_angle_deg"],
            "Below Threshold Std Angle (degrees)": below_stats["std_angle_deg"],
            "Below Threshold Rayleigh p-value": below_stats["rayleigh_p_value"],
            "Below Threshold MRL": below_stats["mrl"],
            "Below Threshold Status": status_below,
            "Above Threshold Mean Angle (degrees)": above_stats["mean_angle_deg"],
            "Above Threshold Std Angle (degrees)": above_stats["std_angle_deg"],
            "Above Threshold Rayleigh p-value": above_stats["rayleigh_p_value"],
            "Above Threshold MRL": above_stats["mrl"],
            "Above Threshold Status": status_above
        })
    
    def _generate_vector_analysis_plots(self, cx_ch1, cy_ch1, cx_ch2, cy_ch2, mrl_threshold=0.7):
        # Calculate vector magnitudes and angles
        centroids_ch1_x = self.data[cx_ch1]
        centroids_ch1_y = self.data[cy_ch1]
        centroids_ch2_x = self.data[cx_ch2]
        centroids_ch2_y = self.data[cy_ch2]
    
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_angles = np.arctan2(vectors_y, vectors_x)
        vector_magnitudes = np.sqrt(vectors_x**2 + vectors_y**2)
    
        # Calculate threshold using the MRL method
        combined = list(zip(vector_magnitudes, vector_angles))
        sort_order = self.sorting_order.get()
        combined_sorted = sorted(combined, key=lambda x: x[0], reverse=(sort_order == 'descending'))
    
        threshold = 0
        mrl_values = []
        cumulative_mrl = []
    
        for i in range(1, len(combined_sorted) + 1):
            current_subset = combined_sorted[:i]
            current_angles = np.array([x[1] for x in current_subset])
            mrl = self.calculate_mrl(current_angles)
            mrl_values.append(mrl)
            cumulative_mrl.append((i, mrl))
            if mrl < mrl_threshold and threshold == 0:  # Use the provided mrl_threshold instead of hardcoded 0.7
                threshold = combined_sorted[i-1][0]
    
        if threshold == 0 and len(combined_sorted) > 0:
            threshold = combined_sorted[-1][0]
    
        # Create masks for below and above threshold
        below_threshold_mask = vector_magnitudes <= threshold
        above_threshold_mask = vector_magnitudes > threshold
    
        # Calculate statistics for plotting
        total_vectors = len(vector_magnitudes)
        num_below_threshold = np.sum(below_threshold_mask)
        percent_below_threshold = (num_below_threshold / total_vectors) * 100 if total_vectors > 0 else 0
        num_above_threshold = np.sum(above_threshold_mask)
        percent_above_threshold = (num_above_threshold / total_vectors) * 100 if total_vectors > 0 else 0
    
        # Calculate statistics for all vectors
        mrl_all = self.calculate_mrl(vector_angles)
        rayleigh_p_value_all = self.rayleigh_test_p_value(vector_angles)
    
        # Calculate corrected p-value for all vectors
        corrected_p_value_all = min(rayleigh_p_value_all * 3, 1.0)  # Bonferroni correction for 3 tests
        formatted_p_value_all = f"{corrected_p_value_all:.5f}" if corrected_p_value_all >= 0.0001 else "p < 0.0001"
    
        p_value_text_all = (
            f'Rayleigh corrected p-value: {formatted_p_value_all}\n'
            f'MRL: {mrl_all:.3f}\n'
            f'Directionality: {"Significant" if mrl_all > mrl_threshold else "Not Significant"}'
        )
    
        # For below threshold vectors
        if np.sum(below_threshold_mask) > 0:
            below_threshold_angles = vector_angles[below_threshold_mask]
            rayleigh_p_value_below = self.rayleigh_test_p_value(below_threshold_angles)
            mrl_below = self.calculate_mrl(below_threshold_angles)
    
            # Calculate corrected p-value
            corrected_p_value_below = min(rayleigh_p_value_below * 3, 1.0)
            formatted_p_value_below = f"{corrected_p_value_below:.5f}" if corrected_p_value_below >= 0.0001 else "p < 0.0001"
    
            p_value_text_below = (
                f'Rayleigh corrected p-value: {formatted_p_value_below}\n'
                f'MRL: {mrl_below:.3f}\n'
                f'Directionality: {"Significant" if mrl_below and mrl_below > mrl_threshold else "Not Significant"}'
            )
    
        # For above threshold vectors
        if np.sum(above_threshold_mask) > 0:
            above_threshold_angles = vector_angles[above_threshold_mask]
            rayleigh_p_value_above = self.rayleigh_test_p_value(above_threshold_angles)
            mrl_above = self.calculate_mrl(above_threshold_angles)
    
            # Calculate corrected p-value
            corrected_p_value_above = min(rayleigh_p_value_above * 3, 1.0)
            formatted_p_value_above = f"{corrected_p_value_above:.5f}" if corrected_p_value_above >= 0.0001 else "p < 0.0001"
    
            p_value_text_above = (
                f'Rayleigh corrected p-value: {formatted_p_value_above}\n'
                f'MRL: {mrl_above:.3f}\n'
                f'Directionality: {"Significant" if mrl_above and mrl_above > mrl_threshold else "Not Significant"}'
            )

        # Create figure with 4 subplots: 3 polar plots and 1 MRL progression plot
        fig = plt.figure(figsize=(20, 15))
    
        # 1. All vectors polar plot
        ax1 = plt.subplot(2, 2, 1, polar=True)
        n_bins = 36
        bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        counts, _ = np.histogram(vector_angles, bins=bin_edges)
        bars = ax1.bar(bin_centers, counts, width=(2*np.pi/n_bins), bottom=0.0,
                      color='green', alpha=0.5, edgecolor='k', linewidth=0.5)
    
        p_value_text_all = (
            f'Rayleigh p-value: {rayleigh_p_value_all:.5f}\n'
            f'MRL: {mrl_all:.3f}\n'
            f'Directionality: {"Significant" if mrl_all > mrl_threshold else "Not Significant"}'
        )
        ax1.text(0.5, 0.1, p_value_text_all, transform=ax1.transAxes, fontsize=10,
                horizontalalignment='center', verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
        ax1.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
        ax1.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
        ax1.set_rlabel_position(30)
        ax1.grid(True)
        ax1.set_title('All Vectors', y=1.1)
    
        # 2. Below threshold vectors polar plot
        ax2 = plt.subplot(2, 2, 2, polar=True)
        if np.sum(below_threshold_mask) > 0:
            below_threshold_angles = vector_angles[below_threshold_mask]
            counts_below, _ = np.histogram(below_threshold_angles, bins=bin_edges)
            bars_below = ax2.bar(bin_centers, counts_below, width=(2*np.pi/n_bins), bottom=0.0,
                                color='blue', alpha=0.5, edgecolor='k', linewidth=0.5)
    
            p_value_text_below = (
                f'Rayleigh p-value: {rayleigh_p_value_below:.5f}\n'
                f'MRL: {mrl_below:.3f}\n'
                f'Directionality: {"Significant" if mrl_below and mrl_below > mrl_threshold else "Not Significant"}'
            )
            ax2.text(0.5, 0.1, p_value_text_below, transform=ax2.transAxes, fontsize=10,
                    horizontalalignment='center', verticalalignment='bottom',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
            ax2.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
            ax2.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
            ax2.set_rlabel_position(30)
            ax2.grid(True)
        else:
            ax2.set_title('No vectors below threshold', y=1.1)
        ax2.set_title(f'Below Threshold ({num_below_threshold} vectors, {percent_below_threshold:.1f}%)', y=1.1)
    
        # 3. Above threshold vectors polar plot
        ax3 = plt.subplot(2, 2, 3, polar=True)
        if np.sum(above_threshold_mask) > 0:
            above_threshold_angles = vector_angles[above_threshold_mask]
            counts_above, _ = np.histogram(above_threshold_angles, bins=bin_edges)
            bars_above = ax3.bar(bin_centers, counts_above, width=(2*np.pi/n_bins), bottom=0.0,
                                color='red', alpha=0.5, edgecolor='k', linewidth=0.5)
    
            p_value_text_above = (
                f'Rayleigh p-value: {rayleigh_p_value_above:.5f}\n'
                f'MRL: {mrl_above:.3f}\n'
                f'Directionality: {"Significant" if mrl_above and mrl_above > mrl_threshold else "Not Significant"}'
            )
            ax3.text(0.5, 0.1, p_value_text_above, transform=ax3.transAxes, fontsize=10,
                    horizontalalignment='center', verticalalignment='bottom',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
            ax3.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
            ax3.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
            ax3.set_rlabel_position(30)
            ax3.grid(True)
        else:
            ax3.set_title('No vectors above threshold', y=1.1)
        ax3.set_title(f'Above Threshold ({num_above_threshold} vectors, {percent_above_threshold:.1f}%)', y=1.1)
    
        # 4. Plot showing how MRL changes as we add more vectors (sorted by magnitude)
        ax4 = plt.subplot(2, 2, 4)
        if len(cumulative_mrl) > 0:
            x_vals, y_vals = zip(*cumulative_mrl)
            ax4.plot(x_vals, y_vals, 'b-', linewidth=2)
            ax4.axhline(y=mrl_threshold, color='r', linestyle='--', label=f'MRL Threshold ({mrl_threshold})')
    
            # Mark the threshold point if it exists
            if threshold > 0:
                # Find the index where we set the threshold
                threshold_idx = None
                for i, (count, mrl) in enumerate(cumulative_mrl):
                    if mrl < mrl_threshold and (i == 0 or cumulative_mrl[i-1][1] >= mrl_threshold):
                        threshold_idx = i
                        break
    
                if threshold_idx is not None:
                    ax4.scatter(x_vals[threshold_idx], y_vals[threshold_idx],
                               color='red', s=100, label=f'Threshold at vector {x_vals[threshold_idx]}')
                    ax4.annotate(f'Threshold\n({x_vals[threshold_idx]} vectors)',
                                xy=(x_vals[threshold_idx], y_vals[threshold_idx]),
                                xytext=(x_vals[threshold_idx]+5, y_vals[threshold_idx]+0.1),
                                arrowprops=dict(facecolor='black', shrink=0.05))
    
            ax4.set_xlabel('Number of vectors (sorted by magnitude)')
            ax4.set_ylabel('MRL')
            ax4.set_title('MRL Progression as Vectors are Added (Sorted by Magnitude)')
            ax4.legend()
            ax4.grid(True)
    
        plt.tight_layout()
    
        # Plot magnitude distribution with threshold
        fig_magnitudes = plt.figure(figsize=(12, 6))
        ax_magnitudes = fig_magnitudes.add_subplot(111)
        if len(vector_magnitudes) > 0:
            sns.histplot(vector_magnitudes, bins=50, stat='density',
                        color='blue', alpha=0.3, ax=ax_magnitudes)
            ax_magnitudes.axvline(threshold, color='purple', linestyle=':',
                                linewidth=2, label=f'Threshold: {threshold:.2f}')
            ax_magnitudes.set_xlabel('Vector Magnitude')
            ax_magnitudes.set_ylabel('Density')
            ax_magnitudes.set_title('Distribution of Vector Magnitudes with MRL Threshold')
            ax_magnitudes.legend()
    
        return fig, fig_magnitudes, threshold
    
    def vector_analysis_mrl_plot(self, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        mrl_threshold = self.mrl_threshold.get()
        fig, fig_magnitudes, threshold = self._generate_vector_analysis_plots(cx_ch1, cy_ch1, cx_ch2, cy_ch2, mrl_threshold=mrl_threshold)
        self.display_plot_in_new_window(fig, "Vector Analysis (MRL Threshold)")
        self.display_plot_in_new_window(fig_magnitudes, "Vector Magnitudes with MRL Threshold")
        plt.close(fig)
        plt.close(fig_magnitudes)
    
    def export_vector_analysis_mrl_plot(self, directory, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        mrl_threshold = self.mrl_threshold.get()
        fig, fig_magnitudes, threshold = self._generate_vector_analysis_plots(cx_ch1, cy_ch1, cx_ch2, cy_ch2, mrl_threshold=mrl_threshold)
        filename_polar = os.path.join(directory, f"Vector_Analysis_MRL_Polar_{cx_ch1}_to_{cx_ch2}.{self.export_format.get()}")
        fig.savefig(filename_polar, format=self.export_format.get())
        plt.close(fig)
    
        filename_magnitudes = os.path.join(directory, f"Vector_Analysis_MRL_Magnitudes_{cx_ch1}_to_{cx_ch2}.{self.export_format.get()}")
        fig_magnitudes.savefig(filename_magnitudes, format=self.export_format.get())
        plt.close(fig_magnitudes)

    def collect_subcluster_correlation_analysis(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values.reshape(-1, 1)
        ch2_intensity = self.data[bkgd_corr_ch2].values.reshape(-1, 1)
        combined_intensities = np.hstack((ch1_intensity, ch2_intensity))
        kmeans = KMeans(n_clusters=3, random_state=0).fit(combined_intensities)
        labels = kmeans.labels_
        for i in range(3):
            cluster_mask = labels == i
            ch1_cluster = ch1_intensity[cluster_mask]
            ch2_cluster = ch2_intensity[cluster_mask]
            if len(ch1_cluster) > 1 and len(ch2_cluster) > 1:
                corr_coeff, p_value = pearsonr(ch1_cluster.flatten(), ch2_cluster.flatten())
                p_value_rounded = f"{p_value:.5f}" if p_value >= 0.0001 else "p < 0.0001"
                status = "Correlated" if (p_value < 0.05 and abs(corr_coeff) > 0.5) else "Not Correlated"
                self.analysis_results.append({
                    "Analysis Type": "Sub-cluster Correlation",
                    "Cluster ID": i + 1,
                    "Correlation Coefficient": corr_coeff,
                    "P-value": p_value_rounded,
                    "Status": status
                })
    
    def _create_subcluster_correlation_plot(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values.reshape(-1, 1)
        ch2_intensity = self.data[bkgd_corr_ch2].values.reshape(-1, 1)
        combined_intensities = np.hstack((ch1_intensity, ch2_intensity))
        kmeans = KMeans(n_clusters=3, random_state=0).fit(combined_intensities)
        labels = kmeans.labels_
        fig, axs = plt.subplots(1, 3, figsize=(18, 5))
        for i in range(3):
            cluster_mask = labels == i
            ch1_cluster = ch1_intensity[cluster_mask]
            ch2_cluster = ch2_intensity[cluster_mask]
            if len(ch1_cluster) > 1 and len(ch2_cluster) > 1:
                corr_coeff, p_value = pearsonr(ch1_cluster.flatten(), ch2_cluster.flatten())
                p_value_rounded = f"{p_value:.5f}" if p_value >= 0.0001 else "p < 0.0001"
                status = "Correlated" if (p_value < 0.05 and abs(corr_coeff) > 0.5) else "Not Correlated"
                self._setup_symlog_scale(axs[i], ch1_cluster.flatten(), ch2_cluster.flatten())
                axs[i].scatter(ch1_cluster, ch2_cluster)
                axs[i].set_xlabel(bkgd_corr_ch1)
                axs[i].set_ylabel(bkgd_corr_ch2)
                axs[i].set_title(f'Sub-cluster {i+1} Correlation Analysis')
    
                # Create text content
                text_content = (
                    f'Correlation Coefficient: {corr_coeff:.2f}\n'
                    f'p-value: {p_value_rounded}\n'
                    f'Status: {status}\n'
                    f'Criteria: p < 0.05 and |r| > 0.5'
                )
    
                # Place text box in the bottom-left corner
                axs[i].text(
                    0.02, 0.02, text_content,
                    transform=axs[i].transAxes,
                    fontsize=9,
                    verticalalignment='bottom',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )
    
        plt.tight_layout()
        return fig

    
    def subcluster_correlation_analysis_plot(self, bkgd_corr_ch1, bkgd_corr_ch2):
        fig = self._create_subcluster_correlation_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        self.display_plot_in_new_window(fig, "Sub-cluster Correlation Analysis")
    
    def export_subcluster_correlation_analysis_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2):
        fig = self._create_subcluster_correlation_plot(bkgd_corr_ch1, bkgd_corr_ch2)
        filename = os.path.join(directory, f"Subcluster_Correlation_Analysis_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)


    def collect_background_comparison(self, bkgd_ch1, bkgd_ch2):
        ch1_background = self.data[bkgd_ch1]
        ch2_background = self.data[bkgd_ch2]
        t_stat, t_p_value = ttest_ind(ch1_background, ch2_background)
        u_stat, u_p_value = mannwhitneyu(ch1_background, ch2_background)
        t_p_value_rounded = f"{t_p_value:.5f}" if t_p_value >= 0.0001 else "p < 0.0001"
        u_p_value_rounded = f"{u_p_value:.5f}" if u_p_value >= 0.0001 else "p < 0.0001"
        self.analysis_results.append({
            "Analysis Type": "Background Comparison",
            "T-test Statistic": t_stat,
            "T-test P-value": t_p_value_rounded,
            "U Statistic": u_stat,
            "U P-value": u_p_value_rounded
        })

    def _create_background_comparison_plot(self, bkgd_ch1, bkgd_ch2):
        ch1_background = self.data[bkgd_ch1]
        ch2_background = self.data[bkgd_ch2]
        t_stat, t_p_value = ttest_ind(ch1_background, ch2_background)
        u_stat, u_p_value = mannwhitneyu(ch1_background, ch2_background)
        t_p_value_rounded = f"{t_p_value:.5f}" if t_p_value >= 0.0001 else "p < 0.0001"
        u_p_value_rounded = f"{u_p_value:.5f}" if u_p_value >= 0.0001 else "p < 0.0001"
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.boxplot([ch1_background, ch2_background], tick_labels=[bkgd_ch1, bkgd_ch2])
        ax.set_title('Background Comparison')
        text_content = f'T-test p-value: {t_p_value_rounded}\n' \
                      f'U-test p-value: {u_p_value_rounded}'
        self.add_plot_text(ax, text_content, position=(0.5, 0.95))
        plt.tight_layout()
        return fig
    
    def background_comparison_plot(self, bkgd_ch1, bkgd_ch2):
        fig = self._create_background_comparison_plot(bkgd_ch1, bkgd_ch2)
        self.display_plot_in_new_window(fig, "Background Comparison")
    
    def export_background_comparison_plot(self, directory, bkgd_ch1, bkgd_ch2):
        fig = self._create_background_comparison_plot(bkgd_ch1, bkgd_ch2)
        filename = os.path.join(directory, f"Background_Comparison_{bkgd_ch1}_vs_{bkgd_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)

    def collect_gmm_thresholding(self, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values.reshape(-1, 1)
        ch2_intensity = self.data[bkgd_corr_ch2].values.reshape(-1, 1)
        ch1_background = self.data[bkgd_ch1].values.reshape(-1, 1)
        ch2_background = self.data[bkgd_ch2].values.reshape(-1, 1)
        combined_data_ch1 = np.hstack((ch1_intensity, ch1_background))
        combined_data_ch2 = np.hstack((ch2_intensity, ch2_background))
        def find_optimal_clusters_gmm(data, max_clusters=10):
            bics = []
            for k in range(1, max_clusters + 1):
                gmm = GaussianMixture(n_components=k, random_state=0).fit(data)
                bics.append(gmm.bic(data))
            optimal_k = np.argmin(bics) + 1  # BIC is minimized for the optimal model
            return optimal_k
        optimal_k_ch1 = find_optimal_clusters_gmm(combined_data_ch1)
        optimal_k_ch2 = find_optimal_clusters_gmm(combined_data_ch2)
        gmm_ch1 = GaussianMixture(n_components=optimal_k_ch1, random_state=0).fit(combined_data_ch1)
        means_ch1 = np.sort(gmm_ch1.means_[:, 0])
        gmm_thresholds_ch1 = [(means_ch1[i] + means_ch1[i + 1]) / 2 for i in range(len(means_ch1) - 1)]
        threshold_ch1 = gmm_thresholds_ch1[0] if gmm_thresholds_ch1 else None
        gmm_ch2 = GaussianMixture(n_components=optimal_k_ch2, random_state=0).fit(combined_data_ch2)
        means_ch2 = np.sort(gmm_ch2.means_[:, 0])
        gmm_thresholds_ch2 = [(means_ch2[i] + means_ch2[i + 1]) / 2 for i in range(len(means_ch2) - 1)]
        threshold_ch2 = gmm_thresholds_ch2[0] if gmm_thresholds_ch2 else None
        self.analysis_results.append({
            "Analysis Type": "GMM Thresholding",
            "Channel 1 Optimal Clusters": optimal_k_ch1,
            "Channel 1 Threshold": threshold_ch1,
            "Channel 2 Optimal Clusters": optimal_k_ch2,
            "Channel 2 Threshold": threshold_ch2
        })


    def collect_triangle_thresholding(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
        triangle_threshold_ch1 = self.triangle_threshold(ch1_intensity)
        triangle_threshold_ch2 = self.triangle_threshold(ch2_intensity)
        self.analysis_results.append({
            "Analysis Type": "Triangle Thresholding",
            "Channel 1 Threshold": triangle_threshold_ch1,
            "Channel 2 Threshold": triangle_threshold_ch2
        })

    def _generate_polar_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2):
        # Get the intensities and coordinates
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
        centroids_ch1_x = self.data[cx_ch1].values
        centroids_ch1_y = self.data[cy_ch1].values
        centroids_ch2_x = self.data[cx_ch2].values
        centroids_ch2_y = self.data[cy_ch2].values
    
        # Calculate vectors
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_angles = np.arctan2(vectors_y, vectors_x)
        vector_magnitudes = np.sqrt(vectors_x**2 + vectors_y**2)
    
        # Classify data points
        positive_both = (ch1_intensity > threshold_ch1) & (ch2_intensity > threshold_ch2)
    
        # Create a figure for the polar plot
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, polar=True)
    
        if np.sum(positive_both) > 0:
            angles_both = vector_angles[positive_both]
            n_bins = 36
            bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            counts, _ = np.histogram(angles_both, bins=bin_edges)
            bars = ax.bar(bin_centers, counts, width=(2*np.pi/n_bins), bottom=0.0, color='green', alpha=0.5, edgecolor='k', linewidth=0.5)
    
            # Calculate Rayleigh p-value and MRL
            mean_angle = circmean(angles_both)
            R = np.sqrt((np.sum(np.cos(angles_both)))**2 + (np.sum(np.sin(angles_both)))**2)
            n = len(angles_both)
            r = R / n if n > 0 else 0
            z = 2 * n * r**2 if n > 0 else 0
            rayleigh_p_value = chi2.sf(z, df=2) if n > 0 else 1.0
            mrl = r
    
            # Determine directionality status
            directionality_status = "Significant" if (rayleigh_p_value < 0.05 and mrl > 0.7) else "Not Significant"
    
            # Display Rayleigh p-value, MRL, and directionality status
            p_value_text = f'Rayleigh p-value: {rayleigh_p_value:.5f}' if rayleigh_p_value >= 0.0001 else 'Rayleigh p-value: < 0.0001'
            mrl_text = f'MRL: {mrl:.3f}'
            status_text = f'Directionality: {directionality_status}'
            combined_text = f'{p_value_text}\n{mrl_text}\n{status_text}'
            ax.text(0.5, 0.1, combined_text, transform=ax.transAxes, fontsize=10, horizontalalignment='center', verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
            # Plot mean angle if directionality is significant
            if rayleigh_p_value < 0.05 and mrl > 0.7 and mean_angle is not None:
                ax.plot([0, mean_angle], [0, np.max(counts)], 'r-', linewidth=2)
                ax.text(0.5, 0.5, f'Mean Angle: {np.degrees(mean_angle):.2f}°', transform=ax.transAxes, fontsize=10, horizontalalignment='center', verticalalignment='center', color='red', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
            ax.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
            ax.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
            ax.set_rlabel_position(30)
            ax.grid(True)
            ax.set_title(f'Both Channels Positive ({np.sum(positive_both)} vectors)', y=1.1)
        else:
            ax.set_title('No vectors for Both Channels Positive', y=1.1)
    
        plt.tight_layout()
        return fig
    
    def generate_single_polar_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2):
        fig = self._generate_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2)
        self.display_plot_in_new_window(fig, "Polar Plot for Both Channels Above Threshold")
    
    def export_single_polar_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2, prefix):
        fig = self._generate_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2)
        filename = os.path.join(directory, f"{prefix}_Single_Polar_Plot_Both_Above_Threshold.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)


    def collect_kmeans_thresholding(self, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values.reshape(-1, 1)
        ch2_intensity = self.data[bkgd_corr_ch2].values.reshape(-1, 1)
        ch1_background = self.data[bkgd_ch1].values.reshape(-1, 1)
        ch2_background = self.data[bkgd_ch2].values.reshape(-1, 1)
        combined_data_ch1 = np.hstack((ch1_intensity, ch1_background))
        combined_data_ch2 = np.hstack((ch2_intensity, ch2_background))

        def find_optimal_clusters(data, max_clusters=10):
            distortions = []
            for k in range(2, max_clusters + 1):
                kmeans = KMeans(n_clusters=k, random_state=0).fit(data)
                distortions.append(kmeans.inertia_)
            deltas = np.diff(distortions)
            optimal_k = np.argmin(deltas / deltas[0]) + 2
            return optimal_k

        optimal_k_ch1 = find_optimal_clusters(combined_data_ch1)
        optimal_k_ch2 = find_optimal_clusters(combined_data_ch2)
        kmeans_ch1 = KMeans(n_clusters=optimal_k_ch1, random_state=0).fit(combined_data_ch1)
        centers_ch1 = np.sort(kmeans_ch1.cluster_centers_[:, 0])
        kmeans_thresholds_ch1 = [(centers_ch1[i] + centers_ch1[i + 1]) / 2 for i in range(len(centers_ch1) - 1)]
        kmeans_ch2 = KMeans(n_clusters=optimal_k_ch2, random_state=0).fit(combined_data_ch2)
        centers_ch2 = np.sort(kmeans_ch2.cluster_centers_[:, 0])
        kmeans_thresholds_ch2 = [(centers_ch2[i] + centers_ch2[i + 1]) / 2 for i in range(len(centers_ch2) - 1)]
        threshold_ch1 = kmeans_thresholds_ch1[0] if kmeans_thresholds_ch1 else None
        threshold_ch2 = kmeans_thresholds_ch2[0] if kmeans_thresholds_ch2 else None
        self.analysis_results.append({
            "Analysis Type": "K-means Thresholding",
            "Channel 1 Optimal Clusters": optimal_k_ch1,
            "Channel 1 Threshold": threshold_ch1,
            "Channel 2 Optimal Clusters": optimal_k_ch2,
            "Channel 2 Threshold": threshold_ch2
        })

    def _create_thresholding_plot(self, ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2):
        fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    
        # Plot channel 1 histogram
        ax[0].hist(ch1_intensity, bins=256, alpha=0.6, color='blue')
        if threshold_ch1:
            ax[0].axvline(threshold_ch1, color='purple', linestyle='dashed', linewidth=1)
            text_content = f'Threshold: {threshold_ch1:.2f}'
            self.add_plot_text(ax[0], text_content)
        ch1_name = self.channel1_var.get()
        ax[0].set_title(f'{ch1_name} Intensity')
    
        # Plot channel 2 histogram
        ax[1].hist(ch2_intensity, bins=256, alpha=0.6, color='green')
        if threshold_ch2:
            ax[1].axvline(threshold_ch2, color='purple', linestyle='dashed', linewidth=1)
            text_content = f'Threshold: {threshold_ch2:.2f}'
            self.add_plot_text(ax[1], text_content)
        ch2_name = self.channel2_var.get()
        ax[1].set_title(f'{ch2_name} Intensity')
    
        # Create the scatter plot with symlog scale
        ax[2].scatter(ch1_intensity, ch2_intensity, alpha=0.3)
        self._setup_symlog_scale(ax[2], ch1_intensity, ch2_intensity)
    
        # Add threshold lines
        if threshold_ch1 and threshold_ch2:
            ax[2].axvline(threshold_ch1, color='red', linestyle='dashed', linewidth=1, label=f'Channel 1 Threshold: {threshold_ch1:.2f}')
            ax[2].axhline(threshold_ch2, color='green', linestyle='dashed', linewidth=1, label=f'Channel 2 Threshold: {threshold_ch2:.2f}')
    
        ax[2].set_xlabel(f'{ch1_name} Intensity')
        ax[2].set_ylabel(f'{ch2_name} Intensity')
        ax[2].set_title('Channel Intensities Scatter Plot (Symlog Scale)')
        ax[2].legend()
    
        plt.tight_layout()
    
        return fig, threshold_ch1, threshold_ch2
    
    def _prepare_data(self, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values.reshape(-1, 1)
        ch2_intensity = self.data[bkgd_corr_ch2].values.reshape(-1, 1)
        ch1_background = self.data[bkgd_ch1].values.reshape(-1, 1)
        ch2_background = self.data[bkgd_ch2].values.reshape(-1, 1)
        combined_data_ch1 = np.hstack((ch1_intensity, ch1_background))
        combined_data_ch2 = np.hstack((ch2_intensity, ch2_background))
        return ch1_intensity, ch2_intensity, combined_data_ch1, combined_data_ch2
    
    def _compute_kmeans_thresholds(self, combined_data_ch1, combined_data_ch2):
        def find_optimal_clusters(data, max_clusters=10):
            distortions = []
            for k in range(2, max_clusters + 1):
                kmeans = KMeans(n_clusters=k, random_state=0).fit(data)
                distortions.append(kmeans.inertia_)
            deltas = np.diff(distortions)
            optimal_k = np.argmin(deltas / deltas[0]) + 2
            return optimal_k
    
        optimal_k_ch1 = find_optimal_clusters(combined_data_ch1)
        optimal_k_ch2 = find_optimal_clusters(combined_data_ch2)
    
        kmeans_ch1 = KMeans(n_clusters=optimal_k_ch1, random_state=0).fit(combined_data_ch1)
        centers_ch1 = np.sort(kmeans_ch1.cluster_centers_[:, 0])
        kmeans_thresholds_ch1 = [(centers_ch1[i] + centers_ch1[i + 1]) / 2 for i in range(len(centers_ch1) - 1)]
        kmeans_ch2 = KMeans(n_clusters=optimal_k_ch2, random_state=0).fit(combined_data_ch2)
        centers_ch2 = np.sort(kmeans_ch2.cluster_centers_[:, 0])
        kmeans_thresholds_ch2 = [(centers_ch2[i] + centers_ch2[i + 1]) / 2 for i in range(len(centers_ch2) - 1)]
    
        threshold_ch1 = kmeans_thresholds_ch1[0] if kmeans_thresholds_ch1 else None
        threshold_ch2 = kmeans_thresholds_ch2[0] if kmeans_thresholds_ch2 else None
    
        return threshold_ch1, threshold_ch2
    
    def kmeans_thresholding_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity, ch2_intensity, combined_data_ch1, combined_data_ch2 = self._prepare_data(bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        threshold_ch1, threshold_ch2 = self._compute_kmeans_thresholds(combined_data_ch1, combined_data_ch2)
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2)
        self.display_plot_in_new_window(fig, "K-means Thresholding")
    
        # Generate single polar plot for both channels above threshold
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        if threshold_ch1 is not None and threshold_ch2 is not None:
            self.generate_single_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2)
    
    def export_kmeans_thresholding_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity, ch2_intensity, combined_data_ch1, combined_data_ch2 = self._prepare_data(bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        threshold_ch1, threshold_ch2 = self._compute_kmeans_thresholds(combined_data_ch1, combined_data_ch2)
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2)
    
        filename = os.path.join(directory, f"Kmeans_Thresholding_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
    
        # Export single polar plot for both channels above threshold
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        if threshold_ch1 is not None and threshold_ch2 is not None:
            self.export_single_polar_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2, "Kmeans_Thresholding")
    
    def _compute_triangle_thresholds(self, ch1_intensity, ch2_intensity):
        triangle_threshold_ch1 = self.triangle_threshold(ch1_intensity)
        triangle_threshold_ch2 = self.triangle_threshold(ch2_intensity)
        return triangle_threshold_ch1, triangle_threshold_ch2
    
    def triangle_thresholding_plot(self, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
    
        threshold_ch1, threshold_ch2 = self._compute_triangle_thresholds(ch1_intensity, ch2_intensity)
    
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2)
        self.display_plot_in_new_window(fig, "Triangle Thresholding")
    
        # Generate single polar plot for both channels above threshold
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        self.generate_single_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2)
    
    def export_triangle_thresholding_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2):
        ch1_intensity = self.data[bkgd_corr_ch1].values
        ch2_intensity = self.data[bkgd_corr_ch2].values
    
        threshold_ch1, threshold_ch2 = self._compute_triangle_thresholds(ch1_intensity, ch2_intensity)
    
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2)
    
        filename = os.path.join(directory, f"Triangle_Thresholding_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
    
        # Export single polar plot for both channels above threshold
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        self.export_single_polar_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2, "Triangle_Thresholding")
    
    def _compute_gmm_thresholds(self, combined_data_ch1, combined_data_ch2):
        def find_optimal_clusters_gmm(data, max_clusters=10):
            bics = []
            for k in range(1, max_clusters + 1):
                gmm = GaussianMixture(n_components=k, random_state=0).fit(data)
                bics.append(gmm.bic(data))
            optimal_k = np.argmin(bics) + 1  # BIC is minimized for the optimal model
            return optimal_k
    
        optimal_k_ch1 = find_optimal_clusters_gmm(combined_data_ch1)
        optimal_k_ch2 = find_optimal_clusters_gmm(combined_data_ch2)
    
        gmm_ch1 = GaussianMixture(n_components=optimal_k_ch1, random_state=0).fit(combined_data_ch1)
        means_ch1 = np.sort(gmm_ch1.means_[:, 0])
        gmm_thresholds_ch1 = [(means_ch1[i] + means_ch1[i + 1]) / 2 for i in range(len(means_ch1) - 1)]
        threshold_ch1 = gmm_thresholds_ch1[0] if gmm_thresholds_ch1 else None
    
        gmm_ch2 = GaussianMixture(n_components=optimal_k_ch2, random_state=0).fit(combined_data_ch2)
        means_ch2 = np.sort(gmm_ch2.means_[:, 0])
        gmm_thresholds_ch2 = [(means_ch2[i] + means_ch2[i + 1]) / 2 for i in range(len(means_ch2) - 1)]
        threshold_ch2 = gmm_thresholds_ch2[0] if gmm_thresholds_ch2 else None
    
        return threshold_ch1, threshold_ch2
    
    def gmm_thresholding_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity, ch2_intensity, combined_data_ch1, combined_data_ch2 = self._prepare_data(bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        threshold_ch1, threshold_ch2 = self._compute_gmm_thresholds(combined_data_ch1, combined_data_ch2)
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2)
        self.display_plot_in_new_window(fig, "GMM Thresholding")
    
        # Generate single polar plot for both channels above threshold
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        if threshold_ch1 is not None and threshold_ch2 is not None:
            self.generate_single_polar_plot(bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2)
    
    def export_gmm_thresholding_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2):
        ch1_intensity, ch2_intensity, combined_data_ch1, combined_data_ch2 = self._prepare_data(bkgd_corr_ch1, bkgd_corr_ch2, bkgd_ch1, bkgd_ch2)
        threshold_ch1, threshold_ch2 = self._compute_gmm_thresholds(combined_data_ch1, combined_data_ch2)
        fig, threshold_ch1, threshold_ch2 = self._create_thresholding_plot(ch1_intensity, ch2_intensity, threshold_ch1, threshold_ch2)
    
        filename = os.path.join(directory, f"GMM_Thresholding_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get())
        plt.close(fig)
    
        # Export single polar plot for both channels above threshold
        cx_ch1 = self.column_vars[2].get()
        cy_ch1 = self.column_vars[3].get()
        cx_ch2 = self.column_vars[6].get()
        cy_ch2 = self.column_vars[7].get()
        if threshold_ch1 is not None and threshold_ch2 is not None:
            self.export_single_polar_plot(directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, threshold_ch1, threshold_ch2, "GMM_Thresholding")

    
    def rayleigh_test(self, angles):
        n = len(angles)
        if n < 2:
            return 1.0
        cos_mean = np.mean(np.cos(angles))
        sin_mean = np.mean(np.sin(angles))
        R = np.sqrt(cos_mean**2 + sin_mean**2)
        test_statistic = 2 * n * R**2
        p_value = chi2.sf(test_statistic, df=2)
        return p_value

    def determine_thresholds_based_on_rayleigh(self, vector_magnitudes, vector_angles, alpha=0.05):
        combined = list(zip(vector_magnitudes, vector_angles))
        sort_order = self.sorting_order.get()
        combined_sorted = sorted(combined, key=lambda x: x[0], reverse=(sort_order == 'descending'))
        if not combined_sorted:
            return 0
        threshold = None
        for i in range(1, len(combined_sorted) + 1):
            current_subset = combined_sorted[:i]
            current_angles = np.array([x[1] for x in current_subset])
            if len(current_angles) < 2:
                continue
            p_value = self.rayleigh_test_p_value(current_angles)
            adjusted_alpha = alpha / i  # Bonferroni correction
            if p_value < adjusted_alpha and threshold is None:
                threshold = current_subset[-1][0]
                break
        if threshold is None:
            if len(combined_sorted) > 0:
                threshold = combined_sorted[0][0]
            else:
                threshold = 0
        return threshold
    
    def determine_intensity_threshold(self, intensities, vector_angles, alpha, threshold_method='mean', mrl_threshold=0.3):
        """
        Determine an intensity threshold that separates vectors with significant directionality.
    
        Parameters:
        - intensities: array of intensity values
        - vector_angles: array of vector angles in radians
        - alpha: significance level for hypothesis testing
        - threshold_method: 'mean', 'median', or 'min'
        - mrl_threshold: minimum MRL value for above-threshold group
    
        Returns:
        - intensity_threshold: the determined threshold value
        - adjusted_p_value: the adjusted p-value for the threshold
        """
        # Combine intensities and angles into a list of tuples
        combined = list(zip(intensities, vector_angles))
    
        # Sort by intensity (ascending)
        sorted_data = sorted(combined, key=lambda x: x[0])
    
        # If no data, return default values
        if len(sorted_data) == 0:
            return 0, 1.0
    
        # Calculate overall statistic based on threshold method
        if threshold_method == 'mean':
            overall_stat = np.mean(intensities)
        elif threshold_method == 'median':
            overall_stat = np.median(intensities)
        else:  # min
            overall_stat = np.min(intensities)
    
        # If we have very few data points, just return the overall statistic
        if len(sorted_data) <= 20:
            return overall_stat, 1.0
    
        # Initialize variables to track the best threshold
        best_threshold = overall_stat
        best_score = 0
        best_p_above = 1.0
        best_mrl_above = 0
    
        # We'll look at potential thresholds at percentiles across the range
        percentiles = np.linspace(0, 100, min(100, len(sorted_data))) if len(sorted_data) > 0 else []
        potential_thresholds = []
    
        for p in percentiles:
            if len(sorted_data) > 0:
                idx = int(p/100 * (len(sorted_data)-1))
                potential_thresholds.append(sorted_data[idx][0])
    
        potential_thresholds = sorted(list(set(potential_thresholds)))
    
        # For each potential threshold
        for threshold in potential_thresholds:
            # Split into below and above threshold
            below_threshold_angles = []
            above_threshold_angles = []
    
            for intensity, angle in sorted_data:
                if intensity <= threshold:
                    below_threshold_angles.append(angle)
                else:
                    above_threshold_angles.append(angle)
    
            # Calculate Rayleigh test p-values and MRL for both groups
            p_below = 1.0
            mrl_below = 0
            if len(below_threshold_angles) > 1:
                _, p_below, mrl_below = self.calculate_rayleigh_p_value(below_threshold_angles)
    
            p_above = 1.0
            mrl_above = 0
            if len(above_threshold_angles) > 1:
                _, p_above, mrl_above = self.calculate_rayleigh_p_value(above_threshold_angles)
    
            # We want to find where above becomes significant and below is not
            # and where the above group has stronger directionality and meets MRL threshold
            if p_above < alpha and (p_below >= alpha or len(below_threshold_angles) <= 1) and mrl_above >= mrl_threshold:
                # Calculate a score that considers both p-value difference and MRL
                current_score = (alpha - p_above) * mrl_above  # Higher when p_above is smaller and MRL is larger
    
                if current_score > best_score:
                    best_score = current_score
                    best_threshold = threshold
                    best_p_above = p_above
                    best_mrl_above = mrl_above
    
        # If we didn't find a suitable threshold, use the overall statistic
        if best_score == 0:
            best_threshold = overall_stat
            best_p_above = 1.0
    
        # Calculate adjusted p-value using Bonferroni correction
        num_tests = len(potential_thresholds) if potential_thresholds else 1
        adjusted_p_value = min(best_p_above * num_tests, 1.0)
    
        return best_threshold, adjusted_p_value

    def collect_vector_analysis(self, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        # Get the alpha value from the UI
        alpha = self.alpha_value.get()
    
        centroids_ch1_x = self.data[cx_ch1]
        centroids_ch1_y = self.data[cy_ch1]
        centroids_ch2_x = self.data[cx_ch2]
        centroids_ch2_y = self.data[cy_ch2]
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_angles = np.arctan2(vectors_y, vectors_x)
        vector_angles_deg = np.degrees(vector_angles)
        vector_magnitudes = np.sqrt(vectors_x**2 + vectors_y**2)
    
        # Use the alpha value from the UI in the threshold determination
        threshold = self.determine_thresholds_based_on_rayleigh(vector_magnitudes, vector_angles, alpha=alpha)
    
        below_threshold_angles = vector_angles[vector_magnitudes <= threshold]
        above_threshold_angles = vector_angles[vector_magnitudes > threshold]
    
        mean_angle = circmean(vector_angles)
        std_angle = circstd(vector_angles)
        R = np.sqrt((np.sum(np.cos(vector_angles)))**2 + (np.sum(np.sin(vector_angles)))**2)
        n = len(vector_angles)
        rayleigh_statistic = R / n if n > 0 else 0
        z = 2 * n * rayleigh_statistic**2 if n > 0 else 0
        rayleigh_p_value = self.rayleigh_test_p_value(vector_angles) if n > 0 else 1.0
        mean_angle_deg = np.degrees(mean_angle)
        std_angle_deg = np.degrees(std_angle)
        mrl = self.calculate_mrl(vector_angles)
    
        def calculate_stats(angles_subset):
            if len(angles_subset) > 0:
                mean_angle_subset = circmean(angles_subset)
                std_angle_subset = circstd(angles_subset)
                R_subset = np.sqrt((np.sum(np.cos(angles_subset)))**2 + (np.sum(np.sin(angles_subset)))**2)
                n_subset = len(angles_subset)
                r_subset = R_subset / n_subset if n_subset > 0 else 0
                z_subset = 2 * n_subset * r_subset**2 if n_subset > 0 else 0
                rayleigh_p_value_subset = self.rayleigh_test_p_value(angles_subset) if n_subset > 0 else 1.0
                mrl_subset = self.calculate_mrl(angles_subset)
                mean_angle_deg_subset = np.degrees(mean_angle_subset)
                std_angle_deg_subset = np.degrees(std_angle_subset)
                return {
                    "mean_angle_deg": mean_angle_deg_subset,
                    "std_angle_deg": std_angle_deg_subset,
                    "rayleigh_p_value": rayleigh_p_value_subset,
                    "mrl": mrl_subset
                }
            else:
                return {
                    "mean_angle_deg": None,
                    "std_angle_deg": None,
                    "rayleigh_p_value": None,
                    "mrl": None
                }
    
        below_stats = calculate_stats(below_threshold_angles)
        above_stats = calculate_stats(above_threshold_angles)
    
        # Apply Bonferroni correction for multiple testing
        corrected_alpha = alpha / 3  # Divide by 3 because we have three tests
    
        # Determine significance based on corrected alpha and MRL
        status_all = "Significant" if rayleigh_p_value < corrected_alpha and mrl > 0.7 else "Not Significant"
        status_below = "Significant" if below_stats["rayleigh_p_value"] is not None and below_stats["rayleigh_p_value"] < corrected_alpha and below_stats["mrl"] > 0.7 else "Not Significant"
        status_above = "Significant" if above_stats["rayleigh_p_value"] is not None and above_stats["rayleigh_p_value"] < corrected_alpha and above_stats["mrl"] > 0.7 else "Not Significant"
    
        self.analysis_results.append({
            "Analysis Type": "Vector Analysis",
            "Mean Angle (radians)": mean_angle,
            "Standard Deviation of Angles (radians)": std_angle,
            "Mean Angle (degrees)": mean_angle_deg,
            "Standard Deviation of Angles (degrees)": std_angle_deg,
            "Angle Range (degrees)": f"{vector_angles_deg.min():.2f} to {vector_angles_deg.max():.2f}",
            "Rayleigh Test p-value": rayleigh_p_value,
            "Mean Resultant Length (MRL)": mrl,
            "Rayleigh Test Status": status_all,
            "Alpha Value Used": alpha,
            "Threshold for Magnitudes": threshold,
            "Below Threshold Count": len(below_threshold_angles),
            "Above Threshold Count": len(above_threshold_angles),
            "Below Threshold Mean Angle (degrees)": below_stats["mean_angle_deg"],
            "Below Threshold Std Angle (degrees)": below_stats["std_angle_deg"],
            "Below Threshold Rayleigh p-value": below_stats["rayleigh_p_value"],
            "Below Threshold MRL": below_stats["mrl"],
            "Below Threshold Rayleigh Status": status_below,
            "Above Threshold Mean Angle (degrees)": above_stats["mean_angle_deg"],
            "Above Threshold Std Angle (degrees)": above_stats["std_angle_deg"],
            "Above Threshold Rayleigh p-value": above_stats["rayleigh_p_value"],
            "Above Threshold MRL": above_stats["mrl"],
            "Above Threshold Rayleigh Status": status_above
        })


    
    def calculate_mrl(self, data):
        """Calculate the Mean Resultant Length (MRL) for circular data."""
        if np.max(data) > 2*np.pi:
            angles_rad = np.deg2rad(data)
        else:
            angles_rad = data
        sum_cos = np.sum(np.cos(angles_rad))
        sum_sin = np.sum(np.sin(angles_rad))
        resultant_length = np.sqrt(sum_cos**2 + sum_sin**2)
        n = len(data)
        mrl = resultant_length / n if n > 0 else 0
        return mrl
    
    def rayleigh_test_p_value(self, data):
        """Calculate the p-value for the Rayleigh test of uniformity."""
        if np.max(data) > 2*np.pi:
            angles_rad = np.deg2rad(data)
        else:
            angles_rad = data
        n = len(angles_rad)
        if n == 0:
            return 1.0
        sum_cos = np.sum(np.cos(angles_rad))
        sum_sin = np.sum(np.sin(angles_rad))
        R = np.sqrt(sum_cos**2 + sum_sin**2)
        R_bar = R / n
        p_value = np.exp(np.sqrt(1 + 4*n + 4*n**2*(1-R_bar)) - (1 + 2*n))
        return p_value
    
    
    def calculate_rayleigh_p_value(self, angles):
        if len(angles) == 0:
            return None, None, None  # p_value, mean_angle, mrl
        p_value = self.rayleigh_test_p_value(angles)
        mean_angle = circmean(angles)
        mrl_value = self.calculate_mrl(angles)
        return p_value, np.degrees(mean_angle), mrl_value

    def _calculate_vector_data(self, cx_ch1, cy_ch1, cx_ch2, cy_ch2, alpha=None):
        if alpha is None:
            alpha = self.alpha_value.get()
    
        centroids_ch1_x = self.data[cx_ch1]
        centroids_ch1_y = self.data[cy_ch1]
        centroids_ch2_x = self.data[cx_ch2]
        centroids_ch2_y = self.data[cy_ch2]
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_angles = np.arctan2(vectors_y, vectors_x)
        vector_magnitudes = np.sqrt(vectors_x**2 + vectors_y**2)
        threshold = self.determine_thresholds_based_on_rayleigh(vector_magnitudes, vector_angles, alpha=alpha)
        below_threshold_mask = vector_magnitudes <= threshold
        above_threshold_mask = vector_magnitudes > threshold
        return {
            'vector_angles': vector_angles,
            'vector_magnitudes': vector_magnitudes,
            'threshold': threshold,
            'below_threshold_mask': below_threshold_mask,
            'above_threshold_mask': above_threshold_mask,
            'vectors_x': vectors_x,
            'vectors_y': vectors_y,
        }
    
    def _calculate_statistics(self, vector_angles, below_threshold_mask, above_threshold_mask):
        stats = {}
    
        # All vectors
        rayleigh_p_value_all, mean_angle_all, mrl_all = self.calculate_rayleigh_p_value(vector_angles)
        stats['all'] = {
            'rayleigh_p_value_rounded': f"{rayleigh_p_value_all:.5f}" if rayleigh_p_value_all is not None and rayleigh_p_value_all >= 0.0001 else "p < 0.0001" if rayleigh_p_value_all is not None else "N/A",
            'mrl_rounded': f"{mrl_all:.3f}" if mrl_all is not None else "N/A",
            'directionality_status': "Significant" if (rayleigh_p_value_all is not None and rayleigh_p_value_all < 0.05 and mrl_all is not None and mrl_all > 0.7) else "Not Significant",
            'rayleigh_p_value': rayleigh_p_value_all,
            'mean_angle': mean_angle_all,
            'mrl': mrl_all,
        }
    
        # Below threshold
        if np.sum(below_threshold_mask) > 0:
            below_threshold_angles = vector_angles[below_threshold_mask]
            rayleigh_p_value_below, mean_angle_below, mrl_below = self.calculate_rayleigh_p_value(below_threshold_angles)
            stats['below_threshold'] = {
                'rayleigh_p_value': rayleigh_p_value_below,
                'mean_angle': mean_angle_below,
                'mrl': mrl_below,
                'count': np.sum(below_threshold_mask),
            }
        else:
            stats['below_threshold'] = {
                'count': 0,
            }
    
        # Above threshold
        if np.sum(above_threshold_mask) > 0:
            above_threshold_angles = vector_angles[above_threshold_mask]
            rayleigh_p_value_above, mean_angle_above, mrl_above = self.calculate_rayleigh_p_value(above_threshold_angles)
            stats['above_threshold'] = {
                'rayleigh_p_value': rayleigh_p_value_above,
                'mean_angle': mean_angle_above,
                'mrl': mrl_above,
                'count': np.sum(above_threshold_mask),
            }
        else:
            stats['above_threshold'] = {
                'count': 0,
            }
    
        return stats
    

    def _create_main_figure(self, vector_data, stats):
        fig = plt.figure(figsize=(18, 12))
    
        # 1. All vectors polar plot
        ax1 = plt.subplot(2, 2, 1, polar=True)
        n_bins = 36
        bin_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        counts, _ = np.histogram(vector_data['vector_angles'], bins=bin_edges)
        bars = ax1.bar(bin_centers, counts, width=(2*np.pi/n_bins), bottom=0.0,
                      color='green', alpha=0.5, edgecolor='k', linewidth=0.5)
    
        # Calculate corrected p-value for all vectors
        if stats["all"]["rayleigh_p_value"] is not None:
            corrected_p_value_all = min(stats["all"]["rayleigh_p_value"] * 3, 1.0)
            formatted_p_value_all = f"{corrected_p_value_all:.5f}" if corrected_p_value_all >= 0.0001 else "p < 0.0001"
        else:
            formatted_p_value_all = "N/A"
        p_value_text_all = (
            f'Rayleigh corrected p-value: {formatted_p_value_all}\n'
            f'MRL: {stats["all"]["mrl_rounded"]}\n'
            f'Directionality: {stats["all"]["directionality_status"]}'
        )
        ax1.text(0.5, 0.1, p_value_text_all, transform=ax1.transAxes, fontsize=10,
                horizontalalignment='center', verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
        if stats["all"]["rayleigh_p_value"] is not None and stats["all"]["rayleigh_p_value"] < 0.05 and stats["all"]["mrl"] is not None and stats["all"]["mrl"] > 0.7 and stats["all"]["mean_angle"] is not None:
            ax1.plot([0, stats["all"]["mean_angle"]], [0, np.max(counts)], 'r-', linewidth=2)
            ax1.text(0.5, 0.5, f'Mean Angle: {stats["all"]["mean_angle"]:.2f}°', transform=ax1.transAxes, fontsize=10, horizontalalignment='center',
                    verticalalignment='center', color='red', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
        ax1.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
        ax1.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
        ax1.set_rlabel_position(30)
        ax1.grid(True)
        ax1.set_title('All Vectors', y=1.1)
    
        # 2. Below threshold vectors polar plot
        ax2 = plt.subplot(2, 2, 2, polar=True)
        if stats['below_threshold']['count'] > 0:
            below_threshold_angles = vector_data['vector_angles'][vector_data['below_threshold_mask']]
            counts_below, _ = np.histogram(below_threshold_angles, bins=bin_edges)
            bars_below = ax2.bar(bin_centers, counts_below, width=(2*np.pi/n_bins), bottom=0.0,
                                color='blue', alpha=0.5, edgecolor='k', linewidth=0.5)
    
            # Calculate corrected p-value for below threshold
            if stats['below_threshold']['rayleigh_p_value'] is not None:
                corrected_p_value_below = min(stats['below_threshold']['rayleigh_p_value'] * 3, 1.0)
                formatted_p_value_below = f"{corrected_p_value_below:.5f}" if corrected_p_value_below >= 0.0001 else "p < 0.0001"
            else:
                formatted_p_value_below = "N/A"
            p_value_text_below = (
                f'Rayleigh corrected p-value: {formatted_p_value_below}\n'
                f'MRL: {stats["below_threshold"]["mrl"]:.3f}\n'
                f'Directionality: {"Significant" if stats["below_threshold"]["mrl"] and stats["below_threshold"]["mrl"] > 0.7 else "Not Significant"}'
            )
            ax2.text(0.5, 0.1, p_value_text_below, transform=ax2.transAxes, fontsize=10,
                    horizontalalignment='center', verticalalignment='bottom',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
            ax2.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
            ax2.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
            ax2.set_rlabel_position(30)
            ax2.grid(True)
        else:
            ax2.set_title('No vectors below threshold', y=1.1)
        ax2.set_title(f'Below Threshold ({stats["below_threshold"]["count"]} vectors)', y=1.1)
    
        # 3. Above threshold vectors polar plot
        ax3 = plt.subplot(2, 2, 3, polar=True)
        if stats['above_threshold']['count'] > 0:
            above_threshold_angles = vector_data['vector_angles'][vector_data['above_threshold_mask']]
            counts_above, _ = np.histogram(above_threshold_angles, bins=bin_edges)
            bars_above = ax3.bar(bin_centers, counts_above, width=(2*np.pi/n_bins), bottom=0.0,
                                color='red', alpha=0.5, edgecolor='k', linewidth=0.5)
    
            # Calculate corrected p-value for above threshold
            if stats['above_threshold']['rayleigh_p_value'] is not None:
                corrected_p_value_above = min(stats['above_threshold']['rayleigh_p_value'] * 3, 1.0)
                formatted_p_value_above = f"{corrected_p_value_above:.5f}" if corrected_p_value_above >= 0.0001 else "p < 0.0001"
            else:
                formatted_p_value_above = "N/A"
            p_value_text_above = (
                f'Rayleigh corrected p-value: {formatted_p_value_above}\n'
                f'MRL: {stats["above_threshold"]["mrl"]:.3f}\n'
                f'Directionality: {"Significant" if stats["above_threshold"]["mrl"] and stats["above_threshold"]["mrl"] > 0.7 else "Not Significant"}'
            )
            ax3.text(0.5, 0.1, p_value_text_above, transform=ax3.transAxes, fontsize=10,
                    horizontalalignment='center', verticalalignment='bottom',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
            ax3.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
            ax3.set_xticklabels(['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'])
            ax3.set_rlabel_position(30)
            ax3.grid(True)
        else:
            ax3.set_title('No vectors above threshold', y=1.1)
        ax3.set_title(f'Above Threshold ({stats["above_threshold"]["count"]} vectors)', y=1.1)
    
        # 4. Rayleigh p-value progression plot with improved scaling
        ax4 = plt.subplot(2, 2, 4)
        combined = list(zip(vector_data['vector_magnitudes'], vector_data['vector_angles']))
        combined_sorted_descending = sorted(combined, key=lambda x: x[0], reverse=True)
    
        cumulative_p_values = []
        alpha = self.alpha_value.get()  # Use the alpha value from the UI
    
        for i in range(1, len(combined_sorted_descending) + 1):
            current_subset = combined_sorted_descending[:i]
            current_angles = np.array([x[1] for x in current_subset])
            p_value = self.rayleigh_test_p_value(current_angles)
            # Apply Bonferroni correction
            corrected_p_value = min(p_value * i, 1.0)  # Multiply by i since we're doing i tests
            cumulative_p_values.append((i, corrected_p_value))
    
        if len(cumulative_p_values) > 0:
            x_vals, y_vals = zip(*cumulative_p_values)
    
            # Plot the p-value progression
            ax4.plot(x_vals, y_vals, 'b-', linewidth=2)
    
            # Add significance threshold line
            ax4.axhline(y=alpha, color='r', linestyle='--', label=f'Significance Threshold ({alpha})')
    
            # Set y-axis to logarithmic scale since p-values span orders of magnitude
            ax4.set_yscale('log')
    
            # Calculate reasonable y-axis limits
            min_p = min(y_vals) if y_vals else 1e-10
            max_p = max(y_vals) if y_vals else 1.0
    
            # Set y-axis limits with padding
            ax4.set_ylim(max(1e-10, min_p/10), min(1.0, max_p*10))
    
            # Set x-axis limits with padding (from 0 to max + 5%)
            max_x = max(x_vals) if x_vals else 0
            ax4.set_xlim(0, max_x * 1.05 if max_x > 0 else 10)
    
            # Add minor ticks for better readability on log scale
            ax4.yaxis.set_minor_locator(plt.LogLocator(subs=(0.2, 0.4, 0.6, 0.8)))
            ax4.yaxis.set_minor_formatter(plt.NullFormatter())  # Hide minor tick labels
    
            # Mark the threshold point if it exists
            if vector_data['threshold'] > 0:
                # Find the index where we set the threshold
                threshold_idx = None
                for i, (count, p_val) in enumerate(cumulative_p_values):
                    if p_val < alpha and (i == 0 or cumulative_p_values[i-1][1] >= alpha):
                        threshold_idx = i
                        break
    
                if threshold_idx is not None:
                    ax4.scatter(x_vals[threshold_idx], y_vals[threshold_idx],
                               color='red', s=100, label=f'Threshold at vector {x_vals[threshold_idx]}')
                    ax4.annotate(f'Threshold\n({x_vals[threshold_idx]} vectors)',
                                xy=(x_vals[threshold_idx], y_vals[threshold_idx]),
                                xytext=(x_vals[threshold_idx]+5, y_vals[threshold_idx]*1.5),
                                arrowprops=dict(facecolor='black', shrink=0.05))
    
            ax4.set_xlabel('Number of vectors (sorted by magnitude in descending order)')
            ax4.set_ylabel('Corrected Rayleigh p-value (log scale)')
            ax4.set_title('Rayleigh p-value Progression with Bonferroni Correction')
            ax4.legend()
            ax4.grid(True, which="both", ls="-", alpha=0.5)  # Grid for both major and minor ticks
        else:
            # If there's no data, set some reasonable defaults
            ax4.set_xlim(0, 10)
            ax4.set_ylim(1e-10, 1.05)
            ax4.set_yscale('log')
            ax4.set_xlabel('Number of vectors (sorted by magnitude in descending order)')
            ax4.set_ylabel('Corrected Rayleigh p-value (log scale)')
            ax4.set_title('Rayleigh p-value Progression with Bonferroni Correction')
            ax4.grid(True)
            ax4.text(0.5, 0.5, 'No data available',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax4.transAxes)
    
        plt.tight_layout()
        return fig


    def _create_magnitude_figure(self, vector_magnitudes, threshold):
        fig_magnitudes = plt.figure(figsize=(12, 6))
        ax4_cartesian = fig_magnitudes.add_subplot(111)
        if len(vector_magnitudes) > 0:
            sns.histplot(vector_magnitudes, bins=50, stat='density', color='blue', alpha=0.3, ax=ax4_cartesian)
            ax4_cartesian.axvline(threshold, color='purple', linestyle=':', linewidth=2, label=f'Threshold: {threshold:.2f}')
            ax4_cartesian.set_xlabel('Distance')
            ax4_cartesian.set_ylabel('Density')
            ax4_cartesian.set_title('Distribution of Vector Magnitudes with Threshold')
            ax4_cartesian.legend()
        return fig_magnitudes
    
    def vector_analysis_plot(self, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        alpha = self.alpha_value.get()
        vector_data = self._calculate_vector_data(cx_ch1, cy_ch1, cx_ch2, cy_ch2, alpha=alpha)
        stats = self._calculate_statistics(vector_data['vector_angles'], vector_data['below_threshold_mask'], vector_data['above_threshold_mask'])
    
        fig = self._create_main_figure(vector_data, stats)
        fig_magnitudes = self._create_magnitude_figure(vector_data['vector_magnitudes'], vector_data['threshold'])
    
        self.display_plot_in_new_window(fig, "Vector Analysis - Polar Plots")
        self.display_plot_in_new_window(fig_magnitudes, "Vector Magnitudes with Threshold")
        plt.close(fig)
        plt.close(fig_magnitudes)

    def export_vector_analysis_plot(self, directory, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        vector_data = self._calculate_vector_data(cx_ch1, cy_ch1, cx_ch2, cy_ch2)
        stats = self._calculate_statistics(vector_data['vector_angles'], vector_data['below_threshold_mask'], vector_data['above_threshold_mask'])
    
        fig = self._create_main_figure(vector_data, stats)
        filename_polar = os.path.join(directory, f"Vector_Analysis_Polar_{cx_ch1}_to_{cx_ch2}.{self.export_format.get()}")
        fig.savefig(filename_polar, format=self.export_format.get())
        plt.close(fig)
    
        fig_magnitudes = self._create_magnitude_figure(vector_data['vector_magnitudes'], vector_data['threshold'])
        filename_magnitudes = os.path.join(directory, f"Vector_Magnitudes_with_Threshold_{cx_ch1}_to_{cx_ch2}.{self.export_format.get()}")
        fig_magnitudes.savefig(filename_magnitudes, format=self.export_format.get())
        plt.close(fig_magnitudes)

    def collect_vector_magnitude_intensity_threshold(self, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2, alpha=None, threshold_method=None, mrl_threshold=None):
        """
        Calculate intensity thresholds for channel 1 and channel 2 based on vector magnitudes.
    
        Args:
            bkgd_corr_ch1: Background corrected intensity column for channel 1
            bkgd_corr_ch2: Background corrected intensity column for channel 2
            cx_ch1: Centroid X column for channel 1
            cy_ch1: Centroid Y column for channel 1
            cx_ch2: Centroid X column for channel 2
            cy_ch2: Centroid Y column for channel 2
            alpha: Significance level for Rayleigh test
            threshold_method: Method for determining threshold ('min', 'mean', 'median')
            mrl_threshold: Threshold for mean resultant length
    
        Returns:
            Tuple of (channel 1 intensity threshold, channel 2 intensity threshold,
                     channel 1 adjusted p-value, channel 2 adjusted p-value)
        """
        # Get parameters from UI if not provided
        if alpha is None:
            alpha = self.vector_magnitude_alpha.get()  # Use the dedicated alpha value
        if threshold_method is None:
            threshold_method = self.threshold_method.get()
        if mrl_threshold is None:
            mrl_threshold = self.vector_magnitude_mrl_threshold.get()  # Use the dedicated MRL threshold
    
        # Calculate vector angles
        centroids_ch1_x = self.data[cx_ch1]
        centroids_ch1_y = self.data[cy_ch1]
        centroids_ch2_x = self.data[cx_ch2]
        centroids_ch2_y = self.data[cy_ch2]
    
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_angles = np.arctan2(vectors_y, vectors_x)
    
        # Get intensity values
        ch1_intensities = self.data[bkgd_corr_ch1].values
        ch2_intensities = self.data[bkgd_corr_ch2].values
    
        # Determine thresholds for each channel
        ch1_intensity_threshold, ch1_adjusted_p_value = self.determine_intensity_threshold(
            ch1_intensities, vector_angles, alpha,
            threshold_method, mrl_threshold
        )
        ch2_intensity_threshold, ch2_adjusted_p_value = self.determine_intensity_threshold(
            ch2_intensities, vector_angles, alpha,
            threshold_method, mrl_threshold
        )
    
        # Store results
        self.analysis_results.append({
            "Analysis Type": "Vector Magnitude Intensity Threshold",
            "Channel 1 Intensity Threshold": ch1_intensity_threshold,
            "Channel 2 Intensity Threshold": ch2_intensity_threshold,
            "Channel 1 Adjusted p-value": ch1_adjusted_p_value,
            "Channel 2 Adjusted p-value": ch2_adjusted_p_value,
            "Alpha": alpha,
            "Threshold Method": threshold_method,
            "MRL Threshold": mrl_threshold
        })
    
        return ch1_intensity_threshold, ch2_intensity_threshold, ch1_adjusted_p_value, ch2_adjusted_p_value
    
    def create_vector_mag_intensity_threshold_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        # Calculate thresholds and get additional information
        ch1_intensity_threshold, ch2_intensity_threshold, ch1_p_value, ch2_p_value = self.collect_vector_magnitude_intensity_threshold(
            bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2,
            alpha=self.vector_magnitude_alpha.get(),
            threshold_method=self.threshold_method.get(),
            mrl_threshold=self.vector_magnitude_mrl_threshold.get()
        )
    
        # Calculate vector magnitudes and angles for additional visualization
        centroids_ch1_x = self.data[cx_ch1]
        centroids_ch1_y = self.data[cy_ch1]
        centroids_ch2_x = self.data[cx_ch2]
        centroids_ch2_y = self.data[cy_ch2]
        vectors_x = centroids_ch2_x - centroids_ch1_x
        vectors_y = centroids_ch2_y - centroids_ch1_y
        vector_magnitudes = np.sqrt(vectors_x**2 + vectors_y**2)
        vector_angles = np.arctan2(vectors_y, vectors_x)
    
        # Determine a magnitude threshold
        magnitude_threshold, _ = self.determine_intensity_threshold(
            vector_magnitudes,
            vector_angles,
            alpha=self.vector_magnitude_alpha.get(),
            threshold_method=self.threshold_method.get(),
            mrl_threshold=self.vector_magnitude_mrl_threshold.get()
        )
    
        # Create the figure with four subplots (2 rows, 2 columns)
        fig, ax = plt.subplots(2, 2, figsize=(18, 12))
    
        # Plot Channel 1 intensity histogram (top left)
        ax[0, 0].hist(self.data[bkgd_corr_ch1], bins=256, alpha=0.6, color='blue')
        ax[0, 0].axvline(ch1_intensity_threshold, color='red', linestyle='dashed', linewidth=1)
        text_content = (
            f'Threshold: {ch1_intensity_threshold:.2f}\n'
            f'Adjusted p-value: {ch1_p_value:.4f}'
        )
        self.add_plot_text(ax[0, 0], text_content)
        ch1_name = self.channel1_var.get()
        ax[0, 0].set_title(f'{ch1_name} Intensity')
        ax[0, 0].set_xlabel('Intensity')
        ax[0, 0].set_ylabel('Frequency')
    
        # Plot Channel 2 intensity histogram (top right)
        ax[0, 1].hist(self.data[bkgd_corr_ch2], bins=256, alpha=0.6, color='green')
        ax[0, 1].axvline(ch2_intensity_threshold, color='red', linestyle='dashed', linewidth=1)
        text_content = (
            f'Threshold: {ch2_intensity_threshold:.2f}\n'
            f'Adjusted p-value: {ch2_p_value:.4f}'
        )
        self.add_plot_text(ax[0, 1], text_content)
        ch2_name = self.channel2_var.get()
        ax[0, 1].set_title(f'{ch2_name} Intensity')
        ax[0, 1].set_xlabel('Intensity')
        ax[0, 1].set_ylabel('Frequency')
    
        # Plot vector magnitude distribution (bottom left)
        ax[1, 0].hist(vector_magnitudes, bins=100, alpha=0.7, color='purple')
        ax[1, 0].axvline(magnitude_threshold, color='orange', linestyle='dashed', linewidth=2,
                        label=f'Magnitude Threshold: {magnitude_threshold:.2f}')
        ax[1, 0].set_xlabel('Vector Magnitude')
        ax[1, 0].set_ylabel('Frequency')
        ax[1, 0].set_title('Distribution of Vector Magnitudes')
        ax[1, 0].legend()
    
        # Create scatter plot of intensities with symlog scale (bottom right)
        ch1_intensity = self.data[bkgd_corr_ch1]
        ch2_intensity = self.data[bkgd_corr_ch2]
        scatter = ax[1, 1].scatter(ch1_intensity, ch2_intensity, alpha=0.3)
        self._setup_symlog_scale(ax[1, 1], ch1_intensity, ch2_intensity)
    
        # Add threshold lines to scatter plot
        ax[1, 1].axvline(ch1_intensity_threshold, color='red', linestyle='dashed', linewidth=1,
                        label=f'Channel 1 Threshold: {ch1_intensity_threshold:.2f}')
        ax[1, 1].axhline(ch2_intensity_threshold, color='green', linestyle='dashed', linewidth=1,
                        label=f'Channel 2 Threshold: {ch2_intensity_threshold:.2f}')
    
        ax[1, 1].set_xlabel(f'{ch1_name} Intensity')
        ax[1, 1].set_ylabel(f'{ch2_name} Intensity')
        ax[1, 1].set_title('Channel Intensities Scatter Plot (Symlog Scale)')
        ax[1, 1].legend()
    
        # Add a main title with some additional information
        fig.suptitle(
            f'Vector Magnitude Intensity Threshold Analysis\n'
            f'Threshold Method: {self.threshold_method.get()}, '
            f'MRL Threshold: {self.vector_magnitude_mrl_threshold.get():.2f}, '
            f'Alpha: {self.vector_magnitude_alpha.get():.2f}',
            y=1.02
        )
    
        plt.tight_layout()
        return fig, ch1_intensity_threshold, ch2_intensity_threshold, magnitude_threshold

    
    def vector_magnitude_intensity_threshold_plot(self, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        fig, ch1_threshold, ch2_threshold, mag_threshold = self.create_vector_mag_intensity_threshold_plot(
            bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2
        )
        self.display_plot_in_new_window(fig, "Vector Magnitude Intensity Threshold Analysis")
        plt.close(fig)
    
    def export_vector_magnitude_intensity_threshold_plot(self, directory, bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2):
        fig, ch1_threshold, ch2_threshold, mag_threshold = self.create_vector_mag_intensity_threshold_plot(
            bkgd_corr_ch1, bkgd_corr_ch2, cx_ch1, cy_ch1, cx_ch2, cy_ch2
        )
        filename = os.path.join(directory, f"Vector_Magnitude_Intensity_Threshold_{bkgd_corr_ch1}_vs_{bkgd_corr_ch2}.{self.export_format.get()}")
        fig.savefig(filename, format=self.export_format.get(), bbox_inches='tight', dpi=300)
        plt.close(fig)
        return filename


if __name__ == "__main__":
    root = tk.Tk()
    app = ThresholdingApp(root)
    root.mainloop()
