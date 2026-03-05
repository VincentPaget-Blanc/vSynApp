# Changelog

All notable changes to vSynApp are documented here.
Versions follow the format `V<major>.<minor>` aligned with Zenodo releases.

---

## [V0.1] — 2025

This version consolidates bug fixes, statistical corrections, visualization improvements, and a new multi-file search feature relative to the initial V0 release (0.094051).

### Bug Fixes

**Duplicate method definition removed**
`toggle_tsne_knn_analysis_widgets` was defined twice in the class body. Python silently discards the first definition, making its logic unreachable. The duplicate has been removed, leaving only the correct implementation.

**`export_tsne_clusters_plot` name mismatch fixed**
`export_plots` called `self.export_tsne_clusters_plot(...)`, but the method was defined as `export_plot_tsne_clusters`. This caused an `AttributeError` at runtime whenever t-SNE cluster plots were exported. The method has been renamed to match its call site.

**Array-in-filename crash fixed**
The t-SNE cluster export filename previously embedded the full `X_tsne` and `labels_tsne` NumPy arrays directly into the file path string, generating an unusable path that crashed `fig.savefig`. The filename now uses the cluster count `k` instead (e.g. `t-SNE_Clusters_k3.png`).

**Double-display and crash for Vector Magnitude plot fixed**
`vector_magnitude_intensity_threshold_plot` handles its own display and `plt.close` internally and returns `None`. The caller in `generate_plots` was then trying to pass that `None` to `display_plot_in_new_window` and `plt.close`, causing the window to appear twice and then crash. The redundant calls have been removed.

**`k=1` division-by-zero and axes array crash fixed**
All three cluster plot methods (`_create_tsne_cluster_plot`, `_create_original_cluster_plot`, `_create_individual_clusters_plot`) used `i / (k - 1)` for colormap indexing, which divides by zero when `k=1`. Additionally, `plt.subplots(1, k)` with `k=1` returns a single `Axes` object rather than an array, causing an `IndexError` on `axes[i]`. All three methods now clamp `k` to a minimum of 2, and `_create_individual_clusters_plot` uses `squeeze=False` to ensure `axes` is always indexable.

**Polar plot "empty" titles silently overwritten fixed**
In both `_generate_vector_analysis_plots` and `_create_main_figure`, the "No vectors below/above threshold" warning title was set inside the `else` branch and then unconditionally overwritten by the count-bearing title on the next line. Empty panels now correctly preserve their warning message; the count title is only applied when the count is greater than zero.

---

### Statistical Corrections

**`rayleigh_test_p_value` replaced with exact implementation**
The previous implementation used Mardia's approximation:
```python
p_value = np.exp(np.sqrt(1 + 4*n + 4*n**2*(1-R_bar)) - (1 + 2*n))
```
This approximation can return values outside `[0, 1]` for small `n` or extreme `R_bar` values and diverges from the true p-value at low sample sizes. It has been replaced with the standard exact form:
```python
p_value = float(chi2.sf(2 * n * R_bar**2, df=2))
```
The result is clipped to `[0, 1]`. The now-redundant `rayleigh_test` method (which already used `chi2.sf` correctly but was barely called) delegates to this unified implementation for backward compatibility.

**Hardcoded MRL threshold of `0.7` eliminated**
Six locations across `_generate_polar_plot`, `collect_vector_analysis`, `_calculate_statistics`, and `_create_main_figure` were ignoring the user-configurable MRL threshold and using a hardcoded value of `0.7`. All now read from `self.mrl_threshold.get()`, so the threshold set in the UI is consistently honoured throughout all analysis pathways.

**Sort order inconsistency between MRL threshold and Rayleigh p-value progression plot fixed**
`_generate_vector_analysis_plots` correctly used the UI sorting order for MRL threshold calculation, but `_create_main_figure` hardcoded descending order for the Rayleigh p-value progression plot. The threshold marker drawn on the plot therefore did not correspond to the threshold computed in the analysis. Both now use the same UI-controlled sorting order.

---

### Visualization Improvements

**`add_plot_text` anchor corrected**
The helper used `verticalalignment='top'` with a default anchor position of `(0.05, 0.05)` (bottom-left corner of the axes). This placed the top edge of the text box at 5% height, clipping most of the text below the visible plot area. Changed to `verticalalignment='bottom'` to match the anchor position.

**Intensity histogram axis labels added**
Both per-channel intensity histograms in `_create_thresholding_plot` were missing x-axis and y-axis labels. They now display "Background-Corrected Intensity" and "Count" respectively.

**Background comparison plot upgraded**
The background comparison was previously a bare boxplot, which is poorly suited for the right-skewed distributions typical of immunofluorescence intensity data. It has been replaced with a violin plot overlaid with a subsampled strip plot (capped at 500 points per channel for readability). This communicates the full distribution shape, density, and individual data points simultaneously. The statistical annotations (t-test and Mann-Whitney p-values) are preserved.

**t-SNE result caching added**
`perform_tsne_knn_analysis` was called independently during both `run_analysis` and `generate_plots`, causing t-SNE to be recomputed twice on identical data in a typical session. t-SNE is non-deterministic (even with `random_state=0`) and computationally expensive. The method now caches its result and skips recomputation if called again with the same channel columns and the same `k` value.

---

### New Features

**Multi-File Search panel**
A new **Multi-File Search** panel has been added at the top of the scrollable UI. It is hidden by default and revealed by ticking the "Enable Multi-File Search" checkbox, leaving the existing single-file workflow completely unchanged.

When enabled, the panel provides:

- **Browse… button** — opens a folder picker to set the root search directory. If a filter string is already entered, the search is triggered automatically on folder selection.
- **"Filename contains" entry** — type any substring (case-insensitive) to filter CSV files by name.
- **Search button** — recursively walks the entire folder tree using `os.walk`, collects every `.csv` file whose name contains the filter string, and reports the number of matches (shown in green) or a "none found" message (in red).
- **Matched file dropdown** — populated with the filename only (e.g. `Bassoon_condition1.csv`) for readability. The full path is stored internally and looked up by dropdown index at load time, so duplicate filenames in different subfolders are handled safely.
- **Load Selected File button** — passes the chosen file through the exact same loading pipeline as the original single-file Load Data button. Channel detection, column auto-fill, and all downstream analyses are identical.

The original **Load Data** button remains fully functional and independent of this panel.

---

## [V0.094051] — 2025 (initial release)

Initial public release. See Zenodo record for details:
https://doi.org/10.5281/zenodo.17533244
