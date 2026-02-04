# Determining the Optimal *p* in the Morris Method via Pareto-Based Multi-Objective Selection

MATLAB code for **selecting the Morris grid level `p`** using a **Pareto-based multi-objective** criterion, as described in the accompanying paper:

> *Determining the Optimal p in the Morris Method via Pareto-Based Multi-Objective Selection* (course paper / technical report)

---

## Motivation

In the Morris / Elementary Effects (EE) method, the grid level **`p`** determines the step size

\[
\Delta(p) = \frac{2}{p-1},
\]

which directly affects the estimated EE statistics (e.g., \(\mu^\*\), \(\sigma\)) and can change screening stability and interpretability.  
Instead of a heuristic `p`, this repo treats **choosing `p` as a design problem** with explicit trade-offs.

---

## Method overview (what we optimize)

We evaluate candidate grid levels \(P=\{4,5,\dots,30\}\) under a fixed budget with repeated runs, and compute **three objectives**:

1. **Ranking reproducibility** \(J_{\text{rank}}(p)\) (higher is better)  
   - Median **Kendall’s \(\tau_b\)** between factor rankings across repeated runs.

2. **Cross-\(p\) convergence** \(D_{\text{conv}}(p)\) (lower is better)  
   - Maximum relative change of robust \(\tilde{\mu}^\*(p)\) when moving to the next candidate \(p^+\).

3. **Diagnostic-plane separability** \(J_{\text{sep}}(p)\) (higher is better)  
   - Separation of factor groups in the \((\tilde{\mu}^\*, \tilde{\sigma})\) plane:
     \[
     J_{\text{sep}}(p) = \frac{d_{\min}(p)}{s_{\max}(p)+\varepsilon}.
     \]

We then compute the **Pareto frontier** for \(\big(J_{\text{rank}}(p), J_{\text{sep}}(p), -D_{\text{conv}}(p)\big)\) and choose a single recommended \(p^\star\) using a deterministic rule:

> **Default selection rule (paper):** enforce a convergence threshold \(\tau_{\text{conv}}\), then maximize \(J_{\text{sep}}\), and break ties by \(J_{\text{rank}}\).

---

## Repository layout




```

.
├── src/                   # helper functions used by the pipeline (see paper / scripts)
├── config.m               # experiment configuration (candidate set, budgets, thresholds, paths, RNG, etc.)
├── main_select_p.m        # main entry: evaluate objectives, Pareto front, deterministic selection of p*
├── export_EE_long_csv.m   # export long-format EE / indices to CSV (for plotting or downstream analysis)
├── export_viz_to_csv.m    # export visualization-ready summaries to CSV (Pareto/front/objectives/planes)
└── README.md

````

---

## Requirements

- MATLAB (recommended R2021b+)
- Toolboxes may be required depending on your local functions (commonly Statistics & Machine Learning Toolbox)

---

## Quick start

### 1) Clone

```bash
git clone https://github.com/loyiv/Determining-the-Optimal-p-in-the-Morris-Method-via-Pareto-Based-Multi-Objective-Selection.git
cd Determining-the-Optimal-p-in-the-Morris-Method-via-Pareto-Based-Multi-Objective-Selection
````

### 2) Run in MATLAB

```matlab
addpath(genpath(pwd));

% (optional) edit experiment settings
edit config.m

% run the full pipeline
run main_select_p.m
```

---

## Reproducing the paper setting (default)

In `config.m`, set/verify:

* Candidate set: `P = 4:30`
* Trajectories per candidate: `h = 1200`
* Repeats per candidate: `R = 8`
* Convergence threshold: `tau_conv = 0.020` (paper example)
* Aggregation across repeats: component-wise **median** (robust)

Then:

```matlab
run main_select_p.m
```

Expected high-level outputs:

* Objective curves over `p`: (J_{\text{rank}}(p)), (D_{\text{conv}}(p)), (J_{\text{sep}}(p))
* Pareto-front visualization (3D or projections)
* ((\tilde{\mu}^*, \tilde{\sigma})) diagnostic plots for candidate `p`
* Final chosen `p_star` (paper default budget: `p_star = 15`)
* CSV exports via:

  * `export_EE_long_csv.m`
  * `export_viz_to_csv.m`

> Output paths / filenames are controlled by `config.m`.

---


---

## References

* Morris, M. D. (1991). *Factorial Sampling Plans for Preliminary Computational Experiments*. **Technometrics**, 33(2), 161–174. doi:10.1080/00401706.1991.10484804
* Campolongo, F., Cariboni, J., & Saltelli, A. (2007). *An effective screening design for sensitivity analysis of large models*. **Environmental Modelling & Software**, 22(10), 1509–1518. doi:10.1016/j.envsoft.2006.10.004
* Herman, J., & Usher, W. (2017). *SALib: An open-source Python library for Sensitivity Analysis*. **JOSS**, 2(9), 97. doi:10.21105/joss.00097

---


