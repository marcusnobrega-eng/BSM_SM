# 🌍 BSM-2D: Boussinesq–Soil Moisture 2D Hydrologic Model

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**BSM-2D** is a two-dimensional, physically based hydrologic model that solves the transient Boussinesq equation for unconfined groundwater flow. It is coupled with an explicit unsaturated zone (UZ) reservoir model to represent vertical moisture dynamics. This framework provides a computationally efficient alternative to full Richards Equation simulations at the hillslope scale while preserving physical interpretability.

The application results of this model are presented in the submitted manuscript titled *"Identifying the State Dependence of Effective Material Properties in a Simplified Hydrologic Hillslope Model,"* which is currently under peer review.

This model was built upon the model presented in:

Brandhorst, N., Erdal, D. and Neuweiler, I., 2021. Coupling saturated and unsaturated flow: comparing the iterative and the non-iterative approach. Hydrology and Earth System Sciences, 25(7), pp.4041-4059.

For a full folder with all datasets used in the paper, please refer to the Zenodo file described in the paper.

---


## 🧠 Model Rationale

At intermediate scales (10–100 m), fully solving the Richards equation in 3D is computationally expensive and often poorly constrained by field data. BSM-2D addresses this by:

- Modeling lateral saturated flow via the **2D Boussinesq equation**
- Representing vertical unsaturated dynamics via a **linear reservoir approximation**
- Supporting **state- and parameter-dependent closure relationships** such as variable specific yield
- Enabling **data assimilation** using particle filters (optional)

---

## 🧪 Key Applications

- Hillslope-scale hydrology
- Seepage discharge modeling
- Coupled vadose–phreatic systems
- Data assimilation with sensor observations
- Conceptual catchment testing

---

## 📄 Governing Equations

### 1. Transient 2D Boussinesq Equation

<div align="center">

$$
\frac{\partial}{\partial x} \left( K_x h \frac{\partial h}{\partial x} \right) +
\frac{\partial}{\partial y} \left( K_y h \frac{\partial h}{\partial y} \right) +
N^*(x, y, t) = S_y(x, y, t) \frac{\partial h}{\partial t}
$$

</div>

**Where:**
- $h(x, y, t)$: water table elevation [L]  
- $K_x$, $K_y$: saturated hydraulic conductivity [L/T]  
- $N^*(x, y, t)$: recharge rate from UZ [L/T]  
- $S_y(x, y, t)$: specific yield [–]  

---

### 2. Specific Yield (van Genuchten-based)

<div align="center">

$$
S_y(h) = (\theta_s - \theta_r) \left[1 - \left(1 + \left(\alpha(h - y)\right)^n \right)^{- \frac{n+1}{n}} \right]
$$

</div>

**Where:**
- $\theta_s$, $\theta_r$: saturated and residual moisture content [–]  
- $\alpha$, $n$: van Genuchten parameters [1/L], [–]  
- $y$: bedrock elevation [L]  

---

### 3. Unsaturated Zone Linear Reservoir

<div align="center">

$$
\frac{dS_{UZ}}{dt} = I^* - N^* \quad ; \quad N^* = k \cdot S_{UZ}
$$

</div>

Discretized as:

<div align="center">

$$
S_{UZ}^{t+\Delta t} = S_{UZ}^t + \Delta t \cdot (I^t - k \cdot S_{UZ}^t)
$$

</div>

**Where:**
- $S_{UZ}$: unsaturated zone storage [L]  
- $I$: input irrigation or rainfall [L/T]  
- $k$: drainage/recession coefficient [1/T]  

---

## 🗂️ Input Data

All input data can be conveniently entered directly into the provided Excel spreadsheet. Each variable within the sheet is clearly defined and thoroughly documented to ensure accurate and user-friendly data entry.

## 🗂️ Repository Structure
```bash
📁 BSM_2D_Model/
├── 📂 GW/ # MATLAB groundwater solver codes
├── 📂 inputs/ # Excel spreadsheets for inputs (DEM, params, forcing)
├── 📂 outputs/ # Model outputs (e.g., heads, discharge, Sy)
├── 📂 Particle_Filter/ # Particle filter source codes
├── 📂 Other folders are from paper results and are not directly required to run the model
