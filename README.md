# 🌍 BSM-2D: Boussinesq–Soil Moisture 2D Hydrologic Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**BSM-2D** is a two-dimensional, physically based hydrologic model that solves the transient Boussinesq equation for unconfined groundwater flow. It is coupled with an explicit unsaturated zone (UZ) reservoir model to represent vertical moisture dynamics. This framework provides a computationally efficient alternative to full Richards Equation simulations at the hillslope scale while preserving physical interpretability.

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

## 🧮 Governing Equations

### 1. Transient 2D Boussinesq Equation

\[
\frac{\partial}{\partial x} \left( K_x h \frac{\partial h}{\partial x} \right) + 
\frac{\partial}{\partial y} \left( K_y h \frac{\partial h}{\partial y} \right) + 
N^*(x, y, t) = S_y(x, y, t) \frac{\partial h}{\partial t}
\]

Where:
- \( h(x,y,t) \): water table elevation [L]
- \( K_x, K_y \): saturated hydraulic conductivity [L/T]
- \( N^*(x,y,t) \): recharge rate from UZ [L/T]
- \( S_y(x,y,t) \): specific yield [–]

---

### 2. Specific Yield (van Genuchten-based)

\[
S_y(h) = (\theta_s - \theta_r) \left[1 - \left(1 + (\alpha(h - y))^n \right)^{- \frac{n+1}{n}} \right]
\]

Where:
- \( \theta_s, \theta_r \): saturated and residual moisture content [–]
- \( \alpha, n \): van Genuchten parameters [1/L], [–]
- \( y \): bedrock elevation [L]

---

### 3. Explicit UZ Linear Reservoir Model

\[
\frac{dS_{UZ}}{dt} = I^* - N^* \quad ; \quad N^* = k \cdot S_{UZ}
\]

Discretized as:

\[
S_{UZ}^{t+\Delta t} = S_{UZ}^t + \Delta t (I^t - k S_{UZ}^t)
\]

Where:
- \( S_{UZ} \): unsaturated zone storage [L]
- \( I \): irrigation or rainfall input [L/T]
- \( k \): recession constant [1/T]

---

