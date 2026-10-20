# Integrated modelling for VUV diagnostic for JT-60SA

**Author:** R. Cicioni  (rachele.cicioni@igi.cnr.it)  
**Last update:** 2026-05-29  

---

## Info

This document provides the necessary context to independently use the codes developed to reproduce line-of-sight VUV signals within the VUV project for JT-60SA, coordinated by L. Carraro.
Both this document and the associated Git repository aim to clearly document the workflow and track the development steps in an organized and reproducible way, ensuring that the work can be easily understood, reused, and extended in the future.

## Context

For the JT-60SA tokamak in Japan, a future transition to tungsten (W) wall components is planned. This project aims to provide density and temperature profiles along arbitrary lines of sight for the VUV diagnostic system, for both argon (Ar) and tungsten (W) impurities. The synthetic data have been generated using integrated modelling with the JINTRAC suite of codes, specifically the COCONUT module. This framework simulates the plasma in coupled 2D maps covering both the core and the scrape-off layer (SOL), allowing for a self-consistent description of the entire plasma region.  
The modelling work was supervised by L. Garzotti in his role as scientific coordinator of the Transport and Confinement (TC) group.

---

## Index

1. [JINTRAC simulations](#1-jintrac-simulations)  
   1.1 [What is JINTRAC](#11-what-is-jintrac)  
   1.2 [Impurities modelling in JINTRAC](#12-impurities-modelling-in-jintrac)  
   1.3 [Simulations list](#13-simulations-list)  
2. [VUV LOS code](#2-vuv-los-code)  
   2.1 [Overview](#21-overview)  
   2.2 [Usage](#22-usage)  
   2.3 [Mesh data structure](#23-mesh-data-structure)  
   2.4 [Available quantities](#24-available-quantities)  
   2.5 [Output visualization](#25-output-visualization)  
   2.6 [Return value](#26-return-value)  
   2.7 [EDGE2D available outputs](#27-edge2d-available-outputs)  
3. [Bibliography](#3-bibliography-and-wikipages)

---

# 1. JINTRAC simulations

## 1.1 What is JINTRAC? - Briefly -

The JINTRAC [1] framework allows the simultaneous simulation of core, edge, and SOL regions, accounting for their physical interconnection in a self-consistent manner. At the core of JINTRAC is the COCONUT [2] coupling scheme, which couples the 1.5D core transport code JETTO [3] with the 2D SOL transport code EDGE2D [4]. In the simulations used in this collaboration, an integrated core+edge+SOL modelling framework (COCONUT) was used to evolve the plasma self-consistently across the entire domain (core and SOL).

---

## 1.2 Impurities modelling in JINTRAC

Impurities in JINTRAC are modelled using the SANCO impurity transport solver, which evolves the density of each ionization stage of a given element. In full detail, this would require solving up to Z+1 coupled equations per impurity species. To reduce computational cost in integrated modelling, a charge-state bundling (***superstage***) approach is often used. In this approach, multiple adjacent ionization stages are grouped into a single “superstage". The number of evolved variables is therefore reduced to: N_superstages ≪ Z. Within each superstage, ionization states are assumed to be in internal ***coronal equilibrium***, meaning that fast atomic processes between bundled charge states are not resolved explicitly [5].  
Example: for W (tungsten) it was used the 6-stages boundling. Instead of solving 74 ionization stages, they are grouped in 6 superstages only.  

***ADAS database coupling***  
Atomic data are taken from the ADAS database provided in JDC server:  
`/usr/local/depot/sim/data/adas/adas_v1.2/adf11`  
JINTRAC uses a hybrid ADAS-year structure, meaning that different physical processes may come from different ADAS “years” to ensure both accuracy and numerical stability. The processes considered are:
- ACD: Coll.-diel. effective recombination coefficients  
- SCD: Coll.-diel. effective ionisation coefficients  
- CCD: Coll.-rad. effective charge exchange coefficients  
- PLT: Coll.-rad. line power driven by excitation of dominant ions  
- PRB: Coll.-diel. continuum and line power driven by recombination and Bremsstrahlung of dominant ions  
- PRC: Coll.-rad. effective charge exchange recombination power loss rate coefficients  
- PLS: Coll.-rad. line power from selected transitions of dominant ions (obsolete?)  
- ZCD: Superstate averaged Z  
- YCD: Superstate averaged Z^2  
- ECD: Superstate ionisation potential (eV)  

Therefore, “ADAS YEAR” in JINTRAC is not a single consistent dataset, but a composed model.  
The parameter "ADAS LABEL" defines the bundling scheme applied to the impurity charge states. Bundling files are located with the coefficients for each processes, year, label are in:
`adf11/plt<year>/plt<year>_<element>02<label>.dat`  
Example: plt42_w_02_6.dat cointaints the coefficients for W for the PLT process, year 42, label '02_6' which corresponds to 6 superstages.  
The correspondence between the label and the number of superstages can be found in this page of JINTRAC wiki: https://users.euro-fusion.org/pages/data-cmg/wiki/ADAS_database_for_EDGE2D.html


---

## 1.3 Simulations list

List of simulations with full paths in JDC (JET Data Center) cluster:
- COCONUT, Ar+W case: run_sa_fk_c_p20_a3_ArW_fc

List of simulations with full paths in EFGW (EUROfusion Gateway) cluster:
- COCONUT, case without impurities: run_sa_fk_c_WAr_p20_s19_a3_50ms_std

List of username and reference person in JDC cluster:
- jnv7243: R. Cicioni
- gm7685: S. Gabriellini
- lgarzot: L. Garzotti

List of username and reference person in EFGW (EUROfusion Gateway) cluster:
- g2rcicio: R. Cicioni
---

# 2. VUV LOS code

## 2.1 Overview

The `edge2d_los_profile.py` module performs Line-of-Sight (LoS) reconstruction of any physical quantity (density, temperature, radiation, etc.) from 2D JETTO/JINTRAC simulations.

**Main functionality:**
- Extract arbitrary physical quantities from JINTRAC tran files
- Compute intersections of a line of sight with the EDGE2D 2D mesh 
- Interpolate values along the LoS path
- Visualize 2D mesh with LoS path and 1D profile plot

**Code structure:**
1. Load mesh structure (vertices, quadrilaterals)
2. Load desired physical quantity (temperature, density, etc.)
3. Map quantity values to mesh centroids using nearest-neighbor search
4. Define arbitrary LoS geometry (two points in R-Z space)
5. Compute intersections and extract 1D profile
6. Generate visualization plots

---

## 2.2 Usage

### From command line
In /home/jnv7243/jt60sa_git/vuv_carraro folder:

```
# Extract electron temperature along specific LoS
python edge2d_los_profile.py 

```
<!--
### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tran_file` | str | Required | Path to JINTO TRAN file |
| `quantity_name` | str | 'TEVE' | Signal name to extract (e.g., 'TEVE', 'NEDE', 'PRADR') |
| `point1` | [R, Z] | [4.5, 0.0] | LoS starting point in meters |
| `point2` | [R, Z] | [3.0, 2.0] | LoS ending point in meters |
| `plot` | bool | True | Generate visualization plots |
| `verbose` | bool | True | Print detailed information |

---

## 2.3 Mesh data structure

JETTO/JINTRAC stores mesh vertices and centroids in a **quintuple structure**:
- Every 5 consecutive points: 4 vertices + 1 centroid per quadrilateral cell
- Data arrays: `RVERTP`, `ZVERTP` (R and Z coordinates)
- Vertices stored at indices [0, 1, 2, 3], centroid at index [4]

The code automatically reconstructs the mesh topology and maps 1D profile quantities (defined at mesh points) to 2D quad centroids using a KDTree nearest-neighbor search for efficiency.

---

## 2.4 Available quantities

To list all available quantities in a TRAN file:

```python
import eproc as ep
signals = ep.names("/path/to/tran")
```

Common quantities:
- **TEVE**: Electron temperature (eV)
- **NEDE**: Electron density (m⁻³)
- **NIMP**: Impurity density (m⁻³)
- **PRADR**: Radiated power density (W/m³)
- **TEV** + stage number: Temperature per impurity stage (e.g., TEV1 for W stage 1)
- **NEV** + stage number: Density per impurity stage

Note: Z coordinate in TRAN files is stored negative; the code automatically applies the sign correction.

---

## 2.5 Output visualization

The code generates a single figure with 2 panels:

**Left panel:** 2D mesh colored by quantity value, with LoS path overlay
- Red line: Line of sight trajectory
- Green circle: Start point (point1)
- Magenta circle: End point (point2)
- Red dots: Centers of intersected mesh cells

**Right panel:** 1D profile plot
- Y-axis: Quantity value along the path
- X-axis: Ray parameter `t` (t=0 at point1, t=1 at point2)
- Colored regions: Values in each intersected cell
- Gray dashed lines: Cell boundaries

Both panels include colorbars with quantity units automatically extracted from TRAN metadata.

---

## 2.6 Return value

The `extract_los_profile()` function returns a dictionary with:

```python
{
    'segments': [(t_enter, t_exit, value), ...],  # List of path segments
    'nodes': ndarray,                              # Mesh node coordinates
    'quads': ndarray,                              # Quad connectivity indices
    'values': ndarray,                             # Quantity at each cell
    'quantity_name': str,                          # Full quantity description
    'rmesh': ndarray,                              # R mesh point coordinates
    'zmesh': ndarray,                              # Z mesh point coordinates
}
```

Each segment tuple `(t_enter, t_exit, value)` describes:
- `t_enter`: Ray parameter at cell entry
- `t_exit`: Ray parameter at cell exit
- `value`: Quantity value in that cell



-->
---

## 2.7 EDGE2D available outputs

To list all signals available from an EDGE2D simulation, run:
```bash
$ edge2d_output.py
```

Signals under "Profiles" section can be mapped along arbitrary LoS by changing the PROFILE parameter in edge2d_los_profile.py 

# 3. Bibliography and wikipages

[1] JINTRAC Wiki: https://users.euro-fusion.org/pages/data-cmg/wiki/index.html    
[2] Taroni A et al. 1996 Proceedings of the 16th iaea international conference on fusion energy Proc. 16th IAEA Int. Conf. on Fusion Energy vol 2 (Montreal, Canada: IAEA) p 477  
[3] Cenacchi G and Taroni A 1988 Jetto a free boundary plasma transport code Tech. rep. ENEA   
[4] Simonini R, Corrigan G, Radford G, Spence J and Taroni A 1994 Contributions to Plasma Physics 34 368–373   
[5] Bundled description of heavy impurities in JETTO/SANCO and EDGE2D, Laura Lauro-Taroni    
https://users.euro-fusion.org/pages/data-cmg/wiki/files/JETTO_bundledimpurities.pdf 
