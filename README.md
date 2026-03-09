# E.S.T.H.E.R. 🛰️
**E**nhanced **S**ymbolic **T**ools for **H**armonic **E**quation **R**esolution

An analytical formulation model for the perturbed orbital motion of spacecraft in low-eccentricity orbits, under the influence of the Earth's gravitational zonal harmonics J2 and J3. Developed purely in MATLAB, this project provides explicit closed-form expressions for relative position and velocity, offering an efficient alternative for preliminary mission analysis and implementation on onboard computers (OBC) with limited resources.

## 🚀 Features
* **Closed-Form Cartesian Solutions:** Derives explicit analytical expressions for relative position and velocity in Cartesian coordinates starting from the linearized Clohessy-Wiltshire equations.
* **Advanced Asymptotic Techniques:** Utilizes Taylor series expansions on the eccentricity and trigonometric arguments to accurately model both circular and elliptical low-eccentricity orbits.
* **J2 & J3 Zonal Harmonics Integration:** Accurately calculates the exact perturbative accelerations caused by the Earth's equatorial bulge (J2) and its hemispherical mass asymmetry (J3).
* **Constant Computational Cost:** Offers a lightweight mathematical model that maintains a virtually constant execution time (~0.18 ms), completely bypassing the heavy iterative loads of standard numerical propagators.
* **High Short-Term Accuracy:** Achieves position errors below 500 meters after two full revolutions for polar circular orbits, making it highly reliable for rapid trajectory planning.
* **Comprehensive Numerical Validation:** Fully validated against high-precision numerical orbit propagation (Cowell's method) and cross-checked with classical analytical methods (including Kéchichian's original solutions and Deprit's radial intermediary).
* **Cylindrical Coordinates Formulation:** Explores a novel mathematical approach by reformulating the perturbed equations of motion into a cylindrical coordinate system to preserve full information regarding orbital curvature.

## 🛠️ Technology Stack
This project uses **MATLAB** (R2024b) alongside the **Symbolic Math Toolbox** to process complex mathematical simplifications, resolve coupled differential equations, and evaluate the final analytical functions.

## 📂 Repository Structure
* `/Benito` - Core mathematical derivations, Taylor series expansions, and analytical solvers for J2 and J3.
* `/Clohessy–Wiltshire` - Base linearized equations of relative motion and state-space matrices.
* `/Functions` - Helper functions for matrix rotations (Euler-Hill to ECI) and state transformations.
* `/Kechichian` - Implementation of Kéchichian's original analytical solutions for comparative error analysis.
* `/Melton` - Additional orbital dynamics and relative motion scripts.
* `/zoomPlot` - Visualization utilities for plotting position and velocity errors.
* `UCS_Satellite_Database.txt` - Real-world satellite orbital parameters used to define the 45 matrix study cases.

## ⚙️ Installation & Usage
Since this project is developed entirely in MATLAB, no external compilation or complex dependency management is required. 

1. **Clone the repository:**
   ```bash
   git clone https://github.com/BFG508/ESTHER.git
2. **Open MATLAB** and navigate to the cloned `ESTHER` directory.
3. **Add to Path**: Ensure that all subdirectories (`/Functions`, `/zoomPlot`, etc.) are added to your MATLAB path so the main scripts can access the helper functions. You can do this by right-clicking the `ESTHER` folder in the Current Folder browser and selecting Add to Path > Selected Folders and Subfolders.
4. **Run the Analysis**: Open and run the desired scripts within the `/Benito`, `/Kechichian`, or `/Clohessy–Wiltshire` folders to execute the analytical models or compare them with the numerical propagations. The scripts will automatically read the orbital parameters from `UCS_Satellite_Database.txt` when testing the different study cases.

## 🎓 Academic Context
This repository contains the source code and mathematical tools developed for the Bachelor's Thesis (Trabajo de Fin de Grado) in Aerospace Engineering at Universidad Rey Juan Carlos (URJC), Spain.
* Author: Benito Fernández González
* Tutor: Hodei Urrutxua Cereijo
* Academic Year: 2024/2025