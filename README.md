# Interplanetary Trajectory Optimizatioon
## Introduction
This model primarily generates the data and figures used in the Acta Astronautice paper titled "Trajectories Optimization for Asteroid Kinetic Deflection Missions: Potential Benefits of Eccentric Impacts". It includes all computational models and result processing. Everyone is welcome to use and modify this code! Just remember to reference the paper! If you publish your work, please cite my article! I would greatly appreciate it!

- Lee, Kinthong, Hexi Baoyin, and Zhaokui Wang. "Trajectories optimization for asteroid kinetic deflection missions: Potential benefits of eccentric impacts." Acta Astronautica (2025).
- Lee, Kinthong, Zhengqing Fang, and Zhaokui Wang. "Investigation of the incremental benefits of eccentric collisions in kinetic deflection of potentially hazardous asteroids." Icarus 425 (2025): 116312.
   
This model aims to solve impulse-low-thrust trajectory optimizing problem to deflect potentially hazardous asteroids. The model is coded by Matlab, and its working on both Window & MacOS system.
However for macOS, the first time running would need to approve some system permissions for the mex file(for JPL's SPICE). After that, it should be working fine.

---

## How to Use
- IMPORTANT !!! BEFORE YOU BEGIN!! Make sure you download JPL planetary ephemerides ".bsp" files, specifically de441_part-1.bsp and de441_part-2.bsp (The lastest JPL planetary ephemerides on 26th Sept 2024) from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/ or from BaiduNetdisk: : https://pan.baidu.com/s/1lEcd3QQUUZNuWDM-AmBgjg?pwd=6t5h passcode: 6t5h (for chinese mainland users) and add it to location "interplanetary_trajectory_optimization/code/kernel". Otherwise the high precision orbit propagator will fails. Because both planetary ephemerides are way larger than github's limitation on files upload (100MB), so the planetary ephemerides are leave for you to download separately.

- The backbone of the entire model is in **`main_code.m`**. Simply open it and run. If you want to explore the details or make edits, just follow the logic inside **`main_code.m`**.

---

## Data Handling

- On each run, the code reads data from **`PHA_table.xlsx`**.  
- After the program finishes, it will automatically generate **`temporary_result.xlsx`** to **avoid overwriting** the original **`PHA_table.xlsx`** with incorrect data.  
- The **.mat** files result will be stored in:  
  **`/interplanetary_trajectory_optimization/output_result`**  
  To prevent overwrite with incorrect data. Make sure it is correct, then only manually move it to:  
  **`/interplanetary_trajectory_optimization/output_result/matfile`**

---

## Updating Results

- Once you confirm the results are correct:  
  - **Option 1:** Manually copy the relevant data from **`temporary_result.xlsx`** into **`PHA_table.xlsx`**.  
  - **Option 2:** Delete **`PHA_table.xlsx`** and rename **`temporary_result.xlsx`** to **`PHA_table.xlsx`**.

---

> ⚠️ **Warning**  
> Do not overwrite **`PHA_table.xlsx`** directly before confirming the results.  
> Always check **`temporary_result.xlsx`** first to avoid losing valid data.

---

## Plotting and Animations

- Every code file starting with **`plot`** is used to generate graphs or animations.  
- The outputs will be stored in:  
  - **`/interplanetary_trajectory_optimization/output_result/figure`**  
  - **`/interplanetary_trajectory_optimization/output_result/animation`**

## PSO Search Process Animations

Below are some animations showing the process of PSO searching history:


https://github.com/user-attachments/assets/d1c2cffd-a14f-4da8-8fa6-f66ab9e65d89

https://github.com/user-attachments/assets/4bf1cb29-a44d-48b4-b9f1-c532c3f5f9df

![轨道转移示意图](https://github.com/KinthongLee/interplanetary_trajectory_optimization/blob/main/output_result/figure/detailed_trajectory/Modified/99942Apophis_transfer.png)

---

## Note on Numerical Differences vs. the Paper
## What Changed (Important Implementation Note)

According to the collision-model formulation used in the paper, when resolving the momentum transfer on a candidate impact face, the following three vectors must be **coplanar**:

- **Relative velocity vector**: `v_r`
- **Impact-face unit normal**: `n`
- **Impact-face unit tangential direction**: `t`

Mathematically, this requirement can be written as:

v_r ∈ span{n, t}
⇔
n · (t × v_r) = 0


In the original implementation (used for the published results), the tangential direction `t` was constructed using a **fixed global reference axis**, which does not generally guarantee that `v_r`, `n`, and `t` lie in the same plane. This issue is particularly noticeable in **BIP (eccentric impact)** cases.

The current code fixes this by constructing `t` from the **incoming relative velocity direction**, i.e. the projection of `-v_r` onto the impact-face tangent plane. This enforces coplanarity **by construction**.

### Impact on the Results

This issue only affected the **quantitative values** of the computed deflection.  
The **qualitative conclusions of the paper remain valid**, in particular:

**Eccentric (off-center) impacts can provide additional deflection benefits compared to centric impacts.**

After correcting this geometric inconsistency, the updated results are **more geometrically consistent** and may yield **slightly improved** deflection performance.

### Example: (99942) Apophis

As a concrete example, the asteroid **(99942) Apophis** highlights the numerical impact of the corrected implementation.

Using the original (paper-version) code:

- **COG strategy deflection distance**: 1,306,060.18  
- **BIP strategy deflection distance**: 1,769,910.92  
- **Relative gain (BIP vs. COG)**: **35.52%**

After fixing the coplanarity issue between `v_r`, `n`, and `t`:

- **COG strategy deflection distance**: 1,359,482.65  
- **BIP strategy deflection distance**: 1,957,621.50  
- **Relative gain (BIP vs. COG)**: **44%**

These results show that enforcing the correct geometric relationship does not change the qualitative conclusion of the paper, but leads to **quantitatively larger and more physically consistent deflection estimates**, particularly for the eccentric (BIP) impact strategy.


---

