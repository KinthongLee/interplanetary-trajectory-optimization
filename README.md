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
  Make sure it is correct, then manually move it to:  
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
