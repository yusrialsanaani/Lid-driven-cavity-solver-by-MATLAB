# Lid-driven-cavity-solver-by-MATLAB

A 2D solver to solve lid-driven cavity flow problem has been implemented here using the Finite Difference method at Re = 100. 
The discretization is performed by using a second-order central scheme in space in addition to the fourth-order RK method for time integration. 
The solution is advanced to steady state using different grid refinements starting from nx=ny=33 to 260. 

At each boundary, an appropriate one-sided second-order scheme in lieu of the second-order central scheme was used.
The simulation done using a domain size of L = 1, Reynolds number of Re = 100, Mach number of Ma = 0.1, Prandtl number of 0.71, and ratio of specific heats = 1.4, and grid points nx = ny = 33.

The results of solution at nx=ny=33 are shown in figure 1 (a) - (d).

**Figure 1: Velocity & pressure contours and line plot of velocity components at nx=ny =33.**

**Figure 1 (a):** ![image](https://user-images.githubusercontent.com/89004966/152657677-c2e91f6f-b573-4920-a4d9-470f995c4126.png)

**Figure 1 (b):** ![image](https://user-images.githubusercontent.com/89004966/152657745-0452d52c-776e-4ae5-a5bc-4fbf2504a084.png)

**Figure 1 (c):** ![image](https://user-images.githubusercontent.com/89004966/152657788-7275f569-4821-43c4-84ce-4d528baba985.png)

**Figure 1 (d):** ![image](https://user-images.githubusercontent.com/89004966/152657807-4dff1bc7-1d2e-48a2-8bf6-d1bb5ff41ced.png)



By comparing the line plots of velocity components, this solver gives similar behavior to those results in the reference (Ghia et al. (Journal of Computational Physics, 48, 387-411 (1982))) with significant difference in the magnitude. Also, the contours are not smooth. 


By performing further refinement to the grid size with nx=ny=260, the line plots of velocity components are very similar and close to those results reported in the reference. 
Also, the contours are finer and more organized as shown in figure 2 (a) - (b). 

**Figure 2: Velocity & pressure contours and line plot of velocity components at nx=ny =260.**

**Figure 2 (a):** ![image](https://user-images.githubusercontent.com/89004966/152658498-c5c304d2-2e4c-4a24-8ac2-d0b0dbcbfa48.png)

**Figure 2 (b):** ![image](https://user-images.githubusercontent.com/89004966/152658508-0482f3be-f6f0-4513-865b-a8b1c5d24403.png)

**Figure 2 (c):** ![image](https://user-images.githubusercontent.com/89004966/152658524-cf855c03-0717-4948-8651-c7978f6a619c.png)

**Figure 2 (d):** ![image](https://user-images.githubusercontent.com/89004966/152658536-f6df897d-b049-491b-8127-37d49db44924.png)

#
