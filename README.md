# Lid-driven-cavity-solver-by-MATLAB

A 2D solver to solve lid-driven cavity flow problem has been implemented here. 
The discretization is performed by using a second-order central scheme in space in addition to the fourth-order RK method for time integration. 
The solution is advanced to steady state using different grid refinements starting from nx=ny=33 to 260. 
The results of solution at nx=ny=33 are shown in figure 1.

(a): ![image](https://user-images.githubusercontent.com/89004966/152657677-c2e91f6f-b573-4920-a4d9-470f995c4126.png)

By comparing the line plots of velocity components, this solver gives similar behavior to those results in the reference with significant difference in the magnitude. 

Also, the contours are not smooth. 
By performing further refinement to the grid size, the solution slowly closes to those in reference as shown in figure 2 and figure 3. 
At nx=ny=260, which is the double of 129, the line plots of velocity components are very similar and close to those results reported in the reference. 
Also, the contours are finer and more organized as shown in figure 4. 
