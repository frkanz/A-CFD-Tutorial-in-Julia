# A CFD Tutorial in Julia
A CFD Tutorial in Julia: Introduction to Laminar Boundary-Layer Theory journal's code repository. The codes are developed in the Computational Aerodynamics and Hypersonics Laboratory.

## **Abstract**
The boundary layer is a thin layer near the surface in which the velocity changes from zero at the surface to the free stream value away from the surface. This behavior can be modeled with boundary layer theory which is crucial for many problems in aerodynamics, including wing stall, the skin friction drag on an object, and the heat transfer that occurs in high-speed flight. Understanding complex physics during engineering education may provide students insight into more complex problems, such as boundary-layer transition or flow separation. A computational fluid dynamics (CFD) student must understand the flow physics, discretization methods, and solution algorithms to conduct accurate and successful simulations. Utilizing computer tools are also another aspect of successful simulations. A relatively new programming language, Julia, is becoming one of the most popular programming language. Moreover, it is one of the most user-friendly and fast programming languages. It is crucial to learn and understand the pros and cons of a programming language before starting any project in this environment to validate that it is a convenient language for the project. In this paper, fundamental flow problems such as Blasius flow, Hiemenz flow, Homann flow, and Falkner-Skan flow will be derived from scratch and solved using the Julia language. We used the finite difference scheme to discretize the governing equations and utilized the Thomas algorithm to solve the resulting linear system. In addition, we compared the results with the literature. We released the Julia codes to shorten the learning curve. Additionally, we discussed the pros and cons of the Julia language.


Oz, F.; Kara, K. A CFD Tutorial in Julia: Introduction to Laminar Boundary-Layer Theory. Fluids, 2021, **_Under review_**.

#### **Instructions**

Julia setup files can be downloaded from their website (https://julialang.org/downloads/). The website also includes instructions on how to install Julia on Windows, Linux, and mac operating systems. It is common to use external packages for Julia. In order to do that, Pkg, which is Julia's built-in package manager, can be used. Once Julia is opened, Pkg can be activated with the "]" button in Windows. In Linux, calling "julia" in the terminal will open it. After that "Pkg.add("Pluto")" will trigger the setup process for that package. In here, we used Pluto as an example because, in GitHub, our codes are developed in the Pluto environment. After Pluto is installed. Pluto can be run with "Pluto.run()". This command will open a new tab in the browser which you can run your Julia codes. After that, the "using Pluto" line must be placed to the top of the file. For "Plots" package, the commands will be "Pkg.add("Plots")" and "using Plots". Since the Plots package does not have a GUI, there is not a command called "Plots.run()".

#### **Example Case: Blasius Equation**
Flow over an infinite flat plate can be represented using the ordinary differential equation (ODE) below:

$$f''' + ½ f f'' = 0$$

Boundary Conditions:

f(0) = 0  
f'(0) = 0  
f'(∞) = 1  

#### **Solution Method:**

Let's assume:

p = f(η)  
h = f'(η)  
		
Then, we obtain two ODEs equations and boundary conditions
ODEs:

h'' + p h' = 0  
p' - h = 0  

BCs:  

p(0) = 0  
h(0) = 0  
h(∞) = 1  

The equations given above are solved using Thomas algorithm.

Feel free to ask questions!

*Furkan Oz*  
*Kursat Kara*  

[foz@okstate.edu](foz@okstate.edu)  
[kursat.kara@okstate.edu](kursat.kara@okstate.edu)  
