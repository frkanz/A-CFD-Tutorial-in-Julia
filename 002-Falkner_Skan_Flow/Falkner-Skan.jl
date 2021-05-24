### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ ac3cc660-08e1-11eb-3db6-0d5949a89423
begin
	using Plots
	using Printf
	using CSV
	using DataFrames
end

# ╔═╡ 59909ca0-08e3-11eb-3348-41080f6a2bd8
md"""

## **Falkner-Skan Similarity Solutions** 

#### **Description:**

This notebook computes the Falkner-Skan Similarity Solutions

		falkner_skan(m, ηmax, N, itermax, ϵProfile, ϵBC)

If the arguments are missing, it will use the default values.
    
		falkner_skan(m=0, ηmax=10, N=30, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)
		
		m = 0 : Blasius flow over a flat plate with a sharp leading edge. 
	0 < m < 1 : Flow over a wedge with half angle θ.
		m = 1 : Hiemenz flow toward a plane stagnation point.
	1 < m < 2 : Flow into a corner with θ.
		m > 2 : No corresponding simple flow.


#### **Falkner-Skan Equation**
Stagnation point flow problem towards an infinite flat plate can be represented using the ordinary differential equation (ODE) below:

$$f''' + ½ (m+1) f f'' + m (1-f')² = 0$$

Boundary Conditions:

$$f(0) = 0$$
$$f'(0) = 0$$
$$f'(η→∞) = 1$$

Please see the lecture notes for [more information] on the derivation of the ODE and numerical solution method.

#### **Solution Method:**

Let's assume:

$$p = f(η)$$		
$$h = f'(η)$$ 
		
Then, we obtain two ODEs equations and boundary conditions
ODEs:

$$h'' + ½ (m+1) p h' + m (1- h²) = 0$$
$$p' - h = 0$$

BCs:  

$$p(0) = 0$$
$$h(0) = 0$$
$$h(∞) = 1$$

The equations given above are solved using Thomas algorithm.

Feel free to ask questions!

*Furkan Oz*
*Kursat Kara*

[foz@okstate.edu](foz@okstate.edu)
[kursat.kara@okstate.edu](kursat.kara@okstate.edu)

"""

# ╔═╡ d47c93d0-08e1-11eb-0e92-21ca5523e6c6
"""
Computes the Falkner-Skan Similarity Solutions
This notebook computes the Falkner-Skan Similarity Solutions

		fs(m, ηmax, N, itermax, ϵProfile, ϵBC)

If the arguments are missing, it will use the default values.
    
		fs(m=0, ηmax=10, N=30, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)
		
		m = 0 : Blasius flow over a flat plate with a sharp leading edge. 
	0 < m < 1 : Flow over a wedge with half angle θ
		m = 1 : Hiemenz flow toward a plane stagnation point.
	1 < m < 2 : Flow into a corner with θ.
		m > 2 : No corresponding simple flow.

Furkan Oz,
[foz@okstate.edu](foz@okstate.edu)
    
"""
function falkner_skan(m=0, ηmax=10, N=30, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)
 
	Δη = ηmax/N
    Δη²= Δη^2;

    #itermax = 40
    #ϵProfile = 1e-6
    #ϵBC = 1e-6

    iter = 0
    errorProfile = 1.
    errorBC = 1.;

    # Initialization of the arrays
    A = zeros(N+1)
    B = zeros(N+1)
    C = zeros(N+1)
    D = zeros(N+1)

    G = zeros(N+1)
    H = zeros(N+1)

    p = zeros(N+1)
    h = zeros(N+1)
    η = zeros(N+1)

    η = [(i-1)*Δη for i=1:N+1];

    #BCs
    h[1]   = 0.0   # h(η=0) = 0
    h[N+1] = 1.0   # h(η=∞) = 1
    p[1]   = 0.0   # p(η=0) = 0

    G[1] = 0.0
    H[1] = 0.0;

    #Solution initialization

    # assume a linear profile for h(η)
    h = [(i-1)/N  for i=1:N+1]

    for i = 2:N+1
        p[i] = p[i-1] + (h[i] + h[i-1])*Δη/2
    end

    println("iter     error            convergence of h(η→∞)")
    println("-----------------------------------------------")

    #while ϵProfile<=errorProfile && ϵBC<=errorBC && iter<itermax
    while ϵProfile<=errorProfile && iter<itermax
        A = [ 1/Δη² + (m+1)*p[i]/(4*Δη) for i=1:N+1]
        B = [-2/Δη² - m*h[i] for i=1:N+1]
        C = [ 1/Δη² - (m+1)*p[i]/(4*Δη) for i=1:N+1]
        D = [ m for i=1:N+1]

        for i=2:N
            G[i] = - ( C[i]*G[i-1] + D[i] )/(B[i] + C[i] * H[i-1])
            H[i] = -                 A[i]  /(B[i] + C[i] * H[i-1])
        end

        hold = copy(h)

        for i=N:-1:2
            h[i] = G[i] + H[i] * h[i+1]
        end 

        errorProfile = maximum(abs.(hold-h))

        errorBC = abs(h[N+1]-h[N])

        for i = 2:N
            p[i] = p[i-1] + (h[i] + h[i-1])*Δη/2
        end

        iter += 1

        #println("iter: $iter,  errorProfile: $errorProfile, errorBC: $errorBC")
        @printf("%4.4d %16.6e %16.6e \n", iter, errorProfile, errorBC)
    end

     if errorProfile<=ϵProfile 
        println("")
        println("Solution converged!")
        println("The maximum change between consecutive profiles is less than the error criteria ϵProfile=$ϵProfile.")
     end

     if errorBC<=ϵBC
        println("")
        println("Solution for the boundary condition converged!")
        println("The difference between h(N) and h(N+1) is less than the error criteria ϵBC=$ϵBC.")
     end
    return η,h
end

# ╔═╡ df424c60-08e1-11eb-319f-2fa7a0f0fb08
falkner_skan();

# ╔═╡ 07798e4e-08e2-11eb-3066-0120057d15b2
ηtest, htest = falkner_skan(0,10,30);

# ╔═╡ f0a39bd0-08e1-11eb-3d16-ffb6b9e631e0
begin
	plot(htest,ηtest,
	        title = "Falkner-Skan Similarity Profile(m=0)",
	        label = "f '",
	        legend = :topleft,
	        xlabel = "h = f'",
	        ylabel = "\\eta",
	        linewidth = 2,
	        linecolor = :black,
	        markershape = :circle,
	        markercolor = :red,
		)
end

# ╔═╡ 13eee630-a6e4-11eb-1573-cd1ba529a131
begin
	ηm4, hm4 = falkner_skan(4,10,100);
	ηm13, hm13 = falkner_skan(1/3,10,100);
	ηm19, hm19 = falkner_skan(1/9,10,100);
	ηmn06, hmn06 = falkner_skan(-0.0654,10,100);
end

# ╔═╡ 6ad9b920-a6e4-11eb-3c66-5975c1e2429c
begin
	datam4 = CSV.read("falkner-schlichting-dataset-m-4.csv",DataFrame)
	datam13 = CSV.read("falkner-schlichting-dataset-m-1d3.csv",DataFrame)
	datam19 = CSV.read("falkner-schlichting-dataset-m-1d9.csv",DataFrame)
	datamn06 = CSV.read("falkner-schlichting-dataset-m--00654.csv",DataFrame)
	ηrefm4 = convert(Array,datam4[1:2:end,1])
	urefm4 = convert(Array,datam4[1:2:end,2])
	ηrefm13 = convert(Array,datam13[1:2:end,1])
	urefm13 = convert(Array,datam13[1:2:end,2])
	ηrefm19 = convert(Array,datam19[1:2:end,1])
	urefm19 = convert(Array,datam19[1:2:end,2])
	ηrefmn06 = convert(Array,datamn06[1:2:end,1])
	urefmn06 = convert(Array,datamn06[1:2:end,2])
	plot_ref = plot([hm4,hm13,hm19,hmn06,urefm4,urefm13,urefm19,urefmn06],[ηm4*sqrt((4+1)/2), ηm13*sqrt((0.3333333+1)/2), ηm19*sqrt((0.1111111+1)/2), ηmn06*sqrt((1-0.0654)/2) ,ηrefm4,ηrefm13,ηrefm19,ηrefmn06],
        title = "Falkner-Skan Similarity Profile",
        label = ["m=4" "m=1/3" "m=1/9" "m=-0.0654" "m=4, Schlichting(2016)" "m=1/3, Schlichting(2016)" "m=1/9, Schlichting(2016)" "m=-0.0654, Schlichting(2016)"],
        legend = :topleft,
        xlabel = "h = f '",
        ylabel = "\\eta",
		ylim = [0,5],
        linewidth = [2.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0],
        linecolor = :auto,
        markershape = [:circle :rect :diamond :hexagon :utriangle :pentagon :heptagon :octagon],
		markeralpha = [1.0 1.0 1.0 1.0 0.6 0.6 0.6 0.6],
        markercolor = :auto
	)
	#png(plot_ref, "falknerskanprofile")
end

# ╔═╡ 638b8090-08e2-11eb-2c26-7b7f6caa7d14
md"""
#### The effect of ηmax:
N is selected as 30 and the error limit is 1e-6.
"""

# ╔═╡ 7bbba7b0-08e9-11eb-3d8c-75a72baa41fe
begin
	η5, h5 = falkner_skan(1,5);
	η10, h10 = falkner_skan(1,10);
	η20, h20 = falkner_skan(1,20);
end
	

# ╔═╡ d8a1b2d0-08e9-11eb-1c9a-91c4d297b41e
plot([h5, h10, h20], [η5, η10, η20],
        title = "Hiemenz Flow - Effect of Viscous Layer Height",
        label = ["\\etamax = 5" "\\etamax = 10" "\\etamax = 20"],
        legend = :topleft,
        xlabel = "h",
        ylabel = "\\eta",
        linewidth = 1,
        linecolor = :black,
        markershape = [:x :utriangle :dtriangle],
        markersize = [5 5 5],
        markeralpha = [0.9 0.6 0.6],
        markercolor = [:black :red :green]
    ,)

# ╔═╡ 2a95ec50-08ea-11eb-06a4-050c654d9c1f
md"""
#### The effect of N:
ηmax is selected as 10 and the error limit is 1e-6.
"""

# ╔═╡ 3649fc80-08ea-11eb-28b5-f58a62c010da
begin
	ηN10, hN10 = falkner_skan(1,10, 10);
	ηN20, hN20 = falkner_skan(1,10, 20);
	ηN40, hN40 = falkner_skan(1,10, 40);
	ηN80, hN80 = falkner_skan(1,10, 80);
end

# ╔═╡ 5241b720-08ea-11eb-3fad-5f8164fe40f6
begin
	gr()
plot([hN10, hN20, hN40, hN80], [ηN10, ηN20, ηN40, ηN80],
        seriestype = :scatter,    
        title = "Hiemenz Flow - Effect of Grid Resolution",
        label = ["N = 10" "N = 20" "N = 40" "N = 80"],
        legend = :topleft,
        xlabel = "h",
        ylabel = "\\eta",
        markershape = [:circle :utriangle :dtriangle :square],
        markersize = [10 8 6 3],
        markeralpha = [0.6 0.6 0.6 0.6],
        markercolor = [:purple :red :green :blue]
    ,)
end

# ╔═╡ Cell order:
# ╟─59909ca0-08e3-11eb-3348-41080f6a2bd8
# ╠═ac3cc660-08e1-11eb-3db6-0d5949a89423
# ╠═d47c93d0-08e1-11eb-0e92-21ca5523e6c6
# ╠═df424c60-08e1-11eb-319f-2fa7a0f0fb08
# ╠═07798e4e-08e2-11eb-3066-0120057d15b2
# ╠═f0a39bd0-08e1-11eb-3d16-ffb6b9e631e0
# ╠═13eee630-a6e4-11eb-1573-cd1ba529a131
# ╠═6ad9b920-a6e4-11eb-3c66-5975c1e2429c
# ╟─638b8090-08e2-11eb-2c26-7b7f6caa7d14
# ╠═7bbba7b0-08e9-11eb-3d8c-75a72baa41fe
# ╟─d8a1b2d0-08e9-11eb-1c9a-91c4d297b41e
# ╟─2a95ec50-08ea-11eb-06a4-050c654d9c1f
# ╠═3649fc80-08ea-11eb-28b5-f58a62c010da
# ╟─5241b720-08ea-11eb-3fad-5f8164fe40f6
