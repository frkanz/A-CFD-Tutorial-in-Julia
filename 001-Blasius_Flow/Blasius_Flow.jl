### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ bb15be80-a6e8-11eb-0275-d1636e2a6445
begin
	using Plots
	using Printf
	using CSV
	using DataFrames
end

# ╔═╡ 4c484400-a6e8-11eb-343d-e50f20f9e747
md"""

## **_Blasius Flow_**

#### **Description:**

This notebook computes the Blasius Flow (two-dimensional flat plate flow)

		blasius(ηmax, N, itermax, ϵProfile, ϵBC)

If the arguments are missing, it will use the default values.
    
		blasius(ηmax=10, N=30, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)


#### **Blasius Equation**
Flow over an infinite flat plate can be represented using the ordinary differential equation (ODE) below:

$$f''' + ½ f f'' = 0$$

Boundary Conditions:

$$f(0) = 0$$
$$f'(0) = 0$$
$$f'(∞) = 1$$

#### **Solution Method:**

Let's assume:

$$p = f(η)$$		
$$h = f'(η)$$ 
		
Then, we obtain two ODEs equations and boundary conditions
ODEs:

$$h'' + p h' = 0$$
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

# ╔═╡ c46503b0-a6e8-11eb-2f66-b7b20fe36bf8
"""
Computes the Blasius Flow (Two-dimensional flow over a flat plate)

    blasius(ηmax, N, itermax, ϵProfile, ϵBC)

If the arguments are missing, it will use the default values.
    
    blasius(ηmax=10, N=30, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)

Furkan Oz,
[foz@okstate.edu](foz@okstate.edu), 
    
"""
function blasius(ηmax=10, N=30, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)
 
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
        A = [ 2/Δη² + p[i]/(2*Δη) for i=1:N+1]
        B = [-4/Δη²  for i=1:N+1]
        C = [ 2/Δη² - p[i]/(2*Δη) for i=1:N+1]
        D = [0 for i=1:N+1]

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

# ╔═╡ 0a40028e-a6e9-11eb-0fa1-31a1e582b549
blasius();

# ╔═╡ 0d602e50-a6e9-11eb-3e38-9770e40c5747
ηtest, htest = blasius(10,40);

# ╔═╡ 18d364f2-a6e9-11eb-2177-0ba14d02a5b4
begin
	data = CSV.read("blasius-schlichting-dataset.csv",DataFrame)
	ηref = convert(Array,data[1:3:end,1])
	uref = convert(Array,data[1:3:end,2])
	plot_ref = plot([htest,uref],[ηtest*sqrt(1/2),ηref],
        title = "Blasius Similarity Profile",
        label = ["Thomas Algorithm" "Schlichting(2016)"],
        legend = :topleft,
        xlabel = "h = f '",
        ylabel = "\\eta",
		ylim = [0,5],
        linewidth = [2.0 0.000],
        linecolor = :auto,
        markershape = :auto,
		markeralpha = [1.0 0.6],
        markercolor = :auto,
	)
	#png(plot_ref, "blasiusprofile")
end

# ╔═╡ 38f42fd0-a6e9-11eb-26fc-8f7fbe77dbb0
begin
	η5, h5 = blasius(5);
	η10, h10 = blasius(10);
	η20, h20 = blasius(20);
end

# ╔═╡ 3ecfd20e-a6e9-11eb-02d9-9dbb4b44a2be
plot([h5, h10, h20], [η5, η10, η20],
        title = "Blasius Flow - Effect of Viscous Layer Height",
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
    )

# ╔═╡ 46bd8b70-a6e9-11eb-226c-39e68021b054
begin
	ηN10, hN10 = blasius(10, 10);
	ηN20, hN20 = blasius(10, 20);
	ηN40, hN40 = blasius(10, 40);
	ηN80, hN80 = blasius(10, 80);
end

# ╔═╡ 4b765892-a6e9-11eb-0d51-f332409372a6
begin
	gr()
plot([hN10, hN20, hN40, hN80], [ηN10, ηN20, ηN40, ηN80],
        seriestype = :scatter,    
        title = "Blasius Flow - Effect of Grid Resolution",
        label = ["N = 10" "N = 20" "N = 40" "N = 80"],
        legend = :topleft,
        xlabel = "h",
        ylabel = "\\eta",
        markershape = [:circle :utriangle :dtriangle :square],
        markersize = [10 8 6 3],
        markeralpha = [0.6 0.6 0.6 0.6],
        markercolor = [:purple :red :green :blue]
    )
end

# ╔═╡ Cell order:
# ╟─4c484400-a6e8-11eb-343d-e50f20f9e747
# ╠═bb15be80-a6e8-11eb-0275-d1636e2a6445
# ╠═c46503b0-a6e8-11eb-2f66-b7b20fe36bf8
# ╠═0a40028e-a6e9-11eb-0fa1-31a1e582b549
# ╠═0d602e50-a6e9-11eb-3e38-9770e40c5747
# ╠═18d364f2-a6e9-11eb-2177-0ba14d02a5b4
# ╠═38f42fd0-a6e9-11eb-26fc-8f7fbe77dbb0
# ╠═3ecfd20e-a6e9-11eb-02d9-9dbb4b44a2be
# ╠═46bd8b70-a6e9-11eb-226c-39e68021b054
# ╠═4b765892-a6e9-11eb-0d51-f332409372a6
