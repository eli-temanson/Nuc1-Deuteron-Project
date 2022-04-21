
using LinearAlgebra
using Plots
using NumericalIntegration
using Optim
using SpecialFunctions

#Constants
ħc = 197.327054 #MeV fm
mp = 938.27231  #Proton Mass
mn = 939.56563  #Neutron Mass
μ = mp*mn/(mp+mn)
Rmax = 25 #fm
N = 1000
dx = 1/N
r = LinRange(0.0, Rmax, N)
dr = Rmax/(2*N)
c1 = 2*μ/ħc^2
a = 137.03599

μp = 2.7928473443 # per nuclear magneton
μn = -1.91304273 # per nuclear magneton

function V18(r)
	vpw = zeros(Float64, 2, 2)
	ccall((:av18pw_,"./av18pot.so"),Cvoid,
		(Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
			Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Array{Float64,2}}),
		1,0,1,1,0,1,-1,r,vpw)
	return vpw
end

function G(r,E) #Array of multipliers for operators in V14
	pot = V18(r)
	g00 = c1*(E + pot[1])
    g01 = c1*pot[2]
    g10 = g01
    g11 = 6/r^2 + c1*(E + pot[4])      

    return [g00 g01; g10 g11] #Generate and return (2x2)
end


f(r,E) = I + G(r,E) * dr^2 / 12 

function Numerov(U,E,dx,Rmax)
    #Compute first step
    c0 = 2*I + 5/6 * dx^2 * G(dx,E)
    cp = I - dx^2 / 12 * G(2*dx,E)

    U0 = cp \ (c0*U)
    Um = U

    # Data Storage of Wavefunctions and r
    r = [0, dx, 2*dx]
    u = [0, U[1], U0[1]]
    w = [0, U[2], U0[2]]

    # Numerov Loop: Start Computing for r = 3*dx
    for i in 2:(Rmax*(1/dx) + 1) # radius array
        rm = (i-1)*dx
        r0 = i*dx
        rp = (i+1)*dx
        append!(r,rp)

        # Coefficent Matrcies 2x2
        cm = G(rm,E)*(dx^2 / 12) - I
        c0 = 2*I + (5/6 * dx^2)*G(r0,E)
        cp = I - (dx^2 / 12)*G(rp,E)

        # Compute Next Step (Ax=b -> x = A\b)
        Up = cp \ (c0*U0 + cm*Um)

        # Proceed
        Um = U0
        U0 = Up

        append!(u,U0[1])
        append!(w,U0[2])
    end

    return r, u, w
end

# Long Range Solutions
fs(r,E) = exp(-√(c1*E)*r) # s-wave l=0
fs′(r,E)= -sqrt(c1*E)*fs(r,E)
fd(r,E) =  0.0256 * fs(r,E) * (1 + 3/(√(c1*E)*r) + 3/(√(c1*E)*r)^2) # d-wave l=2
fd′(r,E)= -fs(r,E) * 0.0256 * (√(c1*E) + 3/r + 6 / (√(c1*E) * r^2) + 6 / (√(c1*E)^2 * r^3))

function Solve(E)
    # Run Numerov
    r,u0,w0 = Numerov([1e-6, 0],E,dx,Rmax)
    r,u2,w2 = Numerov([0, 1e-6],E,dx,Rmax)

    # Compute Derivatives -- Central Difference
    l = length(u0)
 
    ϕ =[u0[l-1], 
        w0[l-1], 
        (u0[l]-u0[l-2])/(2*dx), 
        (w0[l]-w0[l-2])/(2*dx)] # column vector
    
    ψ =[u2[l-1], 
        w2[l-1], 
        (u2[l]-u2[l-2])/(2*dx), 
        (w2[l]-w2[l-2])/(2*dx)] # column vector

    FS=[fs(Rmax,E), 
        0, 
        fs′(Rmax,E), 
        0] # column vector
    
    FD=[0, 
        fd(Rmax,E), 
        0, 
        fd(Rmax,E)] # column vector

    return abs(det([ϕ ψ FS FD])) # return determinant
end


# ==================================================
# Plot: Determinate Minimization
e = LinRange(1,3,150)
plot(e, e->Solve(e), line=(2, :black, :solid), 
    yaxis=:log, 
    xlabel="Binding Energy (MeV)", 
    ylabel="det|M|",label=false,
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 18,
    yguidefontsize = 18,
    legendfontsize = 18)

savefig("det-min.png")


# ==================================================
res = optimize(Solve, 2.21, 2.23)
println("Binding Energy: -",round(Optim.minimizer(res); digits=8), ", MeV")

E = Optim.minimizer(res)
g = sqrt(c1*E)

# Run Numerov
r,u0,w0 = Numerov([1e-6, 0],E,dx,Rmax)
r,u2,w2 = Numerov([0, 1e-6],E,dx,Rmax)

# Compute Derivatives
l = length(u0)
ϕ = [u0[l-1], w0[l-1], (u0[l]-u0[l-2])/(2*dx), (w0[l]-w0[l-2])/(2*dx)]
ψ = [u2[l-1], w2[l-1], (u2[l]-u2[l-2])/(2*dx), (w2[l]-w2[l-2])/(2*dx)]

ϕc, ϕt, ϕc′, ϕt′ = ϕ
ψc, ψt, ψc′, ψt′ = ψ

b = (fs(Rmax,E) * ϕc′ - fs′(Rmax,E) * ϕc) / (ψc * ϕc′ - ψc′ * ϕc)
a = b * (fd(Rmax,E) * ψt′ - ψt * fd′(Rmax,E)) / (ϕt * fd′(Rmax,E) - ϕt′ * fd(Rmax,E))

u = a*u0 + b*u2 # s-wave Solution
w = a*w0 + b*w2 # d-wave Solution

A = 1 / √(integrate(r, u.^2 + w.^2))
u = u .* A
w = w .* A

# ======================================================
# Plot: Plotting wave function and long range solutions
plot(r, [u w], 
    xlabel = "Radius (fm)",
    ylabel = "Radial W.F.",
    linewidth = 2, 
    xtickfontsize = 14,
    ytickfontsize = 14,
    xguidefontsize = 18,
    yguidefontsize = 18,
    legendfontsize = 18,
    lab = ["u(r)" "w(r)"])

r_asym = 2:0.1:Rmax
plot!(r_asym, r_asym->fs(r_asym,E), linewidth=2, lab="fs(r)")
plot!(r_asym, r_asym->fd(r_asym,E), linewidth=2, lab="fd(r)")
savefig("wf.png")



# ==================================================
println("Matter Radius: ",√(0.25*integrate(r, r.^2 .* (u.^2 + w.^2))), ", fm")

Wd = integrate(r, w.^2)
println("l=2 probability: ", Wd)
println("Magnetic Moment: ", (μp + μn) - 3/2*(μp + μn - 1/2)*Wd, ", μN")

Qd = integrate(r, w .* r.^2 .* (u .- w/√(8)))/√(50)
println("Quadrapole Moment: ", Qd, ", e fm^2")

# ==================================================
