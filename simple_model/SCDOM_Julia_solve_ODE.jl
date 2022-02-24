using DifferentialEquations
using Plots

function B₀(γ);
    Bbar, η, γbar = (3.42, 0.5, 1.22);
    Bbar + η * (γbar - γ)
end

function sm(du,u,p,t)
    q₀, q₁, b₀, τ, γ = p; 
    du[1] = -(u[2] - B₀(γ)) - abs(q₀ + q₁*(u[1] - b₀))*(u[1] - b₀);
    du[2] = (u[1] - γ)/τ;
end

u0 = [0.0, 0.0];
p = (-9, 12, 0.625, 0.9, 1.2);
q₀, q₁, b₀, τ, γ = p; 
tspan = (0.0,60.0);

# use Euler-Maruyama (EM)
sm_prob = ODEProblem(sm, u0, tspan, p);
sm_sol = solve(sm_prob, EM(), dt = 1e-1, adaptive=false);

y_sm = [u[1] for u in sm_sol.u];
x_sm = [u[2] for u in sm_sol.u];

tc = 166.66 # use this characteristic time to get the real time
plot(sm_sol.t, x_sm,
    linewidth=2,
    label="x")
plot!(sm_sol.t, y_sm,
    linewidth=2,
    label="y")
