push!(LOAD_PATH, "/home/aaron/projects/LevyProcesses/src")

using LevyProcesses

Δs = [(:put, 10.), (:put, 25.), (:put, 35.), (:put, 45.), (:call, 50.), (:call, 45.), (:call, 35.), (:call, 25.), (:call, 10.)]
Ts = [30, 60, 90, 120, 150, 180, 252, 360]

S = 50.0
r = 0.02

model = HestonModel(S, r)

volsurface = volatility_surface(Δs, Ts, model)

using Plots

pyplot()

# surf(surface)
# savefig(pwd()*"/plots/surface.png")

plot(volsurface)
savefig(pwd()*"/plots/surface_slices.png")