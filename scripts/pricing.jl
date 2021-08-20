push!(LOAD_PATH, "/Users/aaron/projects/LevyProcesses/src")

using LevyProcesses

model = HestonModel(50.0, 0.02)

option = VanillaCallOption(50, 30, model)

print(delta(option))


