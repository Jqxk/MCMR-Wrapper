using MCMRSimulator;
using MRIBuilder;
using Plots;
using CairoMakie;
using GLMakie;

geometry = [];

radii = [0.6, 0.6, 0.6, 0.6, 0.6];
repeats = [2., 2.1, 2.25, 2.75, 4.];
susceptibility_iso = [0.03, 0.031, 0.033, 0.04, 0.06]
off_resonance_inside = [0.01, 0.012, 0.02, 0.08, 0.15]
R1_surface = [0.2, 0.21, 0.24, 0.35, 0.5]
R2_surface = [0.02, 0.021, 0.025, 0.035, 0.6]
g_ratio = [0.6, 0.625, 0.65, 0.7, 0.8]

for i in 1:length(radii)

	Cylinder_group = Cylinders(
		radius = radii[i],
		repeats = [repeats[i], repeats[i]],
		R1_surface = R1_surface[i],
		R1_inside = 0.,
		R2_surface = R2_surface[i],
		R2_inside = 0.,
		off_resonance_surface = 0.,
		off_resonance_inside = off_resonance_inside[i],
		position = [0., 0.],
		g_ratio = g_ratio[i],
		susceptibility_iso = susceptibility_iso[i],
		susceptibility_aniso = 0.,
		use_boundingbox = true)

	push!(geometry, Cylinder_group)

end
return geometry

b_scales = [];
G_scales = [];
orientations = [];

function Get_SO(iterations = 1, orientation = [0.0, 1.0, 0.0], mode = 1)

	empty!(b_scales)
	empty!(G_scales)
	empty!(orientations)

	for i in 1:iterations

		b_scale = ((i - 1)/(iterations - 1))
		G_scale = sqrt(b_scale)

		push!(b_scales, b_scale)
		push!(G_scales, G_scale)
		push!(orientations, orientation)

	end

return b_scales, G_scales, orientations
end

b_scales, G_scales, orientations = Get_SO(21, [1., 0., 0.], 1)

#Sequence Hyperparameters & Scanner

b_val_max = 20.
b_scales_plot = b_val_max*b_scales

TE = 100
TR = 500

PCS = Scanner(;B0 = 20.0,
		gradient = 500.,
		slew_rate = 10000.,
		units=:kHz);

sequence = DWI(bval = b_val_max, TE = TE, TR = TR, scanner = PCS)

adj_sequence = adjust(sequence, diffusion=(orientation = orientations, scale = G_scales))

plot_sequence(adj_sequence)

#Simulation definition part

diffusivity = [2.5, 2.45, 2.2, 1.5, 1.0]
permeability = [0, 0.01, 0.05, 0.15, 0.3]

simulations = [];

for i in 1:length(geometry)

	Simulation_group = Simulation(
		adj_sequence,
		geometry = geometry[i],
		R1 = 1/780,
		R2 = 1/90,
		diffusivity = diffusivity[i],
		off_resonance = 0.01,
		permeability = permeability[i])

	push!(simulations, Simulation_group)

end
return simulations

#Readout definition part

sim_time = 0:0.1:(length(G_scales) * TR)

readouts = [];
transverse = [];

for i in 1:length(simulations)

	Readout_group = readout(
		1000,
		simulations[i],
		sim_time,
		skip_TR = 1,
		return_snapshot = false,
		subset = [Subset(inside=true), Subset(inside=false)])

		#Note: Intracellular and Extracellular spin signal detection

	push!(readouts, Readout_group)

end
return readouts

#For every Array in readouts[i] is a readout for cylinders of different radii.
#readouts[i][:, 1] denotes Intracellular signals [Subset(inside=true)];
#readouts[i][:, 2] denotes Extracellular signals [Subset(inside=false)];
#both are vectors of SpinOrientationSum objects.
#SpinOrientationSum.orients.transverse extracts the transverse signal from every SOS in each vector.

transverse = [];

function Transverse(I_E = 1)

	transverse = Float64[]

	for i in 1:length(radii)
		for j in 1:length(b_scales_plot)

			k = 1+ TE*10 + TR*10*(j-1)

			if I_E == 1

				Transverse_group = getproperty.(getproperty.(readouts[i][:, 1], :orient), :transverse)
				Transverse_readout = Transverse_group[k]
				push!(transverse, Transverse_readout)

			elseif I_E == 2

				Transverse_group = getproperty.(getproperty.(readouts[i][:, 2], :orient), :transverse)
				Transverse_readout = Transverse_group[k]
				push!(transverse, Transverse_readout)

			else

				error("Insert 1 for Intracellular transverse signals OR Insert 2 for Extracellular 				transverse signals.")

			end
		end
	end

	return reshape(transverse, length(b_scales_plot), length(radii))
end

transverse = Transverse(1)

#Plotting definition part

als_labels = ["Healthy", "Mild ALS", "Moderate ALS", "Advanced ALS", "Severe ALS"]

fig = Figure()

ax = Axis(fig[1, 1], 
	title="Transverse Magnetization vs b-value",
	xlabel="b-value scale (relative)", 
	ylabel="Transverse Magnetization")

for i in 1:length(radii)

	lines!(ax,
	b_scales_plot,                
           transverse[:, i],                   
           label = als_labels[i])

end

axislegend(ax)
display(fig)
