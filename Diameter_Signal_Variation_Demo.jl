using MCMRSimulator;
using MRIBuilder;
using Plots;
using CairoMakie;
using GLMakie;
using MAT;
using Meshes;
using GeometryBasics;

geometry = [];
radii = [5.];

for i in 1:length(radii)

	radius = radii[i]

	Cylinder_group = Cylinders(
		radius = radii[i],
		repeats = [radii[i]*2.5, radii[i]*2.5],
		R1_surface = 0.,
		R1_inside = 0.,
		R2_surface = 0.,
		R2_inside = 0.,
		off_resonance_surface = 0.,
		off_resonance_inside = 0.,
		position = [0., 0.],
		susceptibility_iso = 0.,
		susceptibility_aniso = 0.,
		use_boundingbox = true)

	push!(geometry, Cylinder_group)

end
return geometry

#Optionally, define your own mesh

#geometry = load_mesh("Users/jqxk/Documents/MCMR/Axon_045.ply")
#geometry = [geometry]

#radii = [0.45];

#Sequence definition part

b_scales = [];
G_scales = [];
orientations = [];

function Get_SO(iterations = 1, orientation = [0.0, 1.0, 0.0], mode = 1)

	empty!(b_scales)
	empty!(G_scales)
	empty!(orientations)

	g_r = (1 + sqrt(5))/2

	if mode == 1

		if !(orientation isa Vector{Float64}) || length(orientation) != 3

			error("Input must be a 3-dimensional vector of Float64 types.")

		end
	
		for i in 1:iterations
	
			b_scale = ((i - 1)/(iterations - 1))
			G_scale = sqrt(b_scale)

			push!(b_scales, b_scale)
			push!(G_scales, G_scale)
			push!(orientations, orientation)

		end

	elseif mode == 2

		for i in 1:iterations

			b_scale = 1
			G_scale = 1
		
			z = 1 - (2*(i - 1) + 1)/iterations
			r = sqrt(1 - z^2)
			theta = 2π*(i-1)/g_r

			x = r*cos(theta)
			y = r*sin(theta)

			push!(b_scales, b_scale)
			push!(G_scales, G_scale)
			push!(orientations, [x, y, z])

		end

	else

		error("Mode must either be 1 (b-sweep diameter variation) or 2 (3D HARDI).")

	end

return b_scales, G_scales, orientations
end

#MODE 1 EXAMPLE:
#b_scales, G_scales, orientations = Get_SO(21, [1., 0., 0.], 1)

#MODE 2 EXAMPLE: 
b_scales, G_scales, orientations = Get_SO(100, [1., 0., 0.,], 2)

#Sequence Hyperparameters & Scanner

b_val_max = 2.5
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

simulations = [];

for i in 1:length(geometry)

	Simulation_group = Simulation(
		adj_sequence,
		geometry = geometry[i],
		R1 = 1/780,
		R2 = 1/90,
		diffusivity = 1.70,
		off_resonance = 0.,
		permeability = 0.)

	push!(simulations, Simulation_group)

end
return simulations
		
#Readout definition part

sim_time = 0:0.1:(length(G_scales) * TR)

readouts = [];
transverse = [];

for i in 1:length(simulations)

	Readout_group = readout(
		5000,
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

				error("Insert 1 for Intracellular transverse signals OR Insert 2 for Extracellular transverse signals.")

			end
		end
	end

	return reshape(transverse, length(b_scales_plot), length(radii))
end

transverse = Transverse(1)

S0 = transverse[1, :]  

transverse .= transverse./S0'


#HARDI Plotting definition (ONLY FOR HARDI VISUALISATION)

function HARDI_plot(transverse, orientations, mode=1)

	points = Vector{Vector{Float64}}()
	points_orientations = Vector{Vector{Float64}}()

	if mode == 1

        		points = [transverse[i] * orientations[i] for i in 1:length(transverse)]
        		vertices = Point3f0.(points)

        		lat_count = 10
        		lon_count = 20

        		faces = GeometryBasics.TriangleFace{Int}[]

        		for lat in 1:(lat_count - 1)
                		for lon in 1:(lon_count - 1)

                        		i = (lat - 1) * lon_count + lon
                        		i_right = i + 1
                        		i_down = i + lon_count
                       		i_diag = i + lon_count + 1

                    			if i_diag ≤ length(vertices)
                        			push!(faces, TriangleFace(i, i_down, i_diag))
                        			push!(faces, TriangleFace(i, i_diag, i_right))
                    			end
            		end
		end

        		mesh = GeometryBasics.Mesh(vertices, faces)
        
        		fig = Figure()

        		ax = Axis3(fig[1, 1], title = "HARDI Mesh")
        		mesh!(ax, mesh, color = [v[3] for v in vertices], colormap = :coolwarm;
            	transparency = true)

	elseif mode == 2

		for i in 1:length(transverse)

			push!(points, transverse[i]*orientations[i])
			push!(points_orientations, orientations[i])

		end

		x_points = [p[1] for p in points]
		y_points = [p[2] for p in points]
		z_points = [p[3] for p in points]

		x_o = [p[1] for p in points_orientations]
		y_o = [p[2] for p in points_orientations]
		z_o = [p[3] for p in points_orientations]

		fig = Figure()

		ax1 = Axis3(fig[1, 1], 
		title="$(length(orientations)) sampled orientations")

		GLMakie.scatter!(ax1, x_o, y_o, z_o, markersize = 5)

		ax2 = Axis3(fig[1, 2], 
		title="HARDI plot for $(length(orientations)) spherically sampled points")

      		GLMakie.scatter!(ax2, x_points, y_points, z_points, markersize = 5)

	end

display(fig)
return points
end

points = HARDI_plot(transverse, orientations, 1)


	
#Plotting definition part


fig = Figure()

ax = Axis(fig[1, 1], 
	title="Transverse Magnetization vs b-value",
	xlabel="b-value scale (relative)", 
	ylabel="Transverse Magnetization")

for i in 1:length(radii)

	lines!(ax,
	b_scales_plot,                
           transverse[:, i],                   
           label = "Radius $(radii[i]) µm")

end

axislegend(ax)
display(fig)

#Export to MATLAB

matopen("transverse_signals.mat", "w") do file
	write(file, "transverse", transverse)
	write(file, "b_values", b_scales_plot)
	write(file, "radii", radii)
end
