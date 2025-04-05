using CairoMakie; 
using GLMakie;
using MCMRSimulator; 
using MRIBuilder; 
using GeometryBasics;
using Meshes;
using FileIO;
using LinearAlgebra;
using IterTools

function NoColSpheres(size::Float64, td::Float64, mean::Float64, var::Float64, max_iter::Int64, min_r::Float64, max_r::Float64)

	positions, radii = random_positions_radii(
		(size, size, size), #box_size
		td, #target_density
		3, #n_dimensions for sphere
		mean = mean,
		variance = var,
		max_iter = max_iter,
		min_radius = min_r,
		max_radius = max_r
	)

	geometry = Spheres(
		radius = radii,
		position = positions,
		R1_surface = 0.1,
		R1_inside = 0.1,
		R2_surface = 0.1,
		R2_inside = 0.1,
		off_resonance_surface = 0.0,
		off_resonance_inside = 0.0,
		use_boundingbox = false,
		size_scale = 0.5
	)

	return geometry
end

function check_collision(c1::MCMRSimulator.Geometries.User.Obstructions.ObstructionGroups.ObstructionGroup{:BendyCylinder}, geometry::Vector{MCMRSimulator.Geometries.User.Obstructions.ObstructionGroups.ObstructionGroup{:BendyCylinder}})

	for c2 in geometry #Loop over all cylinders

		for p1 in c1.control_point.value #point on 'spline' cylinder 1

			for p2 in c2.control_point.value #point on 'spline' cylinder 2

				if abs(p1[3] - p2[3]) > 0.1
					continue
				end

				distance = norm(p1 .- p2)
				threshold = c1.radius.value + c2.radius.value 
				#println("dist = $distance, threshold = $threshold")
		
				if distance < threshold #collision occurs if dist < their radii sum
					return true
				end
			end
		end
	end
	return false
end

function CylGen(z_step::Int64, z_gap::Float64, tortuosity::Float64, spline::Int64, min_r::Float64, max_r::Float64, xy_bound::Float64)

	Start = [rand(-xy_bound:0.01:xy_bound), rand(-xy_bound:0.01:xy_bound), 0.0]

	c_p = [Start .+ [rand(-0.05:0.005:0.05)*tortuosity, rand(-0.05:0.005:0.05)*tortuosity, z_gap*j] for j in 0:z_step]

	geometry = BendyCylinder(
		control_point = c_p,
		radius = rand(min_r:0.005:max_r),
		spline_order = spline,
		myelin = true,
		use_boundingbox = true
		#use default rotation
	)
	
	return geometry
end

function NoColCylinder(N::Int64, max_attempts::Int64, z_step::Int64, z_gap::Float64, tortuosity::Float64, spline::Int, min_r::Float64, max_r::Float64, xy_bound::Float64)

	geometry = MCMRSimulator.Geometries.User.Obstructions.ObstructionGroups.ObstructionGroup{:BendyCylinder}[] #init vector

	for _ in 1:N

		attempts = 0

		while attempts < max_attempts

			new_cyl = CylGen(z_step, z_gap, tortuosity, spline, min_r, max_r, xy_bound) #generate cylinders for the nth time

			if check_collision(new_cyl, geometry) == true

				attempts += 1

			else

				push!(geometry, new_cyl) #add new cylinder
				println("Valid solution found at attempt: $attempts")
				break

			end
		end

		if attempts == max_attempts
			error("Too many attempts without finding a valid solution.")

		end
	end

	return geometry
end

function MeshGen(geo_params)

	geometry = load_mesh(geo_params)

	return geometry
end

function B_Iter(b_step, k, t_echo, t_rep, grad)

	if !(grad isa AbstractVector{<:AbstractVector})
		error("grad must be of type ::Vector{Vector{Float64}}.")

	end

	sequence = [Vector{Any}(undef, length(grad)) for _ in 1:k] #Initialization, holds upcoming values
	TR = t_rep
	bvals = []

	for i in 1:k

		bval = i*b_step
		push!(bvals, bval)

		for j in 1:length(grad)

			sequence[i][j] = DWI(
				bval = bval, 
				TE = t_echo, 
				TR = t_rep, 
				scanner = Siemens_Terra,
				gradient=(type=:instant, orientation = grad[j])
				)
		
		end
	end
	
	sequence = reduce(vcat, sequence)

	return sequence, TR, bvals

end

function B_IterSE(b_step, k, delay, t_echo, t_rep)

	sequence = Vector{Any}(undef, k)

	TR = t_rep
	bvals = []
	
	for i in 1:k

		bval = i*b_step
		push!(bvals, bval)

		sequence[i] = SpinEcho(
			TE = t_echo,               # Echo time (TE) in ms
			TR = t_rep,
			delay = delay,                    # Delay between readout and spin echo peak in ms
			bval = bval,                     # Diffusion weighting (b-value)
			slice_thickness = 2.0,          # Slice thickness in mm
			#excitation = excitation,        # Defined excitation pulse
			#refocus = refocus,              # Defined refocus pulse
			scanner = Siemens_Terra            # Use the predefined scanner
		    )
        end

        return sequence, TR, bvals

end

function Sim(sequence, geometry, T1::Int64, T2::Int64, diffusivity::Float64, end_time::Float64, n_spins::Int64, TR)

	simulation = Simulation(
		sequence, 
		R1 = 1/T1, 
		R2 = 1/T2, 
		diffusivity = diffusivity, 
		off_resonance = 0.1, 
		geometry = geometry)

	sim_time = 0:0.1:end_time

	if end_time > TR
		error("End time cannot exceed TR.")

	end

	return simulation, n_spins, sim_time

end

function ReadIndex(simulation, n_spins::Int64, sim_time, skip_TR::Int64, nTR::Int64, sim_b::Int, sim_g::Int, nTR_i::Int, l_t::Int, grad)

	avg_signals = readout(
		n_spins,       # Number of spins to be simulated
         	simulation, 
         	sim_time,
        		skip_TR = skip_TR, # Dummy cycles to stabilize re-magnetization
        		nTR = nTR,
        		return_snapshot = true
    		)

	sim_b = sim_b
	sim_g = sim_g
	sim_i = length(grad)*(sim_b - 1) + sim_g

	if sim_i > size(avg_signals, 1)
        		error("sim_i exceeds the number of simulations in avg_signals.")
    	end

	if nTR_i > nTR
       		error("nTR_i exceeds the number of TR cycles in avg_signals.")
    	end

   	if l_t == 1

        		signal_type = "Transverse"
        		signals = transverse.(avg_signals[sim_i, :, nTR_i])


    	elseif l_t == 2

        		signal_type = "Longitudinal"
        		signals = longitudinal.(avg_signals[sim_i, :, nTR_i])


    	else 

        		error("Invalid l_t value. Use 1 for Transverse or 2 for Longitudinal signals.")

    	end
    	return signals, signal_type, sim_b, sim_g, avg_signals

end

function Plot_signals(signals, sim_time)

	f = Figure(resolution = (1000, 750))
	
	ax = Axis(f[1, 1],
		title = "Signal against simulation time",
		xlabel = "Time (ms)",
		ylabel = "Signal strength"
	)

	lines!(ax, sim_time, signals, label = "Signal")

	display(f)

	return f
end

function Direction(grad)

	origin = Point3f0(0.0, 0.0, 0.0)

	f = Figure()
	ax = LScene(f[1, 1], scenekw=(show_axis=true, ))

	for (i, g) in enumerate(grad)

		if length(g) != 3

			error("All vectors in grad must be 3-dimensional.")

		end

		u_x, u_y, u_z = g #The i-th 3D vector
		arrows!(ax, [origin], [Point3f0(u_x, u_y, u_z)], arrowsize = 0.2, label = "Gradient $i")

	end

	return f

end

function GradGen(m, n)

	grad = Vector{Vector{Float64}}()

	for i in 1:m
	
		theta = i*(2*pi)/m
		grad_add = [0, 0, 1]' * [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]
		grad_temp = normalize(vec(grad_add'))

		push!(grad, grad_temp)

		for j in 1:n

			phi = j*(2*pi)/n
			grad_add = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)]*grad_temp
			grad = push!(grad, normalize(grad_add))

		end
	end

	grad = [round.(normalize(g), digits = 1) for g in grad]

	return grad
end


#Wrapper function as below


function MCMRSim(
               geo_mode::Int64, # Defines which geometry
               geo_params, # Defines what parameters
	    seq_mode::Int64, #Defines which mode of seq gen
               seq_params::Tuple, # Defines sequence parameters
               sim_params::NamedTuple, # Defines simulation parameters
               out_params::NamedTuple, # Defines readout and indexing parameters
	    optimize::Bool = false
               )

	geometry = nothing

	if geo_mode == 1

		geometry = NoColSpheres(geo_params...)

    	elseif geo_mode == 2

        		geometry = NoColCylinder(geo_params...)

	elseif geo_mode == 3
		
		geometry = MeshGen(geo_params)

   	else
        		error("Invalid geo_mode. Enter 1 for Spheres or 2 for Bendy Cylinders.")

    	end

	sequence = nothing 
	TR = 0 
	grad = []
	bvals = []

	if seq_mode == 1 

    		grad = seq_params[5]  # Extract gradients for DWI
  		sequence, TR, bvals = B_Iter(seq_params[1], seq_params[2], seq_params[3], seq_params[4], grad)

	elseif seq_mode == 2

   		sequence, TR, bvals = B_IterSE(seq_params[1], seq_params[2], seq_params[3], seq_params[4], seq_params[5])

	else

   		error("Invalid seq_mode. Use 1 for DWI mode or 2 for Spin Echo mode.")

	end


	sim_i = 0 

	if seq_mode == 1 

		sim_i = length(grad) * (out_params.sim_b - 1) + out_params.sim_g 
		#FOR 6 ORIENTATIONS AND 5 GRADS, [1] [2] [3] [4] [5] [6] IS FIXED B

	elseif seq_mode == 2

		sim_i = out_params.sim_b

	end

    	sim_params_up = NamedTuple{(:sequence, :geometry, :T1, :T2, :diffusivity, :end_time, :n_spins, :TR)}(
        		(sequence, geometry, sim_params.T1, sim_params.T2, sim_params.diffusivity, sim_params.end_time, sim_params.n_spins, TR)
    	)

    	simulation, n_spins, sim_time = Sim(sim_params_up...)

    	out_params_up = NamedTuple{(:simulation, :n_spins, :end_time, :skip_TR, :nTR, :sim_b, :sim_g, :nTR_i, :l_t)}(
        		(simulation, n_spins, sim_time, out_params.skip_TR, out_params.nTR, out_params.sim_b, out_params.sim_g, out_params.nTR_i, 	out_params.l_t)
    	)

	signals, signal_type, _, _, avg_signals = ReadIndex( 
		out_params_up.simulation, 
		out_params_up.n_spins, 
		sim_time, 
		out_params_up.skip_TR, 
		out_params_up.nTR, 
		out_params_up.sim_b, 
		out_params_up.sim_g, 
		out_params_up.nTR_i, 
		out_params_up.l_t, 
		grad
		)

	return geometry, simulation, signals, signal_type, sim_time, bvals, avg_signals

end

grad = GradGen(100, 0)

geo_mode = 3
#geo_params = (2.0, 0.6, 1., 1., 2000, 0.01, 0.2) 
geo_params = ("Users/jqxk/Documents/MCMR/Test_Neuron2.ply")
seq_mode = 1
seq_params = (0.05, 50, 80, 300, grad) 
#seq_params = (0.05, 1, 5.0, 50, 250)
sim_params = (T1 = 1000, T2 = 100, diffusivity = 2.0, end_time = 299., n_spins = 10) 
out_params = (skip_TR = 0, nTR = 1, sim_b = 1, sim_g = 1, nTR_i = 1, l_t = 1) 

geometry, simulation, signals, signal_type, sim_time, bvals, avg_signals = MCMRSim(geo_mode, geo_params, seq_mode, seq_params, sim_params, out_params)

#println("Signals: ", signals)
#println("Signal Type: ", signal_type)

signals_extract, signals_plot = HARDI_signals(15, 700, 1, grad, avg_signals)
HARDI_plot(signals_plot)


#UNDER DEVELOPMENT FUNCTIONS AND TO-BE-IMPLEMENTED FEATURES

function HARDI_signals(b_index::Int64, time_index::Int64, nTR_i::Int64, grad, avg_signals)

	signals_extract = []
	signals_plot = []


	for i in 1:length(grad)

		sim_index = i + (b_index - 1) * length(grad)
		
		signal = avg_signals[sim_index, time_index, nTR_i]

		push!(signals_extract, transverse(signal))

		#Extract info from fixed bvals, fixed time, for all gradients, fit them in something else

	end

	grad_red = vcat(grad...)

	for k in 1:length(signals_extract)
	
		j = (k - 1)*3 + 1

		grad_dir = grad_red[j:j+2]
		signal_mag = signals_extract[k] * grad_dir

		push!(signals_plot, signal_mag)

	end

	return signals_extract, signals_plot
end

function HARDI_plot(signals_plot)

	origin = Point3f0(0.0, 0.0, 0.0)

	f = Figure()

	ax = LScene(f[1, 1], scenekw=(show_axis=true, ))

	points = [Point3f0(sp...) for sp in signals_plot]
	lines!(ax, points, label = "HARDI signals")
	#scatter!(ax, points, markersize = 10, label = "HARDI signals")

	update_cam!(ax.scene, Vec3f(1, 0, 0), Vec3f(0, 0, 0), Vec3f(0, 5, 0))

	display(f)
	return f

end

