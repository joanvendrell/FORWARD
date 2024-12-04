"""
######################## CREATE RANDOM SMALL-WORLD GRAPH #########################
Functions to define a random Watts-Strogatz small-world graph in a .dss, being:
    N - Number of loads
    K - Each node is connected to K nearest neighbors in ring topology
    β - Rewiring probability
    G - Number of generators
    file_path - Output file name
"""
module GraphCreator

using Graphs
using Random

function power_small_world_graph_dss(N::Int, K::Int, β::Float64, G::Int, filename::String)
    # Create a small-world network using the Watts-Strogatz model
    total_nodes = N + G
    ws_graph = watts_strogatz(total_nodes, K, β)
    
    # Open file for writing
    path = "../data/datasets/"*filename*"/data/"
    file_path = joinpath(path, filename*".dss")
    open(file_path, "w") do file
        # Header and circuit definition
        write(file, "Clear\n")
        write(file, "Set DefaultBaseFrequency=60\n")
        total_nodes = N+G
        write(file, "new circuit.LANL$(total_nodes)Nodeckt\n")
        write(file, "~ Angle=45\n")
        write(file, "~ basekv=0.4 pu=1 phases=1 bus1=101.1 enabled=n ! VOLTAGE SOURCE DEFINITION\n\n")       
        # Define line codes
        edge_count = 1
        for u in 1:total_nodes
            for v in neighbors(ws_graph, u)
                if u < v  # Avoid duplicating edges
                    r_value = rand()
                    x_value = rand()
                    write(file, "New linecode.mtx$(lpad(edge_count, 3, '0')) nphases=1 Units=m Rmatrix=[$(round(r_value, digits=2))] Xmatrix=[$(round(x_value, digits=2))]\n")
                    edge_count += 1
                end
            end
        end

        # Define loads
        idx = shuffle(102:100 + total_nodes)
        LoadList = Int[]
        GenList = [101]

        write(file, "\n!LOAD DEFINITIONS -- kV is vm_nom and kW is pd_nom\n")
        total_load_power = 0
        load_powers = Float64[]
        
        for i in 1:N
            load_power = rand(100:1000)
            push!(load_powers, load_power)
            total_load_power += load_power
            write(file, "New Load.$(idx[i]) Bus1=$(idx[i]).1 Phases=1 Conn=POWER Model=1 kW=$load_power kvar=$load_power vminpu=0 vmaxpu=12\n")
            push!(LoadList, idx[i])
        end

        # Define generators to meet total load power
        write(file, "\n!ADD GENERATORS\n")
        GenList = [GenList; setdiff(idx, LoadList)]
        generator_power = total_load_power / G  # Divide total power among generators
        
        for i in 1:length(GenList)
            gen_power = round(generator_power, digits=2)
            write(file, "New Generator.$(GenList[i]) Phases=1 Bus1=$(GenList[i]).1 kw=$gen_power kvar=$gen_power\n")
        end

        # Ensure generator power meets total load power
        if generator_power * G < total_load_power
            error("Total generator power must be equal or greater than the total load power.")
        end

        # Define switches based on adjacency matrix
        write(file, "\n!SWITCH DEFINITIONS\n")
        edge_count = 1
        for u in 1:total_nodes
            for v in neighbors(ws_graph, u)
                if u < v  # Avoid duplicating edges
                    write(file, "New Line.$(u + 100)$(v + 100) Phases=1 Bus1=$(u + 100).1 Bus2=$(v + 100).1 Switch=y Enabled=y LineCode=mtx$(lpad(edge_count, 3, '0'))\n")
                    edge_count += 1
                end
            end
        end

        # Final setup
        write(file, "\nSet voltagebases=[0.4]\n")
        write(file, "Set tolerance=0.000001\n")
        write(file, "Calcvoltagebases\n")
        write(file, "Solve\n")
    end
end


end