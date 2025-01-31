"""
######################## Islander function ########################
This function partitions the graph and returns a bi-connected 
sub-graph. Note that Islander receives a PreProcessed graph
G_p where the number of nodes may differ from its real id.
"""
module Islander

using Graphs
using MetaGraphs
using JuMP
using LinearAlgebra

include("Auxiliary.jl")
using .Auxiliary

# Main function
function islander(G::MetaDiGraph{Int64, Float64}, 
                  V_g::Vector{Int}, 
                  optimizer,
                  outlev = 0)
    # Turn the graph into a undirected graph to apply articulation
    G_undir = MetaGraph(G)
    # Find articulation nodes in V_g 
    articulationSources = intersect(articulation(G_undir), [u for u in vertices(G) if get_prop(G, u, :id) in V_g])
    # Find connected components - This would change the correlation node - :id!
    for point in articulationSources
        rem_vertex!(G_undir, point)
    end
    components = connected_components(G_undir)
    # Split the graph on the articulation sources
    subgraphs = []; powerInComps = []
    for comp in components
        # create an empty sub-graph for each partition
        partition = MetaDiGraph(length(comp)) 
        nodesInComp = [get_prop(G_undir, v, :id) for v in comp]
        powerInComp = sum(get_prop(G_undir, v, :p) for v in comp if !(get_prop(G_undir, v, :id) in V_g))
        # operate over each node in the partition
        for (idx, u) in enumerate(comp)
            # add the properties in the new sub-graph
            set_prop!(partition, idx, :p,  get_prop(G_undir, u, :p))
            set_prop!(partition, idx, :id, get_prop(G_undir, u, :id))
            set_prop!(partition, idx, :h, get_prop(G_undir, u, :h))
            # find properties of u in the different reference systems
            u_real = get_prop(G_undir, u, :id)                                      # real u id
            u_idx = findfirst(x -> x == u_real, nodesInComp)                        # u idx in component
            u_idx_in_G = findfirst(x -> get_prop(G, x, :id) == u_real, vertices(G)) # u idx in G
            # look for connections
            if !isempty(inneighbors(G,u_idx_in_G) ∪ outneighbors(G,u_idx_in_G))  
                for v_idx_in_G in inneighbors(G,u_idx_in_G) ∪ outneighbors(G,u_idx_in_G)
                    println()
                    v_real = get_prop(G, v_idx_in_G, :id) 
                    if !(v_real in nodesInComp)                                  # v is in V_g
                        add_vertex!(partition); push!(nodesInComp,v_real)
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :p, get_prop(G, v_idx_in_G, :p))
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :id, v_real)
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :h, get_prop(G, v_idx_in_G, :h))
                    end
                    v_idx = findfirst(x -> x == v_real, nodesInComp)
                    if has_edge(G,u_idx_in_G,v_idx_in_G)
                        add_edge!(partition, Edge(u_idx,v_idx))
                        weight = get_prop(G, Edge(u_idx_in_G,v_idx_in_G), :w)
                        set_prop!(partition, Edge(u_idx,v_idx), :w, weight)
                        try 
                            capacity = get_prop(G, Edge(u_idx_in_G,v_idx_in_G), :c)
                            set_prop!(partition, Edge(u_idx,v_idx), :c, capacity)
                        catch e
                            nothing
                        end
                    elseif has_edge(G,v_idx_in_G,u_idx_in_G)
                        add_edge!(partition, Edge(v_idx,u_idx))
                        weight = get_prop(G, Edge(v_idx_in_G,u_idx_in_G), :w)
                        set_prop!(partition, Edge(v_idx,u_idx), :w, weight)
                        try
                            capacity = get_prop(G, Edge(v_idx_in_G,u_idx_in_G), :c)
                            set_prop!(partition, Edge(v_idx,u_idx), :c, capacity)
                        catch e
                            nothing
                        end
                    end
                end
            end
        end
        # add connections between articulation nodes
        articulationNodesInComp = intersect(nodesInComp,articulationSources)
        pairs = collect(Iterators.product(articulationNodesInComp, articulationNodesInComp))
        for pair in pairs
            if has_edge(G, pair[1], pair[2])
                add_edge!(partition, findfirst(x -> x == pair[1], nodesInComp), findfirst(x -> x == pair[2], nodesInComp))
                edge = Edge(findfirst(x -> x == pair[1], nodesInComp), findfirst(x -> x == pair[2], nodesInComp))
                weight = get_prop(G, pair[1], pair[2], :w)
                set_prop!(partition, edge, :w, weight)
                try
                    capacity = get_prop(G, pair[1], pair[2], :c)
                    set_prop!(partition, edge, :c, capacity)
                catch e
                    nothing
                end
            end
        end
        # save the partition
        push!(subgraphs, partition) 
        push!(powerInComps, powerInComp)
    end
    # Split the flow for each articulation node in each partition
    A = splitFlow(G, V_g, subgraphs, powerInComps, optimizer, outlev)
    for i=1:length(subgraphs)
        nodes = [get_prop(subgraphs[i], v, :id) for v in vertices(subgraphs[i])]
        for j=1:length(V_g)
            if V_g[j] in nodes
                idx = findfirst(x -> x == V_g[j], nodes)
                set_prop!(subgraphs[i], idx, :p, get_prop(subgraphs[i], idx, :p)*A[i,j])
            end
        end
    end
    return subgraphs
end

# Function to split the flow
function splitFlow(G::MetaDiGraph{Int64, Float64}, 
                   V_g::Vector{Int}, subgraphs::Vector{Any}, 
                   flow::Vector{Any},
                   optimizer,
                   outlev = 0)
    # Define problem b-Ax = 0
    model = Model(optimizer)
    outlev == 0 && set_optimizer_attribute(model, "print_level", outlev)
    # Known variables
    x = [get_prop(G, findfirst(x -> get_prop(G,x,:id) == v, vertices(G)), :p) for v in V_g]
    b = -flow
    # Unkown variables
    rows = length(flow)
    columns = length(V_g)
    @variable(model, A[1:rows, 1:columns]) # rows : sub_graph, columns : sources
    for i=1:length(subgraphs)
        nodesInSubg = [get_prop(subgraphs[i], v, :id) for v in vertices(subgraphs[i])]
        sourcesInSubg = intersect(V_g,nodesInSubg)
        sourcesInSubgIdx = [findfirst(==(b), V_g) for b in sourcesInSubg]
        sourcesNotInSubg = setdiff(V_g,nodesInSubg)
        sourcesNotInSubgIdx = [findfirst(==(b), V_g) for b in sourcesNotInSubg]
        # check if the subgraph has positive or negative flow (not all sources are delete)
        if flow[i]>=0
            for source in sourcesInSubgIdx
                flow[i] = flow[i] - get_prop(subgraphs[i], i, :p)
                @constraint(model, A[i,source] == 1)
            end
        end
        for j=1:columns
            if j in sourcesNotInSubgIdx
                @constraint(model, A[i,j] == 0)
            end
            if i==1
                @constraint(model, sum(A[:,j]) <= 1) 
                @constraint(model, -sum(A[:,j]) <= 0)
            end
        end
    end
    # Add constraints to avoid cycles:
    pairs = collect(combinations(collect(1:rows), 2))
    for pair in pairs
        @constraint(model, prod(A[pairs[1],:])*prod(A[pairs[2],:]) == 0)
    end
    @constraint(model, sum(A[:,:])>= 1e-6)
    # Standarize constants
    x = Float64.(x)
    b = Float64.(b)
    # Solve the problem
    @objective(model, Min, sum((b - A * x) .^ 2))
    optimize!(model)
    A_optimal = value.(A)
    return A_optimal
end

end
