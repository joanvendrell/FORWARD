"""
######################## Islander function ########################
This function partitions the graph and returns a bi-connected 
sub-graph.
"""
module Islander

using Graphs
using MetaGraphs
using JuMP

# Main function
function islander(G::MetaDiGraph{Int64, Float64}, V_g::Vector{Int}, optimizer)
    # Turn the graph into a undirected graph to apply articulation
    G_undir = MetaGraph(G)
    # Find articulation nodes in V_g
    articulationSources = intersect(articulation(G_undir),V_g)
    # Find connected components
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
        powerInComp = sum([get_prop(G_undir, v, :p) for v in comp])
        # operate over each node in the partition
        for (idx, u) in enumerate(comp)
            # add the properties in the new sub-graph
            set_prop!(partition, idx, :p,  get_prop(G_undir, u, :p))
            set_prop!(partition, idx, :id, get_prop(G_undir, u, :id))
            u_real = get_prop(G_undir, u, :id)
            # look for parents of u
            if !isempty(inneighbors(G,u_real))
                for v_real in inneighbors(G,u_real)
                    if !(v_real in nodesInComp)                                  # v is in V_g
                        add_vertex!(partition); push!(nodesInComp,v_real)
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :p, get_prop(G, v_real, :p))
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :id, get_prop(G, v_real, :id))
                    end
                    add_edge!(partition, findfirst(x -> x == v_real, nodesInComp), findfirst(x -> x == u_real, nodesInComp))
                    edge = Edge(findfirst(x -> x == v_real, nodesInComp), findfirst(x -> x == u_real, nodesInComp))
                    weight = get_prop(G, v_real, u_real, :w)
                    set_prop!(partition, edge, :w, weight)
                end
            end
            # look for children of v
            if !isempty(outneighbors(G,u_real))
                for v_real in outneighbors(G,u_real)
                    if !(v_real in nodesInComp)                                  # v is in V_g
                        add_vertex!(partition); push!(nodesInComp,v_real)
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :p, get_prop(G, v_real, :p))
                        set_prop!(partition, findfirst(x -> x == v_real, nodesInComp), :id, get_prop(G, v_real, :id))
                    end
                    add_edge!(partition, findfirst(x -> x == u_real, nodesInComp), findfirst(x -> x == v_real, nodesInComp))
                    edge = Edge(findfirst(x -> x == u_real, nodesInComp), findfirst(x -> x == v_real, nodesInComp))
                    weight = get_prop(G, u_real, v_real, :w)
                    set_prop!(partition, edge, :w, weight)
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
            end
        end
        # save the partition
        push!(subgraphs, partition) 
        push!(powerInComps, powerInComp)
    end
    # Split the flow for each articulation node in each partition
    A = splitFlow(G, V_g, subgraphs, powerInComps, optimizer)
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
function splitFlow(G::MetaDiGraph{Int64, Float64}, V_g::Vector{Int}, subgraphs::Vector{Any}, flow::Vector{Any}, optimizer)
    # Define problem b-Ax = 0
    model = Model(optimizer)
    # Known variables
    x = [get_prop(G, v, :p) for v in V_g]
    b = -flow 
    # Unkown variables
    rows = length(flow)
    columns = length(V_g)
    @variable(model, A[1:rows, 1:columns])
    for i=1:length(subgraphs)
        sources = [get_prop(subgraphs[i], v, :id) for v in vertices(subgraphs[i])]
        sourcesNotInSubg = setdiff(V_g,sources)
        idxs = [findfirst(==(b), V_g) for b in sourcesNotInSubg]
        for j in idxs
            @constraint(model, A[i,j] == 0) 
        end
    end
    for j=1:columns
        @constraint(model, sum(A[:,j]) == 1) 
    end
    # Solve the problem
    @objective(model, Min, sum((b - A * x) .^ 2))
    optimize!(model)
    A_optimal = value.(A)
    return A_optimal
end

end