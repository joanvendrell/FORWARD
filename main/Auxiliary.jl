"""
######################## Auxiliary function ########################
Functions to get information from the graph
"""
module Auxiliary

using Graphs
using GraphPlot
using MetaGraphs

# 1- function to find parents
function find_parents(G::MetaDiGraph{Int64, Float64}, node)
    parents = []
    for v in vertices(G)
        if has_edge(G, v, node)  # Check if there's an edge from v to node
            push!(parents, v)    # Add v to the parent list
        end
    end
    return parents
end

#2- function to visualize information
function show_node(G::MetaDiGraph{Int64, Float64}, var::Symbol)
    for node in vertices(G)
        weight = get_prop(G, node, var)
        u = get_prop(G, node, :id)
        println("Node $u has $var: $weight")
    end
end

#3- function to turn a graph into a bidirected graph
function make_bidirected!(G::MetaDiGraph)
    for e in edges(G)
        u, v = src(e), dst(e)
        if !has_edge(G, v, u)
            add_edge!(G, v, u)
            set_prop!(G, v, u, :w, get_prop(G, u, v, :w))
        end
    end
    return G
end

#4- function to turn sample a sub-graph from a subset of nodes
function get_subgraph(G::MetaDiGraph{Int64, Float64},list::Vector{Int64})
    sub_G = MetaDiGraph(length(list))
    for (i,node) in enumerate(list)
        set_prop!(sub_G,i,:superNode,get_prop(G,node,:superNode))
        set_prop!(sub_G,i,:p,get_prop(G,node,:p))
        set_prop!(sub_G,i,:h,get_prop(G,node,:h))
        set_prop!(sub_G,i,:id,get_prop(G,node,:id))
    end
    nodes = [get_prop(sub_G,i,:id) for i in vertices(sub_G)]
    [add_edge!(sub_G, idx, jdx) for edge in edges(G) for idx in [findfirst(x -> src(edge) == get_prop(sub_G, x, :id), vertices(sub_G))] 
                                                     for jdx in [findfirst(x -> dst(edge) == get_prop(sub_G, x, :id), vertices(sub_G))] 
              if idx !== nothing && jdx !== nothing]
    return sub_G
end


end