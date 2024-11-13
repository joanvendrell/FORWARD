"""
######################## Pre-Processor function ########################
This function compacts pendant nodes and returns a bi-connected graph.
Graph = 
   { nodes: Name, p ;
     edges: EndPoints, Weight, Capacity 
   }
"""
module PreProcessor

using Graphs
using MetaGraphs
include("Auxiliary.jl")
using .Auxiliary

# Main function
function preprocessor(G::MetaDiGraph{Int64, Float64})
    G_p = deepcopy(G)
    # Define solution
    S = MetaDiGraph(nv(G_p))
    # Find one-degree edges
    pendantNodes = [v for v in vertices(G_p) if degree(G_p, v) == 1]
    parents = [Auxiliary.find_parents(G_p, node) for node in pendantNodes]
    # Find the IDs
    pendantNodesID = [get_prop(G_p, v, :id) for v in pendantNodes]
    parentsID = [get_prop(G_p, v[1], :id) for v in parents]
    # Start a loop
    while ~isempty(pendantNodes)
        for idx in length(pendantNodesID)
            u = findfirst(x -> get_prop(G_p, x, :id) == pendantNodesID[idx], vertices(G_p))
            v = findfirst(x -> get_prop(G_p, x, :id) == parentsID[idx], vertices(G_p))
            # update parent the data
            set_prop!(G_p, v, :p, get_prop(G_p, u, :p) + get_prop(G_p, v, :p))
            # sample in solution set
            add_edge!(S, parentsID[idx], pendantNodesID[idx]) 
            # delete pendant nodes
            rem_vertex!(G_p, u) 
        end
        # re-compute degree
        pendantNodes = [v for v in vertices(G_p) if degree(G_p, v) == 1] 
        parents = [Auxiliary.find_parents(G_p, node) for node in pendantNodes]
        # re-compute real id
        pendantNodesID = [get_prop(G_p, v, :id) for v in pendantNodes]
        parentsID = [get_prop(G_p, v[1], :id) for v in parents]
    end
    return G_p,S
end

end