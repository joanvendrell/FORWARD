"""
######################## Sampler function #########################
This function samples an edge from the graph.
"""
module Sampler

include("Auxiliary.jl")
include("NetConcad.jl")
using .Auxiliary
using .NetConcad
using Graphs
using MetaGraphs
using StatsBase
using LinearAlgebra
using DataStructures

# 1- Sampler function
function sampler(G::MetaDiGraph{Int64, Float64},
                 G_cond::MetaDiGraph{Int64, Float64},
                 S::MetaDiGraph{Int64, Float64}, 
                 strategy::String,
                 outlev = 0)
    verticesInS = unique(v for e in edges(S) for v in (src(e), dst(e))); n = 1
    termination_status = 1; f_S = 0
    while length(verticesInS) < length(vertices(S))
        # define candidates
        edgeCandidates = [Edge(v,u)  for v in vertices(G) 
                                     for u in outneighbors(G, v)  if get_prop(G, v, :p) >= 0 &&
                                                                     #get_prop(G, v, :p) >= abs(get_prop(G, u, :p)) &&
                                                                     get_prop(G, v, :superNode) !== get_prop(G, u, :superNode) && 
                                                                     !Sampler.edge_delete(G, S, Edge(v, u))] 
        # check if there are candidates
        if isempty(edgeCandidates)
            if sum([get_prop(G,i,:p) for i in vertices(G)])>=0 && outlev>0
                println("FORWARD stopped in a feasible solution")
            else
                outlev>0 && println("FORWARD stopped in an unfeasible solution")
            end
            outlev>0 && Auxiliary.show_node(G,:p)
            termination_status = 0
            return S, termination_status, f_S
        end
        # refine candidates
        edgeCandidates, nodeType = Sampler.feasibility_guidance(G, G_cond, edgeCandidates)
        # define weight
        probabilities = Sampler.weight_definition(G, edgeCandidates)
        # select an edge
        edge, weight = Sampler.edge_selector(edgeCandidates,probabilities,nodeType,strategy)
        outlev>0 && println("The sampled edge sampled is ",get_prop(G,src(edge),:id)," - ",get_prop(G,dst(edge),:id)," with weight ",weight)
        # update the information
        add_edge!(S, edge)
        f_S = f_S + weight
        verticesInS = unique(v for e in edges(S) for v in (src(e), dst(e))); n+=1
        NetConcad.update_graph(G,G_cond,S,edge,weight)
    end
    return S, termination_status, f_S
end
# 2- Weight definition
function weight_definition(G::MetaDiGraph{Int64, Float64}, 
                           edgeCandidates::Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1})
    weight = [get_prop(G, edge, :w)+get_prop(G, src(edge), :h) for edge in edgeCandidates]
    probabilities = normalize(weight,1) 
    return probabilities
end
# 3- Edge Delete procedure: delete and edge if it is already connected and if it has positive flow
function edge_delete(G::MetaDiGraph{Int64, Float64}, 
                     S::MetaDiGraph{Int64, Float64}, 
                     edge::Graphs.SimpleGraphs.SimpleEdge{Int64})
    return has_edge(S, edge) || has_edge(S, Edge(dst(edge),src(edge))) || get_prop(G, dst(edge), :p) >= 0
end
# 4- Step candidates guidance
function feasibility_guidance(G::MetaDiGraph{Int64, Float64}, 
                              G_cond::MetaDiGraph{Int64, Float64},
                              edgeCandidates::Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1})
    # look for pendant superSourceNodes
    pendantNodes = [min(get_prop(G_cond, findfirst(x -> get_prop(G, src(edge), :superNode) == x, vertices(G_cond)), :deg), 
                        get_prop(G_cond, findfirst(x -> get_prop(G, dst(edge), :superNode) == x, vertices(G_cond)), :deg)) == 1
                    for edge in edgeCandidates]
    # find flow independent superSourceNodes
    independentNodes = [get_prop(G_cond, findfirst(x -> get_prop(G, src(edge), :superNode) == x, vertices(G_cond)), :independence) 
                        for edge in edgeCandidates]
    # re-order the priority queue
    nodeType = pendantNodes*3 + independentNodes*2 + ones(Int, length(pendantNodes)) 
    perm = sortperm(nodeType, rev=true)
    return edgeCandidates[perm], nodeType[perm]
end
# 5- Edge selector
function edge_selector(edgeCandidates::Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}, 
                       probabilities::Array{Float64, 1},
                       nodeType::Array{Int, 1}, 
                       strategy::String)
    # filter by priority
    edgeCandidates = [e for (i, e) in enumerate(edgeCandidates) if nodeType[i] == maximum(nodeType)]
    probabilities  = [w for (i, w) in enumerate(probabilities)  if nodeType[i] == maximum(nodeType)]
    # sample and edge considering the choosen strategy
    if strategy == "greedy_max"
        idx = argmax(probabilities)
        edge = edgeCandidates[idx]; weight = probabilities[idx]
    elseif strategy == "greedy_min"
        idx = argmin(probabilities)
        edge = edgeCandidates[idx]; weight = probabilities[idx]
    elseif strategy == "random"
        idx = sample(1:length(edgeCandidates), Weights(probabilities))
        edge = edgeCandidates[idx]; weight = probabilities[idx]
    end
    return edge, weight
end
end
