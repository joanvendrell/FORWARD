"""
######################## FORWARD #########################
Main algorithm implementation.
"""
module Forward

include("PreProcessor.jl")
include("Islander.jl")
include("NetConcad.jl")
include("Sampler.jl")
include("Auxiliary.jl")
using .PreProcessor
using .Islander
using .NetConcad
using .Sampler
using .Auxiliary
using Graphs
using MetaGraphs
using GraphPlot
using Plots
using JuMP
using Ipopt
using Combinatorics

function forward(G::MetaDiGraph{Int64, Float64},
                 V_g::Vector{Int}, 
                 optimizer,
                 strategy::String)
    G_p,S = PreProcessor.preprocessor(G)
    subgraphs = Islander.islander(G,V_g,optimizer)
    Auxiliary.make_bidirected!(G)
    G_cond = NetConcad.netconcad(G,V_g)
    S = Sampler.sampler(G,G_cond,S,strategy)
    return S
end


end