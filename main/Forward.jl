"""
######################## FORWARD #########################
Main algorithm implementation. Returns:
  - S: radial configuration
  - termination_status: 1 if completed, 0 if uncompleted
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

function forward(G::MetaDiGraph{Int64, Float64},
                 optimizer,
                 strategy::String,
                 outlev = 0)
    outlev > 1 && println("Running Forward...")
    S = []; termination_status = 1; f_S = 0
    G_p,S_p = PreProcessor.preprocessor(G)
    outlev > 1 && println("PreProcessor applied!")
    V_g_p = [get_prop(G_p,i,:id) for i in vertices(G_p) if get_prop(G_p,i,:p)>0]
    push!(S, S_p)
    subgraphs = Islander.islander(G_p,V_g_p,optimizer, outlev)
    outlev > 1 && println("Islander applied!")
    for id=1:length(subgraphs)
        S_id = MetaDiGraph(nv(subgraphs[id]))
        [set_prop!(S_id,i,:id,get_prop(subgraphs[id],i,:id)) for i in vertices(S_id)]
        Auxiliary.make_bidirected!(subgraphs[id])
        G_cond = NetConcad.netconcad(subgraphs[id],V_g_p)
        outlev > 1 && println("NetConcad applied for parition ", id, " of ", length(subgraphs))
        S_id,termination_status_i,f_S_id = Sampler.sampler(subgraphs[id],G_cond,S_id,strategy,outlev)
        outlev > 1 && println("Sampler applied for parition ", id, " of ", length(subgraphs))
        push!(S, S_id)
        f_S = f_S + f_S_id
        termination_status = min(termination_status, termination_status_i)
    end
    S = Auxiliary.join_subsets(G,S)
    outlev > 1 && println("Ploytrees joined!")
    return S, termination_status, f_S
end


end