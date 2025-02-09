"""
######################## GRAPH INITIALIZATION #########################
Functions to turn a .dss file into a graph.
"""
module GraphInit

using PowerModelsDistribution
const _PMD = PowerModelsDistribution
using GraphRecipes, Plots
using Graphs
using MetaGraphs

function dc_graph_initialization(path::String)
    case_eng = _PMD.parse_file(path; data_model=ENGINEERING)
    loads = collect(keys(case_eng["load"])); sources = collect(keys(case_eng["generator"]))
    loads = vcat(loads,sources)
    lines = collect(keys(case_eng["switch"]))
    G = MetaDiGraph(length(loads))
    lookup_dict = Dict{String,Int}()
    for line in lines
        node1 = case_eng["switch"][line]["f_bus"]
        node2 = case_eng["switch"][line]["t_bus"]
        node1_idx = 0; node2_idx = 0
        # find the indices
        for idx in 1:length(loads)
            if loads[idx]==node1
                node1_idx = idx
            elseif loads[idx]==node2
                node2_idx = idx
            end
            # create a lookup reference dictionary --> NEW!!!!
            lookup_dict[loads[idx]] = idx
            # add node information
            if loads[idx] in sources
                set_prop!(G,idx,:p,case_eng["generator"][loads[idx]]["pg"][1])
            else
                set_prop!(G,idx,:p,-case_eng["load"][loads[idx]]["pd_nom"][1])
            end
            set_prop!(G,idx,:id,idx)
            set_prop!(G,idx,:h,0)
        end
        # define the edge
        add_edge!(G, node1_idx, node2_idx)
        add_edge!(G, node2_idx, node1_idx)
        # find the resistance
        line_code = case_eng["switch"][line]["linecode"]
        set_prop!(G,Edge(node1_idx,node2_idx),:w,case_eng["linecode"][line_code]["rs"][1])
        set_prop!(G,Edge(node2_idx,node1_idx),:w,case_eng["linecode"][line_code]["rs"][1])
    end
    return G,lookup_dict
end

end
