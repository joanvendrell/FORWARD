"""
######################## Custom Plot #########################
This file implements a custom plot functions.
"""
module CustomPlot

using GraphRecipes, Plots
using PowerModelsDistribution
const _PMD = PowerModelsDistribution

function custom_PMD_solution_plot(case_math::Dict{String, Any},
                                  case_eng::Dict{String, Any},
                                  solution::Dict{String, Any})
    # Create an empty adjacency matrix and fill it
    loads = collect(keys(solution["solution"]["load"]))
    sources = collect(keys(solution["solution"]["gen"]))
    lines = collect(keys(solution["solution"]["switch"]))
    load_name = []
    for idx in 1:length(loads)
        load = case_math["load"][loads[idx]]["name"]
        push!(load_name,load)
    end
    gen_name = []
    for id in sources
        push!(gen_name,case_math["gen"][id]["name"])
    end
    load_name = vcat(load_name,gen_name)
    println("Nodes: ", load_name)
    adjacency = zeros(length(load_name),length(load_name))
    for line in lines
        name = case_math["switch"][line]["name"]
        node1 = case_eng["switch"][case_math["switch"][line]["name"]]["f_bus"]#[1:3]
        node2 = case_eng["switch"][case_math["switch"][line]["name"]]["t_bus"]#[4:6]
        println("Lines: from ",node1," to ",node2)
        # find the indices
        node1_idx = findfirst(x -> x == node1, load_name)
        node2_idx = findfirst(x -> x == node2, load_name)
        # find the direction of the edge
        if solution["solution"]["switch"][line]["state"]>=0.1
            weight = solution["solution"]["switch"][line]["pt"][1]
            idx = 1+(node1_idx-1)*length(load_name) + (node2_idx-1)
            jdx = 1+(node2_idx-1)*length(load_name) + (node1_idx-1)
            if weight<=0
                adjacency[jdx] = weight
            else
                adjacency[idx] = -1*weight
            end
        end
    end
    # Define colors in function of the type of the node
    normal_color = RGB(0, 0.5, 1); gen_color = RGB(0, 0.7, 0)
    colors = fill(normal_color, length(load_name))
    for idx in 1:length(load_name)
        load = load_name[idx]
        if in(load,gen_name)
            for id in sources#gen
                if case_math["gen"][id]["name"] == load
                    if solution["solution"]["gen"][id]["gen_status"]>=0.1
                        colors[idx] = gen_color
                    end
                end
            end
        end
    end
    # Plot
    println("Adjacency: ",adjacency)
    plt = graphplot(adjacency,
              names=load_name,
              markercolor = colors,
              markersize = 0.1,
              fontsize = 7,
              linecolor = :darkgrey)
    plot!(plt, size=(800, 800))
end

end