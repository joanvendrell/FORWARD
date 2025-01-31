"""
######################## NetConcad function #########################
This function condenses the graph and returns a quasi-bipartite graph.
"""

module NetConcad

include("Auxiliary.jl")
using .Auxiliary
using Graphs
using MetaGraphs
using Combinatorics

# Main function
function netconcad(G::MetaDiGraph{Int64, Float64}, 
                   V_g::Vector{Int})
    # Split the graph into G_pos and G_neg
    G_pos = MetaGraph(G); G_neg = MetaGraph(G)
    nodes = [get_prop(G, v, :id) for v in vertices(G)]
    inputOutput = [get_prop(G, v, :p) for v in vertices(G)]
    for i=length(nodes):-1:1  # delete in reverse order to avoid problems with indexing
        if inputOutput[i]>=0
            rem_vertex!(G_neg, i)
        else
            rem_vertex!(G_pos, i)
        end
    end
    # In G_pos, consider V_g as independent trees if not needed - this will never allow connections
    [rem_edge!(G_pos, u, v)&rem_edge!(G_pos, v, u) 
        for (u, v) in collect(combinations([v for v in vertices(G_pos)], 2)) 
            if (has_edge(G_pos, u, v) || has_edge(G_pos, v, u))]
    # Find connected components
    superSourceNodes = connected_components(G_pos) 
    superSinkNodes = connected_components(G_neg)
    # Turn the values of vertices into real id
    superSourceNodesID = map(sublist -> map(v -> get_prop(G_pos, v, :id), sublist), superSourceNodes)
    superSinkNodesID = map(sublist -> map(v -> get_prop(G_neg, v, :id), sublist), superSinkNodes)
    # Compact power in each super node
    superSourceNodesPower = map(sublist -> sum(map(v -> get_prop(G_pos, v, :p), sublist)), superSourceNodes)
    superSinkNodesPower = map(sublist -> sum(map(v -> get_prop(G_neg, v, :p), sublist)), superSinkNodes)
    # Unify the super nodes
    superNodes = superSourceNodesID ∪ superSinkNodesID
    superNodesPower = vcat(superSourceNodesPower, superSinkNodesPower)
    # Create the condensed graph - first k nodes are sources k+1 -> n are sinks
    G_cond = MetaDiGraph(length(superNodes))
    for edge in edges(G)
        u = src(edge); v = dst(edge)
        idx = findfirst(x -> get_prop(G, u, :id) in x, superNodes); jdx = findfirst(x -> get_prop(G, v, :id) in x, superNodes) # Mapping
        set_prop!(G,u,:superNode,idx);set_prop!(G,v,:superNode,jdx)                        # Add superNode information in G
        if idx!=jdx                                                                        # Add edge to condensed graph G_cond
            add_edge!(G_cond,idx,jdx)
            set_prop!(G_cond,Edge(idx,jdx),:w,get_prop(G,Edge(u,v),:w))
            set_prop!(G_cond,idx,:p,superNodesPower[idx])
            set_prop!(G_cond,jdx,:p,superNodesPower[jdx])
            set_prop!(G_cond,idx,:id,get_prop(G, u, :id)); set_prop!(G_cond,jdx,:id,get_prop(G, v, :id))
        end
    end
    # Add degree information
    [set_prop!(G_cond, v, :deg, length(inneighbors(G_cond, v) ∪ outneighbors(G_cond, v))) for v in vertices(G_cond)]
    # Add power inependence information
    [(get_prop(G_cond, v, :p)>=0)&&(sum([1 for u in outneighbors(G_cond, v) 
         if abs(get_prop(G_cond, u, :p)) <= abs(get_prop(G_cond, v, :p))]) >= 1) ? set_prop!(G_cond, v, :independence, 1) : 
             set_prop!(G_cond, v, :independence, 0) for v in vertices(G_cond)]
    return G_cond
end

# Update function
function update_graph(G::MetaDiGraph{Int64, Float64},
                      G_cond::MetaDiGraph{Int64, Float64},
                      S::MetaDiGraph{Int64, Float64},
                      edge::Graphs.SimpleGraphs.SimpleEdge{Int64},
                      weight::Float64)
    # Define the number of nodes previous to optimization as reference
    max_superNode = 0
    for r in vertices(G)
        if get_prop(G, r, :superNode)>=max_superNode
            max_superNode = get_prop(G, r, :superNode)
        end
    end
    modSuperNodeID = get_prop(G,dst(edge),:superNode) # This superNode will be modified
    # Move sampled node into its parent superNode and
    set_prop!(G, dst(edge), :superNode, get_prop(G, src(edge), :superNode))
    for node in vertices(G)
        if (get_prop(G, node, :superNode) == modSuperNodeID)&&(!isempty(inneighbors(S, node))||!isempty(outneighbors(S, node)))
            set_prop!(G, node, :superNode, get_prop(G, src(edge), :superNode))
        end
    end
    # Update the power information on the superNodes dual graph
    set_prop!(G_cond,modSuperNodeID,:p,get_prop(G_cond, modSuperNodeID,:p)-get_prop(G,dst(edge),:p))
    set_prop!(G_cond,get_prop(G,src(edge),:superNode),:p,get_prop(G_cond, get_prop(G,src(edge),:superNode),:p)+get_prop(G,dst(edge),:p))
    # Re-configure the modified superNode, note that after losing one vertex, a superNode may be internally 
    # disconnected into more than one superNode
    superNodeGraph = [v for v in vertices(G) if get_prop(G, v, :superNode) == modSuperNodeID]
    # Check if two superNodes have completely got fusionated 
    if isempty(superNodeGraph) 
        for node in vertices(G_cond)
            if has_edge(G_cond,node,modSuperNodeID)
                rem_edge!(G_cond,node,modSuperNodeID)
                add_edge!(G_cond,node,get_prop(G, src(edge), :superNode))
            elseif has_edge(G_cond,modSuperNodeID,node)
                rem_edge!(G_cond,modSuperNodeID,node)
                add_edge!(G_cond,get_prop(G, src(edge), :superNode),node)
            end
        end
    else
        # Get the sub_graph of the modified superNode and check its connectivity
        sub_G = Auxiliary.get_subgraph(G,superNodeGraph)        
        conComps = connected_components(sub_G)
        # Update the number of superNodes if needed 
        for (i, subList) in enumerate(conComps)
            for u in subList
                if i == 1
                    set_prop!(G, findfirst(x -> get_prop(sub_G, u, :id) == get_prop(G, x, :id), vertices(G)), :superNode,  
                                     get_prop(sub_G, u, :superNode) + (i-1) )
                else
                    set_prop!(G, findfirst(x -> get_prop(sub_G, u, :id) == get_prop(G, x, :id), vertices(G)), :superNode,  
                                     max_superNode + (i-1) )
                end
            end
        end
        #[set_prop!(G, findfirst(x -> get_prop(sub_G, u, :id) == get_prop(G, x, :id), vertices(G)), :superNode,  
        #                             get_prop(sub_G, u, :superNode) + (i-1) ) 
        #                             for (i, subList) in enumerate(conComps) for u in subList]
        for i in 1:(length(conComps)-1)
            add_vertex!(G_cond) 
            set_prop!(G_cond,nv(G_cond),:p,0)
        end
        # Update superNodes information in the dual graph
        for node in vertices(G)
            if get_prop(G,node,:superNode)>nv(G_cond)
                set_prop!(G_cond,modSuperNodeID,:p,get_prop(G_cond, modSuperNodeID,:p)-get_prop(G,node,:p))
                set_prop!(G_cond,get_prop(G,node,:superNode),:p,get_prop(G_cond, get_prop(G,node,:superNode),:p)+get_prop(G,node,:p))
            end
        end
        # Change the edges in the dual graph
        [rem_edge!(G_cond,src(edge),dst(edge)) for edge in edges(G_cond) 
                                               if (src(edge) == modSuperNodeID)||(dst(edge) == modSuperNodeID)]
        for edge in edges(G)
            if !has_edge(G_cond,get_prop(G,src(edge),:superNode),get_prop(G,dst(edge),:superNode))
                add_edge!(G_cond,get_prop(G,src(edge),:superNode),get_prop(G,dst(edge),:superNode))
            end
        end
        # Add degree information
        [set_prop!(G_cond, v, :deg, length(inneighbors(G_cond, v) ∪ outneighbors(G_cond, v))) for v in vertices(G_cond)]
        # Add power inependence information
        [(get_prop(G_cond, v, :p)>=0)&&(sum([1 for u in outneighbors(G_cond, v) 
          if abs(get_prop(G_cond, u, :p)) <= abs(get_prop(G_cond, v, :p))]) >= 1) ? 
              set_prop!(G_cond, v, :independence, 1) : set_prop!(G_cond, v, :independence, 0) for v in vertices(G_cond)]     
    end
    # Update flow on the added node
    set_prop!(G, dst(edge), :h, get_prop(G, src(edge), :h)+weight)
    # Update flow on the build polytree
    newFlow = get_prop(G, src(edge), :p) + get_prop(G, dst(edge), :p)
    [set_prop!(G, v, :p, newFlow) for v in vertices(G) if get_prop(G, v, :superNode) == get_prop(G, src(edge), :superNode)]
end

end
