"""
######################## Custom Objective #########################
This file implements a custom objective functions. 
"""
module CustomObjective

using PowerModelsDistribution
import JuMP
const _PMD = PowerModelsDistribution

function custom_objective_transmission_loss_relaxed(pm::_PMD.AbstractUnbalancedPowerModel)
    objective = 0
    for (n, nw_ref) in _PMD.nws(pm)
        for (i,switch) in ref(pm, n, :switch)
            # Power in branches
            t = _PMD.ref(pm, n, :switch, i)["index"]
            variable_z =_PMD.var(pm, n, :z, i)[t]
            objective = objective + variable_z
        end
    end
    return JuMP.@objective(pm.model, Min, objective)
end

end