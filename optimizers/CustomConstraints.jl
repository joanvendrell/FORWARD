"""
######################## Custom Power Constraints #########################
This file implements discrete constraints as an extension of 
PowerModelsDistribution.
"""
module CustomConstraints

using PowerModelsDistribution
import PowerModelsONM
import JuMP
const _ONM = PowerModelsONM
const _PMD = PowerModelsDistribution



function import_constraints_transmission_loss(pm::_PMD.AbstractUnbalancedPowerModel,  k::Int=2)
    # Load constraint
    for id in ids(pm, :load)
        _PMD.constraint_mc_load_power(pm, id)
    end
    # Generators constraint
    for id in _PMD.ids(pm, :gen)
        _PMD.constraint_mc_generator_power(pm, id)
        _PMD.constraint_mc_gen_power_on_off(pm, id) 
    end
    constraint_cardinality(pm, k)
    println("Generators constraints defined")
    # Balance constraints
    for id in _PMD.ids(pm, :bus)
        _PMD.constraint_mc_power_balance(pm, id)  #it joins :p and :pg for balance constraints
    end
    println("Balance constraints defined")
    # Ohms constraints
    for i in ids(pm, :branch)
        _PMD.constraint_mc_ohms_yt_from(pm, i)
        _PMD.constraint_mc_ohms_yt_to(pm, i)
        _PMD.constraint_mc_voltage_angle_difference(pm, i)
        _PMD.constraint_mc_thermal_limit_from(pm, i)
        _PMD.constraint_mc_thermal_limit_to(pm, i)
        _PMD.constraint_mc_ampacity_from(pm, i)
        _PMD.constraint_mc_ampacity_to(pm, i)
    end
    println("Ohms constraints defined") 
    # Radiality constraint
    for id in _PMD.ids(pm, :switch) 
        _ONM.constraint_mc_switch_state_open_close(pm,id)
        _PMD.constraint_mc_switch_thermal_limit(pm, id)
        _PMD.constraint_mc_switch_ampacity(pm, id)
    end
    _ONM.constraint_radial_topology(pm)
    println("Radiality constraints defined")
    for id in _PMD.ids(pm, :switch)
        constraint_y_variable_positiveness(pm,id) # Linearization
        mccormick_relaxation(pm,id)
    end
    println("McCormick linearization defined")
end

# Auxiliary cardinality constraints on sources
function constraint_cardinality(pm::_PMD.AbstractUnbalancedPowerModel, k::Int=2)
    for (n, nw_ref) in _PMD.nws(pm)
        cardinality = 0 
        for (i,gen) in ref(pm, n, :gen)
            z_gen = _PMD.var(pm, n, :z_gen, i)[_PMD.ref(pm, n, :gen, i)["index"]]
            cardinality = cardinality + z_gen
        end
        JuMP.@constraint(pm.model, cardinality == k)
    end
end

# Auxiliary constraint y variable to be positive
function constraint_y_variable_positiveness(pm::_PMD.AbstractUnbalancedPowerModel, i::Int, nw::Int=_PMD.nw_id_default)
    switch = ref(pm, nw, :switch)[i]
    f_bus, t_bus = switch["f_bus"], switch["t_bus"]
    #-------------------------------------
    f_idx = (i, f_bus, t_bus)  #key
    psw_f = var(pm, nw, :psw, f_idx)
    #-------------------------------------
    t_idx = (i, t_bus, f_bus)  #opposite_key
    psw_t = var(pm, nw, :psw, t_idx)
    #-------------------------------------
    y = var(pm, nw, :y, i)
    connections = ref(pm, nw, :switch, i)["f_connections"]
    for (idx, c) in enumerate(connections)
        JuMP.@constraint(pm.model, y[i] >= psw_f[c]+psw_t[c])
        JuMP.@constraint(pm.model, y[i] >= -psw_f[c]-psw_t[c])
    end
end

# Auxiliary McCormick relaxation
function mccormick_relaxation(pm::_PMD.AbstractUnbalancedPowerModel, i::Int, nw::Int=_PMD.nw_id_default)
    x_lb, x_ub = 0,1
    y_lb, y_ub = 0,200000 #100000
    y = var(pm, nw, :y, i)[i]
    x = var(pm, nw, :switch_state, i)[i]
    z = var(pm, nw, :z, i)[i]
    #Contraint
    JuMP.@constraint(pm.model, z >= x_lb*y + y_lb*x - x_lb*y_lb)
    JuMP.@constraint(pm.model, z >= x_ub*y + y_ub*x - x_ub*y_ub)
    JuMP.@constraint(pm.model, z <= x_lb*y + y_ub*x - x_lb*y_ub)
    JuMP.@constraint(pm.model, z <= x_ub*y + y_lb*x - x_ub*y_lb)
end

end