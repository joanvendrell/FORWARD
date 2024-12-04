"""
######################## Custom Warm-Up #########################
This file implements warm up functions for custom problems
"""
module CustomWarmUp

using PowerModelsDistribution
import JuMP
const _PMD = PowerModelsDistribution

function custom_warmstart_transmission_loss(pm::_PMD.AbstractUnbalancedPowerModel, data::Dict)
    for (nw, nw_ref) in _PMD.nws(pm)
        for key in keys(_PMD.var(pm,nw))
            if string(key) == "psw"
                for id in _PMD.ids(pm, :switch)
                    switch = ref(pm, nw, :switch)[id]
                    f_bus, t_bus = switch["f_bus"], switch["t_bus"]
                    psw_f = var(pm, nw, :psw, (id, f_bus, t_bus))[1]
                    psw_t = var(pm, nw, :psw, (id, t_bus, f_bus))[1]
                    JuMP.set_start_value(psw_f,data["solution"]["switch"][string(id)]["pf"][1])
                    JuMP.set_start_value(psw_t,data["solution"]["switch"][string(id)]["pt"][1])
                end
            elseif string(key) == "qsw"
                for id in _PMD.ids(pm, :switch)
                    switch = ref(pm, nw, :switch)[id]
                    f_bus, t_bus = switch["f_bus"], switch["t_bus"]
                    qsw_f = var(pm, nw, :qsw, (id, f_bus, t_bus))[1]
                    qsw_t = var(pm, nw, :qsw, (id, t_bus, f_bus))[1]
                    JuMP.set_start_value(qsw_f,data["solution"]["switch"][string(id)]["qf"][1])
                    JuMP.set_start_value(qsw_t,data["solution"]["switch"][string(id)]["qt"][1])
                end
            elseif string(key) == "p"
                for id in _PMD.ids(pm, :branch)
                    branch = ref(pm, nw, :branch)[id]
                    f_bus, t_bus = branch["f_bus"], branch["t_bus"]
                    p_f = var(pm, nw, :p, (id, f_bus, t_bus))[1]
                    p_t = var(pm, nw, :p, (id, t_bus, f_bus))[1]
                    JuMP.set_start_value(p_f,data["solution"]["branch"][string(id)]["pf"][1])
                    JuMP.set_start_value(p_t,data["solution"]["branch"][string(id)]["pt"][1])
                end
            elseif string(key) == "q"
                for id in _PMD.ids(pm, :branch)
                    branch = ref(pm, nw, :branch)[id]
                    f_bus, t_bus = branch["f_bus"], branch["t_bus"]
                    q_f = var(pm, nw, :q, (id, f_bus, t_bus))[1]
                    q_t = var(pm, nw, :q, (id, t_bus, f_bus))[1]
                    JuMP.set_start_value(q_f,data["solution"]["branch"][string(id)]["qf"][1])
                    JuMP.set_start_value(q_t,data["solution"]["branch"][string(id)]["qt"][1])
                end
            elseif string(key) == "qd_bus"
                continue 
            elseif string(key) == "pd_bus"
                continue 
            elseif string(key) == "z_gen"
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        m = match(r"\[(\d+)\]", string(variable))
                        JuMP.set_start_value(variable, 
                            Int(round(data["solution"]["gen"][m.match[2:length(m.match)-1]]["gen_status"][1])))
                    end
                end
            elseif string(key) == "z_demand"
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        m = match(r"\[(\d+)\]", string(variable))
                        JuMP.set_start_value(variable, data["solution"]["load"][m.match[2:length(m.match)-1]]["status"][1])
                    end
                end
            elseif string(key) == "switch_state"
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        m = match(r"\[(\d+)\]", string(variable))
                        switch = ref(pm, nw, :switch)[parse(Int64, m.match[2:length(m.match)-1])]
                        f_bus, t_bus = switch["f_bus"], switch["t_bus"]
                        JuMP.set_start_value(variable, 
                            Int(round(data["solution"]["switch"][m.match[2:length(m.match)-1]]["state"][1])))
                    end
                end
            elseif haskey(data["solution"]["gen"]["1"],string(key))
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        JuMP.set_start_value(variable, data["solution"]["gen"][string(id)][string(key)][1])
                    end
                end
            elseif haskey(data["solution"]["switch"]["1"],string(key))
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        JuMP.set_start_value(variable, data["solution"]["switch"][string(id)][string(key)][1])
                    end
                end
            elseif haskey(data["solution"]["load"]["1"],string(key))
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        JuMP.set_start_value(variable, data["solution"]["load"][string(id)][string(key)][1])
                    end
                end
            elseif haskey(data["solution"]["bus"]["1"],string(key))
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        JuMP.set_start_value(variable, data["solution"]["bus"][string(id)][string(key)][1])
                    end
                end
            else
                for id in keys(_PMD.var(pm,nw,key))
                    for variable in _PMD.var(pm,nw,key,id)
                        JuMP.set_start_value(variable, 0)
                    end
                end
            end
        end
    end
end

# Auxiliary function to print initial values
function check_start_values(pm::_PMD.AbstractUnbalancedPowerModel)
    for (nw, nw_ref) in _PMD.nws(pm)
        for key in keys(_PMD.var(pm,nw))
            for id in keys(_PMD.var(pm,nw,key))
                for variable in _PMD.var(pm,nw,key,id)
                    println("Using ",variable," -> ",JuMP.start_value(variable))
                end
            end
        end
    end
end

end