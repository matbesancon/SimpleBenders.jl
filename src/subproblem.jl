
"""
Data structure regrouping the information for a Benders subproblem
"""
struct SubProblemData
    b::Vector{Float64}
    D::Matrix{Float64}
    A::Matrix{Float64}
    c::Vector{Float64}
end

"""
    DualSubProblem

Benders dual subproblem:
max (b - Dŷ)ᵀ α
s.t. Aᵀα ⩽ c
     α ⩾ 0
"""
struct DualSubProblem
    data::SubProblemData
    α::Vector{VariableRef}
    m::Model
end

function DualSubProblem(d::SubProblemData, m::Model)
    α = @variable(m, α[i = 1:size(d.A, 1)] >= 0)
    @constraint(m, d.A' * α .<= d.c)
    return DualSubProblem(d, α, m)
end

"""
    JuMP.optimize!(sp::DualSubProblem, ŷ)

Solve the dual subproblem, returns a pair `(t::Symbol, α::Vector{Float64})`.
Where `t` is a cut type, either `:OptimalityCut` or `FeasibilityCut`.
"""
function JuMP.optimize!(sp::DualSubProblem, ŷ)
    obj = sp.data.b .- sp.data.D * ŷ
    @objective(sp.m, Max, obj' * sp.α)
    optimize!(sp.m)
    st = termination_status(sp.m)
    if st == MOI.OPTIMAL
        α = JuMP.value.(sp.α)
        return (:OptimalityCut, α)
    elseif (st == MOI.DUAL_INFEASIBLE) || (st==MOI.INFEASIBLE_OR_UNBOUNDED)
        # retrieve extreme ray
        α = MOI.get.(sp.m, MOI.VariablePrimal(), sp.α)
        return (:FeasibilityCut, α)
    else
        error("DualSubProblem error: status $st")
    end
end
