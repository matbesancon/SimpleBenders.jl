
"""
Data structure grouping the information for a Benders subproblem.
"""
struct SubProblemData{BT<:AbstractVector, DT<:AbstractMatrix, AT<:AbstractMatrix, CT<:AbstractVector}
    b::BT
    D::DT
    A::AT
    c::CT

    SubProblemData(b::BT, D::DT, A::AT, c::CT) where {BT, DT, AT, CT} = new{BT, DT, AT, CT}(b, D, A, c)
end

"""
    DualSubProblem

Benders dual subproblem:
max (b - Dŷ)ᵀ α
s.t. Aᵀα ⩽ c
     α ⩾ 0
"""
struct DualSubProblem{SPD <: SubProblemData}
    data::SPD
    α::Vector{VariableRef}
    m::Model

    DualSubProblem(data::SPD, α, m) where {SPD} = new{SPD}(data, α, m)
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
    elseif st == MOI.DUAL_INFEASIBLE
        return (:FeasibilityCut, α)
    else
        error("DualSubProblem error: status $status")
    end
end
