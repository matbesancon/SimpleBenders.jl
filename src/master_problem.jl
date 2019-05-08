
function benders_optimize!(m::Model, y::Vector{VariableRef}, sd::SubProblemData, sp_optimizer, f::Union{Function,Type})
    subproblem = Model(with_optimizer(sp_optimizer))
    dsp = DualSubProblem(sd, subproblem)
    @variable(m, η)
    @objective(m, Min, f(y) + η)
    optimize!(m)
    st = MOI.get(m, MOI.TerminationStatus())
    # restricted master has a solution or is unbounded
    nopt_cons, nfeas_cons = (0, 0)
    @info "Initial status $st"
    cuts = Tuple{Symbol, Vector{Float64}}[]
    while (st == MOI.DUAL_INFEASIBLE) || (st == MOI.OPTIMAL)
        optimize!(m)
        st = MOI.get(m, MOI.TerminationStatus())
        ŷ = JuMP.value.(y)
        η0 = JuMP.value(η)
        (res, α) = optimize!(dsp, ŷ)
        if res == :OptimalityCut
            @info "Optimality cut found"
            if η0 ≥ α' * (dsp.data.b - dsp.data.D * ŷ)
                break
            else
                nopt_cons += 1
                @constraint(m, η ≥ α' * (dsp.data.b - dsp.data.D * y))
            end
        else
            @info "Feasibility cut found"
            nfeas_cons += 1
            @constraint(m, 0 ≥ α' * (dsp.data.b - dsp.data.D * y))
        end
        push!(cuts, (res, α))
    end
    return (m, y, cuts, nopt_cons, nfeas_cons)
end
