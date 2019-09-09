using Test

import SimpleBenders
using JuMP

import Clp
import Cbc

# test from http://www.iems.ucf.edu/qzheng/grpmbr/seminar/Yuping_Intro_to_BendersDecomp.pdf

function test_data()
    c = [2., 3.]
    A = [1 2;2 -1]
    D = zeros(2, 1) .+ [1, 3]
    b = [3, 4]
    return SimpleBenders.SubProblemData(b, D, A, c)
end

function test_result()
    d = test_data()
    m = Model(with_optimizer(Cbc.Optimizer, LogLevel = 0))
    @variable(m, x[1:2] >= 0)
    @variable(m, y[1:1] >= 0)
    @objective(m, Min, d.c'*x + 2y[1])
    @constraint(m, d.A * x + d.D * y .>= d.b)
    optimize!(m)
    return (JuMP.value.(x), JuMP.value.(y), JuMP.objective_value(m))
end

@testset "Basic test" begin
    data = test_data()
    f(v) = 2v[1]
    m = Model(with_optimizer(Cbc.Optimizer, LogLevel = 0))
    @variable(m, y[j=1:1] >= 0)
    (m, y, cuts, nopt_cons, nfeas_cons) = SimpleBenders.benders_optimize!(m, y, data, () -> Clp.Optimizer(LogLevel = 0), f)
    (xref, yref, objref) = test_result()
    @test yref[1] ≈ JuMP.value(y[1])
    @test objref ≈ JuMP.objective_value(m)
end
