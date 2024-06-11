# Import the package and define the problem to optimize
using Optimization
obj(x, c) = c[1]*x[1]^2 + c[2]*x[2]^2
x0 = [7.0, 2.0]
c = [1.0, 0.5]

prob = OptimizationProblem(obj, x0, c)

# Import a solver package and solve the optimization problem
using OptimizationOptimJL
sol = solve(prob, NelderMead())

# # Import a different solver package and solve the optimization problem a different way
# using OptimizationBBO
# prob = OptimizationProblem(rosenbrock, u0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
# sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())