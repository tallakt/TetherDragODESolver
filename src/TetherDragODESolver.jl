module TetherDragODESolver

export solve_tether, solve_tether_for_non_trpt, calculate_non_drag_trpt_moment,
calculate_non_drag_lambda, solved_trpt_moment_at_radius1,
solved_trpt_moment_at_radius0, solved_lambda, solved_moment_loss,
solved_efficiency_ratio, solved_drag_coefficient_multiplier

using DifferentialEquations
using LinearAlgebra
using Optim
using Rotations


struct TetherDragODESolution
  solution::SciMLBase.AbstractODESolution
  optimization::Optim.OptimizationResults
  rotational_speed::Float64
  tension::Float64
  radius0::Float64
  radius1::Float64
  length::Float64
  twist::Float64
  tether_diameter::Float64
  tether_drag_coefficient::Float64
  mu::Float64
  rho::Float64
end


"""
    solve_tether(rotational_speed, tension, radius0, radius1, length, twist, tether_diameter)

Performs two stages calculation of the accurate tether shape as defined by the
differential equations related to tether drag and centrifugal forces.

Optional parameters:

- `tether_drag_coefficient`: drag coefficient of the tether, when moving sideways through air. Default value `1.0`
- `mu`: mass per length of tether. Defgault value calculated from tether diameter
- `rho`: density of air. Default value `1.225`
- `reltolode`: relative tolerance of the ODE solver. Default value `1e-10`

The algorithm will first guess the initial derivatives of the tether position
curve vector function, then calculate the shape of the tether. Based on where
the other end of the tether ends up and where it should be according to
`radius1` and `twist`, the initial derivaties are refined by an optimizer.
After the optimizer has completed, the tether shape is recalculated for the
optimal derivatives.

Note the `rotational_speed` and `twist` values are given in radians per second
and radians.
"""
function solve_tether(rotational_speed::Float64, tension::Float64, radius0::Float64, radius1::Float64, length::Float64, twist::Float64, tether_diameter::Float64; tether_drag_coefficient::Float64 = 1.0, mu::Float64 = 558.8 * tether_diameter^2, rho::Float64 = 1.225, reltolode::Float64 = 1e-10)::TetherDragODESolution
  # rotational_speed - rotational speed of rig in rad/s
  # tension - shaft tension along centerline
  # radius0 - radius at cartwheel
  # radius1 - radius at kite
  # length - tether length
  # twist - shaft twist in radians
  # tether_diameter - tether diameter
  # tether_drag_coefficient - drag coefficient of tether
  # mu - mass per meter of tether

  tmp = mu * rotational_speed^2 / tension

  function tether_belly_diff!(du, u, p, t)
    x, dx, y, dy, z = u
    dz = -sqrt(max(0.0, 1 - dx^2 - dy^2))

    centerline_to_kite = [x, y, 0]
    centerline_to_kite_n = normalize(centerline_to_kite)
    rotate_270_deg_z = [0 1 0; -1 0 0; 0 0 1]
    centrifugal_force = mu * rotational_speed^2 * centerline_to_kite
    drag_direction_n = rotate_270_deg_z * centerline_to_kite_n
    drag_force = 1 // 2 * rho * rotational_speed^2 * norm(centerline_to_kite)^2 * tether_drag_coefficient * tether_diameter * norm(cross([dx, dy, dz], drag_direction_n)) * drag_direction_n
    forces = drag_force + centrifugal_force
    forces_normal_to_tether = cross([dx, dy, dz], cross([dx, dy, dz], forces))

    # first equation
    A0 = cross([dx, dy, dz], forces)'
    b0 = 0

    # second equation
    A1 = forces_normal_to_tether'
    b1 = 1 / tension * abs(dz) * dot(forces_normal_to_tether, forces_normal_to_tether)

    # third equation
    A2 = [dx dy dz]
    b2 = 0

    # Solve set of equations
    A = [A0; A1; A2]
    b = [b0; b1; b2]

    double_deriv_position = inv(A) * b;


    du[1] = dx
    du[2] = double_deriv_position[1]
    du[3] = dy
    du[4] = double_deriv_position[2]
    du[5] = dz
    # we ignore double deriv position in z direction, it is not useful
  end

  x0 = float(max(radius1 / 10_000.0, radius0))
  y0 = 0.0
  x1 = radius1 * cos(twist)
  y1 = radius1 * sin(twist)
  z0 = 0.0

  solve_with_initial_derivatives = function(dx0, dy0)
    u0 = [x0, dx0, y0, dy0, z0]
    prob = ODEProblem(tether_belly_diff!, u0, (0, length))
    solve(prob, AutoTsit5(Rosenbrock23()), reltol = reltolode)
  end

  minimizer = function(x)
    (dx0, dy0) = x
    sol = solve_with_initial_derivatives(dx0, dy0)
    (x1s, _, y1s, _, _) = sol.u[end]
    (x1 - x1s)^2 + (y1 - y1s)^2
  end

  opt = Optim.optimize(minimizer, normalize([x1 - x0, y1 - 0.0]))
  (dx0, dy0) = opt.minimizer

  # solve_with_initial_derivatives(dx0, dy0)
  (solve_with_initial_derivatives(dx0, dy0), opt)

  TetherDragODESolution(
                        solve_with_initial_derivatives(dx0, dy0)
                        , opt
                        , rotational_speed
                        , tension
                        , radius0
                        , radius1
                        , length
                        , twist
                        , tether_diameter
                        , tether_drag_coefficient
                        , mu
                        , rho
                       )
end


"""
    solve_tether_for_non_trpt(rotational_speed, tension, radius, length, tether_diameter)

This is a wrapper function for `solve_tether` used for non-TRPT rigs. The
`radius` parameter maps to `radius1` while `radius0` and `ttwist` are assigned
zero values.
"""
function solve_tether_for_non_trpt(rotational_speed::Float64, tension::Float64, radius::Float64, length::Float64, tether_diameter::Float64; tether_drag_coefficient::Float64 = 1.0, mu::Float64 = 558.8 * tether_diameter^2, rho::Float64 = 1.225, reltolode::Float64 = 1e-10)::TetherDragODESolution
  solve_tether(rotational_speed, tension, 0.0, radius, length, 0.0, tether_diameter, tether_drag_coefficient = tether_drag_coefficient, mu = mu, rho = rho, reltolode = reltolode)
end


"""
    calculate_non_drag_trpt_moment(tension, radius0, radius1, length, twist)

Calculated the moment that the TRPT shaft can transfer if the tether was
completely straight, not considering drag or centrifucal forces.
"""
function calculate_non_drag_trpt_moment(tension::Float64, radius0::Float64, radius1::Float64, length::Float64, twist::Float64)
  calculate_non_drag_lambda(tension, radius0, radius1, length, twist) * tension * radius1
end


"""
    calculate_non_drag_lambda(tension, radius0, radius1, length, twist)

Calculated the Lambda value that the TRPT shaft can transfer if the tether was
completely straight, not considering drag or centrifucal forces.

The Lambda value is defines as moment per tension and radius, or the force -
force ratio.
"""
function calculate_non_drag_lambda(tension::Float64, radius0::Float64, radius1::Float64, length::Float64, twist::Float64)
  # radius0 may be zero for yoyo, so we use radius1 as reference, even though the moment
  # at either end must be the same
  x0 = [radius0, 0.0, 0.0]
  x1_2d = radius1 .* [cos(twist), sin(twist)]
  x1 = vcat(x1_2d, x0[3] - sqrt(length^2 - x0[1]^2 - x0[2]^2))
  dot(RotZ(pi / 2 + twist) * [1.0, 0.0, 0.0], normalize(x1 .- x0))
end


"""
    solved_trpt_moment_at_radius1(ode_sol)

Calculates the moment applied to the kite by the tether, based on a solved
tether shape.

The `ode_sol` is found using either `solve_tether` for a TRPT style plant, or
`solve_tether_for_non_trpt` for rigs typically not using a torque transfer to
transfer power.
"""
function solved_trpt_moment_at_radius1(ode_sol::TetherDragODESolution)
  (_, dx, _, dy, _) = ode_sol.solution.u[end]
  ode_sol.tension * ode_sol.radius1 * dot(RotZ(pi / 2 + ode_sol.twist) * [1.0, 0.0, 0.0], [dx, dy, 0])
end


"""
    solved_trpt_moment_at_radius0(ode_sol)

Calculates the moment applied to the ground cartwheel by the tether, based on a
solved tether shape.

The `ode_sol` is found using either `solve_tether` for a TRPT style plant, or
`solve_tether_for_non_trpt` for rigs typically not using a torque transfer to
transfer power.
"""
function solved_trpt_moment_at_radius0(ode_sol::TetherDragODESolution)
  (_, dx0, _, dy0, _) = ode_sol.solution.u[1]
  ode_sol.tension * ode_sol.radius0 * dot(RotZ(pi / 2) * [1.0, 0.0, 0.0], [dx0, dy0, 0])
end


"""
    solved_lambda(ode_sol)

Calculates the usable Lambda value of the shaft, considering a curved, solved
tether shape. The Lambda value represents the moment per tether and tension or
the force - force ratio. It is calculated by looking at the moment applied to
the ground cartwheel. Any drag seen by the kite is assumed to represent
additional drag at the kite, normally.

The `ode_sol` is found using either `solve_tether` for a TRPT style plant, or
`solve_tether_for_non_trpt` for rigs typically not using a torque transfer to
transfer power.
"""
function solved_lambda(ode_sol::TetherDragODESolution)
  (_, dx0, _, dy0, _) = ode_sol.solution.u[1]
  ode_sol.radius0 / ode_sol.radius1 * dot(RotZ(pi / 2) * [1.0, 0.0, 0.0], [dx0, dy0, 0])
end


"""
    solved_moment_loss(ode_sol)

Calculates the moment lost to tether drag between the kite the ground cartwheel. 

The `ode_sol` is found using either `solve_tether` for a TRPT style plant, or
`solve_tether_for_non_trpt` for rigs typically not using a torque transfer to
transfer power.
"""
function solved_moment_loss(ode_sol::TetherDragODESolution)
  solved_trpt_moment_at_radius1(ode_sol) - solved_trpt_moment_at_radius0(ode_sol)
end


"""
    solved_efficiency_ratio(ode_sol)

Calculates the efficiency of the rotary shaft, being the ratio of moment at the
ground cartwheel to the moment transferred at the ground cartwheel. This also
maps directly to power efficiency given a rotational speed of the shaft.

The `ode_sol` is found using either `solve_tether` for a TRPT style plant, or
`solve_tether_for_non_trpt` for rigs typically not using a torque transfer to
transfer power.
"""
function solved_efficiency_ratio(ode_sol::TetherDragODESolution)
  m0 = solved_trpt_moment_at_radius0(ode_sol)
  m1 = solved_trpt_moment_at_radius1(ode_sol)
  m0 / m1 
end


"""
    solved_drag_coefficient_multiplier(ode_sol)

Computes a ratio between the drag of a tether as if the tether was all moving
at the speed of the kite, and the actual tether drag considering a more exact
ODE solution to the tether path as defined by the tether drag and tether
centrifugal forces.

The `ode_sol` is found using either `solve_tether` for a TRPT style plant, or
`solve_tether_for_non_trpt` for rigs typically not using a torque transfer to
transfer power.

For the latter, a value of 1/4 should be expected unless tether tension is very
low.

This function is useful to calculate how much tether drag should be assignes to
the kite to make a simplified model.
"""
function solved_drag_coefficient_multiplier(ode_sol::TetherDragODESolution)
  m0 = solved_trpt_moment_at_radius0(ode_sol)
  m1 = solved_trpt_moment_at_radius1(ode_sol)
  drag_moment = m1 - m0
  drag_force = drag_moment / ode_sol.radius1
  airspeed = ode_sol.rotational_speed * ode_sol.radius1
  drag_100_percent = 1 // 2 * ode_sol.rho * airspeed^2 * ode_sol.tether_drag_coefficient * ode_sol.length * ode_sol.tether_diameter
  drag_force / drag_100_percent
end

end
