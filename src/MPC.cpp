#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
const double dt = 0.1;
const double Lf = 2.67;
const double v_ref = 120;

size_t t_x = 0;
size_t t_y = t_x + N;
size_t t_psi = t_y + N;
size_t t_v = t_psi + N;
size_t t_cte = t_v + N;
size_t t_epsi = t_cte + N;
size_t t_delta = t_epsi + N;
size_t t_a = t_delta + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
      fg[0] = 0;
      for (int i = 0; i < N; i++) {
          fg[0] += 3000 * CppAD::pow(vars[t_cte + i], 2);
          fg[0] += 3000 * CppAD::pow(vars[t_epsi + i], 2);
          fg[0] += 1 * CppAD::pow(vars[t_v + i] - v_ref, 2);
      }
      
      for (int i = 0; i < N - 1; i++) {
          fg[0] += 5 * CppAD::pow(vars[t_delta + i], 2);
          fg[0] += 5 * CppAD::pow(vars[t_a + i], 2);
          fg[0] += 700*CppAD::pow(vars[t_delta + i] * vars[t_v + i], 2);
      }
      
      for (int i = 0; i < N - 2; i++) {
          fg[0] += 200 * CppAD::pow(vars[t_delta + i + 1] - vars[t_delta + i], 2);
          fg[0] += 10 * CppAD::pow(vars[t_a + i + 1] - vars[t_a + i], 2);
      }
      
      // Initial constraints
      fg[1 + t_x] = vars[t_x];
      fg[1 + t_y] = vars[t_y];
      fg[1 + t_psi] = vars[t_psi];
      fg[1 + t_v] = vars[t_v];
      fg[1 + t_cte] = vars[t_cte];
      fg[1 + t_epsi] = vars[t_epsi];
      
      // The rest of the constraints
      for (int i = 1; i < N; i++) {
          // The state at time t+1 .
          AD<double> x1 = vars[t_x + i];
          AD<double> y1 = vars[t_y + i];
          AD<double> psi1 = vars[t_psi + i];
          AD<double> v1 = vars[t_v + i];
          AD<double> cte1 = vars[t_cte + i];
          AD<double> epsi1 = vars[t_epsi + i];
          
          // The state at time t.
          AD<double> x0 = vars[t_x + i - 1];
          AD<double> y0 = vars[t_y + i - 1];
          AD<double> psi0 = vars[t_psi + i - 1];
          AD<double> v0 = vars[t_v + i - 1];
          AD<double> cte0 = vars[t_cte + i - 1];
          AD<double> epsi0 = vars[t_epsi + i - 1];
          
          // Only consider the actuation at time t.
          AD<double> delta0 = vars[t_delta + i - 1];
          AD<double> a0 = vars[t_a + i - 1];
          
          if (i > 1) {   // use previous actuations (to account for latency)
              a0 = vars[t_a + i - 2];
              delta0 = vars[t_delta + i - 2];
          }
          AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
          AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));
          
          
          fg[1 + t_x + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
          fg[1 + t_y + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
          fg[1 + t_psi + i] = psi1 - (psi0 + v0 / Lf * delta0 * dt);
          fg[1 + t_v + i] = v1 - (v0 + a0 * dt);
          fg[1 + t_cte + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
          fg[1 + t_epsi + i] = epsi1 - ((psi0 - psides0) - v0 / Lf * delta0 * dt);
      }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
    bool ok = true;
    size_t i;
    typedef CPPAD_TESTVECTOR(double) Dvector;
    const double x = state[0];
    const double y = state[1];
    const double psi = state[2];
    const double v = state[3];
    const double cte = state[4];
    const double epsi = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
    // Set the initial variable values
    vars[t_x] = x;
    vars[t_y] = y;
    vars[t_psi] = psi;
    vars[t_v] = v;
    vars[t_cte] = cte;
    vars[t_epsi] = epsi;
    
    // Non Actuator limits
    for (int i = 0; i < t_delta; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }
    // The upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians).
    for (int i = t_delta; i < t_a; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }
    // Acceleration/decceleration upper and lower limits
    for (int i = t_a; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
    constraints_lowerbound[t_x] = x;
    constraints_lowerbound[t_y] = y;
    constraints_lowerbound[t_psi] = psi;
    constraints_lowerbound[t_v] = v;
    constraints_lowerbound[t_cte] = cte;
    constraints_lowerbound[t_epsi] = epsi;
    constraints_upperbound[t_x] = x;
    constraints_upperbound[t_y] = y;
    constraints_upperbound[t_psi] = psi;
    constraints_upperbound[t_v] = v;
    constraints_upperbound[t_cte] = cte;
    constraints_upperbound[t_epsi] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
    /*this->mpc_x = {};
    this->mpc_y = {};
    for (int i = 0; i < N; i++) {
        this->mpc_x.push_back(solution.x[t_x + i]);
        this->mpc_y.push_back(solution.x[t_y + i]);
    }
    vector<double> result;
    result.push_back(solution.x[t_delta]);
    result.push_back(solution.x[t_a]);
    return result;*/
    vector<double> result;
    
    result.push_back(solution.x[t_delta]);
    result.push_back(solution.x[t_a]);
    
    for (int i = 0; i < N; i++) {
        result.push_back(solution.x[t_x + i]);
        result.push_back(solution.x[t_y + i]);
    }
    
    return result;
}
