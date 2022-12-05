// (C) Copyright 2014 by Autodesk, Inc.

#include <Base/Utils/BaseErrorInc.hh>

DEFINE_ERROR(UNSPECIFIED_COMISO_FAILURE, "Unspecified CoMISo Failure")
DEFINE_ERROR(UNSPECIFIED_EIGEN_FAILURE, "Unspecified Eigen Failure")

// Misc 3P Solver issues
DEFINE_ERROR(UNSPECIFIED_GUROBI_EXCEPTION, "Unspecified Gurobi Exception")
DEFINE_ERROR(GUROBI_LICENCE_ABSENT, "Gurobi license absent")
DEFINE_ERROR(GUROBI_LICENCE_MODEL_TOO_LARGE, "Model too large for current Gurobi license")

DEFINE_ERROR(
    QP_INITIALIZATION_FAILED, "Quadratic program initialization failed")
DEFINE_ERROR(QP_OPTIMIZATION_FAILED, "Quadratic program optimization failed")
DEFINE_ERROR(QP_MAXIMUM_ITERATIONS_EXCEEDED,
    "Quadratic program maximum iterations exceeded")

// TODO: Obsolete: remove when we can
DEFINE_ERROR(IPOPT_STEP_FAILURE, "IPOPT step failure") 
DEFINE_ERROR(IPOPT_UNSPECIFIED_FAILURE, "IPOPT Unspecified failure")
DEFINE_ERROR(UNSPECIFIED_CBC_EXCEPTION, "Unspecified CBC Exception")

DEFINE_ERROR(MIPS_NO_SOLUTION, "Mixed integer problem solver cannot find a solution")
DEFINE_ERROR(MIPS_SOLUTION_INCORRECT, "Mixed integer problem solution incorrect")

// DOCloud related message (TODO: Obsolete, remove when we can)
DEFINE_ERROR(DOCLOUD_REQUEST_INIT_FAILED, 
  "DO cloud request intialization failure [defect]")
DEFINE_ERROR(DOCLOUD_REQUEST_EXEC_FAILED, 
  "DO cloud request execution failed, [network problem]")
DEFINE_ERROR(DOCLOUD_CONFIG_SET_VALUE_INVALID, 
  "DO cloud configuration set value is invalid")
DEFINE_ERROR(DOCLOUD_JOB_DATA_INVALID, 
  "DO cloud job data invalid [defect]")
DEFINE_ERROR(DOCLOUD_SUBSCRIPTION_LIMIT, 
  "DO cloud job subscription limit exceeded, check account")
DEFINE_ERROR(DOCLOUD_JOB_NOT_FOUND, "DO cloud job not found [defect]")
DEFINE_ERROR(DOCLOUD_JOB_UNRECOGNIZED_FAILURE, 
  "DO cloud job unrecognized failure [defect]")
DEFINE_ERROR(DOCLOUD_JOB_HTTP_CODE_NOT_FOUND, 
  "DO cloud job http status code not found [defect | network problem]")
DEFINE_ERROR(DOCLOUD_JOB_LOCATION_NOT_FOUND, 
  "DO cloud job location not found [defect | network problem]")
DEFINE_ERROR(DOCLOUD_CPLEX_SOLUTION_MISMATCH, 
  "DO cloud CPLEX solution mismatch [defect | network problem]")

// LSQ constrained solver related messages
DEFINE_ERROR(LSQC_UNEXPECTED_VARIABLE, "A variable is not in the expected set")
DEFINE_ERROR(LSQC_SINGULAR, "LSQ constrained problem resulted in a singular system")
DEFINE_ERROR(LSQC_INFEASIBLE, "LSQ problem does not have solution")

// DEFINE_ERROR(DOCLOUD_ROOT_URL_INVALID, "DO cloud root URL invalid")
//DEFINE_ERROR(DOCLOUD_API_KEY_INVALID, "DO cloud API key invalid")
//DEFINE_ERROR(DOCLOUD_INFEASIBLE_TIMEOUT_INVALID, 
//  "DO cloud infeasible timeout invalid")
//DEFINE_ERROR(DOCLOUD_FEASIBLE_TIMEOUT_INVALID, 
//  "DO cloud feasible timeout invalid")
//
