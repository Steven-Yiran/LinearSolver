
#ifndef FLUX_POISSON_H_
#define FLUX_POISSON_H_

#include "spmat.h"

namespace flux {

  void get_model_problem(int n , spmat<double>& A , vecd<double>& f, int IC_TYPE);
  vecd<double> solve_model_problem(int n , double tol, int max_iter, int IC_TYPE);
  vecd<double> solve_model_problem_cg(int n, int IC_TYPE);
  
} // flux

#endif


