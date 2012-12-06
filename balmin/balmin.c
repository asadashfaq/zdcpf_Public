/* Produced by CVXGEN, 2012-03-15 05:13:46 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

// modified by Sarah Becker, University of Frankfurt (Main)

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;
int N = 27; // number of nodes
int L = 44; // number of links

double balmin(double delta[N], double k[L*2], double hmns[L],
  double hpls[L]) {
  int num_iters;
  //printf("%f, %f\n",hmns[5],hpls[5]);
  set_defaults();
  setup_indexing();
  //settings.eps=1.e-2;
  //settings.resid_tol=1.e-2;
  //settings.max_iters=100;
  //settings.kkt_reg=1.e-3;
  load_data(delta, k, hmns, hpls);

  /* Solve problem instance for the record. */
  settings.verbose = 0;
  num_iters = solve();
  if (work.converged != 1)
    printf("Balancing minimization failed to converge!\n");
  return work.optval;
}

void load_data(double delta[N], double k[L*2], double hmns[L],
  double hpls[L]) {
  int i;
  for (i=0; i<N; i++)
    params.Delta[i] = delta[i];
  for (i=0; i<2*L; i++)
    params.K[i] = k[i];
  for (i=0; i<L; i++)
    params.h_mns[i] = hmns[i];
  for (i=0; i<L; i++)
    params.h_pls[i] = hpls[i];
}
