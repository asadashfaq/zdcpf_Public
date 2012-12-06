/* Produced by CVXGEN, 2012-03-15 08:24:15 -0700.  */
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

int flowmin(double delta[N], double k[L*2], double hmns[L],
		double hpls[L], double bmin, double flow[L]) {
  int num_iters;
  set_defaults();
  setup_indexing();
  load_data(delta, k, hmns, hpls, bmin);

  /* Solve problem instance for the record. */
  settings.verbose = 0;
  num_iters = solve();
  if (work.converged != 1)
    printf("Flow minimization failed to converge!\n");

  int i=0;
  for (i=0; i<L; i++)
    flow[i]=vars.F[i]; // return flow by reference
  return 0;
}

void load_data(double delta[N], double k[L*2], double hmns[L],
  double hpls[L], double bmin) {
  int i;
  for (i=0; i<N; i++)
    params.Delta[i] = delta[i];
  for (i=0; i<2*L; i++)
    params.K[i] = k[i];
  for (i=0; i<L; i++)
    params.h_mns[i] = hmns[i];
  for (i=0; i<L; i++)
    params.h_pls[i] = hpls[i];
  params.B_min[0] = bmin;
}
