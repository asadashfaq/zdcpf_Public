/* Produced by CVXGEN, 2012-03-15 08:24:15 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: csolve.c. */
/* Description: mex-able file for running cvxgen solver. */

#include "mex.h"
#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int i, j;
  mxArray *xm, *cell, *xm_cell;
  double *src;
  double *dest;
  double *dest_cell;
  int valid_vars;
  int steps;
  int this_var_errors;
  int warned_diags;
  int prepare_for_c = 0;
  int extra_solves;
  const char *status_names[] = {"optval", "gap", "steps", "converged"};
  mwSize dims1x1of1[1] = {1};
  mwSize dims[1];

  const char *var_names[] = {"F"};
  const int num_var_names = 1;

  /* Avoid compiler warnings of unused variables by using a dummy assignment. */
  warned_diags = j = 0;
  extra_solves = 0;

  set_defaults();

  /* Check we got the right number of arguments. */
  if (nrhs == 0)
    mexErrMsgTxt("Not enough arguments: You need to specify at least the parameters.\n");

  if (nrhs > 1) {
    /* Assume that the second argument is the settings. */
    if (mxGetField(prhs[1], 0, "eps") != NULL)
      settings.eps = *mxGetPr(mxGetField(prhs[1], 0, "eps"));

    if (mxGetField(prhs[1], 0, "max_iters") != NULL)
      settings.max_iters = *mxGetPr(mxGetField(prhs[1], 0, "max_iters"));

    if (mxGetField(prhs[1], 0, "refine_steps") != NULL)
      settings.refine_steps = *mxGetPr(mxGetField(prhs[1], 0, "refine_steps"));

    if (mxGetField(prhs[1], 0, "verbose") != NULL)
      settings.verbose = *mxGetPr(mxGetField(prhs[1], 0, "verbose"));

    if (mxGetField(prhs[1], 0, "better_start") != NULL)
      settings.better_start = *mxGetPr(mxGetField(prhs[1], 0, "better_start"));

    if (mxGetField(prhs[1], 0, "verbose_refinement") != NULL)
      settings.verbose_refinement = *mxGetPr(mxGetField(prhs[1], 0,
            "verbose_refinement"));

    if (mxGetField(prhs[1], 0, "debug") != NULL)
      settings.debug = *mxGetPr(mxGetField(prhs[1], 0, "debug"));

    if (mxGetField(prhs[1], 0, "kkt_reg") != NULL)
      settings.kkt_reg = *mxGetPr(mxGetField(prhs[1], 0, "kkt_reg"));

    if (mxGetField(prhs[1], 0, "s_init") != NULL)
      settings.s_init = *mxGetPr(mxGetField(prhs[1], 0, "s_init"));

    if (mxGetField(prhs[1], 0, "z_init") != NULL)
      settings.z_init = *mxGetPr(mxGetField(prhs[1], 0, "z_init"));

    if (mxGetField(prhs[1], 0, "resid_tol") != NULL)
      settings.resid_tol = *mxGetPr(mxGetField(prhs[1], 0, "resid_tol"));

    if (mxGetField(prhs[1], 0, "extra_solves") != NULL)
      extra_solves = *mxGetPr(mxGetField(prhs[1], 0, "extra_solves"));
    else
      extra_solves = 0;

    if (mxGetField(prhs[1], 0, "prepare_for_c") != NULL)
      prepare_for_c = *mxGetPr(mxGetField(prhs[1], 0, "prepare_for_c"));
  }

  valid_vars = 0;

  this_var_errors = 0;

  xm = mxGetField(prhs[0], 0, "B_min");

  if (xm == NULL) {
    printf("could not find params.B_min.\n");
  } else {
    if (!((mxGetM(xm) == 1) && (mxGetN(xm) == 1))) {
      printf("B_min must be size (1,1), not (%d,%d).\n", mxGetM(xm), mxGetN(xm));
      this_var_errors++;
    }

    if (mxIsComplex(xm)) {
      printf("parameter B_min must be real.\n");
      this_var_errors++;
    }

    if (!mxIsClass(xm, "double")) {
      printf("parameter B_min must be a full matrix of doubles.\n");
      this_var_errors++;
    }

    if (mxIsSparse(xm)) {
      printf("parameter B_min must be a full matrix.\n");
      this_var_errors++;
    }

    if (this_var_errors == 0) {
      dest = params.B_min;
      src = mxGetPr(xm);

      for (i = 0; i < 1; i++)
        *dest++ = *src++;

      valid_vars++;
    }
  }

  this_var_errors = 0;

  xm = mxGetField(prhs[0], 0, "Delta");

  if (xm == NULL) {
    printf("could not find params.Delta.\n");
  } else {
    if (!((mxGetM(xm) == 27) && (mxGetN(xm) == 1))) {
      printf("Delta must be size (27,1), not (%d,%d).\n", mxGetM(xm), mxGetN(xm));
      this_var_errors++;
    }

    if (mxIsComplex(xm)) {
      printf("parameter Delta must be real.\n");
      this_var_errors++;
    }

    if (!mxIsClass(xm, "double")) {
      printf("parameter Delta must be a full matrix of doubles.\n");
      this_var_errors++;
    }

    if (mxIsSparse(xm)) {
      printf("parameter Delta must be a full matrix.\n");
      this_var_errors++;
    }

    if (this_var_errors == 0) {
      dest = params.Delta;
      src = mxGetPr(xm);

      for (i = 0; i < 27; i++)
        *dest++ = *src++;

      valid_vars++;
    }
  }

  this_var_errors = 0;

  xm = mxGetField(prhs[0], 0, "K");

  if (xm == NULL) {
    printf("could not find params.K.\n");
  } else {
    if (!((mxGetM(xm) == 27) && (mxGetN(xm) == 44))) {
      printf("K must be size (27,44), not (%d,%d).\n", mxGetM(xm), mxGetN(xm));
      this_var_errors++;
    }

    if (mxIsComplex(xm)) {
      printf("parameter K must be real.\n");
      this_var_errors++;
    }

    if (!mxIsClass(xm, "double")) {
      printf("parameter K must be a full matrix of doubles.\n");
      this_var_errors++;
    }

    if (mxIsSparse(xm)) {
      printf("parameter K must be a full matrix.\n");
      this_var_errors++;
    }

    if (this_var_errors == 0) {
      dest = params.K;
      src = mxGetPr(xm);

      dest[0] = src[0];  /* (1,1) entry. */
      dest[1] = src[12];  /* (13,1) entry. */
      dest[2] = src[27];  /* (1,2) entry. */
      dest[3] = src[42];  /* (16,2) entry. */
      dest[4] = src[54];  /* (1,3) entry. */
      dest[5] = src[70];  /* (17,3) entry. */
      dest[6] = src[81];  /* (1,4) entry. */
      dest[7] = src[99];  /* (19,4) entry. */
      dest[8] = src[108];  /* (1,5) entry. */
      dest[9] = src[130];  /* (23,5) entry. */
      dest[10] = src[135];  /* (1,6) entry. */
      dest[11] = src[158];  /* (24,6) entry. */
      dest[12] = src[163];  /* (2,7) entry. */
      dest[13] = src[182];  /* (21,7) entry. */
      dest[14] = src[191];  /* (3,8) entry. */
      dest[15] = src[195];  /* (7,8) entry. */
      dest[16] = src[218];  /* (3,9) entry. */
      dest[17] = src[234];  /* (19,9) entry. */
      dest[18] = src[246];  /* (4,10) entry. */
      dest[19] = src[256];  /* (14,10) entry. */
      dest[20] = src[273];  /* (4,11) entry. */
      dest[21] = src[287];  /* (18,11) entry. */
      dest[22] = src[301];  /* (5,12) entry. */
      dest[23] = src[303];  /* (7,12) entry. */
      dest[24] = src[328];  /* (5,13) entry. */
      dest[25] = src[331];  /* (8,13) entry. */
      dest[26] = src[355];  /* (5,14) entry. */
      dest[27] = src[363];  /* (13,14) entry. */
      dest[28] = src[382];  /* (5,15) entry. */
      dest[29] = src[396];  /* (19,15) entry. */
      dest[30] = src[409];  /* (5,16) entry. */
      dest[31] = src[427];  /* (23,16) entry. */
      dest[32] = src[436];  /* (5,17) entry. */
      dest[33] = src[456];  /* (25,17) entry. */
      dest[34] = src[464];  /* (6,18) entry. */
      dest[35] = src[479];  /* (21,18) entry. */
      dest[36] = src[491];  /* (6,19) entry. */
      dest[37] = src[507];  /* (22,19) entry. */
      dest[38] = src[520];  /* (8,20) entry. */
      dest[39] = src[532];  /* (20,20) entry. */
      dest[40] = src[548];  /* (9,21) entry. */
      dest[41] = src[555];  /* (16,21) entry. */
      dest[42] = src[575];  /* (9,22) entry. */
      dest[43] = src[585];  /* (19,22) entry. */
      dest[44] = src[602];  /* (9,23) entry. */
      dest[45] = src[620];  /* (27,23) entry. */
      dest[46] = src[630];  /* (10,24) entry. */
      dest[47] = src[631];  /* (11,24) entry. */
      dest[48] = src[657];  /* (10,25) entry. */
      dest[49] = src[662];  /* (15,25) entry. */
      dest[50] = src[684];  /* (10,26) entry. */
      dest[51] = src[692];  /* (18,26) entry. */
      dest[52] = src[712];  /* (11,27) entry. */
      dest[53] = src[724];  /* (23,27) entry. */
      dest[54] = src[740];  /* (12,28) entry. */
      dest[55] = src[753];  /* (25,28) entry. */
      dest[56] = src[768];  /* (13,29) entry. */
      dest[57] = src[774];  /* (19,29) entry. */
      dest[58] = src[795];  /* (13,30) entry. */
      dest[59] = src[805];  /* (23,30) entry. */
      dest[60] = src[823];  /* (14,31) entry. */
      dest[61] = src[826];  /* (17,31) entry. */
      dest[62] = src[850];  /* (14,32) entry. */
      dest[63] = src[854];  /* (18,32) entry. */
      dest[64] = src[877];  /* (14,33) entry. */
      dest[65] = src[887];  /* (24,33) entry. */
      dest[66] = src[905];  /* (15,34) entry. */
      dest[67] = src[907];  /* (17,34) entry. */
      dest[68] = src[932];  /* (15,35) entry. */
      dest[69] = src[935];  /* (18,35) entry. */
      dest[70] = src[960];  /* (16,36) entry. */
      dest[71] = src[963];  /* (19,36) entry. */
      dest[72] = src[987];  /* (16,37) entry. */
      dest[73] = src[998];  /* (27,37) entry. */
      dest[74] = src[1015];  /* (17,38) entry. */
      dest[75] = src[1016];  /* (18,38) entry. */
      dest[76] = src[1042];  /* (17,39) entry. */
      dest[77] = src[1052];  /* (27,39) entry. */
      dest[78] = src[1071];  /* (19,40) entry. */
      dest[79] = src[1073];  /* (21,40) entry. */
      dest[80] = src[1098];  /* (19,41) entry. */
      dest[81] = src[1101];  /* (22,41) entry. */
      dest[82] = src[1125];  /* (19,42) entry. */
      dest[83] = src[1132];  /* (26,42) entry. */
      dest[84] = src[1154];  /* (21,43) entry. */
      dest[85] = src[1155];  /* (22,43) entry. */
      dest[86] = src[1183];  /* (23,44) entry. */
      dest[87] = src[1184];  /* (24,44) entry. */
      valid_vars++;
    }
  }

  this_var_errors = 0;

  xm = mxGetField(prhs[0], 0, "h_mns");

  if (xm == NULL) {
    printf("could not find params.h_mns.\n");
  } else {
    if (!((mxGetM(xm) == 44) && (mxGetN(xm) == 1))) {
      printf("h_mns must be size (44,1), not (%d,%d).\n", mxGetM(xm), mxGetN(xm));
      this_var_errors++;
    }

    if (mxIsComplex(xm)) {
      printf("parameter h_mns must be real.\n");
      this_var_errors++;
    }

    if (!mxIsClass(xm, "double")) {
      printf("parameter h_mns must be a full matrix of doubles.\n");
      this_var_errors++;
    }

    if (mxIsSparse(xm)) {
      printf("parameter h_mns must be a full matrix.\n");
      this_var_errors++;
    }

    if (this_var_errors == 0) {
      dest = params.h_mns;
      src = mxGetPr(xm);

      for (i = 0; i < 44; i++)
        *dest++ = *src++;

      valid_vars++;
    }
  }

  this_var_errors = 0;

  xm = mxGetField(prhs[0], 0, "h_pls");

  if (xm == NULL) {
    printf("could not find params.h_pls.\n");
  } else {
    if (!((mxGetM(xm) == 44) && (mxGetN(xm) == 1))) {
      printf("h_pls must be size (44,1), not (%d,%d).\n", mxGetM(xm), mxGetN(xm));
      this_var_errors++;
    }

    if (mxIsComplex(xm)) {
      printf("parameter h_pls must be real.\n");
      this_var_errors++;
    }

    if (!mxIsClass(xm, "double")) {
      printf("parameter h_pls must be a full matrix of doubles.\n");
      this_var_errors++;
    }

    if (mxIsSparse(xm)) {
      printf("parameter h_pls must be a full matrix.\n");
      this_var_errors++;
    }

    if (this_var_errors == 0) {
      dest = params.h_pls;
      src = mxGetPr(xm);

      for (i = 0; i < 44; i++)
        *dest++ = *src++;

      valid_vars++;
    }
  }

  if (valid_vars != 5) {
    printf("Error: %d parameters are invalid.\n", 5 - valid_vars);
    mexErrMsgTxt("invalid parameters found.");
  }

  if (prepare_for_c) {
    printf("settings.prepare_for_c == 1. thus, outputting for C.\n");
    for (i = 0; i < 27; i++)
      printf("  params.Delta[%d] = %.6g;\n", i, params.Delta[i]);
    for (i = 0; i < 88; i++)
      printf("  params.K[%d] = %.6g;\n", i, params.K[i]);
    for (i = 0; i < 1; i++)
      printf("  params.B_min[%d] = %.6g;\n", i, params.B_min[i]);
    for (i = 0; i < 44; i++)
      printf("  params.h_mns[%d] = %.6g;\n", i, params.h_mns[i]);
    for (i = 0; i < 44; i++)
      printf("  params.h_pls[%d] = %.6g;\n", i, params.h_pls[i]);
  }

  /* Perform the actual solve in here. */
  steps = solve();

  /* For profiling purposes, allow extra silent solves if desired. */
  settings.verbose = 0;
  for (i = 0; i < extra_solves; i++)
    solve();

  /* Update the status variables. */
  plhs[1] = mxCreateStructArray(1, dims1x1of1, 4, status_names);

  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[1], 0, "optval", xm);
  *mxGetPr(xm) = work.optval;

  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[1], 0, "gap", xm);
  *mxGetPr(xm) = work.gap;

  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[1], 0, "steps", xm);
  *mxGetPr(xm) = steps;

  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[1], 0, "converged", xm);
  *mxGetPr(xm) = work.converged;

  /* Extract variable values. */
  plhs[0] = mxCreateStructArray(1, dims1x1of1, num_var_names, var_names);

  xm = mxCreateDoubleMatrix(44, 1, mxREAL);
  mxSetField(plhs[0], 0, "F", xm);
  dest = mxGetPr(xm);
  src = vars.F;
  for (i = 0; i < 44; i++) {
    *dest++ = *src++;
  }
}
