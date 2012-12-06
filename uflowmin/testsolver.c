/* Produced by CVXGEN, 2012-03-18 11:07:47 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 0

int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif

  set_defaults();
  setup_indexing();
  load_default_data();

  /* Solve problem instance for the record. */
  settings.verbose = 1;
  num_iters = solve();

#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;

  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();

  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;

  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif

  return 0;
}

void load_default_data(void) {
  params.Delta[0] = 0.203191610298302;
  params.Delta[1] = 0.832591290472419;
  params.Delta[2] = -0.836381044348223;
  params.Delta[3] = 0.0433104207906521;
  params.Delta[4] = 1.57178781739062;
  params.Delta[5] = 1.58517235573375;
  params.Delta[6] = -1.49765875814465;
  params.Delta[7] = -1.17102848744725;
  params.Delta[8] = -1.79413118679668;
  params.Delta[9] = -0.236760625397454;
  params.Delta[10] = -1.88049515648573;
  params.Delta[11] = -0.172667102421156;
  params.Delta[12] = 0.596576190459043;
  params.Delta[13] = -0.886050869408099;
  params.Delta[14] = 0.705019607920525;
  params.Delta[15] = 0.363451269665403;
  params.Delta[16] = -1.90407247049134;
  params.Delta[17] = 0.235416351963528;
  params.Delta[18] = -0.962990212370138;
  params.Delta[19] = -0.339595211959721;
  params.Delta[20] = -0.865899672914725;
  params.Delta[21] = 0.772551673251985;
  params.Delta[22] = -0.238185129317042;
  params.Delta[23] = -1.37252904610015;
  params.Delta[24] = 0.178596072127379;
  params.Delta[25] = 1.12125905804547;
  params.Delta[26] = -0.774545870495281;
  params.K[0] = -1.11216846427127;
  params.K[1] = -0.448114969777405;
  params.K[2] = 1.74553459944172;
  params.K[3] = 1.90398168989174;
  params.K[4] = 0.689534703651255;
  params.K[5] = 1.61133643415359;
  params.K[6] = 1.38300348517272;
  params.K[7] = -0.488023834684443;
  params.K[8] = -1.6311319645131;
  params.K[9] = 0.613643610094145;
  params.K[10] = 0.231363049553804;
  params.K[11] = -0.553740947749688;
  params.K[12] = -1.09978198064067;
  params.K[13] = -0.373920334495006;
  params.K[14] = -0.124239005203324;
  params.K[15] = -0.923057686995755;
  params.K[16] = -0.83282890309827;
  params.K[17] = -0.169254402708088;
  params.K[18] = 1.44213565178771;
  params.K[19] = 0.345011617871286;
  params.K[20] = -0.866048550271161;
  params.K[21] = -0.888089973505595;
  params.K[22] = -0.181511697912213;
  params.K[23] = -1.17835862158005;
  params.K[24] = -1.19448515582771;
  params.K[25] = 0.0561402392697676;
  params.K[26] = -1.65108252487678;
  params.K[27] = -0.0656578705936539;
  params.K[28] = -0.551295150448667;
  params.K[29] = 0.830746487262684;
  params.K[30] = 0.986984892408018;
  params.K[31] = 0.764371687423057;
  params.K[32] = 0.756721655019656;
  params.K[33] = -0.505599503404287;
  params.K[34] = 0.67253921894107;
  params.K[35] = -0.640605344172728;
  params.K[36] = 0.2911754794755;
  params.K[37] = -0.696771367740502;
  params.K[38] = -0.219419802945872;
  params.K[39] = -1.75388427668024;
  params.K[40] = -1.02929831126265;
  params.K[41] = 1.88641042469427;
  params.K[42] = -1.0776631825797;
  params.K[43] = 0.765910043789321;
  params.K[44] = 0.601907432854958;
  params.K[45] = 0.895756557749928;
  params.K[46] = -0.0996455574622748;
  params.K[47] = 0.386655098407451;
  params.K[48] = -1.73212230426869;
  params.K[49] = -1.70975144871107;
  params.K[50] = -1.20409589481169;
  params.K[51] = -1.39255601196584;
  params.K[52] = -1.59958262167422;
  params.K[53] = -1.48282454156458;
  params.K[54] = 0.213110927230614;
  params.K[55] = -1.24874070030449;
  params.K[56] = 1.80840497212483;
  params.K[57] = 0.726447115229707;
  params.K[58] = 0.164078693439085;
  params.K[59] = 0.828722403231591;
  params.K[60] = -0.944453316189946;
  params.K[61] = 1.70690273701491;
  params.K[62] = 1.35677223119988;
  params.K[63] = 0.905277993712149;
  params.K[64] = -0.0790401756583599;
  params.K[65] = 1.36841274350659;
  params.K[66] = 0.979009293697437;
  params.K[67] = 0.64130362559845;
  params.K[68] = 1.65590106802375;
  params.K[69] = 0.534662255150299;
  params.K[70] = -0.536237660589562;
  params.K[71] = 0.211378292601782;
  params.K[72] = -1.21447769319945;
  params.K[73] = -1.23171081442559;
  params.K[74] = 0.902678495731283;
  params.K[75] = 1.13974681372452;
  params.K[76] = 1.88839345473506;
  params.K[77] = 1.40388566816601;
  params.K[78] = 0.174377306383291;
  params.K[79] = -1.64083652190774;
  params.K[80] = -0.0445070215355488;
  params.K[81] = 1.7117453902485;
  params.K[82] = 1.15047279801391;
  params.K[83] = -0.0596230957836474;
  params.K[84] = -0.178882554076455;
  params.K[85] = -1.12805692636259;
  params.K[86] = -1.29114647679271;
  params.K[87] = -1.70550532312257;
  params.B_min[0] = 1.56957275034837;
}
