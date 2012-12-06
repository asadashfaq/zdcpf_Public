/* Produced by CVXGEN, 2012-03-18 10:53:46 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */

#include "solver.h"

void multbymA(double *lhs, double *rhs) {
}

void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
  lhs[10] = 0;
  lhs[11] = 0;
  lhs[12] = 0;
  lhs[13] = 0;
  lhs[14] = 0;
  lhs[15] = 0;
  lhs[16] = 0;
  lhs[17] = 0;
  lhs[18] = 0;
  lhs[19] = 0;
  lhs[20] = 0;
  lhs[21] = 0;
  lhs[22] = 0;
  lhs[23] = 0;
  lhs[24] = 0;
  lhs[25] = 0;
  lhs[26] = 0;
  lhs[27] = 0;
  lhs[28] = 0;
  lhs[29] = 0;
  lhs[30] = 0;
  lhs[31] = 0;
  lhs[32] = 0;
  lhs[33] = 0;
  lhs[34] = 0;
  lhs[35] = 0;
  lhs[36] = 0;
  lhs[37] = 0;
  lhs[38] = 0;
  lhs[39] = 0;
  lhs[40] = 0;
  lhs[41] = 0;
  lhs[42] = 0;
  lhs[43] = 0;
  lhs[44] = 0;
  lhs[45] = 0;
  lhs[46] = 0;
  lhs[47] = 0;
  lhs[48] = 0;
  lhs[49] = 0;
  lhs[50] = 0;
  lhs[51] = 0;
  lhs[52] = 0;
  lhs[53] = 0;
  lhs[54] = 0;
  lhs[55] = 0;
  lhs[56] = 0;
  lhs[57] = 0;
  lhs[58] = 0;
  lhs[59] = 0;
  lhs[60] = 0;
  lhs[61] = 0;
  lhs[62] = 0;
  lhs[63] = 0;
  lhs[64] = 0;
  lhs[65] = 0;
  lhs[66] = 0;
  lhs[67] = 0;
  lhs[68] = 0;
  lhs[69] = 0;
  lhs[70] = 0;
}

void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(params.K[0])-rhs[1]*(params.K[2])-rhs[2]*(params.K[4])-rhs[3]*(params.K[6])-rhs[4]*(params.K[8])-rhs[5]*(params.K[10])-rhs[44]*(-1);
  lhs[1] = -rhs[6]*(params.K[12])-rhs[45]*(-1);
  lhs[2] = -rhs[7]*(params.K[14])-rhs[8]*(params.K[16])-rhs[46]*(-1);
  lhs[3] = -rhs[9]*(params.K[18])-rhs[10]*(params.K[20])-rhs[47]*(-1);
  lhs[4] = -rhs[11]*(params.K[22])-rhs[12]*(params.K[24])-rhs[13]*(params.K[26])-rhs[14]*(params.K[28])-rhs[15]*(params.K[30])-rhs[16]*(params.K[32])-rhs[48]*(-1);
  lhs[5] = -rhs[17]*(params.K[34])-rhs[18]*(params.K[36])-rhs[49]*(-1);
  lhs[6] = -rhs[7]*(params.K[15])-rhs[11]*(params.K[23])-rhs[50]*(-1);
  lhs[7] = -rhs[12]*(params.K[25])-rhs[19]*(params.K[38])-rhs[51]*(-1);
  lhs[8] = -rhs[20]*(params.K[40])-rhs[21]*(params.K[42])-rhs[22]*(params.K[44])-rhs[52]*(-1);
  lhs[9] = -rhs[23]*(params.K[46])-rhs[24]*(params.K[48])-rhs[25]*(params.K[50])-rhs[53]*(-1);
  lhs[10] = -rhs[23]*(params.K[47])-rhs[26]*(params.K[52])-rhs[54]*(-1);
  lhs[11] = -rhs[27]*(params.K[54])-rhs[55]*(-1);
  lhs[12] = -rhs[0]*(params.K[1])-rhs[13]*(params.K[27])-rhs[28]*(params.K[56])-rhs[29]*(params.K[58])-rhs[56]*(-1);
  lhs[13] = -rhs[9]*(params.K[19])-rhs[30]*(params.K[60])-rhs[31]*(params.K[62])-rhs[32]*(params.K[64])-rhs[57]*(-1);
  lhs[14] = -rhs[24]*(params.K[49])-rhs[33]*(params.K[66])-rhs[34]*(params.K[68])-rhs[58]*(-1);
  lhs[15] = -rhs[1]*(params.K[3])-rhs[20]*(params.K[41])-rhs[35]*(params.K[70])-rhs[36]*(params.K[72])-rhs[59]*(-1);
  lhs[16] = -rhs[2]*(params.K[5])-rhs[30]*(params.K[61])-rhs[33]*(params.K[67])-rhs[37]*(params.K[74])-rhs[38]*(params.K[76])-rhs[60]*(-1);
  lhs[17] = -rhs[10]*(params.K[21])-rhs[25]*(params.K[51])-rhs[31]*(params.K[63])-rhs[34]*(params.K[69])-rhs[37]*(params.K[75])-rhs[61]*(-1);
  lhs[18] = -rhs[3]*(params.K[7])-rhs[8]*(params.K[17])-rhs[14]*(params.K[29])-rhs[21]*(params.K[43])-rhs[28]*(params.K[57])-rhs[35]*(params.K[71])-rhs[39]*(params.K[78])-rhs[40]*(params.K[80])-rhs[41]*(params.K[82])-rhs[62]*(-1);
  lhs[19] = -rhs[19]*(params.K[39])-rhs[63]*(-1);
  lhs[20] = -rhs[6]*(params.K[13])-rhs[17]*(params.K[35])-rhs[39]*(params.K[79])-rhs[42]*(params.K[84])-rhs[64]*(-1);
  lhs[21] = -rhs[18]*(params.K[37])-rhs[40]*(params.K[81])-rhs[42]*(params.K[85])-rhs[65]*(-1);
  lhs[22] = -rhs[4]*(params.K[9])-rhs[15]*(params.K[31])-rhs[26]*(params.K[53])-rhs[29]*(params.K[59])-rhs[43]*(params.K[86])-rhs[66]*(-1);
  lhs[23] = -rhs[5]*(params.K[11])-rhs[32]*(params.K[65])-rhs[43]*(params.K[87])-rhs[67]*(-1);
  lhs[24] = -rhs[16]*(params.K[33])-rhs[27]*(params.K[55])-rhs[68]*(-1);
  lhs[25] = -rhs[41]*(params.K[83])-rhs[69]*(-1);
  lhs[26] = -rhs[22]*(params.K[45])-rhs[36]*(params.K[73])-rhs[38]*(params.K[77])-rhs[70]*(-1);
  lhs[27] = -rhs[44]*(-1);
  lhs[28] = -rhs[45]*(-1);
  lhs[29] = -rhs[46]*(-1);
  lhs[30] = -rhs[47]*(-1);
  lhs[31] = -rhs[48]*(-1);
  lhs[32] = -rhs[49]*(-1);
  lhs[33] = -rhs[50]*(-1);
  lhs[34] = -rhs[51]*(-1);
  lhs[35] = -rhs[52]*(-1);
  lhs[36] = -rhs[53]*(-1);
  lhs[37] = -rhs[54]*(-1);
  lhs[38] = -rhs[55]*(-1);
  lhs[39] = -rhs[56]*(-1);
  lhs[40] = -rhs[57]*(-1);
  lhs[41] = -rhs[58]*(-1);
  lhs[42] = -rhs[59]*(-1);
  lhs[43] = -rhs[60]*(-1);
  lhs[44] = -rhs[61]*(-1);
  lhs[45] = -rhs[62]*(-1);
  lhs[46] = -rhs[63]*(-1);
  lhs[47] = -rhs[64]*(-1);
  lhs[48] = -rhs[65]*(-1);
  lhs[49] = -rhs[66]*(-1);
  lhs[50] = -rhs[67]*(-1);
  lhs[51] = -rhs[68]*(-1);
  lhs[52] = -rhs[69]*(-1);
  lhs[53] = -rhs[70]*(-1);
}

void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(params.K[0])-rhs[12]*(params.K[1]);
  lhs[1] = -rhs[0]*(params.K[2])-rhs[15]*(params.K[3]);
  lhs[2] = -rhs[0]*(params.K[4])-rhs[16]*(params.K[5]);
  lhs[3] = -rhs[0]*(params.K[6])-rhs[18]*(params.K[7]);
  lhs[4] = -rhs[0]*(params.K[8])-rhs[22]*(params.K[9]);
  lhs[5] = -rhs[0]*(params.K[10])-rhs[23]*(params.K[11]);
  lhs[6] = -rhs[1]*(params.K[12])-rhs[20]*(params.K[13]);
  lhs[7] = -rhs[2]*(params.K[14])-rhs[6]*(params.K[15]);
  lhs[8] = -rhs[2]*(params.K[16])-rhs[18]*(params.K[17]);
  lhs[9] = -rhs[3]*(params.K[18])-rhs[13]*(params.K[19]);
  lhs[10] = -rhs[3]*(params.K[20])-rhs[17]*(params.K[21]);
  lhs[11] = -rhs[4]*(params.K[22])-rhs[6]*(params.K[23]);
  lhs[12] = -rhs[4]*(params.K[24])-rhs[7]*(params.K[25]);
  lhs[13] = -rhs[4]*(params.K[26])-rhs[12]*(params.K[27]);
  lhs[14] = -rhs[4]*(params.K[28])-rhs[18]*(params.K[29]);
  lhs[15] = -rhs[4]*(params.K[30])-rhs[22]*(params.K[31]);
  lhs[16] = -rhs[4]*(params.K[32])-rhs[24]*(params.K[33]);
  lhs[17] = -rhs[5]*(params.K[34])-rhs[20]*(params.K[35]);
  lhs[18] = -rhs[5]*(params.K[36])-rhs[21]*(params.K[37]);
  lhs[19] = -rhs[7]*(params.K[38])-rhs[19]*(params.K[39]);
  lhs[20] = -rhs[8]*(params.K[40])-rhs[15]*(params.K[41]);
  lhs[21] = -rhs[8]*(params.K[42])-rhs[18]*(params.K[43]);
  lhs[22] = -rhs[8]*(params.K[44])-rhs[26]*(params.K[45]);
  lhs[23] = -rhs[9]*(params.K[46])-rhs[10]*(params.K[47]);
  lhs[24] = -rhs[9]*(params.K[48])-rhs[14]*(params.K[49]);
  lhs[25] = -rhs[9]*(params.K[50])-rhs[17]*(params.K[51]);
  lhs[26] = -rhs[10]*(params.K[52])-rhs[22]*(params.K[53]);
  lhs[27] = -rhs[11]*(params.K[54])-rhs[24]*(params.K[55]);
  lhs[28] = -rhs[12]*(params.K[56])-rhs[18]*(params.K[57]);
  lhs[29] = -rhs[12]*(params.K[58])-rhs[22]*(params.K[59]);
  lhs[30] = -rhs[13]*(params.K[60])-rhs[16]*(params.K[61]);
  lhs[31] = -rhs[13]*(params.K[62])-rhs[17]*(params.K[63]);
  lhs[32] = -rhs[13]*(params.K[64])-rhs[23]*(params.K[65]);
  lhs[33] = -rhs[14]*(params.K[66])-rhs[16]*(params.K[67]);
  lhs[34] = -rhs[14]*(params.K[68])-rhs[17]*(params.K[69]);
  lhs[35] = -rhs[15]*(params.K[70])-rhs[18]*(params.K[71]);
  lhs[36] = -rhs[15]*(params.K[72])-rhs[26]*(params.K[73]);
  lhs[37] = -rhs[16]*(params.K[74])-rhs[17]*(params.K[75]);
  lhs[38] = -rhs[16]*(params.K[76])-rhs[26]*(params.K[77]);
  lhs[39] = -rhs[18]*(params.K[78])-rhs[20]*(params.K[79]);
  lhs[40] = -rhs[18]*(params.K[80])-rhs[21]*(params.K[81]);
  lhs[41] = -rhs[18]*(params.K[82])-rhs[25]*(params.K[83]);
  lhs[42] = -rhs[20]*(params.K[84])-rhs[21]*(params.K[85]);
  lhs[43] = -rhs[22]*(params.K[86])-rhs[23]*(params.K[87]);
  lhs[44] = -rhs[0]*(-1)-rhs[27]*(-1);
  lhs[45] = -rhs[1]*(-1)-rhs[28]*(-1);
  lhs[46] = -rhs[2]*(-1)-rhs[29]*(-1);
  lhs[47] = -rhs[3]*(-1)-rhs[30]*(-1);
  lhs[48] = -rhs[4]*(-1)-rhs[31]*(-1);
  lhs[49] = -rhs[5]*(-1)-rhs[32]*(-1);
  lhs[50] = -rhs[6]*(-1)-rhs[33]*(-1);
  lhs[51] = -rhs[7]*(-1)-rhs[34]*(-1);
  lhs[52] = -rhs[8]*(-1)-rhs[35]*(-1);
  lhs[53] = -rhs[9]*(-1)-rhs[36]*(-1);
  lhs[54] = -rhs[10]*(-1)-rhs[37]*(-1);
  lhs[55] = -rhs[11]*(-1)-rhs[38]*(-1);
  lhs[56] = -rhs[12]*(-1)-rhs[39]*(-1);
  lhs[57] = -rhs[13]*(-1)-rhs[40]*(-1);
  lhs[58] = -rhs[14]*(-1)-rhs[41]*(-1);
  lhs[59] = -rhs[15]*(-1)-rhs[42]*(-1);
  lhs[60] = -rhs[16]*(-1)-rhs[43]*(-1);
  lhs[61] = -rhs[17]*(-1)-rhs[44]*(-1);
  lhs[62] = -rhs[18]*(-1)-rhs[45]*(-1);
  lhs[63] = -rhs[19]*(-1)-rhs[46]*(-1);
  lhs[64] = -rhs[20]*(-1)-rhs[47]*(-1);
  lhs[65] = -rhs[21]*(-1)-rhs[48]*(-1);
  lhs[66] = -rhs[22]*(-1)-rhs[49]*(-1);
  lhs[67] = -rhs[23]*(-1)-rhs[50]*(-1);
  lhs[68] = -rhs[24]*(-1)-rhs[51]*(-1);
  lhs[69] = -rhs[25]*(-1)-rhs[52]*(-1);
  lhs[70] = -rhs[26]*(-1)-rhs[53]*(-1);
}

void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
  lhs[10] = 0;
  lhs[11] = 0;
  lhs[12] = 0;
  lhs[13] = 0;
  lhs[14] = 0;
  lhs[15] = 0;
  lhs[16] = 0;
  lhs[17] = 0;
  lhs[18] = 0;
  lhs[19] = 0;
  lhs[20] = 0;
  lhs[21] = 0;
  lhs[22] = 0;
  lhs[23] = 0;
  lhs[24] = 0;
  lhs[25] = 0;
  lhs[26] = 0;
  lhs[27] = 0;
  lhs[28] = 0;
  lhs[29] = 0;
  lhs[30] = 0;
  lhs[31] = 0;
  lhs[32] = 0;
  lhs[33] = 0;
  lhs[34] = 0;
  lhs[35] = 0;
  lhs[36] = 0;
  lhs[37] = 0;
  lhs[38] = 0;
  lhs[39] = 0;
  lhs[40] = 0;
  lhs[41] = 0;
  lhs[42] = 0;
  lhs[43] = 0;
  lhs[44] = 0;
  lhs[45] = 0;
  lhs[46] = 0;
  lhs[47] = 0;
  lhs[48] = 0;
  lhs[49] = 0;
  lhs[50] = 0;
  lhs[51] = 0;
  lhs[52] = 0;
  lhs[53] = 0;
  lhs[54] = 0;
  lhs[55] = 0;
  lhs[56] = 0;
  lhs[57] = 0;
  lhs[58] = 0;
  lhs[59] = 0;
  lhs[60] = 0;
  lhs[61] = 0;
  lhs[62] = 0;
  lhs[63] = 0;
  lhs[64] = 0;
  lhs[65] = 0;
  lhs[66] = 0;
  lhs[67] = 0;
  lhs[68] = 0;
  lhs[69] = 0;
  lhs[70] = 0;
}

void fillq(void) {
  work.q[0] = 0;
  work.q[1] = 0;
  work.q[2] = 0;
  work.q[3] = 0;
  work.q[4] = 0;
  work.q[5] = 0;
  work.q[6] = 0;
  work.q[7] = 0;
  work.q[8] = 0;
  work.q[9] = 0;
  work.q[10] = 0;
  work.q[11] = 0;
  work.q[12] = 0;
  work.q[13] = 0;
  work.q[14] = 0;
  work.q[15] = 0;
  work.q[16] = 0;
  work.q[17] = 0;
  work.q[18] = 0;
  work.q[19] = 0;
  work.q[20] = 0;
  work.q[21] = 0;
  work.q[22] = 0;
  work.q[23] = 0;
  work.q[24] = 0;
  work.q[25] = 0;
  work.q[26] = 0;
  work.q[27] = 0;
  work.q[28] = 0;
  work.q[29] = 0;
  work.q[30] = 0;
  work.q[31] = 0;
  work.q[32] = 0;
  work.q[33] = 0;
  work.q[34] = 0;
  work.q[35] = 0;
  work.q[36] = 0;
  work.q[37] = 0;
  work.q[38] = 0;
  work.q[39] = 0;
  work.q[40] = 0;
  work.q[41] = 0;
  work.q[42] = 0;
  work.q[43] = 0;
  work.q[44] = 1;
  work.q[45] = 1;
  work.q[46] = 1;
  work.q[47] = 1;
  work.q[48] = 1;
  work.q[49] = 1;
  work.q[50] = 1;
  work.q[51] = 1;
  work.q[52] = 1;
  work.q[53] = 1;
  work.q[54] = 1;
  work.q[55] = 1;
  work.q[56] = 1;
  work.q[57] = 1;
  work.q[58] = 1;
  work.q[59] = 1;
  work.q[60] = 1;
  work.q[61] = 1;
  work.q[62] = 1;
  work.q[63] = 1;
  work.q[64] = 1;
  work.q[65] = 1;
  work.q[66] = 1;
  work.q[67] = 1;
  work.q[68] = 1;
  work.q[69] = 1;
  work.q[70] = 1;
}

void fillh(void) {
  work.h[0] = params.Delta[0];
  work.h[1] = params.Delta[1];
  work.h[2] = params.Delta[2];
  work.h[3] = params.Delta[3];
  work.h[4] = params.Delta[4];
  work.h[5] = params.Delta[5];
  work.h[6] = params.Delta[6];
  work.h[7] = params.Delta[7];
  work.h[8] = params.Delta[8];
  work.h[9] = params.Delta[9];
  work.h[10] = params.Delta[10];
  work.h[11] = params.Delta[11];
  work.h[12] = params.Delta[12];
  work.h[13] = params.Delta[13];
  work.h[14] = params.Delta[14];
  work.h[15] = params.Delta[15];
  work.h[16] = params.Delta[16];
  work.h[17] = params.Delta[17];
  work.h[18] = params.Delta[18];
  work.h[19] = params.Delta[19];
  work.h[20] = params.Delta[20];
  work.h[21] = params.Delta[21];
  work.h[22] = params.Delta[22];
  work.h[23] = params.Delta[23];
  work.h[24] = params.Delta[24];
  work.h[25] = params.Delta[25];
  work.h[26] = params.Delta[26];
  work.h[27] = 0;
  work.h[28] = 0;
  work.h[29] = 0;
  work.h[30] = 0;
  work.h[31] = 0;
  work.h[32] = 0;
  work.h[33] = 0;
  work.h[34] = 0;
  work.h[35] = 0;
  work.h[36] = 0;
  work.h[37] = 0;
  work.h[38] = 0;
  work.h[39] = 0;
  work.h[40] = 0;
  work.h[41] = 0;
  work.h[42] = 0;
  work.h[43] = 0;
  work.h[44] = 0;
  work.h[45] = 0;
  work.h[46] = 0;
  work.h[47] = 0;
  work.h[48] = 0;
  work.h[49] = 0;
  work.h[50] = 0;
  work.h[51] = 0;
  work.h[52] = 0;
  work.h[53] = 0;
}

void fillb(void) {
}

void pre_ops(void) {
}
