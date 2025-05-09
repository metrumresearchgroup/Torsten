functions {
  vector f_0_arg(real t, vector z) {
    return z;
  }
  vector f_1_arg(real t, vector z, real a) {
    return z;
  }
  vector f_2_arg(real t, vector z, int b, real a) {
    return z;
  }
  vector f_3_arg(real t, vector z, array[] real c, int b, real a) {
    return z;
  }
  vector f_4_arg(real t, vector z, array[] int d, array[] real c, int b,
                 real a) {
    return z;
  }
  vector f_5_arg(real t, vector z, vector e, array[] int d, array[] real c,
                 int b, real a) {
    return z;
  }
  vector f_6_arg(real t, vector z, row_vector f, vector e, array[] int d,
                 array[] real c, int b, real a) {
    return z;
  }
  vector f_7_arg(real t, vector z, matrix g, row_vector f, vector e,
                 array[] int d, array[] real c, int b, real a) {
    return z;
  }
  vector f_8_arg(real t, vector z, array[,] real h, matrix g, row_vector f,
                 vector e, array[] int d, array[] real c, int b, real a) {
    return z;
  }
  vector f_9_arg(real t, vector z, array[,] int i, array[,] real h, matrix g,
                 row_vector f, vector e, array[] int d, array[] real c,
                 int b, real a) {
    return z;
  }
  vector f_10_arg(real t, vector z, array[] vector j, array[,] int i,
                  array[,] real h, matrix g, row_vector f, vector e,
                  array[] int d, array[] real c, int b, real a) {
    return z;
  }
  vector f_11_arg(real t, vector z, array[] row_vector k, array[] vector j,
                  array[,] int i, array[,] real h, matrix g, row_vector f,
                  vector e, array[] int d, array[] real c, int b, real a) {
    return z;
  }
  vector f_12_arg(real t, vector z, array[] matrix l, array[] row_vector k,
                  array[] vector j, array[,] int i, array[,] real h,
                  matrix g, row_vector f, vector e, array[] int d,
                  array[] real c, int b, real a) {
    return z;
  }
}
data {
  int N;
  int id;
  real rd;
  array[N] real rad;
  array[N] int iad;
  vector[N] vd;
  row_vector[N] rvd;
  matrix[N, N] md;
  array[N, N] real raad;
  array[N, N] int iaad;
  array[N] vector[N] vad;
  array[N] row_vector[N] rvad;
  array[N] matrix[N, N] mad;
}
transformed data {
  // ODE
  array[N] vector[N] zd = pmx_ode_bdf(f_12_arg, vd, rd, rad, mad, rvad, vad,
                                    iaad, raad, md, rvd, vd, iad, rad, id,
                                    rd);
  zd = pmx_ode_bdf(f_11_arg, vd, rd, rad, rvad, vad, iaad, raad, md, rvd, vd,
                 iad, rad, id, rd);
  zd = pmx_ode_bdf(f_10_arg, vd, rd, rad, vad, iaad, raad, md, rvd, vd, iad,
                 rad, id, rd);
  zd = pmx_ode_bdf(f_9_arg, vd, rd, rad, iaad, raad, md, rvd, vd, iad, rad, id,
                 rd);
  zd = pmx_ode_bdf(f_8_arg, vd, rd, rad, raad, md, rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf(f_7_arg, vd, rd, rad, md, rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf(f_6_arg, vd, rd, rad, rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf(f_5_arg, vd, rd, rad, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf(f_4_arg, vd, rd, rad, iad, rad, id, rd);
  zd = pmx_ode_bdf(f_3_arg, vd, rd, rad, rad, id, rd);
  zd = pmx_ode_bdf(f_2_arg, vd, rd, rad, id, rd);
  zd = pmx_ode_bdf(f_1_arg, vd, rd, rad, rd);
  zd = pmx_ode_bdf(f_0_arg, vd, rd, rad);
  // ODE with control
  zd = pmx_ode_bdf_ctrl(f_12_arg, vd, rd, rad, 1e-6, 1e-6,
                                        100, mad, rvad, vad, iaad, raad, md,
                                        rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_11_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvad, vad, iaad,
                     raad, md, rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_10_arg, vd, rd, rad, 1e-6, 1e-6, 100, vad, iaad, raad,
                     md, rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_9_arg, vd, rd, rad, 1e-6, 1e-6, 100, iaad, raad, md,
                     rvd, vd, iad, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_8_arg, vd, rd, rad, 1e-6, 1e-6, 100, raad, md, rvd,
                     vd, iad, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_7_arg, vd, rd, rad, 1e-6, 1e-6, 100, md, rvd, vd, iad,
                     rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_6_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvd, vd, iad,
                     rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_5_arg, vd, rd, rad, 1e-6, 1e-6, 100, vd, iad, rad, id,
                     rd);
  zd = pmx_ode_bdf_ctrl(f_4_arg, vd, rd, rad, 1e-6, 1e-6, 100, iad, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_3_arg, vd, rd, rad, 1e-6, 1e-6, 100, rad, id, rd);
  zd = pmx_ode_bdf_ctrl(f_2_arg, vd, rd, rad, 1e-6, 1e-6, 100, id, rd);
  zd = pmx_ode_bdf_ctrl(f_1_arg, vd, rd, rad, 1e-6, 1e-6, 100, rd);
  zd = pmx_ode_bdf_ctrl(f_0_arg, vd, rd, rad, 1e-6, 1e-6, 100);
}
parameters {
  real r;
  array[N] real ra;
  vector[N] v;
  row_vector[N] rv;
  matrix[N, N] m;
  array[N, N] real raa;
  array[N] vector[N] va;
  array[N] row_vector[N] rva;
  array[N] matrix[N, N] ma;
}
transformed parameters {
  // ODE
  array[N] vector[N] z = pmx_ode_bdf(f_12_arg, vd, rd, rad, mad, rvad, vad,
                                   iaad, raad, md, rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf(f_11_arg, vd, rd, rad, rvad, vad, iaad, raad, md, rvd, vd,
                iad, rad, id, rd);
  z = pmx_ode_bdf(f_10_arg, vd, rd, rad, vad, iaad, raad, md, rvd, vd, iad,
                rad, id, rd);
  z = pmx_ode_bdf(f_9_arg, vd, rd, rad, iaad, raad, md, rvd, vd, iad, rad, id,
                rd);
  z = pmx_ode_bdf(f_8_arg, vd, rd, rad, raad, md, rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf(f_7_arg, vd, rd, rad, md, rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf(f_6_arg, vd, rd, rad, rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf(f_5_arg, vd, rd, rad, vd, iad, rad, id, rd);
  z = pmx_ode_bdf(f_4_arg, vd, rd, rad, iad, rad, id, rd);
  z = pmx_ode_bdf(f_3_arg, vd, rd, rad, rad, id, rd);
  z = pmx_ode_bdf(f_2_arg, vd, rd, rad, id, rd);
  z = pmx_ode_bdf(f_1_arg, vd, rd, rad, rd);
  z = pmx_ode_bdf(f_0_arg, vd, rd, rad);
  z = pmx_ode_bdf(f_12_arg, v, rd, ra, ma, rva, va, iaad, raa, m, rv, v, iad,
                ra, id, r);
  z = pmx_ode_bdf(f_11_arg, v, rd, ra, rva, va, iaad, raa, m, rv, v, iad, ra,
                id, r);
  z = pmx_ode_bdf(f_10_arg, v, rd, ra, va, iaad, raa, m, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf(f_9_arg, v, rd, ra, iaad, raa, m, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf(f_8_arg, v, rd, ra, raa, m, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf(f_7_arg, v, rd, ra, m, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf(f_6_arg, v, rd, ra, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf(f_5_arg, v, rd, ra, v, iad, ra, id, r);
  z = pmx_ode_bdf(f_4_arg, v, rd, ra, iad, ra, id, r);
  z = pmx_ode_bdf(f_3_arg, v, rd, ra, ra, id, r);
  z = pmx_ode_bdf(f_2_arg, v, rd, ra, id, r);
  z = pmx_ode_bdf(f_1_arg, v, rd, ra, r);
  z = pmx_ode_bdf(f_0_arg, v, rd, ra);
//   z = pmx_ode_bdf(f_12_arg, v, r, ra, ma, rva, va, iaad, raa, m, rv, v, iad,
//                 ra, id, r);
//   z = pmx_ode_bdf(f_11_arg, v, r, ra, rva, va, iaad, raa, m, rv, v, iad, ra,
//                 id, r);
//   z = pmx_ode_bdf(f_10_arg, v, r, ra, va, iaad, raa, m, rv, v, iad, ra, id, r);
//   z = pmx_ode_bdf(f_9_arg, v, r, ra, iaad, raa, m, rv, v, iad, ra, id, r);
//   z = pmx_ode_bdf(f_8_arg, v, r, ra, raa, m, rv, v, iad, ra, id, r);
//   z = pmx_ode_bdf(f_7_arg, v, r, ra, m, rv, v, iad, ra, id, r);
//   z = pmx_ode_bdf(f_6_arg, v, r, ra, rv, v, iad, ra, id, r);
//   z = pmx_ode_bdf(f_5_arg, v, r, ra, v, iad, ra, id, r);
//   z = pmx_ode_bdf(f_4_arg, v, r, ra, iad, ra, id, r);
//   z = pmx_ode_bdf(f_3_arg, v, r, ra, ra, id, r);
//   z = pmx_ode_bdf(f_2_arg, v, r, ra, id, r);
//   z = pmx_ode_bdf(f_1_arg, v, r, ra, r);
//   z = pmx_ode_bdf(f_0_arg, v, r, ra);
  // ODE with control
  z = pmx_ode_bdf_ctrl(f_12_arg, vd, rd, rad, 1e-6, 1e-6,
                                       100, mad, rvad, vad, iaad, raad, md,
                                       rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_11_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvad, vad, iaad,
                    raad, md, rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_10_arg, vd, rd, rad, 1e-6, 1e-6, 100, vad, iaad, raad,
                    md, rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_9_arg, vd, rd, rad, 1e-6, 1e-6, 100, iaad, raad, md,
                    rvd, vd, iad, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_8_arg, vd, rd, rad, 1e-6, 1e-6, 100, raad, md, rvd, vd,
                    iad, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_7_arg, vd, rd, rad, 1e-6, 1e-6, 100, md, rvd, vd, iad,
                    rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_6_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvd, vd, iad, rad,
                    id, rd);
  z = pmx_ode_bdf_ctrl(f_5_arg, vd, rd, rad, 1e-6, 1e-6, 100, vd, iad, rad, id,
                    rd);
  z = pmx_ode_bdf_ctrl(f_4_arg, vd, rd, rad, 1e-6, 1e-6, 100, iad, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_3_arg, vd, rd, rad, 1e-6, 1e-6, 100, rad, id, rd);
  z = pmx_ode_bdf_ctrl(f_2_arg, vd, rd, rad, 1e-6, 1e-6, 100, id, rd);
  z = pmx_ode_bdf_ctrl(f_1_arg, vd, rd, rad, 1e-6, 1e-6, 100, rd);
  z = pmx_ode_bdf_ctrl(f_0_arg, vd, rd, rad, 1e-6, 1e-6, 100);
  z = pmx_ode_bdf_ctrl(f_12_arg, v, rd, ra, 1e-6, 1e-6, 100, ma, rva, va, iaad,
                    raa, m, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_11_arg, v, rd, ra, 1e-6, 1e-6, 100, rva, va, iaad, raa,
                    m, rv, v, iad, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_10_arg, v, rd, ra, 1e-6, 1e-6, 100, va, iaad, raa, m,
                    rv, v, iad, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_9_arg, v, rd, ra, 1e-6, 1e-6, 100, iaad, raa, m, rv, v,
                    iad, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_8_arg, v, rd, ra, 1e-6, 1e-6, 100, raa, m, rv, v, iad,
                    ra, id, r);
  z = pmx_ode_bdf_ctrl(f_7_arg, v, rd, ra, 1e-6, 1e-6, 100, m, rv, v, iad, ra,
                    id, r);
  z = pmx_ode_bdf_ctrl(f_6_arg, v, rd, ra, 1e-6, 1e-6, 100, rv, v, iad, ra, id,
                    r);
  z = pmx_ode_bdf_ctrl(f_5_arg, v, rd, ra, 1e-6, 1e-6, 100, v, iad, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_4_arg, v, rd, ra, 1e-6, 1e-6, 100, iad, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_3_arg, v, rd, ra, 1e-6, 1e-6, 100, ra, id, r);
  z = pmx_ode_bdf_ctrl(f_2_arg, v, rd, ra, 1e-6, 1e-6, 100, id, r);
  z = pmx_ode_bdf_ctrl(f_1_arg, v, rd, ra, 1e-6, 1e-6, 100, r);
  z = pmx_ode_bdf_ctrl(f_0_arg, v, rd, ra, 1e-6, 1e-6, 100);
}
model {
  // ODE
  array[N] vector[N] zm = pmx_ode_bdf(f_12_arg, vd, rd, rad, mad, rvad, vad,
                                    iaad, raad, md, rvd, vd, iad, rad, id,
                                    rd);
  zm = pmx_ode_bdf(f_11_arg, vd, rd, rad, rvad, vad, iaad, raad, md, rvd, vd,
                 iad, rad, id, rd);
  zm = pmx_ode_bdf(f_10_arg, vd, rd, rad, vad, iaad, raad, md, rvd, vd, iad,
                 rad, id, rd);
  zm = pmx_ode_bdf(f_9_arg, vd, rd, rad, iaad, raad, md, rvd, vd, iad, rad, id,
                 rd);
  zm = pmx_ode_bdf(f_8_arg, vd, rd, rad, raad, md, rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf(f_7_arg, vd, rd, rad, md, rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf(f_6_arg, vd, rd, rad, rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf(f_5_arg, vd, rd, rad, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf(f_4_arg, vd, rd, rad, iad, rad, id, rd);
  zm = pmx_ode_bdf(f_3_arg, vd, rd, rad, rad, id, rd);
  zm = pmx_ode_bdf(f_2_arg, vd, rd, rad, id, rd);
  zm = pmx_ode_bdf(f_1_arg, vd, rd, rad, rd);
  zm = pmx_ode_bdf(f_0_arg, vd, rd, rad);
  zm = pmx_ode_bdf(f_12_arg, v, rd, ra, ma, rva, va, iaad, raa, m, rv, v, iad,
                 ra, id, r);
  zm = pmx_ode_bdf(f_11_arg, v, rd, ra, rva, va, iaad, raa, m, rv, v, iad, ra,
                 id, r);
  zm = pmx_ode_bdf(f_10_arg, v, rd, ra, va, iaad, raa, m, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf(f_9_arg, v, rd, ra, iaad, raa, m, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf(f_8_arg, v, rd, ra, raa, m, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf(f_7_arg, v, rd, ra, m, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf(f_6_arg, v, rd, ra, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf(f_5_arg, v, rd, ra, v, iad, ra, id, r);
  zm = pmx_ode_bdf(f_4_arg, v, rd, ra, iad, ra, id, r);
  zm = pmx_ode_bdf(f_3_arg, v, rd, ra, ra, id, r);
  zm = pmx_ode_bdf(f_2_arg, v, rd, ra, id, r);
  zm = pmx_ode_bdf(f_1_arg, v, rd, ra, r);
  zm = pmx_ode_bdf(f_0_arg, v, rd, ra);
//   zm = pmx_ode_bdf(f_12_arg, v, r, ra, ma, rva, va, iaad, raa, m, rv, v, iad,
//                  ra, id, r);
//   zm = pmx_ode_bdf(f_11_arg, v, r, ra, rva, va, iaad, raa, m, rv, v, iad, ra,
//                  id, r);
//   zm = pmx_ode_bdf(f_10_arg, v, r, ra, va, iaad, raa, m, rv, v, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_9_arg, v, r, ra, iaad, raa, m, rv, v, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_8_arg, v, r, ra, raa, m, rv, v, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_7_arg, v, r, ra, m, rv, v, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_6_arg, v, r, ra, rv, v, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_5_arg, v, r, ra, v, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_4_arg, v, r, ra, iad, ra, id, r);
//   zm = pmx_ode_bdf(f_3_arg, v, r, ra, ra, id, r);
//   zm = pmx_ode_bdf(f_2_arg, v, r, ra, id, r);
//   zm = pmx_ode_bdf(f_1_arg, v, r, ra, r);
//   zm = pmx_ode_bdf(f_0_arg, v, r, ra);
  // ODE with control
  zm = pmx_ode_bdf_ctrl(f_12_arg, vd, rd, rad, 1e-6, 1e-6,
                                        100, mad, rvad, vad, iaad, raad, md,
                                        rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_11_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvad, vad, iaad,
                     raad, md, rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_10_arg, vd, rd, rad, 1e-6, 1e-6, 100, vad, iaad, raad,
                     md, rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_9_arg, vd, rd, rad, 1e-6, 1e-6, 100, iaad, raad, md,
                     rvd, vd, iad, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_8_arg, vd, rd, rad, 1e-6, 1e-6, 100, raad, md, rvd,
                     vd, iad, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_7_arg, vd, rd, rad, 1e-6, 1e-6, 100, md, rvd, vd, iad,
                     rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_6_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvd, vd, iad,
                     rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_5_arg, vd, rd, rad, 1e-6, 1e-6, 100, vd, iad, rad, id,
                     rd);
  zm = pmx_ode_bdf_ctrl(f_4_arg, vd, rd, rad, 1e-6, 1e-6, 100, iad, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_3_arg, vd, rd, rad, 1e-6, 1e-6, 100, rad, id, rd);
  zm = pmx_ode_bdf_ctrl(f_2_arg, vd, rd, rad, 1e-6, 1e-6, 100, id, rd);
  zm = pmx_ode_bdf_ctrl(f_1_arg, vd, rd, rad, 1e-6, 1e-6, 100, rd);
  zm = pmx_ode_bdf_ctrl(f_0_arg, vd, rd, rad, 1e-6, 1e-6, 100);
  zm = pmx_ode_bdf_ctrl(f_12_arg, v, rd, ra, 1e-6, 1e-6, 100, ma, rva, va, iaad,
                     raa, m, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_11_arg, v, rd, ra, 1e-6, 1e-6, 100, rva, va, iaad, raa,
                     m, rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_10_arg, v, rd, ra, 1e-6, 1e-6, 100, va, iaad, raa, m,
                     rv, v, iad, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_9_arg, v, rd, ra, 1e-6, 1e-6, 100, iaad, raa, m, rv, v,
                     iad, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_8_arg, v, rd, ra, 1e-6, 1e-6, 100, raa, m, rv, v, iad,
                     ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_7_arg, v, rd, ra, 1e-6, 1e-6, 100, m, rv, v, iad, ra,
                     id, r);
  zm = pmx_ode_bdf_ctrl(f_6_arg, v, rd, ra, 1e-6, 1e-6, 100, rv, v, iad, ra, id,
                     r);
  zm = pmx_ode_bdf_ctrl(f_5_arg, v, rd, ra, 1e-6, 1e-6, 100, v, iad, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_4_arg, v, rd, ra, 1e-6, 1e-6, 100, iad, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_3_arg, v, rd, ra, 1e-6, 1e-6, 100, ra, id, r);
  zm = pmx_ode_bdf_ctrl(f_2_arg, v, rd, ra, 1e-6, 1e-6, 100, id, r);
  zm = pmx_ode_bdf_ctrl(f_1_arg, v, rd, ra, 1e-6, 1e-6, 100, r);
  zm = pmx_ode_bdf_ctrl(f_0_arg, v, rd, ra, 1e-6, 1e-6, 100);
  r ~ normal(0, 1);
}
generated quantities {
  // ODE
  array[N] vector[N] zg = pmx_ode_bdf(f_12_arg, vd, rd, rad, mad, rvad, vad,
                                    iaad, raad, md, rvd, vd, iad, rad, id,
                                    rd);
  zg = pmx_ode_bdf(f_11_arg, vd, rd, rad, rvad, vad, iaad, raad, md, rvd, vd,
                 iad, rad, id, rd);
  zg = pmx_ode_bdf(f_10_arg, vd, rd, rad, vad, iaad, raad, md, rvd, vd, iad,
                 rad, id, rd);
  zg = pmx_ode_bdf(f_9_arg, vd, rd, rad, iaad, raad, md, rvd, vd, iad, rad, id,
                 rd);
  zg = pmx_ode_bdf(f_8_arg, vd, rd, rad, raad, md, rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf(f_7_arg, vd, rd, rad, md, rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf(f_6_arg, vd, rd, rad, rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf(f_5_arg, vd, rd, rad, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf(f_4_arg, vd, rd, rad, iad, rad, id, rd);
  zg = pmx_ode_bdf(f_3_arg, vd, rd, rad, rad, id, rd);
  zg = pmx_ode_bdf(f_2_arg, vd, rd, rad, id, rd);
  zg = pmx_ode_bdf(f_1_arg, vd, rd, rad, rd);
  zg = pmx_ode_bdf(f_0_arg, vd, rd, rad);
  zg = pmx_ode_bdf(f_12_arg, v, rd, ra, ma, rva, va, iaad, raa, m, rv, v, iad,
                 ra, id, r);
  zg = pmx_ode_bdf(f_11_arg, v, rd, ra, rva, va, iaad, raa, m, rv, v, iad, ra,
                 id, r);
  zg = pmx_ode_bdf(f_10_arg, v, rd, ra, va, iaad, raa, m, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf(f_9_arg, v, rd, ra, iaad, raa, m, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf(f_8_arg, v, rd, ra, raa, m, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf(f_7_arg, v, rd, ra, m, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf(f_6_arg, v, rd, ra, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf(f_5_arg, v, rd, ra, v, iad, ra, id, r);
  zg = pmx_ode_bdf(f_4_arg, v, rd, ra, iad, ra, id, r);
  zg = pmx_ode_bdf(f_3_arg, v, rd, ra, ra, id, r);
  zg = pmx_ode_bdf(f_2_arg, v, rd, ra, id, r);
  zg = pmx_ode_bdf(f_1_arg, v, rd, ra, r);
  zg = pmx_ode_bdf(f_0_arg, v, rd, ra);
  // zg = pmx_ode_bdf(f_12_arg, v, r, ra, ma, rva, va, iaad, raa, m, rv, v, iad,
  //                ra, id, r);
  // zg = pmx_ode_bdf(f_11_arg, v, r, ra, rva, va, iaad, raa, m, rv, v, iad, ra,
  //                id, r);
  // zg = pmx_ode_bdf(f_10_arg, v, r, ra, va, iaad, raa, m, rv, v, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_9_arg, v, r, ra, iaad, raa, m, rv, v, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_8_arg, v, r, ra, raa, m, rv, v, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_7_arg, v, r, ra, m, rv, v, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_6_arg, v, r, ra, rv, v, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_5_arg, v, r, ra, v, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_4_arg, v, r, ra, iad, ra, id, r);
  // zg = pmx_ode_bdf(f_3_arg, v, r, ra, ra, id, r);
  // zg = pmx_ode_bdf(f_2_arg, v, r, ra, id, r);
  // zg = pmx_ode_bdf(f_1_arg, v, r, ra, r);
  // zg = pmx_ode_bdf(f_0_arg, v, r, ra);
  // ODE with control
  zg = pmx_ode_bdf_ctrl(f_12_arg, vd, rd, rad, 1e-6, 1e-6,
                                        100, mad, rvad, vad, iaad, raad, md,
                                        rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_11_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvad, vad, iaad,
                     raad, md, rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_10_arg, vd, rd, rad, 1e-6, 1e-6, 100, vad, iaad, raad,
                     md, rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_9_arg, vd, rd, rad, 1e-6, 1e-6, 100, iaad, raad, md,
                     rvd, vd, iad, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_8_arg, vd, rd, rad, 1e-6, 1e-6, 100, raad, md, rvd,
                     vd, iad, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_7_arg, vd, rd, rad, 1e-6, 1e-6, 100, md, rvd, vd, iad,
                     rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_6_arg, vd, rd, rad, 1e-6, 1e-6, 100, rvd, vd, iad,
                     rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_5_arg, vd, rd, rad, 1e-6, 1e-6, 100, vd, iad, rad, id,
                     rd);
  zg = pmx_ode_bdf_ctrl(f_4_arg, vd, rd, rad, 1e-6, 1e-6, 100, iad, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_3_arg, vd, rd, rad, 1e-6, 1e-6, 100, rad, id, rd);
  zg = pmx_ode_bdf_ctrl(f_2_arg, vd, rd, rad, 1e-6, 1e-6, 100, id, rd);
  zg = pmx_ode_bdf_ctrl(f_1_arg, vd, rd, rad, 1e-6, 1e-6, 100, rd);
  zg = pmx_ode_bdf_ctrl(f_0_arg, vd, rd, rad, 1e-6, 1e-6, 100);
  zg = pmx_ode_bdf_ctrl(f_12_arg, v, rd, ra, 1e-6, 1e-6, 100, ma, rva, va, iaad,
                     raa, m, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_11_arg, v, rd, ra, 1e-6, 1e-6, 100, rva, va, iaad, raa,
                     m, rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_10_arg, v, rd, ra, 1e-6, 1e-6, 100, va, iaad, raa, m,
                     rv, v, iad, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_9_arg, v, rd, ra, 1e-6, 1e-6, 100, iaad, raa, m, rv, v,
                     iad, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_8_arg, v, rd, ra, 1e-6, 1e-6, 100, raa, m, rv, v, iad,
                     ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_7_arg, v, rd, ra, 1e-6, 1e-6, 100, m, rv, v, iad, ra,
                     id, r);
  zg = pmx_ode_bdf_ctrl(f_6_arg, v, rd, ra, 1e-6, 1e-6, 100, rv, v, iad, ra, id,
                     r);
  zg = pmx_ode_bdf_ctrl(f_5_arg, v, rd, ra, 1e-6, 1e-6, 100, v, iad, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_4_arg, v, rd, ra, 1e-6, 1e-6, 100, iad, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_3_arg, v, rd, ra, 1e-6, 1e-6, 100, ra, id, r);
  zg = pmx_ode_bdf_ctrl(f_2_arg, v, rd, ra, 1e-6, 1e-6, 100, id, r);
  zg = pmx_ode_bdf_ctrl(f_1_arg, v, rd, ra, 1e-6, 1e-6, 100, r);
  zg = pmx_ode_bdf_ctrl(f_0_arg, v, rd, ra, 1e-6, 1e-6, 100);
}
