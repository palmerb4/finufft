#include "defs.h"
#include "dirft.h"
#include <iostream>

// This is basically a port of dirft2d.f from CMCL package, except with
// the 1/nj prefactors for type-1 removed.

void dirft5d1(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu,
              BIGINT mv, BIGINT mw, CPX *f)
/* Direct computation of 5D type-1 nonuniform FFT. Interface same as finufft5d1.
c                        nj-1
c     f[k1,k2,k3,k4] =   SUM  c[j] exp(+-i (k1 x[j] + k2 y[j] + k3 z[j] + k4 p[j] + k5 q[j]))
c                        j=0
c
c     for -ms/2 <= k1 <= (ms-1)/2,  -mt/2 <= k2 <= (mt-1)/2,
          -mu/2 <= k3 <= (mu-1)/2,  -mv/2 <= k4 <= (mv-1)/2,
          -mw/2 <= k5 <= (mw-1)/2
c     The output array is in increasing k1 ordering (fast), then increasing
      k2 ordering (medium), then increasing k3 (fast), then increasing kf (fast). If iflag>0 the + sign is
c     used, otherwise the - sign is used, in the exponential.
*  Uses C++ complex type and winding trick.  Barnett 2/1/17
*/
{
  BIGINT k1min = -(ms / 2), k2min = -(mt / 2), k3min = -(mu / 2), k4min = -(mv / 2),
         k5min = -(mw / 2);          // integer divide
  BIGINT N = ms * mt * mu * mv * mw; // total # output modes
  for (BIGINT m = 0; m < N; ++m)
    f[m] = CPX(0, 0);               // it knows f is complex type
  for (BIGINT j = 0; j < nj; ++j) { // src pts
    CPX a1 = (iflag > 0) ? exp(IMA * x[j]) : exp(-IMA * x[j]);
    CPX a2 = (iflag > 0) ? exp(IMA * y[j]) : exp(-IMA * y[j]);
    CPX a3 = (iflag > 0) ? exp(IMA * z[j]) : exp(-IMA * z[j]);
    CPX a4 = (iflag > 0) ? exp(IMA * p[j]) : exp(-IMA * p[j]);
    CPX a5 = (iflag > 0) ? exp(IMA * q[j]) : exp(-IMA * q[j]);
    CPX sp1 = pow(a1, (FLT)k1min); // starting phase for most neg k1 freq
    CPX sp2 = pow(a2, (FLT)k2min);
    CPX sp3 = pow(a3, (FLT)k3min);
    CPX sp4 = pow(a4, (FLT)k4min);
    CPX p5 = pow(a5, (FLT)k5min);
    CPX cc = c[j]; // no 1/nj norm
    BIGINT m = 0;  // output pointer
    for (BIGINT m5 = 0; m5 < mw; ++m5) {
      CPX p4 = sp4;
      for (BIGINT m4 = 0; m4 < mv; ++m4) {
        CPX p3 = sp3;
        for (BIGINT m3 = 0; m3 < mu; ++m3) {
          CPX p2 = sp2;
          for (BIGINT m2 = 0; m2 < mt; ++m2) {
            CPX p1 = sp1;
            for (BIGINT m1 = 0; m1 < ms; ++m1) {
              f[m++] += cc * p1 * p2 * p3 * p4 * p5;
              p1 *= a1;
            }
            p2 *= a2;
          }
          p3 *= a3;
        }
        p4 *= a4;
      }
      p5 *= a5;
    }
  }
}

void dirft5d2(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu,
              BIGINT mv, BIGINT mw, CPX *f)
/* Direct computation of 5D type-2 nonuniform FFT. Interface same as finufft5d2

     c[j] =   SUM    f[k1,k2,k3,k3,k4] exp(+-i (k1 x[j] + k2 y[j] + k3 z[j] + k4 p[j] + k5 q[j]))
         k1,k2,k3,k4,k5
                            for j = 0,...,nj-1
    where sum is over -ms/2 <= k1 <= (ms-1)/2,  -mt/2 <= k2 <= (mt-1)/2,
                      -mu/2 <= k3 <= (mu-1)/2,  -mv/2 <= k4 <= (mv-1)/2,
                      -mw/2 <= k5 <= (mw-1)/2

    The input array is in increasing k1 ordering (fast), then increasing
    k2 ordering (medium), then increasing k3 (fast), then increasing k4 (fast), then increasing k5 (fast).
    If iflag>0 the + sign is used, otherwise the - sign is used, in the
    exponential.
    Uses C++ complex type and winding trick.  Barnett 2/1/17
*/
{
  BIGINT k1min = -(ms / 2), k2min = -(mt / 2), k3min = -(mu / 2), k4min = -(mv / 2),
         k5min = -(mw / 2); // integer divide
  for (BIGINT j = 0; j < nj; ++j) {
    CPX a1 = (iflag > 0) ? exp(IMA * x[j]) : exp(-IMA * x[j]);
    CPX a2 = (iflag > 0) ? exp(IMA * y[j]) : exp(-IMA * y[j]);
    CPX a3 = (iflag > 0) ? exp(IMA * z[j]) : exp(-IMA * z[j]);
    CPX a4 = (iflag > 0) ? exp(IMA * p[j]) : exp(-IMA * p[j]);
    CPX a5 = (iflag > 0) ? exp(IMA * q[j]) : exp(-IMA * q[j]);
    CPX sp1 = pow(a1, (FLT)k1min);
    CPX sp2 = pow(a2, (FLT)k2min);
    CPX sp3 = pow(a3, (FLT)k3min);
    CPX sp4 = pow(a4, (FLT)k4min);
    CPX p5 = pow(a5, (FLT)k5min);
    CPX cc = CPX(0, 0);
    BIGINT m = 0; // input pointer
    for (BIGINT m4 = 0; m4 < mv; ++m4) {
      CPX p4 = sp4;
      for (BIGINT m4 = 0; m4 < mv; ++m4) {
        CPX p3 = sp3;
        for (BIGINT m3 = 0; m3 < mu; ++m3) {
          CPX p2 = sp2;
          for (BIGINT m2 = 0; m2 < mt; ++m2) {
            CPX p1 = sp1;
            for (BIGINT m1 = 0; m1 < ms; ++m1) {
              cc += f[m++] * p1 * p2 * p3 * p4 * p5;
              p1 *= a1;
            }
            p2 *= a2;
          }
          p3 *= a3;
        }
        p4 *= a4;
      }
      p5 *= a5;
    }
    c[j] = cc;
  }
}

void dirft5d3(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int iflag, BIGINT nk, FLT *s, FLT *t, FLT *u,
              FLT *v, FLT *w, CPX *f)
/* Direct computation of 5D type-3 nonuniform FFT. Interface same as finufft5d3
c               nj-1
c     f[k]  =   SUM   c[j] exp(+-i (s[k] x[j] + t[k] y[j] + u[k] z[j] + v[k] p[j] + w[k] q[j]))
c               j=0
c                    for k = 0, ..., nk-1
c  If iflag>0 the + sign is used, otherwise the - sign is used, in the
c  exponential. Uses C++ complex type. Simple brute force.  Barnett 2/1/17
*/
{
  for (BIGINT k = 0; k < nk; ++k) {
    CPX ss = (iflag > 0) ? IMA * s[k] : -IMA * s[k];
    CPX tt = (iflag > 0) ? IMA * t[k] : -IMA * t[k];
    CPX uu = (iflag > 0) ? IMA * u[k] : -IMA * u[k];
    CPX vv = (iflag > 0) ? IMA * v[k] : -IMA * v[k];
    CPX ww = (iflag > 0) ? IMA * w[k] : -IMA * w[k];
    f[k] = CPX(0, 0);
    for (BIGINT j = 0; j < nj; ++j)
      f[k] += c[j] * exp(ss * x[j] + tt * y[j] + uu * z[j] + vv * p[j] + ww * q[j]);
  }
}
