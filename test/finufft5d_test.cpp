#include <test_defs.h>
// this enforces recompilation, responding to SINGLE...
#include "directft/dirft5d.cpp"
using namespace std;

const char *help[] = {"Tester for FINUFFT in 5d, all 3 types, either precision.",
                      "",
                      "Usage: finufft5d_test Nmodes1 Nmodes2 Nmodes3 Nmodes4 Nmodes5 Nsrc [tol [debug [spread_sort "
                      "[upsampfac [errfail]]]]]",
                      "\teg:\tfinufft5d_test 100 20 20 20 20 1e6 1e-12 0 2 0.0 1e-11",
                      "\tnotes:\tif errfail present, exit code 1 if any error > errfail",
                      NULL};
// Barnett 2/2/17 onwards.

int main(int argc, char *argv[]) {
  BIGINT M, N1, N2, N3, N4, N5; // M = # srcs, N1,N2,N3,N4,N5 = # modes
  double tmp, tol = 1e-6;       // default
  double err, errfail = INFINITY, errmax = 0;
  nufft_opts opts;
  FINUFFT_DEFAULT_OPTS(&opts);
  // opts.fftw = FFTW_MEASURE;  // change from usual FFTW_ESTIMATE
  // opts.spread_max_sp_size = 3e4; // override test
  // opts.spread_nthr_atomic = 15;  // "
  int isign = +1; // choose which exponential sign to test
  if (argc < 7 || argc > 12) {
    for (int i = 0; help[i]; ++i)
      fprintf(stderr, "%s\n", help[i]);
    return 2;
  }
  sscanf(argv[1], "%lf", &tmp);
  N1 = (BIGINT)tmp;
  sscanf(argv[2], "%lf", &tmp);
  N2 = (BIGINT)tmp;
  sscanf(argv[3], "%lf", &tmp);
  N3 = (BIGINT)tmp;
  sscanf(argv[4], "%lf", &tmp);
  N4 = (BIGINT)tmp;
  sscanf(argv[5], "%lf", &tmp);
  N5 = (BIGINT)tmp;
  sscanf(argv[6], "%lf", &tmp);
  M = (BIGINT)tmp;
  if (argc > 7)
    sscanf(argv[7], "%lf", &tol);
  if (argc > 8)
    sscanf(argv[8], "%d", &opts.debug);         // can be 0,1 or 2
  opts.spread_debug = (opts.debug > 1) ? 1 : 0; // see output from spreader
  if (argc > 9)
    sscanf(argv[9], "%d", &opts.spread_sort);
  if (argc > 10) {
    sscanf(argv[10], "%lf", &tmp);
    opts.upsampfac = (FLT)tmp;
  }
  if (argc > 11)
    sscanf(argv[11], "%lf", &errfail);

  cout << scientific << setprecision(15);
  BIGINT N = N1 * N2 * N3 * N4 * N5;

  FLT *x = (FLT *)malloc(sizeof(FLT) * M); // NU pts x coords
  FLT *y = (FLT *)malloc(sizeof(FLT) * M); // NU pts y coords
  FLT *z = (FLT *)malloc(sizeof(FLT) * M); // NU pts z coords
  FLT *p = (FLT *)malloc(sizeof(FLT) * M); // NU pts p coords
  FLT *q = (FLT *)malloc(sizeof(FLT) * M); // NU pts q coords
  CPX *c = (CPX *)malloc(sizeof(CPX) * M); // strengths
  CPX *F = (CPX *)malloc(sizeof(CPX) * N); // mode ampls
#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM(); // needed for parallel random #s
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT j = 0; j < M; ++j) {
      x[j] = M_PI * randm11r(&se);
      y[j] = M_PI * randm11r(&se);
      z[j] = M_PI * randm11r(&se);
      p[j] = M_PI * randm11r(&se);
      q[j] = M_PI * randm11r(&se);
      c[j] = crandm11r(&se);
    }
  }

  printf("test 5d type 1:\n"); // -------------- type 1
  CNTime timer;
  timer.start();
  int ier = FINUFFT5D1(M, x, y, z, p, q, c, isign, tol, N1, N2, N3, N4, N5, F, &opts);
  double ti = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("     %lld NU pts to (%lld,%lld,%lld,%lld,%lld) modes in %.3g s \t%.3g NU pts/s\n", (long long)M,
           (long long)N1, (long long)N2, (long long)N3, (long long)N4, (long long)N5, ti, M / ti);

  BIGINT nt1 = (BIGINT)(0.37 * N1), nt2 = (BIGINT)(0.26 * N2), nt3 = (BIGINT)(-0.39 * N3), nt4 = (BIGINT)(-0.39 * N4),
         nt5 = (BIGINT)(-0.39 * N5); // choose mode to check
  FLT Ftr = 0, Fti = 0;              // crude direct...
#pragma omp parallel for schedule(static, TEST_RANDCHUNK) reduction(+ : Ftr, Fti)
  for (BIGINT j = 0; j < M; ++j) { // Ft += c[j] * exp(J*(nt1*x[j] + nt2*y[j] + nt3*z[j] + nt4*p[j] + nt5*q[j]))
    FLT theta = (FLT)isign * (nt1 * x[j] + nt2 * y[j] + nt3 * z[j] + nt4 * p[j] + nt5 * q[j]), co = cos(theta),
        si = sin(theta);
    Ftr += real(c[j]) * co - imag(c[j]) * si; // cpx arith by hand
    Fti += imag(c[j]) * co + real(c[j]) * si;
  }
  // index in complex F as 1d array...
  BIGINT it = N1 / 2 + nt1 + N1 * (N2 / 2 + nt2) + N1 * N2 * (N3 / 2 + nt3) + N1 * N2 * N3 * (N4 / 2 + nt4) +
              N1 * N2 * N3 * N4 * (N5 / 2 + nt5);
  err = abs(Ftr + IMA * Fti - F[it]) / infnorm(N, F);
  errmax = max(err, errmax);
  printf("\tone mode: rel err in F[%lld,%lld,%lld,%lld,%lld] is %.3g\n", (long long)nt1, (long long)nt2, (long long)nt3,
         (long long)nt4, (long long)nt5, err);
  if ((int64_t)M * N <= TEST_BIGPROB) { // also check vs full direct eval
    CPX *Ft = (CPX *)malloc(sizeof(CPX) * N);
    dirft5d1(M, x, y, z, p, q, c, isign, N1, N2, N3, N4, N5, Ft);
    err = relerrtwonorm(N, Ft, F);
    errmax = max(err, errmax);
    printf("\tdirft5d: rel l2-err of result F is %.3g\n", err);
    free(Ft);
  }

  printf("test 5d type 2:\n"); // -------------- type 2
#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT m = 0; m < N; ++m)
      F[m] = crandm11r(&se);
  }
  timer.restart();
  ier = FINUFFT5D2(M, x, y, z, p, q, c, isign, tol, N1, N2, N3, N4, N5, F, &opts);
  ti = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("     (%lld,%lld,%lld,%lld,%lld) modes to %lld NU pts in %.3g s \t%.3g NU pts/s\n", (long long)N1,
           (long long)N2, (long long)N3, (long long)N4, (long long)N5, (long long)M, ti, M / ti);

  BIGINT jt = M / 2; // check arbitrary choice of one targ pt
  CPX ct = CPX(0, 0);
  BIGINT m = 0;
  for (BIGINT m5 = -(N5 / 2); m5 <= (N5 - 1) / 2; ++m5) // loop in F order
    for (BIGINT m4 = -(N4 / 2); m4 <= (N4 - 1) / 2; ++m4)
      for (BIGINT m3 = -(N3 / 2); m3 <= (N3 - 1) / 2; ++m3)
        for (BIGINT m2 = -(N2 / 2); m2 <= (N2 - 1) / 2; ++m2)
          for (BIGINT m1 = -(N1 / 2); m1 <= (N1 - 1) / 2; ++m1)
            ct += F[m++] * exp(IMA * (FLT)isign * (m1 * x[jt] + m2 * y[jt] + m3 * z[jt] + m4 * p[jt] + m5 * q[jt]));
  err = abs(ct - c[jt]) / infnorm(M, c);
  errmax = max(err, errmax);
  printf("\tone targ: rel err in c[%lld] is %.3g\n", (long long)jt, err);
  if ((int64_t)M * N <= TEST_BIGPROB) { // also full direct eval
    CPX *ct = (CPX *)malloc(sizeof(CPX) * M);
    dirft5d2(M, x, y, z, p, q, ct, isign, N1, N2, N3, N4, N5, F);
    err = relerrtwonorm(M, ct, c);
    errmax = max(err, errmax);
    printf("\tdirft5d: rel l2-err of result c is %.3g\n", err);
    free(ct);
  }

  printf("test 5d type 3:\n"); // -------------- type 3
                               // reuse the strengths c, interpret N as number of targs:
#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT j = 0; j < M; ++j) {
      x[j] = 2.0 + M_PI * randm11r(&se);  // new x_j srcs, offset from origin
      y[j] = -3.0 + M_PI * randm11r(&se); // " y_j
      z[j] = 1.0 + M_PI * randm11r(&se);  // " z_j
      p[j] = 0.0 + M_PI * randm11r(&se);  // " p_j
      q[j] = 0.1 + M_PI * randm11r(&se);  // " q_j
    }
  }
  FLT *s = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (1-cmpt)
  FLT *t = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (2-cmpt)
  FLT *u = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (3-cmpt)
  FLT *v = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (4-cmpt)
  FLT *w = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (5-cmpt)
  FLT S1 = (FLT)N1 / 2;                    // choose freq range sim to type 1
  FLT S2 = (FLT)N2 / 2;
  FLT S3 = (FLT)N3 / 2;
  FLT S4 = (FLT)N4 / 2;
  FLT S5 = (FLT)N5 / 2;
#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT k = 0; k < N; ++k) {
      s[k] = S1 * (1.7 + randm11r(&se)); // S*(1.7 + k/(FLT)N); // offset the freqs
      t[k] = S2 * (-0.5 + randm11r(&se));
      u[k] = S3 * (0.9 + randm11r(&se));
      v[k] = S4 * (0.5 + randm11r(&se));
      w[k] = S5 * (0.1 + randm11r(&se));
    }
  }
  timer.restart();
  ier = FINUFFT5D3(M, x, y, z, p, q, c, isign, tol, N, s, t, u, v, w, F, &opts);
  ti = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("\t%lld NU to %lld NU in %.3g s         \t%.3g tot NU pts/s\n", (long long)M, (long long)N, ti,
           (M + N) / ti);

  BIGINT kt = N / 2; // check arbitrary choice of one targ pt
  Ftr = 0, Fti = 0;  // crude direct...
#pragma omp parallel for schedule(static, TEST_RANDCHUNK) reduction(+ : Ftr, Fti)
  for (BIGINT j = 0; j < M;
       ++j) { // Ft += c[j] * exp(IMA*(FLT)isign*(s[kt]*x[j] + t[kt]*y[j] + u[kt]*z[j] + v[kt]*p[j] + w[kt]*q[j))
    FLT theta = (FLT)isign * (s[kt] * x[j] + t[kt] * y[j] + u[kt] * z[j] + v[kt] * p[j] + w[kt] * q[j]),
        co = cos(theta), si = sin(theta);
    Ftr += real(c[j]) * co - imag(c[j]) * si; // cpx arith by hand
    Fti += imag(c[j]) * co + real(c[j]) * si;
  }
  err = abs(Ftr + IMA * Fti - F[kt]) / infnorm(N, F);
  errmax = max(err, errmax);
  printf("\tone targ: rel err in F[%lld] is %.3g\n", (long long)kt, err);
  if (((int64_t)M) * N <= TEST_BIGPROB) { // also full direct eval
    CPX *Ft = (CPX *)malloc(sizeof(CPX) * N);
    dirft5d3(M, x, y, z, p, q, c, isign, N, s, t, u, v, w, Ft); // writes to F
    err = relerrtwonorm(N, Ft, F);
    errmax = max(err, errmax);
    printf("\tdirft5d: rel l2-err of result F is %.3g\n", err);
    // cout<<"s t u, v, F, Ft, F/Ft:\n"; for (int k=0;k<N;++k) cout<<s[k]<<" "<<t[k]<<" "<<u[k]<<" "<<v[k]<<",
    // "<<F[k]<<",\t"<<Ft[k]<<",\t"<<F[k]/Ft[k]<<endl;
    free(Ft);
  }

  free(x);
  free(y);
  free(z);
  free(p);
  free(q);
  free(c);
  free(F);
  free(s);
  free(t);
  free(u);
  free(v);
  free(w);
  return (errmax > errfail);
}
