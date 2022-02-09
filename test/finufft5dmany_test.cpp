#include <test_defs.h>
// this enforces recompilation, responding to SINGLE...
#include "directft/dirft5d.cpp"
using namespace std;

const char *help[] = {
    "Tester for FINUFFT in 5d, vectorized, all 3 types, either precision.",
    "",
    "Usage: finufft5dmany_test ntrans Nmodes1 Nmodes2 Nmodes3 Nmodes4 Nmodes5 Nsrc [tol [debug [spread_thread "
    "[maxbatchsize [spreadsort [upsampfac [errfail]]]]]]]",
    "\teg:\tfinufft5dmany_test 100 20 20 20 20 20 1e5 1e-3 1 0 0 2 0.0 1e-2",
    "\tnotes:\tif errfail present, exit code 1 if any error > errfail",
    NULL};
// Malleo 2019 based on Shih 2018. Tidied, extra args, Barnett 5/25/20.

int main(int argc, char *argv[]) {
  BIGINT M, N1, N2, N3, N4, N5; // M = # srcs, N1,N2,N3,N4,N5 = # modes
  int ntransf;              // # of vectors for "many" interface
  double tmp, tol = 1e-6;     // default
  double err, errfail = INFINITY, errmax = 0;
  nufft_opts opts;
  FINUFFT_DEFAULT_OPTS(&opts);
  // opts.fftw = FFTW_MEASURE;  // change from usual FFTW_ESTIMATE
  int isign = +1; // choose which exponential sign to test
  if (argc < 8 || argc > 15) {
    for (int i = 0; help[i]; ++i)
      fprintf(stderr, "%s\n", help[i]);
    return 2;
  }
  sscanf(argv[1], "%lf", &tmp);
  ntransf = (int)tmp;
  sscanf(argv[2], "%lf", &tmp);
  N1 = (BIGINT)tmp;
  sscanf(argv[3], "%lf", &tmp);
  N2 = (BIGINT)tmp;
  sscanf(argv[4], "%lf", &tmp);
  N3 = (BIGINT)tmp;
  sscanf(argv[5], "%lf", &tmp);
  N4 = (BIGINT)tmp;
  sscanf(argv[6], "%lf", &tmp);
  N5 = (BIGINT)tmp;
  sscanf(argv[7], "%lf", &tmp);
  M = (BIGINT)tmp;
  if (argc > 8)
    sscanf(argv[8], "%lf", &tol);
  if (argc > 9)
    sscanf(argv[9], "%d", &opts.debug);
  opts.spread_debug = (opts.debug > 1) ? 1 : 0; // see output from spreader
  if (argc > 10)
    sscanf(argv[10], "%d", &opts.spread_thread);
  if (argc > 11)
    sscanf(argv[11], "%d", &opts.maxbatchsize);
  if (argc > 12)
    sscanf(argv[12], "%d", &opts.spread_sort);
  if (argc > 13) {
    sscanf(argv[13], "%lf", &tmp);
    opts.upsampfac = (FLT)tmp;
  }
  if (argc > 14)
    sscanf(argv[14], "%lf", &errfail);

  cout << scientific << setprecision(15);
  BIGINT N = N1 * N2 * N3 * N4 * N5;

  FLT *x = (FLT *)malloc(sizeof(FLT) * M);           // NU pts x coords
  FLT *y = (FLT *)malloc(sizeof(FLT) * M);           // NU pts y coords
  FLT *z = (FLT *)malloc(sizeof(FLT) * M);           // NU pts z coords
  FLT *p = (FLT *)malloc(sizeof(FLT) * M);           // NU pts p coords
  FLT *q = (FLT *)malloc(sizeof(FLT) * M);           // NU pts q coords
  CPX *c = (CPX *)malloc(sizeof(CPX) * M * ntransf); // strengths
  CPX *F = (CPX *)malloc(sizeof(CPX) * N * ntransf); // mode ampls

#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT j = 0; j < M; ++j) {
      x[j] = M_PI * randm11r(&se);
      y[j] = M_PI * randm11r(&se);
      z[j] = M_PI * randm11r(&se);
      p[j] = M_PI * randm11r(&se);
      q[j] = M_PI * randm11r(&se);
    }
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT j = 0; j < ntransf * M; ++j) {
      c[j] = crandm11r(&se);
    }
  }

  printf("test 5d1 many vs repeated single: ------------------------------------\n");
  CNTime timer;
  timer.start();
  int ier = FINUFFT5D1MANY(ntransf, M, x, y, z, p, q, c, isign, tol, N1, N2, N3, N4, N5, F, &opts);
  double ti = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("ntr=%d: %lld NU pts to (%lld,%lld,%lld,%lld,%lld) modes in %.3g s \t%.3g NU pts/s\n", ntransf, (long long)M,
           (long long)N1, (long long)N2, (long long)N3, (long long)N4, (long long)N5, ti, ntransf * M / ti);

  int i = ntransf - 1; // choose a data to check
  BIGINT nt1 = (BIGINT)(0.37 * N1), nt2 = (BIGINT)(0.26 * N2), nt3 = (BIGINT)(-0.39 * N3),
         nt4 = (BIGINT)(-0.39 * N4), nt5 = (BIGINT)(-0.39 * N5); // choose some mode index to check
  CPX Ft = CPX(0, 0), J = IMA * (FLT)isign;
  for (BIGINT j = 0; j < M; ++j)
    Ft += c[j + i * M] * exp(J * (nt1 * x[j] + nt2 * y[j] + nt3 * z[j] + nt4 * p[j] + nt5 * q[j])); // crude direct
  BIGINT it = N1 / 2 + nt1 + N1 * (N2 / 2 + nt2) + N1 * N2 * (N3 / 2 + nt3) + N1 * N2 * N3 * (N4 / 2 + nt4) + N1 * N2 * N3 * N4 * (N5 / 2 + nt5);           // index in complex F as 1d array
  err = abs(Ft - F[it + i * N]) / infnorm(N, F + i * N);
  errmax = max(err, errmax);
  printf("\tone mode: rel err in F[%lld,%lld,%lld,%lld,%lld] of trans#%d is %.3g\n", (long long)nt1, (long long)nt2,
         (long long)nt3, (long long)nt4, (long long)nt5, i, err);

  // compare the result with FINUFFT5D1
  FFTW_FORGET_WISDOM();
  nufft_opts simpleopts = opts;
  simpleopts.debug = 0; // don't output timing for calls of FINUFFT5D1
  simpleopts.spread_debug = 0;

  CPX *cstart;
  CPX *Fstart;
  CPX *F_5d1 = (CPX *)malloc(sizeof(CPX) * N * ntransf);
  timer.restart();
  for (int k = 0; k < ntransf; ++k) {
    cstart = c + k * M;
    Fstart = F_5d1 + k * N;
    ier = FINUFFT5D1(M, x, y, z, p, q, cstart, isign, tol, N1, N2, N3, N4, N5, Fstart, &simpleopts);
  }
  double t = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("%d of: %lld NU pts to (%lld,%lld,%lld,%lld,%lld) modes in %.3g s  \t%.3g NU pts/s\n", ntransf, (long long)M,
           (long long)N1, (long long)N2, (long long)N3, (long long)N4, (long long)N5, t, ntransf * M / t);
  printf("\t\t\tspeedup \t T_FINUFFT5D1 / T_finufft5d1many = %.3g\n", t / ti);

  // Check accuracy (worst over the ntransf)
  double maxerror = 0.0;
  for (int k = 0; k < ntransf; ++k)
    maxerror = max(maxerror, (double)relerrtwonorm(N, F_5d1 + k * N, F + k * N));
  errmax = max(maxerror, errmax);
  printf("\tconsistency check: sup ( ||f_many-f||_2 / ||f||_2 ) =  %.3g\n", maxerror);
  free(F_5d1);

  printf("test 5d2 many vs repeated single: ------------------------------------\n");
#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT m = 0; m < N * ntransf; ++m)
      F[m] = crandm11r(&se);
  }
  FFTW_FORGET_WISDOM();
  timer.restart();
  ier = FINUFFT5D2MANY(ntransf, M, x, y, z, p, q, c, isign, tol, N1, N2, N3, N4, N5, F, &opts);
  ti = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("ntr=%d: (%lld,%lld,%lld,%lld,%lld) modes to %lld NU pts in %.3g s \t%.3g NU pts/s\n", ntransf, (long long)N1,
           (long long)N2, (long long)N3, (long long)N4, (long long)N5, (long long)M, ti, ntransf * M / ti);

  i = ntransf - 1;   // choose a data to check
  BIGINT jt = M / 2; // check arbitrary choice of one targ pt
  CPX ct = CPX(0, 0);
  BIGINT m = 0;
  for (BIGINT m5 = -(N5 / 2); m5 <= (N5 - 1) / 2; ++m5) {
    for (BIGINT m4 = -(N4 / 2); m4 <= (N4 - 1) / 2; ++m4) {
      for (BIGINT m3 = -(N3 / 2); m3 <= (N3 - 1) / 2; ++m3) {
        for (BIGINT m2 = -(N2 / 2); m2 <= (N2 - 1) / 2; ++m2) { // loop in correct order over F
          for (BIGINT m1 = -(N1 / 2); m1 <= (N1 - 1) / 2; ++m1) {
            ct += F[i * N + m++] * exp(J * (m1 * x[jt] + m2 * y[jt] + m3 * z[jt] + m4 * p[jt] + m5 * q[jt])); // crude direct
          }
        }
      }
    }
  }
  err = abs(ct - c[jt + i * M]) / infnorm(M, c + i * M);
  errmax = max(err, errmax);
  printf("\tone targ: rel err in c[%lld] of trans#%d is %.3g\n", (long long)jt, i, err);

  FFTW_FORGET_WISDOM();
  // compare the result with FINUFFT5D2...
  CPX *c_5d2 = (CPX *)malloc(sizeof(CPX) * M * ntransf);
  timer.restart();
  for (int k = 0; k < ntransf; ++k) {
    cstart = c_5d2 + k * M;
    Fstart = F + k * N;
    ier = FINUFFT5D2(M, x, y, z, p, q, cstart, isign, tol, N1, N2, N3, N4, N5, Fstart, &simpleopts);
  }
  t = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("%d of: (%lld,%lld,%lld,%lld,%lld) modes to %lld NU pts in %.3g s \t%.3g NU pts/s\n", ntransf, (long long)N1,
           (long long)N2, (long long)N3, (long long)N4, (long long)N5, (long long)M, t, ntransf * M / t);
  printf("\t\t\tspeedup \t T_FINUFFT5D2 / T_finufft5d2many = %.3g\n", t / ti);

  maxerror = 0.0; // worst error over the ntransf
  for (int k = 0; k < ntransf; ++k)
    maxerror = max(maxerror, (double)relerrtwonorm(M, c_5d2 + k * M, c + k * M));
  errmax = max(maxerror, errmax);
  printf("\tconsistency check: sup ( ||c_many-c||_2 / ||c||_2 ) =  %.3g\n", maxerror);
  free(c_5d2);

  printf("test 5d3 many vs repeated single: ------------------------------------\n");
  FFTW_FORGET_WISDOM();
  // reuse the strengths c, interpret N as number of targs:
#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT j = 0; j < M; ++j) {
      x[j] = 2.0 + M_PI * randm11r(&se);  // new x_j srcs, offset from origin
      y[j] = -3.0 + M_PI * randm11r(&se); // " y_j
      z[j] = 1.0 + M_PI * randm11r(&se);  // " z_j
      p[j] = 0.5 + M_PI * randm11r(&se);  // " p_j
      q[j] = 0.1 + M_PI * randm11r(&se);  // " q_j
    }
  }
  FLT *s_freq = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (1-cmpt)
  FLT *t_freq = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (2-cmpt)
  FLT *u_freq = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (3-cmpt)
  FLT *v_freq = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (4-cmpt)
  FLT *w_freq = (FLT *)malloc(sizeof(FLT) * N); // targ freqs (5-cmpt)
  FLT S1 = (FLT)N1 / 2;                         // choose freq range sim to type 1
  FLT S2 = (FLT)N2 / 2;
  FLT S3 = (FLT)N3 / 2;
  FLT S4 = (FLT)N4 / 2;
  FLT S5 = (FLT)N5 / 2;

#pragma omp parallel
  {
    unsigned int se = MY_OMP_GET_THREAD_NUM();
#pragma omp for schedule(static, TEST_RANDCHUNK)
    for (BIGINT k = 0; k < N; ++k) {
      s_freq[k] = S1 * (1.7 + randm11r(&se)); // S*(1.7 + k/(FLT)N); // offset the freqs
      t_freq[k] = S2 * (-0.5 + randm11r(&se));
      u_freq[k] = S3 * (0.9 + randm11r(&se));
      v_freq[k] = S4 * (0.1 + randm11r(&se));
      w_freq[k] = S5 * (0.4 + randm11r(&se));
    }
  }

  timer.restart();
  ier = FINUFFT5D3MANY(ntransf, M, x, y, z, p, q, c, isign, tol, N, s_freq, t_freq, u_freq, v_freq, w_freq, F, &opts);
  ti = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("ntr=%d: %lld NU to %lld NU in %.3g s \t%.3g tot NU pts/s\n", ntransf, (long long)M, (long long)N, ti,
           ntransf * (M + N) / ti);

  i = ntransf - 1;   // choose a transform to check
  BIGINT kt = N / 4; // check arbitrary choice of one targ pt
  Ft = CPX(0, 0);
  for (BIGINT j = 0; j < M; ++j)
    Ft += c[i * M + j] * exp(J * (s_freq[kt] * x[j] + t_freq[kt] * y[j] + u_freq[kt] * z[j] + v_freq[kt] * p[j] + w_freq[kt] * q[j]));
  err = abs(Ft - F[kt + i * N]) / infnorm(N, F + i * N);
  errmax = max(err, errmax);
  printf("\t one targ: rel err in F[%lld] of trans#%d is %.3g\n", (long long)kt, i, err);

  FFTW_FORGET_WISDOM();
  // compare the result with FINUFFT5D3...
  CPX *f_5d3 = (CPX *)malloc(sizeof(CPX) * N * ntransf);
  timer.restart();
  for (int k = 0; k < ntransf; ++k) {
    Fstart = f_5d3 + k * N;
    cstart = c + k * M;
    ier = FINUFFT5D3(M, x, y, z, p, q, cstart, isign, tol, N, s_freq, t_freq, u_freq, v_freq, w_freq, Fstart, &simpleopts);
  }
  t = timer.elapsedsec();
  if (ier > 1) {
    printf("error (ier=%d)!\n", ier);
    return ier;
  } else
    printf("%d of: %lld NU to %lld NU in %.3g s   \t%.3g tot NU pts/s\n", ntransf, (long long)M, (long long)N, t,
           ntransf * (M + N) / t);
  printf("\t\t\tspeedup \t T_FINUFFT5D3 / T_finufft5d3many = %.3g\n", t / ti);

  maxerror = 0.0; // worst error over the ntransf
  for (int k = 0; k < ntransf; ++k)
    maxerror = max(maxerror, (double)relerrtwonorm(N, f_5d3 + k * N, F + k * N));
  errmax = max(maxerror, errmax);
  printf("\tconsistency check: sup ( ||f_many-f||_2 / ||f||_2 ) =  %.3g\n", maxerror);
  free(f_5d3);

  free(x);
  free(y);
  free(z);
  free(p);
  free(q);
  free(c);
  free(F);
  free(s_freq);
  free(t_freq);
  free(u_freq);
  free(v_freq);
  free(w_freq);
  return (errmax > errfail);
}
