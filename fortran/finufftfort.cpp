#include <dataTypes.h>
#include <finufft_eitherprec.h>
#include <finufft_plan_eitherprec.h>
#include <nufft_opts.h>

/* C++ layer for calling FINUFFT from fortran, in f77 style + derived type
   for the nufft_opts C-struct. The ptr to finufft_plan is passed as an "opaque"
   pointer (as in FFTW3 legacy fortran interface).
   Note the trailing underscore name-mangle which is not needed from fortran.

   Note our typedefs:
   FLT = double (or float, depending on compilation precision)
   CPX = double complex (or float complex, depending on compilation precision)
   BIGINT = int64 (integer*8 in fortran)

   Make sure you call this library with matching fortran types

   For a demo see: examples/simple1d1.f

   Barnett 2/17/17. Single prec 4/5/17. Libin Lu & Alex Barnett, May 2020.
   Garrett Wright dual-prec 6/28/20.
*/

// macros to set distinct names for each precision...
#ifdef SINGLE
#define FINUFFT_MAKEPLAN_ finufftf_makeplan_
#define FINUFFT_SETPTS_ finufftf_setpts_
#define FINUFFT_EXECUTE_ finufftf_execute_
#define FINUFFT_DESTROY_ finufftf_destroy_
#define FINUFFT_DEFAULT_OPTS_ finufftf_default_opts_
#define FINUFFT1D1_ finufftf1d1_
#define FINUFFT1D1MANY_ finufftf1d1many_
#define FINUFFT1D2_ finufftf1d2_
#define FINUFFT1D2MANY_ finufftf1d2many_
#define FINUFFT1D3_ finufftf1d3_
#define FINUFFT1D3MANY_ finufftf1d3many_
#define FINUFFT2D1_ finufftf2d1_
#define FINUFFT2D1MANY_ finufftf2d1many_
#define FINUFFT2D2_ finufftf2d2_
#define FINUFFT2D2MANY_ finufftf2d2many_
#define FINUFFT2D3_ finufftf2d3_
#define FINUFFT2D3MANY_ finufftf2d3many_
#define FINUFFT3D1_ finufftf3d1_
#define FINUFFT3D1MANY_ finufftf3d1many_
#define FINUFFT3D2_ finufftf3d2_
#define FINUFFT3D2MANY_ finufftf3d2many_
#define FINUFFT3D3_ finufftf3d3_
#define FINUFFT3D3MANY_ finufftf3d3many_
#define FINUFFT4D1_ finufftf4d1_
#define FINUFFT4D1MANY_ finufftf4d1many_
#define FINUFFT4D2_ finufftf4d2_
#define FINUFFT4D2MANY_ finufftf4d2many_
#define FINUFFT4D3_ finufftf4d3_
#define FINUFFT4D3MANY_ finufftf4d3many_
#define FINUFFT5D1_ finufftf5d1_
#define FINUFFT5D1MANY_ finufftf5d1many_
#define FINUFFT5D2_ finufftf5d2_
#define FINUFFT5D2MANY_ finufftf5d2many_
#define FINUFFT5D3_ finufftf5d3_
#define FINUFFT5D3MANY_ finufftf5d3many_
#else
#define FINUFFT_MAKEPLAN_ finufft_makeplan_
#define FINUFFT_SETPTS_ finufft_setpts_
#define FINUFFT_EXECUTE_ finufft_execute_
#define FINUFFT_DESTROY_ finufft_destroy_
#define FINUFFT_DEFAULT_OPTS_ finufft_default_opts_
#define FINUFFT1D1_ finufft1d1_
#define FINUFFT1D1MANY_ finufft1d1many_
#define FINUFFT1D2_ finufft1d2_
#define FINUFFT1D2MANY_ finufft1d2many_
#define FINUFFT1D3_ finufft1d3_
#define FINUFFT1D3MANY_ finufft1d3many_
#define FINUFFT2D1_ finufft2d1_
#define FINUFFT2D1MANY_ finufft2d1many_
#define FINUFFT2D2_ finufft2d2_
#define FINUFFT2D2MANY_ finufft2d2many_
#define FINUFFT2D3_ finufft2d3_
#define FINUFFT2D3MANY_ finufft2d3many_
#define FINUFFT3D1_ finufft3d1_
#define FINUFFT3D1MANY_ finufft3d1many_
#define FINUFFT3D2_ finufft3d2_
#define FINUFFT3D2MANY_ finufft3d2many_
#define FINUFFT3D3_ finufft3d3_
#define FINUFFT3D3MANY_ finufft3d3many_
#define FINUFFT4D1_ finufft4d1_
#define FINUFFT4D1MANY_ finufft4d1many_
#define FINUFFT4D2_ finufft4d2_
#define FINUFFT4D2MANY_ finufft4d2many_
#define FINUFFT4D3_ finufft4d3_
#define FINUFFT4D3MANY_ finufft4d3many_
#define FINUFFT5D1_ finufft5d1_
#define FINUFFT5D1MANY_ finufft5d1many_
#define FINUFFT5D2_ finufft5d2_
#define FINUFFT5D2MANY_ finufft5d2many_
#define FINUFFT5D3_ finufft5d3_
#define FINUFFT5D3MANY_ finufft5d3many_
#endif

#ifdef __cplusplus
extern "C" {
#endif

// --------------------- guru interface from fortran ------------------------
void FINUFFT_MAKEPLAN_(int *type, int *n_dims, BIGINT *n_modes, int *iflag, int *n_transf, FLT *tol, FINUFFT_PLAN *plan,
                       nufft_opts *o, int *ier) {
  if (!plan)
    fprintf(stderr, "%s fortran: plan must be allocated as at least the size of a C pointer (usually 8 bytes)!\n",
            __func__);
  else {
    // pass o whether it's a NULL or pointer to a fortran-allocated nufft_opts:
    *ier = FINUFFT_MAKEPLAN(*type, *n_dims, n_modes, *iflag, *n_transf, *tol, plan, o);
  }
}

void FINUFFT_SETPTS_(FINUFFT_PLAN *plan, BIGINT *M, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, BIGINT *nk, FLT *s,
                     FLT *t, FLT *u, FLT *v, FLT *w, int *ier) {
  if (!*plan) {
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
    return;
  }
  int nk_safe = 0; // catches the case where user passes NULL in
  if (nk)
    nk_safe = *nk;
  *ier = FINUFFT_SETPTS(*plan, *M, xj, yj, zj, pj, qj, nk_safe, s, t, u, v, w);
}

void FINUFFT_EXECUTE_(FINUFFT_PLAN *plan, CPX *weights, CPX *result, int *ier) {
  if (!plan)
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
  else
    *ier = FINUFFT_EXECUTE(*plan, weights, result);
}

void FINUFFT_DESTROY_(FINUFFT_PLAN *plan, int *ier) {
  if (!plan)
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
  else
    *ier = FINUFFT_DESTROY(*plan);
}

// ------------ use FINUFFT to set the default options ---------------------
// (Note the nufft_opts is created in f90-style derived types, not here)
void FINUFFT_DEFAULT_OPTS_(nufft_opts *o) {
  if (!o)
    fprintf(stderr, "%s fortran: opts must be allocated!\n", __func__);
  else
    // o is a ptr to already-allocated fortran nufft_opts derived type...
    FINUFFT_DEFAULT_OPTS(o);
}

// -------------- simple and many-vector interfaces --------------------
// --- 1D ---
void FINUFFT1D1_(BIGINT *nj, FLT *xj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT1D1(*nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

void FINUFFT1D1MANY_(int *ntransf, BIGINT *nj, FLT *xj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, CPX *fk,
                     nufft_opts *o, int *ier) {
  *ier = FINUFFT1D1MANY(*ntransf, *nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

void FINUFFT1D2_(BIGINT *nj, FLT *xj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT1D2(*nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

void FINUFFT1D2MANY_(int *ntransf, BIGINT *nj, FLT *xj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, CPX *fk,
                     nufft_opts *o, int *ier) {
  *ier = FINUFFT1D2MANY(*ntransf, *nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

void FINUFFT1D3_(BIGINT *nj, FLT *x, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s, CPX *f, nufft_opts *o,
                 int *ier) {
  *ier = FINUFFT1D3(*nj, x, c, *iflag, *eps, *nk, s, f, o);
}

void FINUFFT1D3MANY_(int *ntransf, BIGINT *nj, FLT *x, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s, CPX *f,
                     nufft_opts *o, int *ier) {
  *ier = FINUFFT1D3MANY(*ntransf, *nj, x, c, *iflag, *eps, *nk, s, f, o);
}

// --- 2D ---
void FINUFFT2D1_(BIGINT *nj, FLT *xj, FLT *yj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt, CPX *fk,
                 nufft_opts *o, int *ier) {
  *ier = FINUFFT2D1(*nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}
void FINUFFT2D1MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt,
                     CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT2D1MANY(*ntransf, *nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}

void FINUFFT2D2_(BIGINT *nj, FLT *xj, FLT *yj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt, CPX *fk,
                 nufft_opts *o, int *ier) {
  *ier = FINUFFT2D2(*nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}
void FINUFFT2D2MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt,
                     CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT2D2MANY(*ntransf, *nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}

void FINUFFT2D3_(BIGINT *nj, FLT *x, FLT *y, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s, FLT *t, CPX *f,
                 nufft_opts *o, int *ier) {
  *ier = FINUFFT2D3(*nj, x, y, c, *iflag, *eps, *nk, s, t, f, o);
}

void FINUFFT2D3MANY_(int *ntransf, BIGINT *nj, FLT *x, FLT *y, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s, FLT *t,
                     CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT2D3MANY(*ntransf, *nj, x, y, c, *iflag, *eps, *nk, s, t, f, o);
}

// --- 3D ---
void FINUFFT3D1_(BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt,
                 BIGINT *mu, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT3D1(*nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

void FINUFFT3D1MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms,
                     BIGINT *mt, BIGINT *mu, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT3D1MANY(*ntransf, *nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

void FINUFFT3D2_(BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt,
                 BIGINT *mu, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT3D2(*nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

void FINUFFT3D2MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms,
                     BIGINT *mt, BIGINT *mu, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT3D2MANY(*ntransf, *nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

void FINUFFT3D3_(BIGINT *nj, FLT *x, FLT *y, FLT *z, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s, FLT *t, FLT *u,
                 CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT3D3(*nj, x, y, z, c, *iflag, *eps, *nk, s, t, u, f, o);
}

void FINUFFT3D3MANY_(int *ntransf, BIGINT *nj, FLT *x, FLT *y, FLT *z, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s,
                     FLT *t, FLT *u, CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT3D3MANY(*ntransf, *nj, x, y, z, c, *iflag, *eps, *nk, s, t, u, f, o);
}

// --- 4D ---
void FINUFFT4D1_(BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt,
                 BIGINT *mu, BIGINT *mv, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT4D1(*nj, xj, yj, zj, pj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, fk, o);
}

void FINUFFT4D1MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int *iflag, FLT *eps,
                     BIGINT *ms, BIGINT *mt, BIGINT *mu, BIGINT *mv, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT4D1MANY(*ntransf, *nj, xj, yj, zj, pj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, fk, o);
}

void FINUFFT4D2_(BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms, BIGINT *mt,
                 BIGINT *mu, BIGINT *mv, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT4D2(*nj, xj, yj, zj, pj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, fk, o);
}

void FINUFFT4D2MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int *iflag, FLT *eps,
                     BIGINT *ms, BIGINT *mt, BIGINT *mu, BIGINT *mv, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT4D2MANY(*ntransf, *nj, xj, yj, zj, pj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, fk, o);
}

void FINUFFT4D3_(BIGINT *nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s, FLT *t,
                 FLT *u, FLT *v, CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT4D3(*nj, x, y, z, p, c, *iflag, *eps, *nk, s, t, u, v, f, o);
}

void FINUFFT4D3MANY_(int *ntransf, BIGINT *nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *c, int *iflag, FLT *eps, BIGINT *nk,
                     FLT *s, FLT *t, FLT *u, FLT *v, CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT4D3MANY(*ntransf, *nj, x, y, z, p, c, *iflag, *eps, *nk, s, t, u, v, f, o);
}

// --- 5D ---
void FINUFFT5D1_(BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms,
                 BIGINT *mt, BIGINT *mu, BIGINT *mv, BIGINT *mw, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT5D1(*nj, xj, yj, zj, pj, qj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, *mw, fk, o);
}

void FINUFFT5D1MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int *iflag,
                     FLT *eps, BIGINT *ms, BIGINT *mt, BIGINT *mu, BIGINT *mv, BIGINT *mw, CPX *fk, nufft_opts *o,
                     int *ier) {
  *ier = FINUFFT5D1MANY(*ntransf, *nj, xj, yj, zj, pj, qj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, *mw, fk, o);
}

void FINUFFT5D2_(BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int *iflag, FLT *eps, BIGINT *ms,
                 BIGINT *mt, BIGINT *mu, BIGINT *mv, BIGINT *mw, CPX *fk, nufft_opts *o, int *ier) {
  *ier = FINUFFT5D2(*nj, xj, yj, zj, pj, qj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, *mw, fk, o);
}

void FINUFFT5D2MANY_(int *ntransf, BIGINT *nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int *iflag,
                     FLT *eps, BIGINT *ms, BIGINT *mt, BIGINT *mu, BIGINT *mv, BIGINT *mw, CPX *fk, nufft_opts *o,
                     int *ier) {
  *ier = FINUFFT5D2MANY(*ntransf, *nj, xj, yj, zj, pj, qj, cj, *iflag, *eps, *ms, *mt, *mu, *mv, *mw, fk, o);
}

void FINUFFT5D3_(BIGINT *nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int *iflag, FLT *eps, BIGINT *nk, FLT *s,
                 FLT *t, FLT *u, FLT *v, FLT *w, CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT5D3(*nj, x, y, z, p, q, c, *iflag, *eps, *nk, s, t, u, v, w, f, o);
}

void FINUFFT5D3MANY_(int *ntransf, BIGINT *nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int *iflag, FLT *eps,
                     BIGINT *nk, FLT *s, FLT *t, FLT *u, FLT *v, FLT *w, CPX *f, nufft_opts *o, int *ier) {
  *ier = FINUFFT5D3MANY(*ntransf, *nj, x, y, z, p, q, c, *iflag, *eps, *nk, s, t, u, v, w, f, o);
}

#ifdef __cplusplus
}
#endif
