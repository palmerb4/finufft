// Switchable-precision interface template for FINUFFT. Used by finufft.h
// Internal use only: users should link to finufft.h
// Barnett 7/1/20

#if (!defined(FINUFFT_H) && !defined(SINGLE)) || (!defined(FINUFFTF_H) && defined(SINGLE))
// (note we entered one level of conditional until the end of this header)
// Make sure we don't include double and single headers more than once each...
#ifndef SINGLE
#define FINUFFT_H
#else
#define FINUFFTF_H
#endif

// Here just what's needed to describe the headers for what finufft provides
#include <dataTypes.h>
#include <finufft_plan_eitherprec.h>
#include <nufft_opts.h>

// clear the macros so we can define w/o warnings...
#undef FINUFFT_DEFAULT_OPTS
#undef FINUFFT_MAKEPLAN
#undef FINUFFT_SETPTS
#undef FINUFFT_EXECUTE
#undef FINUFFT_DESTROY
#undef FINUFFT1D1
#undef FINUFFT1D1MANY
#undef FINUFFT1D2
#undef FINUFFT1D2MANY
#undef FINUFFT1D3
#undef FINUFFT1D3MANY
#undef FINUFFT2D1
#undef FINUFFT2D1MANY
#undef FINUFFT2D2
#undef FINUFFT2D2MANY
#undef FINUFFT2D3
#undef FINUFFT2D3MANY
#undef FINUFFT3D1
#undef FINUFFT3D1MANY
#undef FINUFFT3D2
#undef FINUFFT3D2MANY
#undef FINUFFT3D3
#undef FINUFFT3D3MANY
#undef FINUFFT4D1
#undef FINUFFT4D1MANY
#undef FINUFFT4D2
#undef FINUFFT4D2MANY
#undef FINUFFT4D3
#undef FINUFFT4D3MANY
#undef FINUFFT5D1
#undef FINUFFT5D1MANY
#undef FINUFFT5D2
#undef FINUFFT5D2MANY
#undef FINUFFT5D3
#undef FINUFFT5D3MANY
// precision-switching macros for interfaces FINUFFT provides to outside world
#ifdef SINGLE
#define FINUFFT_DEFAULT_OPTS finufftf_default_opts
#define FINUFFT_MAKEPLAN finufftf_makeplan
#define FINUFFT_SETPTS finufftf_setpts
#define FINUFFT_EXECUTE finufftf_execute
#define FINUFFT_DESTROY finufftf_destroy
#define FINUFFT1D1 finufftf1d1
#define FINUFFT1D1MANY finufftf1d1many
#define FINUFFT1D2 finufftf1d2
#define FINUFFT1D2MANY finufftf1d2many
#define FINUFFT1D3 finufftf1d3
#define FINUFFT1D3MANY finufftf1d3many
#define FINUFFT2D1 finufftf2d1
#define FINUFFT2D1MANY finufftf2d1many
#define FINUFFT2D2 finufftf2d2
#define FINUFFT2D2MANY finufftf2d2many
#define FINUFFT2D3 finufftf2d3
#define FINUFFT2D3MANY finufftf2d3many
#define FINUFFT3D1 finufftf3d1
#define FINUFFT3D1MANY finufftf3d1many
#define FINUFFT3D2 finufftf3d2
#define FINUFFT3D2MANY finufftf3d2many
#define FINUFFT3D3 finufftf3d3
#define FINUFFT3D3MANY finufftf3d3many
#define FINUFFT4D1 finufftf4d1
#define FINUFFT4D1MANY finufftf4d1many
#define FINUFFT4D2 finufftf4d2
#define FINUFFT4D2MANY finufftf4d2many
#define FINUFFT4D3 finufftf4d3
#define FINUFFT4D3MANY finufftf4d3many
#define FINUFFT5D1 finufftf5d1
#define FINUFFT5D1MANY finufftf5d1many
#define FINUFFT5D2 finufftf5d2
#define FINUFFT5D2MANY finufftf5d2many
#define FINUFFT5D3 finufftf5d3
#define FINUFFT5D3MANY finufftf5d3many
#else
#define FINUFFT_DEFAULT_OPTS finufft_default_opts
#define FINUFFT_MAKEPLAN finufft_makeplan
#define FINUFFT_SETPTS finufft_setpts
#define FINUFFT_EXECUTE finufft_execute
#define FINUFFT_DESTROY finufft_destroy
#define FINUFFT1D1 finufft1d1
#define FINUFFT1D1MANY finufft1d1many
#define FINUFFT1D2 finufft1d2
#define FINUFFT1D2MANY finufft1d2many
#define FINUFFT1D3 finufft1d3
#define FINUFFT1D3MANY finufft1d3many
#define FINUFFT2D1 finufft2d1
#define FINUFFT2D1MANY finufft2d1many
#define FINUFFT2D2 finufft2d2
#define FINUFFT2D2MANY finufft2d2many
#define FINUFFT2D3 finufft2d3
#define FINUFFT2D3MANY finufft2d3many
#define FINUFFT3D1 finufft3d1
#define FINUFFT3D1MANY finufft3d1many
#define FINUFFT3D2 finufft3d2
#define FINUFFT3D2MANY finufft3d2many
#define FINUFFT3D3 finufft3d3
#define FINUFFT3D3MANY finufft3d3many
#define FINUFFT4D1 finufft4d1
#define FINUFFT4D1MANY finufft4d1many
#define FINUFFT4D2 finufft4d2
#define FINUFFT4D2MANY finufft4d2many
#define FINUFFT4D3 finufft4d3
#define FINUFFT4D3MANY finufft4d3many
#define FINUFFT5D1 finufft5d1
#define FINUFFT5D1MANY finufft5d1many
#define FINUFFT5D2 finufft5d2
#define FINUFFT5D2MANY finufft5d2many
#define FINUFFT5D3 finufft5d3
#define FINUFFT5D3MANY finufft5d3many
#endif

// all interfaces are C-style even when used from C++...
#ifdef __cplusplus
extern "C" {
#endif

// ------------------ the guru interface ------------------------------------
// (sources in finufft.cpp)

void FINUFFT_DEFAULT_OPTS(nufft_opts *o);
int FINUFFT_MAKEPLAN(int type, int dim, BIGINT *n_modes, int iflag, int n_transf, FLT tol, FINUFFT_PLAN *plan,
                     nufft_opts *o);
int FINUFFT_SETPTS(FINUFFT_PLAN plan, BIGINT M, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, BIGINT N, FLT *s, FLT *t,
                   FLT *u, FLT *v, FLT *w);
int FINUFFT_EXECUTE(FINUFFT_PLAN plan, CPX *weights, CPX *result);
int FINUFFT_DESTROY(FINUFFT_PLAN plan);

// ----------------- the 30 simple interfaces -------------------------------
// (sources in simpleinterfaces.cpp)

int FINUFFT1D1(BIGINT nj, FLT *xj, CPX *cj, int iflag, FLT eps, BIGINT ms, CPX *fk, nufft_opts *opts);
int FINUFFT1D1MANY(int ntransf, BIGINT nj, FLT *xj, CPX *cj, int iflag, FLT eps, BIGINT ms, CPX *fk, nufft_opts *opts);

int FINUFFT1D2(BIGINT nj, FLT *xj, CPX *cj, int iflag, FLT eps, BIGINT ms, CPX *fk, nufft_opts *opts);
int FINUFFT1D2MANY(int ntransf, BIGINT nj, FLT *xj, CPX *cj, int iflag, FLT eps, BIGINT ms, CPX *fk, nufft_opts *opts);
int FINUFFT1D3(BIGINT nj, FLT *x, CPX *c, int iflag, FLT eps, BIGINT nk, FLT *s, CPX *f, nufft_opts *opts);
int FINUFFT1D3MANY(int ntransf, BIGINT nj, FLT *x, CPX *c, int iflag, FLT eps, BIGINT nk, FLT *s, CPX *f,
                   nufft_opts *opts);

int FINUFFT2D1(BIGINT nj, FLT *xj, FLT *yj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt, CPX *fk,
               nufft_opts *opts);
int FINUFFT2D1MANY(int ndata, BIGINT nj, FLT *xj, FLT *yj, CPX *c, int iflag, FLT eps, BIGINT ms, BIGINT mt, CPX *fk,
                   nufft_opts *opts);
int FINUFFT2D2(BIGINT nj, FLT *xj, FLT *yj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt, CPX *fk,
               nufft_opts *opts);
int FINUFFT2D2MANY(int ndata, BIGINT nj, FLT *xj, FLT *yj, CPX *c, int iflag, FLT eps, BIGINT ms, BIGINT mt, CPX *fk,
                   nufft_opts *opts);
int FINUFFT2D3(BIGINT nj, FLT *x, FLT *y, CPX *cj, int iflag, FLT eps, BIGINT nk, FLT *s, FLT *t, CPX *fk,
               nufft_opts *opts);
int FINUFFT2D3MANY(int ntransf, BIGINT nj, FLT *x, FLT *y, CPX *cj, int iflag, FLT eps, BIGINT nk, FLT *s, FLT *t,
                   CPX *fk, nufft_opts *opts);

int FINUFFT3D1(BIGINT nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt, BIGINT mu,
               CPX *fk, nufft_opts *opts);
int FINUFFT3D1MANY(int ntransfs, BIGINT nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int iflag, FLT eps, BIGINT ms,
                   BIGINT mt, BIGINT mu, CPX *fk, nufft_opts *opts);
int FINUFFT3D2(BIGINT nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt, BIGINT mu,
               CPX *fk, nufft_opts *opts);
int FINUFFT3D2MANY(int ntransf, BIGINT nj, FLT *xj, FLT *yj, FLT *zj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt,
                   BIGINT mu, CPX *fk, nufft_opts *opts);
int FINUFFT3D3(BIGINT nj, FLT *x, FLT *y, FLT *z, CPX *cj, int iflag, FLT eps, BIGINT nk, FLT *s, FLT *t, FLT *u,
               CPX *fk, nufft_opts *opts);
int FINUFFT3D3MANY(int ntransf, BIGINT nj, FLT *x, FLT *y, FLT *z, CPX *cj, int iflag, FLT eps, BIGINT nk, FLT *s,
                   FLT *t, FLT *u, CPX *fk, nufft_opts *opts);

int FINUFFT4D1(BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt,
               BIGINT mu, BIGINT mv, CPX *fk, nufft_opts *opts);
int FINUFFT4D1MANY(int ntransfs, BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int iflag, FLT eps, BIGINT ms,
                   BIGINT mt, BIGINT mu, BIGINT mv, CPX *fk, nufft_opts *opts);
int FINUFFT4D2(BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int iflag, FLT eps, BIGINT ms, BIGINT mt,
               BIGINT mu, BIGINT mv, CPX *fk, nufft_opts *opts);
int FINUFFT4D2MANY(int ntransf, BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, CPX *cj, int iflag, FLT eps, BIGINT ms,
                   BIGINT mt, BIGINT mu, BIGINT mv, CPX *fk, nufft_opts *opts);
int FINUFFT4D3(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *cj, int iflag, FLT eps, BIGINT nk, FLT *s, FLT *t,
               FLT *u, FLT *v, CPX *fk, nufft_opts *opts);
int FINUFFT4D3MANY(int ntransf, BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *cj, int iflag, FLT eps, BIGINT nk,
                   FLT *s, FLT *t, FLT *u, FLT *v, CPX *fk, nufft_opts *opts);

int FINUFFT5D1(BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int iflag, FLT eps, BIGINT ms,
               BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw, CPX *fk, nufft_opts *opts);
int FINUFFT5D1MANY(int ntransfs, BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int iflag, FLT eps,
                   BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw, CPX *fk, nufft_opts *opts);
int FINUFFT5D2(BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int iflag, FLT eps, BIGINT ms,
               BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw, CPX *fk, nufft_opts *opts);
int FINUFFT5D2MANY(int ntransf, BIGINT nj, FLT *xj, FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX *cj, int iflag, FLT eps,
                   BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw, CPX *fk, nufft_opts *opts);
int FINUFFT5D3(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *cj, int iflag, FLT eps, BIGINT nk, FLT *s,
               FLT *t, FLT *u, FLT *v, FLT *w, CPX *fk, nufft_opts *opts);
int FINUFFT5D3MANY(int ntransf, BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *cj, int iflag, FLT eps,
                   BIGINT nk, FLT *s, FLT *t, FLT *u, FLT *v, FLT *w, CPX *fk, nufft_opts *opts);
#ifdef __cplusplus
}
#endif

#endif // FINUFFT_H or FINUFFTF_H
