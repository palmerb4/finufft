#include <finufft_eitherprec.h>
#include <dataTypes.h>

#include <cstdio>
using namespace std;

/* ---------------------------------------------------------------------------
   The 30 simple interfaces (= 5 dims * 3 types * {singlecall,many}) to FINUFFT.
   As of v1.2 these simply invoke the guru interface, through a helper layer.
   See ../docs/usage.rst or http://finufft.readthedocs.io for documentation
   all routines here.
   This compiles in either double or single precision (based on -DSINGLE),
   producing functions finufft?d?{many} or finufftf?1?{many} respectively.

   Authors: Andrea Malleo and Alex Barnett, 2019-2020.
   ---------------------------------------------------------------------------
*/


// Helper layer ...........................................................

int invokeGuruInterface(int n_dims, int type, int n_transf, BIGINT nj, FLT* xj,
                        FLT *yj, FLT *zj, FLT *pj, FLT *qj, CPX* cj,int iflag, FLT eps,
                        BIGINT *n_modes, BIGINT nk, FLT *s, FLT *t, FLT *u, FLT *v, FLT *w,
                        CPX* fk, nufft_opts *popts)
// Helper layer between simple interfaces (with opts) and the guru functions.
// Author: Andrea Malleo, 2019.
{
  FINUFFT_PLAN plan;
  int ier = FINUFFT_MAKEPLAN(type, n_dims, n_modes, iflag, n_transf, eps,
                             &plan, popts);  // popts (ptr to opts) can be NULL
  fprintf(stderr,"%s FINUFFT_MAKEPLAN\n",__func__);

  if (ier>1) {   // since 1 (a warning) still allows proceeding...
    fprintf(stderr, "FINUFFT invokeGuru: plan error (ier=%d)!\n", ier);
    return ier;
  }

  int ier2 = FINUFFT_SETPTS(plan, nj, xj, yj, zj, pj, qj, nk, s, t, u, v, w);
  fprintf(stderr,"%s FINUFFT_SETPTS\n",__func__);
  if (ier2>1) {
    fprintf(stderr,"FINUFFT invokeGuru: setpts error (ier=%d)!\n", ier2);
    return ier2;
  }

  int ier3 = FINUFFT_EXECUTE(plan, cj, fk);
  fprintf(stderr,"%s FINUFFT_EXECUTE\n",__func__);
  if (ier3>1) {
    fprintf(stderr,"FINUFFT invokeGuru: execute error (ier=%d)!\n", ier3);
    return ier3;
  }

  FINUFFT_DESTROY(plan);
  return max(max(ier,ier2),ier3);  // in case any one gave a (positive!) warning
}



// Dimension 1111111111111111111111111111111111111111111111111111111111111111

int FINUFFT1D1(BIGINT nj,FLT* xj,CPX* cj,int iflag,FLT eps,BIGINT ms,
	       CPX* fk, nufft_opts *opts)
//  Type-1 1D complex nonuniform FFT. See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,1,1,1,1};
  int n_dims = 1;
  int n_transf = 1;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, NULL, NULL, NULL, NULL, cj,
			 iflag, eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT1D1MANY(int n_transf, BIGINT nj,FLT* xj,CPX* cj,int iflag,FLT eps,
                   BIGINT ms, CPX* fk, nufft_opts *opts)
// Type-1 1D complex nonuniform FFT for many vectors. See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,1,1,1,1};
  int n_dims = 1;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, NULL, NULL, NULL, NULL, cj,
		      iflag, eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT1D2(BIGINT nj,FLT* xj,CPX* cj,int iflag,FLT eps,BIGINT ms,
	       CPX* fk, nufft_opts *opts)
//  Type-2 1D complex nonuniform FFT. See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,1,1,1,1};
  int n_dims = 1;
  int n_transf = 1;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, NULL, NULL, NULL, NULL, cj,
			  iflag, eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT1D2MANY(int n_transf, BIGINT nj,FLT* xj,CPX* cj,int iflag,FLT eps,BIGINT ms,
	       CPX* fk, nufft_opts *opts)
//  Type-2 1D complex nonuniform FFT, many vectors. See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,1,1,1,1};
  int n_dims = 1;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, NULL, NULL, NULL, NULL, cj,
		      	iflag, eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT1D3(BIGINT nj,FLT* xj,CPX* cj,int iflag, FLT eps, BIGINT nk, FLT* s, CPX* fk, nufft_opts *opts)
// Type-3 1D complex nonuniform FFT. See ../docs/usage.rst
{
  int n_dims = 1;
  int n_transf = 1;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, NULL, NULL, NULL, NULL, cj,
				iflag, eps, NULL, nk, s, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT1D3MANY(int n_transf, BIGINT nj,FLT* xj,CPX* cj,int iflag, FLT eps, BIGINT nk, FLT* s, CPX* fk, nufft_opts *opts)
  // Type-3 1D complex nonuniform FFT, many vectors. See ../docs/usage.rst
{
  int n_dims = 1;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, NULL, NULL, NULL, NULL, cj,
				iflag, eps, NULL, nk, s, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}


// Dimension 22222222222222222222222222222222222222222222222222222222222222222

int FINUFFT2D1(BIGINT nj,FLT* xj,FLT *yj,CPX* cj,int iflag,
	       FLT eps, BIGINT ms, BIGINT mt, CPX* fk, nufft_opts* opts)
//  Type-1 2D complex nonuniform FFT. See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,1,1,1};
  int n_dims = 2;
  int n_transf = 1;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, NULL, NULL, NULL, cj,
                          iflag, eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT2D1MANY(int n_transf, BIGINT nj, FLT* xj, FLT *yj, CPX* c,
		   int iflag, FLT eps, BIGINT ms, BIGINT mt, CPX* fk,
		   nufft_opts *opts)
//  Type-1 2D complex nonuniform FFT, many vectors. See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,1,1,1};
  int n_dims = 2;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, NULL, NULL, NULL, c,
                        iflag, eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT2D2(BIGINT nj,FLT* xj,FLT *yj,CPX* cj,int iflag,FLT eps,
	       BIGINT ms, BIGINT mt, CPX* fk, nufft_opts *opts)
//  Type-2 2D complex nonuniform FFT.  See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,1,1,1};
  int n_dims = 2;
  int n_transf = 1;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, NULL, NULL, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT2D2MANY(int n_transf, BIGINT nj, FLT* xj, FLT *yj, CPX* c, int iflag,
		   FLT eps, BIGINT ms, BIGINT mt, CPX* fk, nufft_opts *opts)
//  Type-2 2D complex nonuniform FFT, many vectors.  See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,1,1,1};
  int n_dims = 2;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, NULL, NULL, NULL, c, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT2D3(BIGINT nj,FLT* xj,FLT* yj,CPX* cj,int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t, CPX* fk, nufft_opts *opts)
// Type-3 2D complex nonuniform FFT.  See ../docs/usage.rst
{
  int n_dims = 2;
  int type = 3;
  int n_transf = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, NULL, NULL, NULL, cj, iflag, 
        eps, NULL, nk, s, t, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT2D3MANY(int n_transf, BIGINT nj,FLT* xj,FLT* yj,CPX* cj,int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t, CPX* fk, nufft_opts *opts)
// Type-3 2D complex nonuniform FFT, many vectors.  See ../docs/usage.rst
{
  int n_dims = 2;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, NULL, NULL, NULL, cj, iflag, 
        eps, NULL, nk, s, t, NULL, NULL, NULL, fk, opts);
  return ier;
}



// Dimension 3333333333333333333333333333333333333333333333333333333333333333

int FINUFFT3D1(BIGINT nj,FLT* xj,FLT *yj,FLT *zj,CPX* cj,int iflag,
	       FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, CPX* fk,
	       nufft_opts *opts)
//  Type-1 3D complex nonuniform FFT.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,1,1};
  int n_dims = 3;
  int n_transf = 1;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, NULL, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}


int FINUFFT3D1MANY(int n_transf, BIGINT nj,FLT* xj,FLT *yj,FLT *zj,CPX* cj,
                   int iflag, FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, CPX* fk,
                   nufft_opts *opts)
// Type-1 3D complex nonuniform FFT, many vectors.  See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,1,1};
  int n_dims = 3;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, NULL, NULL, cj, iflag,
				eps, n_modes, 0,  NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT3D2(BIGINT nj,FLT* xj,FLT *yj,FLT *zj,CPX* cj,
	       int iflag,FLT eps, BIGINT ms, BIGINT mt, BIGINT mu,
	       CPX* fk, nufft_opts *opts)
// Type-2 3D complex nonuniform FFT.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,1,1};
  int n_dims = 3;
  int n_transf = 1;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, NULL, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT3D2MANY(int n_transf, BIGINT nj,FLT* xj,FLT *yj,FLT *zj,CPX* cj,
	       int iflag,FLT eps, BIGINT ms, BIGINT mt, BIGINT mu,
	       CPX* fk, nufft_opts *opts)
// Type-2 3D complex nonuniform FFT, many vectors.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,1,1};
  int n_dims = 3;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, NULL, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT3D3(BIGINT nj,FLT* xj,FLT* yj,FLT *zj, CPX* cj,
	       int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t,
	       FLT *u, CPX* fk, nufft_opts *opts)
//  Type-3 3D complex nonuniform FFT.   See ../docs/usage.rst
{
  int n_dims = 3;
  int n_transf = 1;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, NULL, NULL, cj, iflag,
				eps, NULL, nk, s ,t ,u, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT3D3MANY(int n_transf, BIGINT nj,FLT* xj,FLT* yj,FLT *zj, CPX* cj,
	       int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t,
	       FLT *u, CPX* fk, nufft_opts *opts)
//  Type-3 3D complex nonuniform FFT, many vectors.   See ../docs/usage.rst
{
  int n_dims = 3;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, NULL, NULL, cj, iflag,
				eps, NULL, nk, s ,t ,u, NULL, NULL, fk, opts);
  return ier;
}


// Dimension 4444444444444444444444444444444444444444444444444444444444444444444

int FINUFFT4D1(BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,CPX* cj,int iflag,
	       FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, CPX* fk,
	       nufft_opts *opts)
//  Type-1 4D complex nonuniform FFT.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,1};
  int n_dims = 4;
  int n_transf = 1;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}


int FINUFFT4D1MANY(int n_transf, BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,CPX* cj,
                   int iflag, FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, CPX* fk,
                   nufft_opts *opts)
// Type-1 4D complex nonuniform FFT, many vectors.  See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,1};
  int n_dims = 4;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, NULL, cj, iflag,
				eps, n_modes, 0,  NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT4D2(BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,CPX* cj,
	       int iflag,FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv,
	       CPX* fk, nufft_opts *opts)
// Type-2 4D complex nonuniform FFT.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,1};
  int n_dims = 4;
  int n_transf = 1;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT4D2MANY(int n_transf, BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,CPX* cj,
	       int iflag,FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv,
	       CPX* fk, nufft_opts *opts)
// Type-2 4D complex nonuniform FFT, many vectors.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,1};
  int n_dims = 4;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, NULL, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT4D3(BIGINT nj,FLT* xj,FLT* yj,FLT *zj,FLT *pj, CPX* cj,
	       int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t,
	       FLT *u, FLT *v, CPX* fk, nufft_opts *opts)
//  Type-3 4D complex nonuniform FFT.   See ../docs/usage.rst
{
  int n_dims = 4;
  int n_transf = 1;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, NULL, cj, iflag,
				eps, NULL, nk, s ,t ,u, v, NULL, fk, opts);
  return ier;
}

int FINUFFT4D3MANY(int n_transf, BIGINT nj,FLT* xj,FLT* yj,FLT *zj,FLT *pj, CPX* cj,
	       int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t,
	       FLT *u, FLT *v, CPX* fk, nufft_opts *opts)
//  Type-3 4D complex nonuniform FFT, many vectors.   See ../docs/usage.rst
{
  int n_dims = 4;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, NULL, cj, iflag,
				eps, NULL, nk, s ,t ,u, v, NULL, fk, opts);
  return ier;
}

// Dimension 55555555555555555555555555555555555555555555555555555555555555

int FINUFFT5D1(BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,FLT *qj,CPX* cj,int iflag,
	       FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw, CPX* fk,
	       nufft_opts *opts)
//  Type-1 5D complex nonuniform FFT.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,mw,1};
  int n_dims = 5;
  int n_transf = 1;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, qj, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}


int FINUFFT5D1MANY(int n_transf, BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,FLT *qj,CPX* cj,
                   int iflag, FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw, CPX* fk,
                   nufft_opts *opts)
// Type-1 5D complex nonuniform FFT, many vectors.  See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,mw,1};
  int n_dims = 5;
  int type = 1;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, qj, cj, iflag,
				eps, n_modes, 0,  NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT5D2(BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,FLT *qj,CPX* cj,
	       int iflag,FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw,
	       CPX* fk, nufft_opts *opts)
// Type-2 5D complex nonuniform FFT.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,mw,1};
  int n_dims = 5;
  int n_transf = 1;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, qj, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT5D2MANY(int n_transf, BIGINT nj,FLT* xj,FLT *yj,FLT *zj,FLT *pj,FLT *qj,CPX* cj,
	       int iflag,FLT eps, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv, BIGINT mw,
	       CPX* fk, nufft_opts *opts)
// Type-2 5D complex nonuniform FFT, many vectors.   See ../docs/usage.rst
{
  BIGINT n_modes[]={ms,mt,mu,mv,mw,1};
  int n_dims = 5;
  int type = 2;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, qj, cj, iflag,
				eps, n_modes, 0, NULL, NULL, NULL, NULL, NULL, fk, opts);
  return ier;
}

int FINUFFT5D3(BIGINT nj,FLT* xj,FLT* yj,FLT *zj,FLT *pj,FLT *qj, CPX* cj,
	       int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t,
	       FLT *u, FLT *v, FLT *w, CPX* fk, nufft_opts *opts)
//  Type-3 5D complex nonuniform FFT.   See ../docs/usage.rst
{
  int n_dims = 5;
  int n_transf = 1;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, qj, cj, iflag,
				eps, NULL, nk, s ,t ,u, v, w, fk, opts);
  return ier;
}

int FINUFFT5D3MANY(int n_transf, BIGINT nj,FLT* xj,FLT* yj,FLT *zj,FLT *pj,FLT *qj, CPX* cj,
	       int iflag, FLT eps, BIGINT nk, FLT* s, FLT *t,
	       FLT *u, FLT *v, FLT *w, CPX* fk, nufft_opts *opts)
//  Type-3 5D complex nonuniform FFT, many vectors.   See ../docs/usage.rst
{
  int n_dims = 5;
  int type = 3;
  int ier = invokeGuruInterface(n_dims, type, n_transf, nj, xj, yj, zj, pj, qj, cj, iflag,
				eps, NULL, nk, s ,t ,u, v, w, fk, opts);
  return ier;
}