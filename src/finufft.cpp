#include <finufft.h>
#include <utils.h>
#include <iostream>
#include <common.h>
#include <iomanip>

/* The main guru functions for FINUFFT.

   Original guru interface written by Andrea Malleo, summer 2019, mentored
   by Alex Barnett. Many rewrites in early 2020 by Alex Barnett, Libin Lu.

   As of v1.2 these replace the old hand-coded separate 9 finufft?d?() functions
   and the two finufft2d?many() functions.
   The (now 18) simple C++ interfaces are in simpleinterfaces.cpp

Notes on algorithms taken from old finufft?d?() documentation, Feb-Jun 2017:

   TYPE 1:
     The type 1 NUFFT proceeds in three main steps:
     1) spread data to oversampled regular mesh using kernel.
     2) compute FFT on uniform mesh
     3) deconvolve by division of each Fourier mode independently by the kernel
        Fourier series coeffs (not merely FFT of kernel), shuffle to output.
     The kernel coeffs are precomputed in what is called step 0 in the code.
   Written with FFTW style complex arrays. Step 3a internally uses CPX,
   and Step 3b internally uses real arithmetic and FFTW style complex.

   TYPE 2:
     The type 2 algorithm proceeds in three main steps:
     1) deconvolve (amplify) each Fourier mode, dividing by kernel Fourier coeff
     2) compute inverse FFT on uniform fine grid
     3) spread (dir=2, ie interpolate) data to regular mesh
     The kernel coeffs are precomputed in what is called step 0 in the code.
   Written with FFTW style complex arrays. Step 0 internally uses CPX,
   and Step 1 internally uses real arithmetic and FFTW style complex.

   TYPE 3:
     The type 3 algorithm is basically a type 2 (which is implemented precisely
     as call to type 2) replacing the middle FFT (Step 2) of a type 1.
     Beyond this, the new twists are:
     i) nf1, number of upsampled points for the type-1, depends on the product
       of interval widths containing input and output points (X*S).
     ii) The deconvolve (post-amplify) step is division by the Fourier transform
       of the scaled kernel, evaluated on the *nonuniform* output frequency
       grid; this is done by direct approximation of the Fourier integral
       using quadrature of the kernel function times exponentials.
     iii) Shifts in x (real) and s (Fourier) are done to minimize the interval
       half-widths X and S, hence nf1.
   No references to FFTW are needed here. CPX arithmetic is used.

   MULTIPLE STRENGTH VECTORS FOR THE SAME NONUNIFORM POINTS (n_transf>1):
     maxBatchSize (set to max_num_omp_threads) times the RAM is needed, so
     this is good only for small problems.


Design notes for guru interface implementation:

* Since finufft_plan is C-compatible, we need to use malloc/free for its
  allocatable arrays, keeping it quite low-level. We can't use std::vector
  since that would  only survive in the scope of each function.

*/


int* gridsize_for_fftw(finufft_plan* p){
// helper func returns a new int array of length dim, extracted from
// the finufft plan, that fft_plan_many_dft needs as its 2nd argument.
  int* nf;
  if(p->dim == 1){ 
    nf = new int[1];
    nf[0] = (int)p->nf1;
  }
  else if (p->dim == 2){ 
    nf = new int[2];
    nf[0] = (int)p->nf2;
    nf[1] = (int)p->nf1; 
  }   // fftw enforced row major ordering, ie dims are backwards ordered
  else{ 
    nf = new int[3];
    nf[0] = (int)p->nf3;
    nf[1] = (int)p->nf2;
    nf[2] = (int)p->nf1;
  }
  return nf;
}


// PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
int finufft_makeplan(int type, int dim, BIGINT* n_modes, int iflag,
                     int ntrans, FLT tol, finufft_plan* p, nufft_opts* opts)
// Populates the fields of finufft_plan which is pointed to by "p".
// opts is ptr to a nufft_opts to set options, or NULL to use defaults.
// For types 1,2 allocates memory for internal working arrays,
// evaluates spreading kernel coefficients, and instantiates the fftw_plan
{  
  cout << scientific << setprecision(15);  // for debug outputs

  if((type!=1)&&(type!=2)&&(type!=3)) {
    fprintf(stderr, "Invalid type (%d), should be 1, 2 or 3.",type);
    return ERR_TYPE_NOTVALID;
  }
  if((dim!=1)&&(dim!=2)&&(dim!=3)) {
    fprintf(stderr, "Invalid dim (%d), should be 1, 2 or 3.",dim);
    return ERR_DIM_NOTVALID;
  }
  if (ntrans<1) {
    fprintf(stderr,"ntrans (%d) should be at least 1.\n",ntrans);
    return ERR_NTRANS_NOTVALID;
  }

  if (opts==NULL)                        // use default opts
    finufft_default_opts(&(p->opts));
  else                                   // or read from what's passed in
    p->opts = *opts;    // does deep copy; changing *opts now has no effect
  // write into plan's spread options...
  int ier = setup_spreader_for_nufft(p->spopts, tol, p->opts);
  if (ier)
    return ier;

  // get stuff from args...
  p->type = type;
  p->dim = dim;
  p->ntrans = ntrans;
  p->tol = tol;
  p->fftSign = (iflag>=0) ? 1 : -1;         // clean up flag input

  // choose batchSize for types 1,2 or 3... (uses int ceil(b/a)=1+(b-1)/a trick)
  int nth = min(MY_OMP_GET_MAX_THREADS(), MAX_USEFUL_NTHREADS);   // limit it
  if (p->opts.maxbatchsize==0) {            // logic to auto-set best batchsize
    p->nbatch = 1+(ntrans-1)/nth;           // min # batches poss
    p->batchSize = 1+(ntrans-1)/p->nbatch;  // then cut # thr in each b
  } else {                                  // batchSize override by user
    p->batchSize = min(p->opts.maxbatchsize,ntrans);
    p->nbatch = 1+(ntrans-1)/p->batchSize;  // resulting # batches
  }
  if (p->opts.spread_thread==0)
    p->opts.spread_thread=2;                // the auto choice
  
  // set others as defaults (or unallocated for arrays)...
  p->X = NULL; p->Y = NULL; p->Z = NULL;
  p->phiHat1 = NULL; p->phiHat2 = NULL; p->phiHat3 = NULL; 
  p->nf1 = 1; p->nf2 = 1; p->nf3 = 1;  // crucial to leave as 1 for unused dims
  p->ms = 1; p->mt = 1; p->mu = 1;     // crucial to leave as 1 for unused dims

  //  ------------------------ types 1,2: planning needed ---------------------
  if((type == 1) || (type == 2)) {

    int nth_fft = MY_OMP_GET_MAX_THREADS();   // give FFTW all it has access to
    // *** should limt max # threads here too? or set equal to batchsize?
    // *** put in logic for setting FFTW # thr based on o.spread_thread?
    FFTW_INIT();           // only does anything when OMP=ON for >1 threads
    FFTW_PLAN_TH(nth_fft); // "  (not batchSize since can be 1 but want mul-thr)
    p->spopts.spread_direction = type;
    
    // read user mode array dims then determine fine grid sizes, sanity check...
    p->ms = n_modes[0];
    ier = set_nf_type12(p->ms,p->opts,p->spopts,&(p->nf1));
    if (ier) return ier;    // nf too big; we're done
    p->phiHat1 = (FLT*)malloc(sizeof(FLT)*(p->nf1/2 + 1));
    if (dim > 1) {
      p->mt = n_modes[1];
      ier = set_nf_type12(p->mt, p->opts, p->spopts, &(p->nf2));
      if (ier) return ier;
      p->phiHat2 = (FLT*)malloc(sizeof(FLT)*(p->nf2/2 + 1));
    }
    if (dim > 2) {
      p->mu = n_modes[2];
      ier = set_nf_type12(p->mu, p->opts, p->spopts, &(p->nf3)); 
      if (ier) return ier;
      p->phiHat3 = (FLT*)malloc(sizeof(FLT)*(p->nf3/2 + 1));
    }

    if (p->opts.debug) { // "long long" here is to avoid warnings with printf...
      printf("[finufft_plan] %dd%d: (ms,mt,mu)=(%lld,%lld,%lld) (nf1,nf2,nf3)=(%lld,%lld,%lld)\n               ntrans=%d nth=%d batchSize=%d",
             dim, type, (long long)p->ms,(long long)p->mt,
             (long long) p->mu, (long long)p->nf1,(long long)p->nf2,
             (long long)p->nf3, ntrans, nth, p->batchSize);
      if (p->batchSize==1)          // spread_thread has no effect in this case
        printf("\n");
      else
        printf(" spread_thread=%d\n", p->opts.spread_thread);
    }

    // STEP 0: get Fourier coeffs of spreading kernel along each fine grid dim
    CNTime timer; timer.start();
    onedim_fseries_kernel(p->nf1, p->phiHat1, p->spopts);
    if (dim>1) onedim_fseries_kernel(p->nf2, p->phiHat2, p->spopts);
    if (dim>2) onedim_fseries_kernel(p->nf3, p->phiHat3, p->spopts);
    if (p->opts.debug) printf("[finufft_plan] kernel fser (ns=%d):\t\t%.3g s\n", p->spopts.nspread, timer.elapsedsec());

    timer.restart();
    p->N = p->ms*p->mt*p->mu;          // N = total # modes
    p->nf = p->nf1*p->nf2*p->nf3;      // fine grid total number of points
    p->fwBatch = FFTW_ALLOC_CPX(p->nf * p->batchSize);    // the big workspace
    if (p->opts.debug) printf("[finufft_plan] fwBatch %.2fGB alloc:\t\t%.3g s\n", (double)1E-09*sizeof(CPX)*p->nf*p->batchSize, timer.elapsedsec());
    if(!p->fwBatch) {      // we don't catch all such mallocs, just this big one
      fprintf(stderr, "FFTW malloc failed for fwBatch (working fine grids)!\n");
      free(p->phiHat1); free(p->phiHat2); free(p->phiHat3);
      return ERR_ALLOC; 
    }
   
    timer.restart();            // plan the FFTW
    int *ns = gridsize_for_fftw(p);
    // fftw_plan_many_dft args: rank, gridsize/dim, howmany, in, inembed, istride, idist, ot, onembed, ostride, odist, sign, flags 
    p->fftwPlan = FFTW_PLAN_MANY_DFT(dim, ns, p->batchSize, p->fwBatch,
         NULL, 1, p->nf, p->fwBatch, NULL, 1, p->nf, p->fftSign, p->opts.fftw);
    if (p->opts.debug) printf("[finufft_plan] FFTW plan (mode %d, nth=%d):\t%.3g s\n", p->opts.fftw, nth_fft, timer.elapsedsec());
    delete []ns;
    
  } else {  // -------------------------- type 3 (no planning) ------------

    if (p->opts.debug) printf("[finufft_plan] %dd%d: ntrans=%d\n",dim,type,ntrans);
    // in case destroy occurs before setpts, need safe dummy ptrs/plans...
    p->CpBatch = NULL;
    p->Sp = NULL; p->Tp = NULL; p->Up = NULL;
    p->prephase = NULL;
    p->deconv = NULL;
    p->innerT2plan = NULL;
    // Type 3 will call finufft_makeplan for type 2; no need to init FFTW
    // Note we don't even know nj or nk yet, so can't do anything else!
  }
  return 0;
}


// SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
int finufft_setpts(finufft_plan* p, BIGINT nj, FLT* xj, FLT* yj, FLT* zj,
                   BIGINT nk, FLT* s, FLT* t, FLT* u)
/* For type 1,2: just checks and (possibly) sorts the NU points, in prep for
   spreading.
   For type 3: allocates internal working arrays, scales/centers the NU points
   and NU target freqs, evaluates spreading kernel FT at all target freqs.
*/
{
  int d = p->dim;     // abbrev for spatial dim
  CNTime timer; timer.start();
  p->nj = nj;    // the user only now chooses how many NU (x,y,z) pts

  if (p->type!=3) {  // ------------------ TYPE 1,2 SETPTS -------------------
                     // (all we can do is check and maybe bin-sort the NU pts)
    p->X = xj;       // plan must keep pointers to user's fixed NU pts
    p->Y = yj;
    p->Z = zj;
    int ier = spreadcheck(p->nf1, p->nf2, p->nf3, p->nj, xj, yj, zj, p->spopts);
    if (p->opts.debug>1) printf("[finufft_setpts] spreadcheck (%d):\t%.3g s\n", p->spopts.chkbnds, timer.elapsedsec());
    if (ier)
      return ier;    
    timer.restart();
    p->sortIndices = (BIGINT *)malloc(sizeof(BIGINT)*p->nj);
    if (!p->sortIndices) {
      fprintf(stderr,"[finufft_setpts] failed to allocate sortIndices!\n");
      return ERR_SPREAD_ALLOC;
    }
    p->didSort = indexSort(p->sortIndices, p->nf1, p->nf2, p->nf3, p->nj, xj, yj, zj, p->spopts);
    if (p->opts.debug) printf("[finufft_setpts] sort (didSort=%d):\t\t%.3g s\n", p->didSort, timer.elapsedsec());

    
  } else {   // ------------------------- TYPE 3 SETPTS -----------------------
             // (here we can precompute pre/post-phase factors and plan the t2)

    p->nk = nk;     // user set # targ freq pts
    p->S = s;       // keep pointers to user's input target pts
    p->T = t;
    p->U = u;

    // pick x, s intervals & shifts & # fine grid pts (nf) in each dim...
    FLT S1,S2,S3;       // get half-width X, center C, which contains {x_j}...
    arraywidcen(nj,xj,&(p->t3P.X1),&(p->t3P.C1));
    arraywidcen(nk,s,&S1,&(p->t3P.D1));      // same D, S, but for {s_k}
    set_nhg_type3(S1,p->t3P.X1,p->opts,p->spopts,
           &(p->nf1),&(p->t3P.h1),&(p->t3P.gam1));  // applies twist i)
    p->t3P.C2 = 0.0;        // their defaults if dim 2 unused, etc
    p->t3P.D2 = 0.0;
    if (d>1) {
      arraywidcen(nj,yj,&(p->t3P.X2),&(p->t3P.C2));     // {y_j}
      arraywidcen(nk,t,&S2,&(p->t3P.D2));               // {t_k}
      set_nhg_type3(S2,p->t3P.X2,p->opts,p->spopts,&(p->nf2),
                    &(p->t3P.h2),&(p->t3P.gam2));
    }    
    p->t3P.C3 = 0.0;
    p->t3P.D3 = 0.0;
    if (d>2) {
      arraywidcen(nj,zj,&(p->t3P.X3),&(p->t3P.C3));     // {z_j}
      arraywidcen(nk,u,&S3,&(p->t3P.D3));               // {u_k}
      set_nhg_type3(S3,p->t3P.X3,p->opts,p->spopts,
                    &(p->nf3),&(p->t3P.h3),&(p->t3P.gam3));
    }

    if (p->opts.debug) {  // report on choices of shifts, centers, etc...
      printf("\tM=%lld N=%lld\n",(long long)nj,(long long)nk);
      printf("\tX1=%.3g C1=%.3g S1=%.3g D1=%.3g gam1=%g nf1=%lld\t\n", p->t3P.X1, p->t3P.C1,S1, p->t3P.D1, p->t3P.gam1,(long long) p->nf1);
      if (d>1)
        printf("\tX2=%.3g C2=%.3g S2=%.3g D2=%.3g gam2=%g nf2=%lld\n",p->t3P.X2, p->t3P.C2,S2, p->t3P.D2, p->t3P.gam2,(long long) p->nf2);
      if (d>2)
        printf("\tX3=%.3g C3=%.3g S3=%.3g D3=%.3g gam3=%g nf3=%lld\n", p->t3P.X3, p->t3P.C3,S3, p->t3P.D3, p->t3P.gam3,(long long) p->nf3);
    }
    p->nf = p->nf1*p->nf2*p->nf3;      // fine grid total number of points
    p->fwBatch = FFTW_ALLOC_CPX(p->nf * p->batchSize);    // maybe big workspace
    // (note FFTW_ALLOC is not needed but matches its type)
    p->CpBatch = (CPX*)malloc(sizeof(CPX) * nj*p->batchSize);  // batch c' work
    if (p->opts.debug) printf("[finufft_setpts t3] widcen, batch %.2fGB alloc:\t%.3g s\n", (double)1E-09*sizeof(CPX)*(p->nf+nj)*p->batchSize, timer.elapsedsec());
    if(!p->fwBatch || !p->CpBatch) {
      fprintf(stderr, "finufft_setpts t3 malloc fail for fwBatch or CpBatch!\n");
      return ERR_ALLOC; 
    }

    // alloc rescaled NU src pts x'_j (in X etc), rescaled NU targ pts s'_k ...
    p->X = (FLT*)malloc(sizeof(FLT)*nj);
    p->Sp = (FLT*)malloc(sizeof(FLT)*nk);
    if (d>1) {
      p->Y = (FLT*)malloc(sizeof(FLT)*nj);
      p->Tp = (FLT*)malloc(sizeof(FLT)*nk);
    }
    if (d>2) {
      p->Z = (FLT*)malloc(sizeof(FLT)*nj);
      p->Up = (FLT*)malloc(sizeof(FLT)*nk);
    }

    // always shift as use gam to rescale x_j to x'_j, etc (twist iii)...
    FLT ig1 = 1.0/p->t3P.gam1, ig2=0.0, ig3=0.0;   // "reciprocal-math" optim
    if (d>1)
      ig2 = 1.0/p->t3P.gam2;
    if (d>2)
      ig3 = 1.0/p->t3P.gam3;
#pragma omp parallel for schedule(static)
    for (BIGINT j=0;j<nj;++j) {
      p->X[j] = (xj[j] - p->t3P.C1) * ig1;         // rescale x_j
      if (d>1)        // (ok to do inside loop because of branch predict)
        p->Y[j] = (yj[j]- p->t3P.C2) * ig2;        // rescale y_j
      if (d>2)
        p->Z[j] = (zj[j] - p->t3P.C3) * ig3;       // rescale z_j
    }

    // set up prephase array...
    CPX imasign = (p->fftSign>=0) ? IMA : -IMA;             // +-i
    p->prephase = (CPX*)malloc(sizeof(CPX)*nj);
    if (p->t3P.D1!=0.0 || p->t3P.D2!=0.0 || p->t3P.D3!=0.0) {
#pragma omp parallel for schedule(static)
      for (BIGINT j=0;j<nj;++j) {          // ... loop over src NU locs
        FLT phase = p->t3P.D1*xj[j];
        if (d>1)
          phase += p->t3P.D2*yj[j];
        if (d>2)
          phase += p->t3P.D3*zj[j];
        p->prephase[j] = cos(phase)+imasign*sin(phase);   // Euler e^{+-i.phase}
      }
    } else
      for (BIGINT j=0;j<nj;++j)
        p->prephase[j] = (CPX)1.0;     // *** or keep flag so no mult in exec??
      
    // rescale the target s_k etc to s'_k etc...
#pragma omp parallel for schedule(static)
    for (BIGINT k=0;k<nk;++k) {
      p->Sp[k] = p->t3P.h1*p->t3P.gam1*(s[k]- p->t3P.D1);  // so |s'_k| < pi/R
      if (d>1)
        p->Tp[k] = p->t3P.h2*p->t3P.gam2*(t[k]- p->t3P.D2);  // so |t'_k| < pi/R
      if (d>2)
        p->Up[k] = p->t3P.h3*p->t3P.gam3*(u[k]- p->t3P.D3);  // so |u'_k| < pi/R
    }
    
    // (old STEP 3a) Compute deconvolution post-factors array (per targ pt)...
    // (exploits that FT separates because kernel is prod of 1D funcs)
    p->deconv = (CPX*)malloc(sizeof(CPX)*nk);
    FLT *phiHatk1 = (FLT*)malloc(sizeof(FLT)*nk);  // don't confuse w/ p->phiHat
    onedim_nuft_kernel(nk, p->Sp, phiHatk1, p->spopts);         // fill phiHat1
    FLT *phiHatk2 = NULL, *phiHatk3 = NULL;
    if (d>1) {
      phiHatk2 = (FLT*)malloc(sizeof(FLT)*nk);
      onedim_nuft_kernel(nk, p->Tp, phiHatk2, p->spopts);       // fill phiHat2
    }
    if (d>2) {
      phiHatk3 = (FLT*)malloc(sizeof(FLT)*nk);
      onedim_nuft_kernel(nk, p->Up, phiHatk3, p->spopts);       // fill phiHat3
    }
    // *** check if C1 etc can be nonfinite ?
    int Cnonzero = (p->t3P.C1!=0.0 || p->t3P.C2!=0.0 || p->t3P.C3!=0.0);  // cen
#pragma omp parallel for schedule(static)
    for (BIGINT k=0;k<nk;++k) {         // .... loop over NU targ freqs
      FLT phiHat = phiHatk1[k];
      if (d>1)
        phiHat *= phiHatk2[k];
      if (d>2)
        phiHat *= phiHatk3[k];
      p->deconv[k] = (CPX)(1.0 / phiHat);
      if (Cnonzero) {
        FLT phase = (s[k] - p->t3P.D1) * p->t3P.C1;
        if (d>1)
          phase += (t[k] - p->t3P.D2) * p->t3P.C2;
        if (d>2)
          phase += (u[k] - p->t3P.D3) * p->t3P.C3;
        p->deconv[k] *= cos(phase)+imasign*sin(phase);   // Euler e^{+-i.phase}
      }
    }
    free(phiHatk1); free(phiHatk2); free(phiHatk3);  // done w/ deconv fill
    if (p->opts.debug) printf("[finufft_setpts t3] phase & deconv factors:\t%.3g s\n", timer.elapsedsec());

    // Set up sort for spreading Cp (from primed NU src pts X, Y, Z) to fw...
    timer.restart();
    p->sortIndices = (BIGINT *)malloc(sizeof(BIGINT)*p->nj);
    if (!p->sortIndices) {
      fprintf(stderr,"finufft_setpts t3 failed to allocate sortIndices!\n");
      return ERR_SPREAD_ALLOC;
    }
    p->didSort = indexSort(p->sortIndices, p->nf1, p->nf2, p->nf3, p->nj, p->X, p->Y, p->Z, p->spopts);
    if (p->opts.debug) printf("[finufft_setpts t3] sort (didSort=%d):\t\t%.3g s\n", p->didSort, timer.elapsedsec());
 
    // Plan and setpts once, for the (repeated) inner type 2 finufft call...
    timer.restart();
    p->innerT2plan = new finufft_plan;            // ptr to a plan
    BIGINT t2nmodes[] = {p->nf1,p->nf2,p->nf3};   // t2 input is actually fw
    nufft_opts t2opts = p->opts;                  // deep copy, since not ptrs
    t2opts.debug = max(0,p->opts.debug-1);        // don't print as much detail
    t2opts.spread_debug = max(0,p->opts.spread_debug-1);
    // (...could vary other t2opts here?)
    int ier = finufft_makeplan(2, d, t2nmodes, p->fftSign, p->batchSize, p->tol,
                               p->innerT2plan, &t2opts);
    if (ier) {
      fprintf(stderr,"finufft_setpts t3: inner type 2 plan creation failed!\n");
      return ier;
    }
    ier = finufft_setpts(p->innerT2plan, nk, p->Sp, p->Tp, p->Up, 0, NULL, NULL, NULL);  // note nk = # output points (not nj)
    if (ier) {
      fprintf(stderr,"finufft_setpts t3: inner type 2 setpts failed!\n");
      return ier;
    }
    if (p->opts.debug) printf("[finufft_setpts t3] inner t2 plan & setpts: \t%.3g s\n", timer.elapsedsec());

  }
  return 0;
}
// ............ end setpts ..................................................





// BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

int spreadinterpSortedBatch(int batchSize, finufft_plan* p, CPX* cBatch)
/*
  Spreads (or interpolates) a batch of batchSize strength vectors in cBatch
  to (or from) the batch of fine working grids p->fwBatch, using the same set of
  (index-sorted) NU points p->X,Y,Z for each vector in the batch.
  The direction (spread vs interpolate) is set by p->spopts.spread_direction.
  Returns 0 (no error reporting for now).
  Notes:
  1) cBatch is already assumed to have the correct offset, ie here we
     read from the start of cBatch (unlike Malleo). fwBatch also has zero offset
  2) this routine is a batched version of spreadinterpSorted in spreadinterp.cpp
  Barnett 5/19/20, based on Malleo 2019.
*/
{
  // OMP nesting. 0: any omp-parallelism inside the loop sees only 1 thread;
  // note this doesn't change omp_get_max_nthreads()
  // 1: omp par inside the loop sees all threads.  *** ?
  MY_OMP_SET_NESTED(p->opts.spread_thread!=2);
  int nthr_outer = p->opts.spread_thread==1 ? 1 : batchSize;
  
#pragma omp parallel for num_threads(nthr_outer)
  for (int i=0; i<batchSize; i++) {
    FFTW_CPX *fwi = p->fwBatch + i*p->nf;  // start of i'th fw array in wkspace
    CPX *ci = cBatch + i*p->nj;            // start of i'th c array in cBatch
    spreadinterpSorted(p->sortIndices, p->nf1, p->nf2, p->nf3, (FLT*)fwi, p->nj,
                       p->X, p->Y, p->Z, (FLT*)ci, p->spopts, p->didSort);
  }

  MY_OMP_SET_NESTED(0);                    // back to default
  return 0;
}

int deconvolveBatch(int batchSize, finufft_plan* p, CPX* fkBatch)
/*
  Type 1: deconvolves (amplifies) from each interior fw array in p->fwBatch
  into each output array fk in fkBatch.
  Type 2: deconvolves from user-supplied input fk to 0-padded interior fw,
  again looping over fk in fkBatch and fw in p->fwBatch.
  The direction (spread vs interpolate) is set by p->spopts.spread_direction.
  This is mostly a loop calling deconvolveshuffle?d for the needed dim batchSize
  times.
  Barnett 5/21/20, simplified from Malleo 2019 (eg t3 logic won't be in here)
*/
{
  // since deconvolveshuffle?d are single-thread, omp par seems to help here...
#pragma omp parallel for
  for (int i=0; i<batchSize; i++) {
    FFTW_CPX *fwi = p->fwBatch + i*p->nf;  // start of i'th fw array in wkspace
    CPX *fki = fkBatch + i*p->N;           // start of i'th fk array in fkBatch
    
    // Call routine from common.cpp for the dim; prefactors hardcoded to 1.0...
    if (p->dim == 1)
      deconvolveshuffle1d(p->spopts.spread_direction, 1.0, p->phiHat1,
                          p->ms, (FLT *)fki,
                          p->nf1, fwi, p->opts.modeord);
    else if (p->dim == 2)
      deconvolveshuffle2d(p->spopts.spread_direction,1.0, p->phiHat1,
                          p->phiHat2, p->ms, p->mt, (FLT *)fki,
                          p->nf1, p->nf2, fwi, p->opts.modeord);
    else
      deconvolveshuffle3d(p->spopts.spread_direction, 1.0, p->phiHat1,
                          p->phiHat2, p->phiHat3, p->ms, p->mt, p->mu,
                          (FLT *)fki, p->nf1, p->nf2, p->nf3,
                          fwi, p->opts.modeord);
  }
  return 0;
}



// EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
int finufft_exec(finufft_plan* p, CPX* cj, CPX* fk){
/* For given (batch of) weights cj, performs NUFFTs with existing
   (sorted) NU pts and existing plan.
   Performs spread/interp, pre/post deconvolve, and fftw_exec as appropriate
   for each of the 3 types.
   For cases of ntrans>1, performs work in blocks of size up to batchSize.
   Return value 0, no error reporting yet.
   Barnett 5/20/20 based on Malleo 2019.
*/
  CNTime timer; timer.start();
  
  if (p->type!=3){ // --------------------- TYPE 1,2 EXEC ------------------
  
    double t_sprint = 0.0, t_fft = 0.0, t_deconv = 0.0;  // accumulated timing
    if (p->opts.debug)
      printf("[finufft_exec] start ntrans=%d (%d batches, bsize=%d)...\n", p->ntrans, p->nbatch, p->batchSize);
    
    for (int b=0; b*p->batchSize < p->ntrans; b++) { // .....loop b over batches

      // current batch is either batchSize, or possibly truncated if last one
      int thisBatchSize = min(p->ntrans - b*p->batchSize, p->batchSize);
      int bB = b*p->batchSize;         // index of vector, since batchsizes same
      CPX* cjb = cj + bB*p->nj;        // point to batch of weights
      CPX* fkb = fk + bB*p->N;         // point to batch of mode coeffs
      if (p->opts.debug>1) printf("[finufft_exec] start batch %d (size %d):\n",b,thisBatchSize);
      
      // STEP 1: (varies by type)
      timer.restart();
      if (p->type == 1) {  // type 1: spread NU pts p->X, weights cj, to fw grid
        spreadinterpSortedBatch(thisBatchSize, p, cjb);
        t_sprint += timer.elapsedsec();
      } else {          //  type 2: amplify Fourier coeffs fk into 0-padded fw
        deconvolveBatch(thisBatchSize, p, fkb);
        t_deconv += timer.elapsedsec();
      }
             
      // STEP 2: call the pre-planned FFT on this batch
      timer.restart();
      FFTW_EX(p->fftwPlan);   // if thisBatchSize<batchSize it wastes some flops
      t_fft += timer.elapsedsec();
      if (p->opts.debug>1)
        printf("\tFFTW exec:\t\t%.3g s\n", timer.elapsedsec());
      
      // STEP 3: (varies by type)
      timer.restart();        
      if (p->type == 1) {   // type 1: deconvolve (amplify) fw and shuffle to fk
        deconvolveBatch(thisBatchSize, p, fkb);
        t_deconv += timer.elapsedsec();
      } else {          // type 2: interpolate unif fw grid to NU target pts
        spreadinterpSortedBatch(thisBatchSize, p, cjb);
        t_sprint += timer.elapsedsec(); 
      }
    }                                                   // ........end b loop
    
    if (p->opts.debug) {  // report total times in their natural order...
      if(p->type == 1) {
        printf("[finufft_exec] done. tot spread:\t\t%.3g s\n",t_sprint);
        printf("               tot FFT:\t\t\t\t%.3g s\n", t_fft);
        printf("               tot deconvolve:\t\t\t%.3g s\n", t_deconv);
      } else {
        printf("[finufft_exec] done. tot deconvolve:\t\t%.3g s\n", t_deconv);
        printf("               tot FFT:\t\t\t\t%.3g s\n", t_fft);
        printf("               tot interp:\t\t\t%.3g s\n",t_sprint);
      }
    }
  }

  else {  // ----------------------------- TYPE 3 EXEC ---------------------

    double t_pre=0.0, t_spr=0.0, t_t2=0.0, t_deconv=0.0;  // accumulated timings
    if (p->opts.debug)
      printf("[finufft_exec t3] start ntrans=%d (%d batches, bsize=%d)...\n", p->ntrans, p->nbatch, p->batchSize);

    for (int b=0; b*p->batchSize < p->ntrans; b++) { // .....loop b over batches

      // batching and pointers to this batch, identical to t1,2 above...
      int thisBatchSize = min(p->ntrans - b*p->batchSize, p->batchSize);
      int bB = b*p->batchSize;
      CPX* cjb = cj + bB*p->nj;           // batch of input strengths
      CPX* fkb = fk + bB*p->nk;           // batch of output strengths
      if (p->opts.debug>1) printf("[finufft_exec t3] start batch %d (size %d):\n",b,thisBatchSize);
      
      // STEP 0: pre-phase (possibly) the c_j input strengths into c'_j batch...
      timer.restart();
#pragma omp parallel for
      for (int i=0; i<thisBatchSize; i++) {
        BIGINT ioff = i*p->nj;
        for (BIGINT j=0;j<p->nj;++j)
          p->CpBatch[ioff+j] = p->prephase[j] * cjb[ioff+j];
      }
      t_pre += timer.elapsedsec(); 
      
      // STEP 1: spread c'_j batch (x'_j NU pts) into fw batch grid...
      timer.restart();
      p->spopts.spread_direction = 1;                         // spread
      spreadinterpSortedBatch(thisBatchSize, p, p->CpBatch);  // p->X are primed
      t_spr += timer.elapsedsec();

      // STEP 2: type 2 NUFFT from fw batch to user output fk array batch...
      timer.restart();
      // illegal possible shrink of ntrans *after* plan for smaller last batch:
      p->innerT2plan->ntrans = thisBatchSize;      // do not try this at home!
      /* (alarming that FFTW not shrunk, but safe, because t2's fwBatch array
         still the same size, as Andrea explained; just wastes a few flops) */
      finufft_exec(p->innerT2plan, fkb, (CPX*)(p->fwBatch));
      t_t2 += timer.elapsedsec();

      // STEP 3: apply deconvolve (precomputed 1/phiHat(targ_k), phasing too)...
      timer.restart();
#pragma omp parallel for
      for (int i=0; i<thisBatchSize; i++) {
        BIGINT ioff = i*p->nk;
        for (BIGINT k=0;k<p->nk;++k)
          fkb[ioff+k] *= p->deconv[k];
      }
      t_deconv += timer.elapsedsec();
    }                                                   // ........end b loop

    if (p->opts.debug) {  // report total times in their natural order...
      printf("[finufft_exec t3] done. tot prephase:\t\t%.3g s\n",t_pre);
      printf("                  tot spread:\t\t\t%.3g s\n",t_spr);
      printf("                  tot type 2:\t\t\t%.3g s\n", t_t2);
      printf("                  tot deconvolve:\t\t%.3g s\n", t_deconv);
    }    
  }
  return 0; 
}


// DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
int finufft_destroy(finufft_plan* p)
// Free everything we allocated inside of finufft_plan pointed to by p.
// Also must not crash if called immediately after finufft_makeplan.
{ 
  FFTW_FR(p->fwBatch);   // free the big FFTW (or t3 spread) working array
  if (p->type==1 || p->type==2) {
    FFTW_DE(p->fftwPlan);
    free(p->phiHat1);
    free(p->phiHat2);
    free(p->phiHat3);
    free(p->sortIndices);
  } else {          // free the stuff alloc for type 3
    if (!p->innerT2plan)
      finufft_destroy(p->innerT2plan);
    free(p->CpBatch);
    free(p->Sp); free(p->Tp); free(p->Up);
    free(p->prephase);
    free(p->deconv);
  }
  return 0;
}