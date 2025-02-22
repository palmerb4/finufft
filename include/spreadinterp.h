// Defines interface to spreading/interpolation code.
// Note: see defs.h for definition of MAX_NSPREAD (as of 9/24/18).
// RESCALE macro moved to spreadinterp.cpp, 7/15/20.

#ifndef SPREADINTERP_H
#define SPREADINTERP_H

#include <dataTypes.h>
#include <spread_opts.h>

/* Bitwise debugging timing flag (TF) definitions; see spread_opts.flags.
    This is an unobtrusive way to determine the time contributions of the
    different components of spreading/interp by selectively leaving them out.
    For example, running the following two tests shows the effect of the exp()
    in the kernel evaluation (the last argument is the flag):
    > perftest/spreadtestnd 3 8e6 8e6 1e-6 1 0 0 1 0
    > perftest/spreadtestnd 3 8e6 8e6 1e-6 1 4 0 1 0
    NOTE: non-zero values are for experts only, since
    NUMERICAL OUTPUT MAY BE INCORRECT UNLESS spread_opts.flags=0 !
*/
#define TF_OMIT_WRITE_TO_GRID 1        // don't add subgrids to out grid (dir=1)
#define TF_OMIT_EVALUATE_KERNEL 2      // don't evaluate the kernel at all
#define TF_OMIT_EVALUATE_EXPONENTIAL 4 // omit exp() in kernel (kereval=0 only)
#define TF_OMIT_SPREADING 8            // don't interp/spread (dir=1: to subgrids)

// things external (spreadinterp) interface needs...
int spreadinterp(BIGINT N1, BIGINT N2, BIGINT N3, BIGINT N4, BIGINT N5, FLT *data_uniform, BIGINT M, FLT *kx, FLT *ky,
                 FLT *kz, FLT *kp, FLT *kq, FLT *data_nonuniform, spread_opts opts);
int spreadcheck(BIGINT N1, BIGINT N2, BIGINT N3, BIGINT N4, BIGINT N5, BIGINT M, FLT *kx, FLT *ky, FLT *kz, FLT *kp,
                FLT *kq, spread_opts opts);
int indexSort(BIGINT *sort_indices, BIGINT N1, BIGINT N2, BIGINT N3, BIGINT N4, BIGINT N5, BIGINT M, FLT *kx, FLT *ky,
              FLT *kz, FLT *kp, FLT *kq, spread_opts opts);
int interpSorted(BIGINT *sort_indices, BIGINT N1, BIGINT N2, BIGINT N3, BIGINT N4, BIGINT N5, FLT *data_uniform,
                 BIGINT M, FLT *kx, FLT *ky, FLT *kz, FLT *kp, FLT *kq, FLT *data_nonuniform, spread_opts opts,
                 int did_sort);
int spreadSorted(BIGINT *sort_indices, BIGINT N1, BIGINT N2, BIGINT N3, BIGINT N4, BIGINT N5, FLT *data_uniform,
                 BIGINT M, FLT *kx, FLT *ky, FLT *kz, FLT *kp, FLT *kq, FLT *data_nonuniform, spread_opts opts,
                 int did_sort);
int spreadinterpSorted(BIGINT *sort_indices, BIGINT N1, BIGINT N2, BIGINT N3, BIGINT N4, BIGINT N5, FLT *data_uniform,
                       BIGINT M, FLT *kx, FLT *ky, FLT *kz, FLT *kp, FLT *kq, FLT *data_nonuniform, spread_opts opts,
                       int did_sort);
FLT evaluate_kernel(FLT x, const spread_opts &opts);
FLT evaluate_kernel_noexp(FLT x, const spread_opts &opts);
int setup_spreader(spread_opts &opts, FLT eps, double upsampfac, int kerevalmeth, int debug, int showwarn, int dim);

#endif // SPREADINTERP_H
