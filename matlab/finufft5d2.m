% FINUFFT5D2   5D complex nonuniform FFT of type 2 (uniform to nonuniform).
%
% c = finufft5d2(x,y,z,p,q,isign,eps,f)
% c = finufft5d2(x,y,z,p,q,isign,eps,f,opts)
%
% This computes, to relative precision eps, via a fast algorithm:
%
%    c[j] =   SUM   f[k1,k2,k3,k4,k5] exp(+/-i (k1 x[j] + k2 y[j]
%                                             + k3 z[j] + k4 p[j]
%                                             + k5 q[j]))
%       k1,k2,k3,k4,k5
%                            for j = 1,..,nj
%     where sum is over -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2,
%                       -mu/2 <= k3 <= (mu-1)/2, -mv/2 <= k4 <= (mv-1)/2,
%                       -mw/2 <= k5 <= (mw-1)/2
%
%  Inputs:
%     x,y,z,p,q coordinates of nonuniform targets on the hypercube [-3pi,3pi)^5,
%           each a vector of length nj
%     f     complex Fourier coefficient array, whose size sets (ms,mt,mu,mv).
%           (Mode ordering given by opts.modeord, in each dimension.)
%           If a 6D array, 6th dimension sets ntrans, and each of ntrans
%           5D arrays is transformed with the same nonuniform targets.
%     isign if >=0, uses + sign in exponential, otherwise - sign.
%     eps   relative precision requested (generally between 1e-15 and 1e-1)
%     opts   optional struct with optional fields controlling the following:
%     opts.debug:   0 (silent, default), 1 (timing breakdown), 2 (debug info).
%     opts.spread_debug: spreader: 0 (no text, default), 1 (some), or 2 (lots)
%     opts.spread_sort:  0 (don't sort NU pts), 1 (do), 2 (auto, default)
%     opts.spread_kerevalmeth:  0: exp(sqrt()), 1: Horner ppval (faster)
%     opts.spread_kerpad: (iff kerevalmeth=0)  0: don't pad to mult of 4, 1: do
%     opts.fftw: FFTW plan mode, 64=FFTW_ESTIMATE (default), 0=FFTW_MEASURE, etc
%     opts.upsampfac:   sigma.  2.0 (default), or 1.25 (low RAM, smaller FFT)
%     opts.spread_thread:   for ntrans>1 only. 0:auto, 1:seq multi, 2:par, etc
%     opts.maxbatchsize:  for ntrans>1 only. max blocking size, or 0 for auto.
%     opts.nthreads:   number of threads, or 0: use all available (default)
%     opts.modeord: 0 (CMCL increasing mode ordering, default), 1 (FFT ordering)
%     opts.chkbnds: 0 (don't check NU points valid), 1 (do, default)
%  Outputs:
%     c     complex column vector of nj answers at targets, or,
%           if ntrans>1, matrix of size (nj,ntrans).
%
% Notes:
%  * The vectorized (many vector) interface, ie ntrans>1, can be much faster
%    than repeated calls with the same nonuniform points. Note that here the I/O
%    data ordering is stacked rather than interleaved. See ../docs/matlab.rst
%  * The class of input x (double vs single) controls whether the double or
%    single precision library are called; precisions of all data should match.
%  * For more details about the opts fields, see ../docs/opts.rst
%  * See ERRHANDLER, VALID_* and FINUFFT_PLAN for possible warning/error IDs.
%  * Full documentation is given in ../finufft-manual.pdf and online at
%    http://finufft.readthedocs.io

function c = finufft5d2(x,y,z,p,q,isign,eps,f,o)

if nargin<8, o.dummy=1; end
valid_setpts(2,5,x,y,z,p,q);
o.floatprec=class(x);                      % should be 'double' or 'single'
[ms,mt,mu,mv,mw,n_transf] = size(f);       % if f 5D array, n_transf=1
plan = finufft_plan(2,[ms;mt;mu;mv;mw],isign,n_transf,eps,o);
plan.setpts(x,y,z,p,q);
c = plan.execute(f);
