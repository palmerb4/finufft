#ifndef DIRFT_H
#define DIRFT_H

#include <dataTypes.h>

void dirft1d1(BIGINT nj, FLT *x, CPX *c, int isign, BIGINT ms, CPX *f);
void dirft1d2(BIGINT nj, FLT *x, CPX *c, int iflag, BIGINT ms, CPX *f);
void dirft1d3(BIGINT nj, FLT *x, CPX *c, int iflag, BIGINT nk, FLT *s, CPX *f);

void dirft2d1(BIGINT nj, FLT *x, FLT *y, CPX *c, int iflag, BIGINT ms, BIGINT mt, CPX *f);
void dirft2d2(BIGINT nj, FLT *x, FLT *y, CPX *c, int iflag, BIGINT ms, BIGINT mt, CPX *f);
void dirft2d3(BIGINT nj, FLT *x, FLT *y, CPX *c, int iflag, BIGINT nk, FLT *s, FLT *t, CPX *f);

void dirft3d1(BIGINT nj, FLT *x, FLT *y, FLT *z, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu, CPX *f);
void dirft3d2(BIGINT nj, FLT *x, FLT *y, FLT *z, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu, CPX *f);
void dirft3d3(BIGINT nj, FLT *x, FLT *y, FLT *z, CPX *c, int iflag, BIGINT nk, FLT *s, FLT *t, FLT *u, CPX *f);

void dirft4d1(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv,
              CPX *f);
void dirft4d2(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu, BIGINT mv,
              CPX *f);
void dirft4d3(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, CPX *c, int iflag, BIGINT nk, FLT *s, FLT *t, FLT *u, FLT *v,
              CPX *f);

void dirft5d1(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu,
              BIGINT mv, BIGINT mw, CPX *f);
void dirft5d2(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int iflag, BIGINT ms, BIGINT mt, BIGINT mu,
              BIGINT mv, BIGINT mw, CPX *f);
void dirft5d3(BIGINT nj, FLT *x, FLT *y, FLT *z, FLT *p, FLT *q, CPX *c, int iflag, BIGINT nk, FLT *s, FLT *t, FLT *u,
              FLT *v, FLT *w, CPX *f);
#endif
