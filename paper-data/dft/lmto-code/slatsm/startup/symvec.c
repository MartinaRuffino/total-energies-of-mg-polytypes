/* symvec.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
   NB: This was slightly modified from f2c to use pointers.
   Alternatively, compile with -DNO_LINUX.
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static double c_b19 = 1.;

/* #define unix */
/* #define POINTER */
/* Subroutine */ int addsvv_0_(n__, nam, nelt, ival, iopt, i1, i2, vec, nvar, 
	first, last, ifi, nam_len)
int n__;
char *nam;
integer *nelt, *ival, *iopt, *i1, *i2;
double *vec;
integer *nvar, *first, *last, *ifi;
ftnlen nam_len;
{
    /* Initialized data */

    static integer nnam = 0;

    /* Format strings */
    static char fmt_332[] = "(\002  Vec       Name            Size   Val[1..\
n]\002)";
    static char fmt_333[] = "(i4,4x,a20,i4,2g14.5)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop();
/*     integer s_wsfe(), e_wsfe(), do_fio(), s_cmp(); */
    void s_copy(); integer s_cmp();

    /* Local variables */
    static double *iarr,*symptr[24];
    static integer size[24], i__;
    extern /* Subroutine */ int faloc_();
    static integer io;
    extern /* Subroutine */ int locase_();
    static char tmpnam[16], symnam[16*24];
    extern /* Subroutine */ int ptrcop_();
    static integer i2x;
    static double arr[20];

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 0, 0, fmt_332, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_333, 0 };


/* - Add a symbolic vector to list */
/* ---------------------------------------------------------------- */
/* i Inputs */
/* i   nam:  name of variable */
/* i   nelt: number of elements of the vector */
/* o Outputs */
/* o   ival  index to which variable is declared or accessed */
/* r Remarks */
/* r   addsvv  adds a symbolic name and value to the internal table, */
/* r           and allocates memory for the vector. */
/* r   lodsvv  sets a range of elements of a vector associated with */
/* r           a name or an index, depending on iopt. */
/* r           iopt=0: index associated with name */
/* r           iopt=1: name associated with index */
/* r   getsvv  gets a range of elements of a vector associated with */
/* r           a name or an index, depending on iopt. */
/* r   sizsvv  returns the length of a vector associated with */
/* r           a name or an index, depending on iopt. */
/* r   numsvv  returns the number of variables now declared */
/* r   watsvv  returns name associated with index */
/* r   This implementation works only with fortran compilers allowing */
/* r   pointers. */
/* ---------------------------------------------------------------- */
/*     implicit none */
/* Passed parameters */
/* Local parameters */
/* #ifdefC POINTER */
/*      pointer (iarr, arr) */
/* #elseifC IPOINTER */
/* #endif */
    /* Parameter adjustments */
    if (vec) {
	--vec;
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_lodsvv;
	case 2: goto L_getsvv;
	case 3: goto L_sizsvv;
	case 4: goto L_numsvv;
	case 5: goto L_watsvv;
	case 6: goto L_shosvv;
	}

/* --- Start of addsvv --- */
    ++nnam;
    if (nnam > 24) {rx_("addsvv: too many names", 22L);}

    s_copy(symnam + (nnam - 1 << 4), nam, 16L, nam_len);
/* #ifdef unix */
    locase_(symnam + (nnam - 1 << 4), 16L);
/* #endif */
    *ival = nnam;
    faloc_(&iarr, &c__4, nelt);
    symptr[nnam - 1] = iarr;
    size[nnam - 1] = *nelt;
    return 0;
/* --- lodsvv, getsvv --- */

L_lodsvv:
    io = -1;
    goto L10;

L_getsvv:
    io = 1;
    goto L10;

L_sizsvv:
    io = -2;
    goto L10;
/* --- lodsvv, getsvv --- */

L_numsvv:
    *nvar = nnam;
    return 0;
/* --- watsvv --- */

L_watsvv:
    s_copy(nam, " ", nam_len, 1L);
    if (*ival <= nnam) {
	s_copy(nam, symnam + (*ival - 1 << 4), nam_len, 16L);
    }
    return 0;
/* --- Print out table --- */

L_shosvv:
    io___7.ciunit = *ifi;
#ifndef NO_LINUX
    s_wsfe(&io___7);
    e_wsfe();
    i__1 = min(*last,nnam);
    for (i__ = max(*first,1); i__ <= i__1; ++i__) {
	iarr = symptr[i__ - 1];
/* L60: */
	io___9.ciunit = *ifi;
	s_wsfe(&io___9);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, symnam + (i__ - 1 << 4), 16L);
	do_fio(&c__1, (char *)&size[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iarr[0], (ftnlen)sizeof(double));
	do_fio(&c__1, (char *)&iarr[size[i__ - 1] - 1], (ftnlen)sizeof(
		double));
	e_wsfe();
    }
#endif
    return 0;
/* --- Find an index associated with a name --- */
L10:
/* ... If iopt=0, find the index associated with this name */
    if (*iopt == 0) {
	*ival = 0;
/* #ifndefC POINTER */
/*        return */
/* #endif */
	s_copy(tmpnam, nam, 16L, nam_len);
/* #ifdef unix */
	locase_(tmpnam, 16L);
/* #endif */
	i__1 = nnam;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (s_cmp(tmpnam, symnam + (i__ - 1 << 4), 16L, 16L) != 0) {
		goto L16;
	    }
	    *ival = i__;
	    goto L20;
L16:
	    ;
	}
    }
/* --- Set/Retrieve portions of an array[index], depending on io --- */
L20:
    if (io == 0) {
	return 0;
    }
    if (io == -2) {
	*i1 = size[*ival - 1];
	return 0;
    }
/* ... Return unless ival in range */
    if (*ival <= 0 || *ival > nnam) {
	return 0;
    }
    iarr = symptr[*ival - 1];
/* Computing MIN */
    i__1 = *i2, i__2 = size[*ival - 1];
    i2x = min(i__1,i__2);
    if (i2x < *i1) {
	return 0;
    }
/*      if (io == -1) arr(i1) = -99 */
/*      if (io == 1) vec(1) = -99 */
/*      if (io == -1) print *,'copy to',i1,vec(1),arr(i1) */
/*      if (io == 1) print *,'copy from',i1,vec(1),arr(i1) */
    i__1 = i2x - *i1 + 1;
    ptrcop_(&iarr, &vec[1], &i__1, i1, &c__1, &c_b19, &io);
/*     if (io == -1) call dpscop(vec,arr,i2x-i1+1,1,i1,1d0) */
/*     if (io == 1) call dpscop(arr,vec,i2x-i1+1,i1,1,1d0) */
/*      print *, 'done',vec(1),arr(i1) */
} /* addsvv_ */

/* Subroutine */ int addsvv_(nam, nelt, ival, nam_len)
char *nam;
integer *nelt, *ival;
ftnlen nam_len;
{
    return addsvv_0_(0, nam, nelt, ival, (integer *)0, (integer *)0, (integer 
	    *)0, (double *)0, (integer *)0, (integer *)0, (integer *)0, (
	    integer *)0, nam_len);
    }

/* Subroutine */ int lodsvv_(nam, ival, iopt, i1, i2, vec, nam_len)
char *nam;
integer *ival, *iopt, *i1, *i2;
double *vec;
ftnlen nam_len;
{
    return addsvv_0_(1, nam, (integer *)0, ival, iopt, i1, i2, vec, (integer *
	    )0, (integer *)0, (integer *)0, (integer *)0, nam_len);
    }

/* Subroutine */ int getsvv_(nam, ival, iopt, i1, i2, vec, nam_len)
char *nam;
integer *ival, *iopt, *i1, *i2;
double *vec;
ftnlen nam_len;
{
    return addsvv_0_(2, nam, (integer *)0, ival, iopt, i1, i2, vec, (integer *
	    )0, (integer *)0, (integer *)0, (integer *)0, nam_len);
    }

/* Subroutine */ int sizsvv_(nam, ival, iopt, i1, nam_len)
char *nam;
integer *ival, *iopt, *i1;
ftnlen nam_len;
{
    return addsvv_0_(3, nam, (integer *)0, ival, iopt, i1, (integer *)0, (
	    double *)0, (integer *)0, (integer *)0, (integer *)0, (
	    integer *)0, nam_len);
    }

/* Subroutine */ int numsvv_(nvar)
integer *nvar;
{
    return addsvv_0_(4, (char *)0, (integer *)0, (integer *)0, (integer *)0, (
	    integer *)0, (integer *)0, (double *)0, nvar, (integer *)0, (
	    integer *)0, (integer *)0, (ftnint)0);
    }

/* Subroutine */ int watsvv_(nam, ival, nam_len)
char *nam;
integer *ival;
ftnlen nam_len;
{
    return addsvv_0_(5, nam, (integer *)0, ival, (integer *)0, (integer *)0, (
	    integer *)0, (double *)0, (integer *)0, (integer *)0, (
	    integer *)0, (integer *)0, nam_len);
    }

/* Subroutine */ int shosvv_(first, last, ifi)
integer *first, *last, *ifi;
{
    return addsvv_0_(6, (char *)0, (integer *)0, (integer *)0, (integer *)0, (
	    integer *)0, (integer *)0, (double *)0, (integer *)0, first, 
	    last, ifi, (ftnint)0);
    }

/* Subroutine */ int parsvv_(recrd, recl, indx, mxelt, i1, ip, recrd_len)
char *recrd;
integer *recl, *indx, *mxelt, *i1, *ip;
ftnlen recrd_len;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer nelt;
    extern integer a2vec_();
    static integer i__, k, ix;
    extern /* Subroutine */ int skipbl_(), lodsvv_();
    static double res;

/* - Parses a string for one or more elements of a vector variable */
/*     implicit none */
/* Passed parameters */
/* Local parameters */
    nelt = 0;
    for (i__ = 1; i__ <= 999; ++i__) {
	skipbl_(recrd, recl, ip, 100L);
	if (*ip >= *recl || nelt >= *mxelt) {
	    goto L99;
	}
	k = a2vec_(recrd, recl, ip, &c__4, " ", &c__1, &c__1, &c__1, &ix, &
		res, 100L, 1L);
	if (k == -1) {
	    return 0;
	}
	i__1 = *i1 + nelt;
	i__2 = *i1 + nelt;
	lodsvv_(" ", indx, &c__1, &i__1, &i__2, &res, 1L);
	nelt += k;
/* L33: */
    }
L99:
    ;
} /* parsvv_ */

