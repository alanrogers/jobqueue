#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#ifdef _WIN32
#include <windows.h>
#elif defined(MACOS)
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif
#include <execinfo.h>
#include "misc.h"

/*
 * Describe an option. For use in "usage" functions.
 */
void tellopt(const char *opt, const char *description) {
    fprintf(stderr, "   %s\n      %s\n", opt, description);
    return;
}

/** Convert NULL-terminated string to lower case */
char       *strlowercase(char *s) {
    char       *p = s;

    for(p = s; *p != '\0'; ++p)
        *p = tolower(*p);
    return s;
}

void checkmem( /*@null@ */ void *obj, const char *file, int line) {
    if(obj == NULL)
        die("allocation error", file, line);
    return;
}

void assertFiniteArray(const double *x, size_t dim, const char *file,
                       int line) {
    if(!isfinite(x[cblas_idamax(dim, x, 1)])) {
        fprintf(stderr, "ERR@%s:%d: Array is not finite: [", file, line);
        while(dim > 1) {
            fprintf(stderr, "% g,", *x++);
            --dim;
        }
        if(dim)
            fprintf(stderr, " %g", *x);
        fprintf(stderr, "]\n");
        exit(EXIT_FAILURE);
    }
    return;
}

void printsqrmat(const char *msg, unsigned dim, double m[][dim]) {
    int         i, j;

    if(msg != NULL)
        printf("%s:\n", msg);
    for(i = 0; i < dim; ++i) {
        for(j = 0; j < dim; ++j)
            printf(" %12.3g", m[i][j]);
        putchar('\n');
    }
}

void printgslmat(const char *msg, size_t dim, gsl_matrix * m) {
    size_t      i, j;

    if(msg != NULL)
        printf("%s:\n", msg);
    for(i = 0; i < dim; ++i) {
        for(j = 0; j < dim; ++j)
            printf(" %8.4f", gsl_matrix_get(m, i, j));
        putchar('\n');
    }
}

int matIsFinite(unsigned dim, double m[][dim]) {
    unsigned    i, j;

    for(i = 0; i < dim; ++i)
        for(j = 0; j < dim; ++j)
            if(!isfinite(m[i][j]))
                return 0;
    return 1;
}

/*
 * Calculate the relative absolute difference between two vectors.
 */
double getreldiff(int dim, double x[], double y[], int verbose) {
    int         i;
    double      absdiff, abssum, relerr;

    abssum = absdiff = 0.0;
    for(i = 0; i < dim; ++i) {
        abssum += fabs(y[i]);
        absdiff += fabs(y[i] - x[i]);
        if(verbose)
            printf("%12.4g %12.4g\n", x[i], y[i]);
    }
    relerr = absdiff / abssum;

    if(verbose)
        printf("getreldiff: absdiff=%g abssum=%g relerr=%g\n",
               absdiff, abssum, relerr);

    return relerr;
}

/**
 * Center string "text" in a field of width "width". The centered string
 * is written into the character string "buff", whose size is
 * "buffsize".
 */
char       *strcenter(const char *text, unsigned width,
                      char *buff, unsigned buffsize) {
    int         i, j, lpad = 0, rpad = 0, txtwid;

    txtwid = strlen(text);
    if(txtwid >= buffsize) {
        snprintf(buff, buffsize, "%s", text);
        return buff;
    }
    if(width > txtwid)
        lpad = (width - txtwid) / 2;
    rpad = width - txtwid - lpad;
    for(i = 0; i < lpad; ++i)
        buff[i] = '-';
    for(j = 0; j < txtwid; ++j)
        buff[i + j] = text[j];
    for(j = 0; j < rpad; ++j)
        buff[i + txtwid + j] = '-';
    buff[lpad + txtwid + rpad] = '\0';
    return buff;
}

/*
 * An almost platform-independent function that returns the number of
 * CPU cores on the current machine.
 * Source: http://stackoverflow.com/questions/150355/
 * programmatically-find-the-number-of-cores-on-a-machine. If I
 * understand the webpage correctly, this code was written by Dirk-Jan
 * Kroon.
 */
int getNumCores(void) {
#ifdef WIN32
    SYSTEM_INFO sysinfo;

    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif defined(MACOS)
    int         nm[2];
    size_t      len = 4;
    uint32_t    count;

    nm[0] = CTL_HW;
    nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) {
            count = 1;
        }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

/*
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of first element in sorted array
 * that is >= val.  The function assumes without checking that the
 * input is sorted. If val > vec[len-1], the function returns len.
 */
long long_first_geq(long val, long *v, long len) {
    register long lo, mid, hi;

    myassert(len > 0);
    lo = 0;
    hi = len - 1;
    if(val > v[hi])
        return len;
    while(lo < hi) {
        mid = lo + (hi - lo) / 2;
        if(mid == lo)
            break;
        if(v[mid] < val)
            lo = mid;
        else
            hi = mid;
    }
    if(v[lo] >= val)
        hi = lo;

    myassert(hi >= 0);
    myassert(hi < len);
    myassert(v[hi] >= val);
    myassert(hi == 0 || v[hi - 1] < val);

    return hi;
}

/*
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of last element in sorted array
 * that is <= val.  The function assumes without checking that the
 * input is sorted. If val < vec[0], the function returns -1.
 */
long long_last_leq(long val, long *v, long len) {
    register long lo, mid, hi;

    myassert(len > 0);
    lo = 0;
    hi = len - 1;
    if(val < v[0])
        return -1;
    while(lo < hi) {
        mid = hi - (hi - lo) / 2;
        if(mid == hi)
            break;
        if(v[mid] > val)
            hi = mid;
        else
            lo = mid;
    }
    if(v[hi] <= val)
        lo = hi;

    myassert(lo >= 0);
    myassert(lo < len);
    myassert(v[lo] <= val);
    myassert(lo == len - 1 || v[lo + 1] > val);

    return lo;
}

/* print message and exit */
void die(const char *msg, const char *file, int line) {
    fflush(stdout);
    fprintf(stderr, "ERR@%s:%d: %s \n", file, line, msg);
    exit(EXIT_FAILURE);
}

/* eprintf: print error message and exit */
void eprintf(const char *fmt, ...) {
    va_list     args;

    fflush(stdout);

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    if(fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
        fprintf(stderr, " %s", strerror(errno));
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

/* duplicate memory block */
void       *memdup(const void *p, size_t n) {
    void       *q;

    myassert(p != NULL);
    myassert(n > 0);

    q = malloc(n);
    checkmem(q, __FILE__, __LINE__);
    memcpy(q, p, n);
    return q;
}

/*
 * In string str, count the number of contiguous chunks of characters
 * belonging to set.
 */
int strCountSetChunks(const char *str, const char *sep) {
    int         nchunks = 0, i;

    while(*str != '\0') {
        i = strcspn(str, sep);  /* skip chars not in set */
        str += i;
        i = strspn(str, sep);   /* skip chars in set */
        if(i > 0) {
            ++nchunks;
            str += i;
        }
    }
    return nchunks;
}

/** Return 1 if the first non-white char in string s is '#'; 0 otherwise. */
int strcomment(const char *s) {
    const char *p = s;

    while(isspace(*p))
        ++p;
    if(*p == '#')
        return 1;
    return 0;
}

/** strip comment ('#' to eol) from a string */
char       *stripComment(char *s) {
    char       *p = strchr(s, '#');

    if(p && (*p == '#'))
        *p = '\0';

    return s;
}

/** Return 1 if string contains only whitespace; 0 otherwise. */
int strempty(const char *s) {
    const char *p = s;

    while(isspace(*p))
        ++p;
    if(*p == '\0')
        return 1;
    return 0;
}

/**
 * Test equality of x and y. In this test NaN == NaN, Inf == Inf, and
 * -Inf == -Inf.
 */
#pragma GCC diagnostic ignored "-Wfloat-equal"
int dblEquals(double x, double y) {
    int         type_x, type_y;

    type_x = fpclassify(x);
    type_y = fpclassify(y);

    if(type_y != type_x)
        return 0;               /* x != y because types differ */

    assert(type_x == type_y);

    if(type_y == FP_NAN)
        return 1;

    return (x == y);
}

#pragma GCC diagnostic pop

#define CALLSTACK_SIZE 128
void dostacktrace(const char *file, int line, FILE * ofp) {
    void       *callstack[CALLSTACK_SIZE];
    int         nsymbols = backtrace(callstack, CALLSTACK_SIZE);

    fprintf(ofp, "backtrace returned %d\n", nsymbols);
    fprintf(ofp, "dostacktrace called from %s:%d:\n", file, line);
    backtrace_symbols_fd(callstack, nsymbols, fileno(ofp));
}

/**
 * Fold x back and forth across the boundaries "lo" and "hi" to obtain a value
 * y such that lo <= y <= hi.
 */
double reflect(double x, double lo, double hi) {
    assert(hi > lo);
    x -= lo;
    hi -= lo;

    double      z = fabs(fmod(x, 2.0 * hi));

    /*    printf("initially z=%lg\n", z); */

    if(z > hi)
        z = 2.0 * hi - z;

    z += lo;

    assert(z >= lo && z <= hi + lo);

    return z;
}

void unitTstResult(const char *facility, const char *result) {
    printf("%-26s %s\n", facility, result);
}

#ifdef TEST

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <float.h>

int main(int argc, char **argv) {

    int         verbose = 0;

    long        v[] = { 0, 0, 1, 1, 1, 2, 2 };
    long        len = 7;
    const char *str1 = "abc1cd23efgh4";
    const char *str2 = "999999999abc1cd23efgh4";
    const char *str3 = "abc1cd23efgh";
    const char *str4 = "abccdefgh";
    const char *set = "0123456789";

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xmisc [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xmisc [-v]\n");
    }

    assert(dblEquals(0.0 / 0.0, 0.0 / 0.0));
    assert(!dblEquals(0.0 / 0.0, DBL_MAX));
    assert(dblEquals(1.0 / 0.0, 1.0 / 0.0));
    assert(!dblEquals(-1.0 / 0.0, 1.0 / 0.0));
    assert(dblEquals(1.23, 1.23));
    unitTstResult("dblEquals", "OK");

    assert(long_last_leq(-1, v, len) == -1);
    assert(long_last_leq(0, v, len) == 1);
    assert(long_last_leq(1, v, len) == 4);
    assert(long_last_leq(2, v, len) == 6);
    assert(long_last_leq(3, v, len) == 6);
    assert(long_last_leq(-1, v, 1) == -1);
    assert(long_last_leq(0, v, 1) == 0);
    assert(long_last_leq(1, v, 1) == 0);
    unitTstResult("long_last_leq", "OK");

    assert(long_first_geq(-1, v, len) == 0);
    assert(long_first_geq(0, v, len) == 0);
    assert(long_first_geq(1, v, len) == 2);
    assert(long_first_geq(2, v, len) == 5);
    assert(long_first_geq(3, v, len) == 7);
    assert(long_first_geq(-1, v, 1) == 0);
    assert(long_first_geq(0, v, 1) == 0);
    assert(long_first_geq(1, v, 1) == 1);
    unitTstResult("long_first_geq", "OK");

    assert(strCountSetChunks(str1, set) == 3);
    assert(strCountSetChunks(str2, set) == 4);
    assert(strCountSetChunks(str3, set) == 2);
    assert(strCountSetChunks(str4, set) == 0);
    assert(strCountSetChunks("", set) == 0);
    assert(strCountSetChunks(set, set) == 1);
    unitTstResult("strCountSetChunks", "OK");

    double      x[] = { 1.0, 3.0, 4.3 };
    assertFiniteArray(x, 3, __FILE__, __LINE__);
    unitTstResult("assertFiniteArray", "OK");

#if 0
    x[1] = 1.0 / 0.0;
    printf("next line should die\n");
    assertFiniteArray(x, 3, __FILE__, __LINE__);
#endif

    assert(strcomment("   ab cde") == 0);
    assert(strcomment("   #ab cde") == 1);
    unitTstResult("strcomment", "OK");

    char        buff[30], *s;

    snprintf(buff, sizeof(buff), " asdfaf #comment");

    s = stripComment(buff);
    assert(s == buff);
    assert(strcmp(buff, " asdfaf ") == 0);
    assert(strlen(buff) == strlen(" asdfaf "));
    unitTstResult("stripComment", "OK");

    assert(reflect(0.0, 1.0, 2.0) == 2.0);
    assert(reflect(1.5, 1.0, 2.0) == 1.5);
    assert(reflect(1.0, 1.0, 2.0) == 1.0);
    assert(reflect(2.0, 1.0, 2.0) == 2.0);
    assert(reflect(2.25, 1.0, 2.0) == 1.75);
    assert(reflect(3.25, 1.0, 2.0) == 1.25);
    assert(reflect(4.75, 1.0, 2.0) == 1.25);
    unitTstResult("reflect", "OK");

    assert(encode01('0') == 0);
    assert(encode01('1') == 1);
    assert(encode01('2') == 255);
    assert(encode01('h') == 255);
    unitTstResult("encode01", "OK");

    unsigned char gtype[100];
    unsigned    i, nGtype;
    const char *gtypeString = "1001h1100";

    nGtype = encodeDiploid(gtype, sizeof(gtype), gtypeString);
    if(verbose) {
        printf("diploid input: %s\n", gtypeString);
        printf("nGtype=%u\ngtype:", nGtype);
        for(i = 0; i < nGtype; ++i)
            printf(" %u", (unsigned) gtype[i]);
        putchar('\n');
    }
    assert(nGtype == 5);
    assert(gtype[0] == 2);
    assert(gtype[1] == 1);
    assert(gtype[2] == UNPHASED_HETEROZYGOTE);
    assert(gtype[3] == 3);
    assert(gtype[4] == 0);

    unitTstResult("encodeDiploid", "OK");

    gtypeString = "10011";
    nGtype = encodeHaploid(gtype, sizeof(gtype), gtypeString);
    if(verbose) {
        printf("haploid input: %s\n", gtypeString);
        printf("nGtype=%u\ngtype:", nGtype);
        for(i = 0; i < nGtype; ++i)
            printf(" %u", (unsigned) gtype[i]);
        putchar('\n');
    }
    assert(nGtype == 5);
    assert(gtype[0] == 1);
    assert(gtype[1] == 0);
    assert(gtype[2] == 0);
    assert(gtype[3] == 1);
    assert(gtype[4] == 1);
    unitTstResult("encodeHaploid", "OK");

    return 0;
}
#endif
