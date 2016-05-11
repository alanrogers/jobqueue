#ifndef LDPSIZ_MISC
#define LDPSIZ_MISC

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#define TRUE (1)
#define FALSE (0)
#define UNPHASED_HETEROZYGOTE 4

#define ERR(code, msg) do{\
    fprintf(stderr,"%s:%s:%d: %s %d (%s)\n",\
            __FILE__,__func__,__LINE__,\
            (msg), (code), strerror((code)));   \
    exit(1);\
}while(0)

#include <math.h>
#include <float.h>

static inline int Dbl_near(double x, double y);
void        tellopt(const char *opt, const char *description);
char       *strlowercase(char *s);
void        checkmem( /*@null@ */ void *obj, const char *file, int line);
void        printsqrmat(const char *msg, unsigned dim, double m[][dim]);
void        printgslmat(const char *msg, size_t dim, gsl_matrix * m);
double      getreldiff(int dim, double x[], double y[], int verbose);
int         matIsFinite(unsigned dim, double m[][dim]);
void        assertFiniteArray(const double *x, size_t dim, const char *file,
                              int line);
char       *strcenter(const char *text, unsigned width, char *buff,
                      unsigned buffsize);
int         strCountSetChunks(const char *str, const char *sep);
int         getNumCores(void);
long        long_first_geq(long val, long *v, long len);
long        long_last_leq(long val, long *v, long len);
void        die(const char *msg, const char *file, int line);
void        eprintf(const char *fmt, ...);
void       *memdup(const void *p, size_t n);
int         strempty(const char *s);
int         strcomment(const char *s);
char       *stripComment(char *s);
int         dblEquals(double x, double y);
void        dostacktrace(const char *file, int line, FILE * ofp);
double      reflect(double x, double lo, double hi);
static inline int encodeDiploid(unsigned char *gtype, unsigned gtypeSize,
                                const char *str);
static inline int encodeHaploid(unsigned char *gtype, unsigned gtypeSize,
                                const char *str);
static inline unsigned encode01(char c);
void        unitTstResult(const char *facility, const char *result);

/**
 * Encode a single character, which should equal either '0' or '1' on
 * input.  Returned value is 0 for '0', 1 for '1', and 255 for anything
 * else.
 */
static inline unsigned encode01(char c) {
    if(c == '0')
        return 0;
    if(c == '1')
        return 1;
    return 255;
}

/**
 * Convert diploid genotype data from the character string format used
 * in input to the binary format used internally.
 *
 * In the input string, a genotype may be any of the following: "00",
 * "01", "10", "11", or "h". The four 2-character strings represent
 * phased genotypes. The 1-character "h" is an unphased heterozygote.
 * This function translates these codes into the integers 0, 1, 2, 3,
 * and UNPHASED_HETEROZYGOTE. The latter value is a macro defined
 * elsewhere.  
 *
 * @param[out] gtype an array of unsigned char values into which the
 * binary-ecoded genotype values will be written.
 *
 * @param[in] gtypeSize the size of the gtype array. No more than this
 * number of genotypes will be written into the array.
 *
 * @param[in] str The input, which represents genotypes as a
 * NULL-terminated character string. 
 *
 * @returns The number of genotypes written into array gtype. This is
 * *not* a NULL-terminated string, as the value 0 is valid in the
 * interior of the array.
 */
static inline int encodeDiploid(unsigned char *gtype, unsigned gtypeSize,
                                const char *str) {
    register unsigned char *g = gtype;
    register const char *s = str;
    register unsigned oneBit, twoBit;

    /*
     * The code below loops over the input string
     */
    while(*s != '\0' && g < gtype + gtypeSize) {
        if(*s == 'h') {
            *g = UNPHASED_HETEROZYGOTE;
        } else {
            twoBit = encode01(*s);
            if(twoBit > 1)
                eprintf("ERR@%s:%d: bad value in genotype: \"%c\". str=%s\n",
                        __FILE__, __LINE__, *s, str);

            /*
             * In the final genotype in the string is incorrect,
             * consisting of single character rather than two, then
             * the line below will read the '\0' at the end of the
             * string. This will generate an error but not a
             * segmentation fault.
             */
            ++s;
            if(*s == '\0')
                eprintf("ERR@%s:%d: : malformed input: %s\n",
                        __FILE__, __LINE__, str);
            oneBit = encode01(*s);
            if(oneBit > 1)
                eprintf("ERR@%s:%d: bad value in genotype: \"%c\". str=%s\n",
                        __FILE__, __LINE__, *s, str);

            *g = 2 * twoBit + oneBit;   /* either 0, 1, 2, or 3 */
        }
        ++g;
        ++s;
    }
    return g - gtype;
}

/**
 * Convert haploid genotype data from the character string format used
 * in input to the binary format used internally.
 *
 * In the input buffer a genotype may be either "0" or "1". 
 * This function translates these codes into the integers 0 and 1,
 * using the function encode01.
 *
 * @param[out] gtype an array of unsigned char values into which the
 * binary-ecoded genotype values will be written.
 *
 * @param[in] gtypeSize the size of the gtype array. No more than this
 * number of genotypes will be written into the array.
 *
 * @param[in] str The input, which represents genotypes as a
 * NULL-terminated character string. 
 *
 * @returns The number of genotypes written into array gtype. This is
 * *not* a NULL-terminated string, as the value 0 is valid in the
 * interior of the array.
 */
static inline int encodeHaploid(unsigned char *gtype, unsigned gtypeSize,
                                const char *str) {
    register unsigned char *g = gtype;
    register const char *s = str;

    while(*s != '\0' && g < gtype + gtypeSize) {
        *g = encode01(*s);
        ++g;
        ++s;
    }
    return g - gtype;
}

/**
 * Return 1 if the relative difference between x and y is less than or
 * equal to DBL_EPSILON.
 */
static inline int Dbl_near(double x, double y) {
    return fabs(x - y) <= fmax(fabs(x), fabs(y)) * DBL_EPSILON;
}

#ifndef NDEBUG
#define myassert(x) do { if (!(x)) { dostacktrace(__FILE__,__LINE__,stderr); assert(x); } } while(0)
#else
#define myassert(x)
#endif
#endif

#define REQUIRE(x,file,lineno) do { \
  if (!(x)) { \
    dostacktrace(__FILE__,__LINE__,stderr); \
    eprintf("ERR@%s:%d->%s:%d: Sanity check FAIL\n",\
            (file),(lineno),__FILE__,__LINE__); \
   }\
} while(0)
