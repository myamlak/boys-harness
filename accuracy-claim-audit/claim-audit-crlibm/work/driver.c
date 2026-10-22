/* Driver for the crlibm claim audit.
 * Reads lines "func mode hexin" from stdin, writes "func mode hexout" to stdout.
 * func: exp log log2 log10 log1p expm1 sin cos tan asin acos atan sinh cosh
 *       sinpi cospi tanpi atanpi pow
 * mode: rn ru rd rz
 * hexin: 16 hex digits, big-endian bit pattern of the input double.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "crlibm.h"

static double from_hex(const char *s) {
  uint64_t u = strtoull(s, NULL, 16);
  double d;
  memcpy(&d, &u, 8);
  return d;
}

static void to_hex(double d, char *out) {
  uint64_t u;
  memcpy(&u, &d, 8);
  sprintf(out, "%016llx", (unsigned long long)u);
}

typedef double (*fn1)(double);
typedef int (*fn2)(double, double);

static fn1 lookup(const char *f, const char *m) {
  int rn = strcmp(m, "rn") == 0, ru = strcmp(m, "ru") == 0;
  int rd = strcmp(m, "rd") == 0, rz = strcmp(m, "rz") == 0;
#define P(a, b, c, d, e) \
  if (strcmp(f, a) == 0) { if (rn) return b; if (ru) return c; if (rd) return d; if (rz) return e; }
  P("exp",    exp_rn,    exp_ru,    exp_rd,    exp_rz)
  P("log",    log_rn,    log_ru,    log_rd,    log_rz)
  P("log2",   log2_rn,   log2_ru,   log2_rd,   log2_rz)
  P("log10",  log10_rn,  log10_ru,  log10_rd,  log10_rz)
  P("log1p",  log1p_rn,  log1p_ru,  log1p_rd,  log1p_rz)
  P("expm1",  expm1_rn,  expm1_ru,  expm1_rd,  expm1_rz)
  P("sin",    sin_rn,    sin_ru,    sin_rd,    sin_rz)
  P("cos",    cos_rn,    cos_ru,    cos_rd,    cos_rz)
  P("tan",    tan_rn,    tan_ru,    tan_rd,    tan_rz)
  P("asin",   asin_rn,   asin_ru,   asin_rd,   asin_rz)
  P("acos",   acos_rn,   acos_ru,   acos_rd,   acos_rz)
  P("atan",   atan_rn,   atan_ru,   atan_rd,   atan_rz)
  P("sinh",   sinh_rn,   sinh_ru,   sinh_rd,   sinh_rz)
  P("cosh",   cosh_rn,   cosh_ru,   cosh_rd,   cosh_rz)
  P("sinpi",  sinpi_rn,  sinpi_ru,  sinpi_rd,  sinpi_rz)
  P("cospi",  cospi_rn,  cospi_ru,  cospi_rd,  cospi_rz)
  P("tanpi",  tanpi_rn,  tanpi_ru,  tanpi_rd,  tanpi_rz)
  P("atanpi", atanpi_rn, atanpi_ru, atanpi_rd, atanpi_rz)
  P("sinpi",  sinpi_rn,  sinpi_ru,  sinpi_rd,  sinpi_rz)
  P("cospi",  cospi_rn,  cospi_ru,  cospi_rd,  cospi_rz)
  P("asinpi", asinpi_rn, asinpi_ru, asinpi_rd, asinpi_rz)
  P("acospi", acospi_rn, acospi_ru, acospi_rd, acospi_rz)
#undef P
  return NULL;
}

int main(void) {
  char f[32], m[8], h[32], out[32];
  crlibm_init();
  setvbuf(stdout, NULL, _IOFBF, 1 << 20);
  while (scanf("%31s %7s %31s", f, m, h) == 3) {
    if (strcmp(f, "pow") == 0 && strcmp(m, "rn") == 0) {
      char h2[32];
      if (scanf("%31s", h2) != 1) break;
      double x = from_hex(h), y = from_hex(h2);
      double r = pow_rn(x, y);
      to_hex(r, out);
      printf("pow rn %s %s %s\n", h, h2, out);
      continue;
    }
    fn1 g = lookup(f, m);
    if (!g) { printf("%s %s %s ERR\n", f, m, h); continue; }
    double x = from_hex(h);
    double r = g(x);
    to_hex(r, out);
    printf("%s %s %s %s\n", f, m, h, out);
  }
  fflush(stdout);
  return 0;
}
