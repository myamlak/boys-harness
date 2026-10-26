/* crlibm sinh witness driver.
 *
 * Reads one 64-bit input pattern per line from the file named on the command
 * line and prints, for each input, the bit patterns returned by crlibm's four
 * sinh routines (rn, ru, rd, rz).  Bits are authoritative; the decimal column
 * is printed after the hardware rounding mode is restored to nearest.
 *
 * Build and run instructions are in ../README.md.
 */
#include <fenv.h>
#include <float.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include <crlibm.h>

static double bits_to_double(uint64_t b) {
  union {
    uint64_t u;
    double d;
  } v;
  v.u = b;
  return v.d;
}

static uint64_t double_to_bits(double d) {
  union {
    uint64_t u;
    double d;
  } v;
  v.d = d;
  return v.u;
}

int main(int argc, char **argv) {
  FILE *f;
  char line[256];
  uint64_t b;

  if (argc != 2) {
    fprintf(stderr, "usage: %s inputs.txt\n", argv[0]);
    return 2;
  }
  f = fopen(argv[1], "r");
  if (f == NULL) {
    perror(argv[1]);
    return 2;
  }

  printf("# FLT_EVAL_METHOD=%d sizeof(double)=%zu sizeof(long double)=%zu"
         " fegetround_at_start=%d\n",
         (int)FLT_EVAL_METHOD, sizeof(double), sizeof(long double),
         (int)fegetround());
  printf("input_bits,x,sinh_rn_bits,sinh_rn,sinh_ru_bits,sinh_ru,"
         "sinh_rd_bits,sinh_rd,sinh_rz_bits,sinh_rz\n");

  while (fgets(line, sizeof line, f) != NULL) {
    double x, rn, ru, rd, rz;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
      continue;
    }
    if (sscanf(line, "%" SCNx64, &b) != 1) {
      continue;
    }
    x = bits_to_double(b);
    rn = sinh_rn(x);
    ru = sinh_ru(x);
    rd = sinh_rd(x);
    rz = sinh_rz(x);
    /* crlibm promises to restore the caller's rounding mode; do it anyway so
     * the decimal rendering below is independent of any library slip. */
    fesetround(FE_TONEAREST);
    printf("%016" PRIx64 ",%.17g,%016" PRIx64 ",%.17g,%016" PRIx64 ",%.17g,"
           "%016" PRIx64 ",%.17g,%016" PRIx64 ",%.17g\n",
           b, x, double_to_bits(rn), rn, double_to_bits(ru), ru,
           double_to_bits(rd), rd, double_to_bits(rz), rz);
  }
  fclose(f);
  return 0;
}
