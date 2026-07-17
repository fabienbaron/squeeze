/* philox.h : Philox4x32-10 counter-based random number generator.

   Algorithm: Salmon, Moraes, Dror & Shaw, "Parallel Random Numbers: As Easy as 1, 2, 3",
   SC'11. This is a self-contained minimal implementation of the parts SQUEEZE needs.

   Why counter-based rather than a conventional generator: a stream here is a pure function
   of its id, so chain i always gets stream i no matter which thread reaches it first, and
   there is no shared state to race over. The generator this replaced handed out streams by
   mutating a shared seed, which meant two chains could be given the same stream (measured
   at ~2 runs in 10 with 16 chains) and the chain/stream pairing changed between runs.

   Statistical quality is not the reason for the choice: MRG32k3a, Philox and PCG score
   within 0.1% of each other on TestU01 BigCrush over 1000+ streams (Wartel & Hill 2026).

   Draws are cheap (~5 ns) next to one likelihood evaluation over nuv points, so callers
   should take a fresh variate per decision rather than sharing bits between decisions.
*/
#ifndef PHILOX_H
#define PHILOX_H

#include <stdint.h>

typedef struct {
  uint32_t ctr[4]; /* {draw counter lo, hi, stream id, 0} */
  uint32_t key[2]; /* {seed, Weyl constant} */
  uint32_t buf[4]; /* one block yields four 32-bit outputs */
  int have;        /* outputs still unused in buf */
} Philox;

/* Every SQUEEZE stream needs a distinct id. Chains use their chain index, so these
   reserved ids must stay clear of [0, nchains). */
#define PHILOX_STREAM_INIT 0x7FFF0001u /* initial image */
#define PHILOX_STREAM_FLAT 0x7FFF0002u /* flat-image chi2 benchmark */
#define PHILOX_DEFAULT_SEED 12345u

static inline uint32_t philox_mulhilo(uint32_t a, uint32_t b, uint32_t *hi) {
  uint64_t p = (uint64_t)a * (uint64_t)b;
  *hi = (uint32_t)(p >> 32);
  return (uint32_t)p;
}

static inline void philox_block(const uint32_t ctr_in[4], const uint32_t key_in[2], uint32_t out[4]) {
  uint32_t c0 = ctr_in[0], c1 = ctr_in[1], c2 = ctr_in[2], c3 = ctr_in[3];
  uint32_t k0 = key_in[0], k1 = key_in[1];
  for (int r = 0; r < 10; r++) { /* 10 rounds is the standard, Crush-resistant setting */
    uint32_t hi0, hi1;
    uint32_t lo0 = philox_mulhilo(0xD2511F53u, c0, &hi0);
    uint32_t lo1 = philox_mulhilo(0xCD9E8D57u, c2, &hi1);
    uint32_t n0 = hi1 ^ c1 ^ k0, n1 = lo1, n2 = hi0 ^ c3 ^ k1, n3 = lo0;
    c0 = n0;
    c1 = n1;
    c2 = n2;
    c3 = n3;
    k0 += 0x9E3779B9u;
    k1 += 0xBB67AE85u;
  }
  out[0] = c0;
  out[1] = c1;
  out[2] = c2;
  out[3] = c3;
}

/* A stream is fully determined by (stream_id, seed), which is the property we want here. */
static inline Philox philox_stream(uint32_t stream_id, uint32_t seed) {
  Philox g;
  g.ctr[0] = 0;
  g.ctr[1] = 0;
  g.ctr[2] = stream_id;
  g.ctr[3] = 0;
  g.key[0] = seed;
  g.key[1] = 0x1BD11BDAu;
  g.have = 0;
  return g;
}

static inline uint32_t philox_u32(Philox *g) {
  if (g->have == 0) {
    philox_block(g->ctr, g->key, g->buf);
    if (++g->ctr[0] == 0)
      ++g->ctr[1]; /* 2^66 draws per stream before wrap */
    g->have = 4;
  }
  return g->buf[4 - g->have--];
}

/* Uniform in (0,1), open at both ends so log(u) is always finite. */
static inline double philox_u01(Philox *g) { return ((double)philox_u32(g) + 0.5) * (1.0 / 4294967296.0); }

/* Uniform integer in [i,j]. Rejection, so there is no modulo bias. */
static inline int philox_randint(Philox *g, int i, int j) {
  uint32_t range = (uint32_t)(j - i) + 1u;
  if (range == 0u)
    return i; /* full 32-bit range requested */
  uint32_t lim = UINT32_MAX - (UINT32_MAX % range), v;
  do {
    v = philox_u32(g);
  } while (v >= lim);
  return i + (int)(v % range);
}

#endif
