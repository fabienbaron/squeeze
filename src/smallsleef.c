#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include "nonnumber.h"

#define PI4_A 0.78539816290140151978
#define PI4_B 4.9604678871439933374e-10
#define PI4_C 1.1258708853173288931e-18
#define PI4_D 1.7607799325916000908e-27

#define M_4_PI 1.273239544735162542821171882678754627704620361328125

#define L2U .69314718055966295651160180568695068359375
#define L2L .28235290563031577122588448175013436025525412068e-12
#define R_LN2 1.442695040888963407359924681001892137426645954152985934135449406931

static inline int64_t doubleToRawLongBits(double d) {
  union {
    double f;
    int64_t i;
  } tmp;
  tmp.f = d;
  return tmp.i;
}

static inline double longBitsToDouble(int64_t i) {
  union {
    double f;
    int64_t i;
  } tmp;
  tmp.i = i;
  return tmp.f;
}

static inline double xfabs(double x) {
  return longBitsToDouble(0x7fffffffffffffffLL & doubleToRawLongBits(x));
}

static inline double mulsign(double x, double y) {
  return longBitsToDouble(doubleToRawLongBits(x) ^ (doubleToRawLongBits(y) & (1LL << 63)));
}

static inline double sign(double d) { return mulsign(1, d); }
static inline double mla(double x, double y, double z) { return x * y + z; }
static inline double xrint(double x) { return x < 0 ? (int)(x - 0.5) : (int)(x + 0.5); }

static inline int xisnan(double x) { return x != x; }
static inline int xisinf(double x) { return x == INFINITY || x == -INFINITY; }
static inline int xisminf(double x) { return x == -INFINITY; }
static inline int xispinf(double x) { return x == INFINITY; }

static inline double pow2i(int q) {
  return longBitsToDouble(((int64_t)(q + 0x3ff)) << 52);
}

static inline double ldexpk(double x, int q) {
  double u;
  int m;
  m = q >> 31;
  m = (((m + q) >> 9) - m) << 7;
  q = q - (m << 2);
  m += 0x3ff;
  m = m < 0     ? 0     : m;
  m = m > 0x7ff ? 0x7ff : m;
  u = longBitsToDouble(((int64_t)m) << 52);
  x = x * u * u * u * u;
  u = longBitsToDouble(((int64_t)(q + 0x3ff)) << 52);
  return x * u;
}

static inline int ilogbp1(double d) {
  int m = d < 4.9090934652977266E-91;
  d = m ? 2.037035976334486E90 * d : d;
  int q = (doubleToRawLongBits(d) >> 52) & 0x7ff;
  q = m ? q - (300 + 0x03fe) : q - 0x03fe;
  return q;
}

typedef struct {
  double x, y;
} double2;

#ifndef NDEBUG
static int checkfp(double x) {
  if (xisinf(x) || xisnan(x)) return 1;
  return 0;
}
#endif

static inline double upper(double d) {
  return longBitsToDouble(doubleToRawLongBits(d) & 0xfffffffff8000000LL);
}

static inline double2 dd(double h, double l) {
  double2 ret;
  ret.x = h; ret.y = l;
  return ret;
}

static inline double2 ddnormalize_d2_d2(double2 t) {
  double2 s;

  s.x = t.x + t.y;
  s.y = t.x - s.x + t.y;

  return s;
}

static inline double2 ddscale_d2_d2_d(double2 d, double s) {
  double2 r;

  r.x = d.x * s;
  r.y = d.y * s;

  return r;
}

static inline double2 ddneg_d2_d2(double2 d) {
  double2 r;

  r.x = -d.x;
  r.y = -d.y;

  return r;
}

static inline double2 ddadd_d2_d_d(double x, double y) {
  // |x| >= |y|

  double2 r;

#ifndef NDEBUG
  if (!(checkfp(x) || checkfp(y) || xfabs(x) >= xfabs(y))) fprintf(stderr, "[ddadd_d2_d_d : %g, %g]", x, y);
#endif

  r.x = x + y;
  r.y = x - r.x + y;

  return r;
}

static inline double2 ddadd2_d2_d_d(double x, double y) {
  double2 r;

  r.x = x + y;
  double v = r.x - x;
  r.y = (x - (r.x - v)) + (y - v);

  return r;
}

static inline double2 ddadd_d2_d2_d(double2 x, double y) {
  // |x| >= |y|

  double2 r;

#ifndef NDEBUG
  if (!(checkfp(x.x) || checkfp(y) || xfabs(x.x) >= xfabs(y))) fprintf(stderr, "[ddadd_d2_d2_d : %g %g]", x.x, y);
#endif

  r.x = x.x + y;
  r.y = x.x - r.x + y + x.y;

  return r;
}

static inline double2 ddadd2_d2_d2_d(double2 x, double y) {
  // |x| >= |y|

  double2 r;

  r.x  = x.x + y;
  double v = r.x - x.x;
  r.y = (x.x - (r.x - v)) + (y - v);
  r.y += x.y;

  return r;
}

static inline double2 ddadd_d2_d_d2(double x, double2 y) {
  // |x| >= |y|

  double2 r;

#ifndef NDEBUG
  if (!(checkfp(x) || checkfp(y.x) || xfabs(x) >= xfabs(y.x))) fprintf(stderr, "[ddadd_d2_d_d2 : %g %g]", x, y.x);
#endif

  r.x = x + y.x;
  r.y = x - r.x + y.x + y.y;

  return r;
}

static inline double2 ddadd2_d2_d_d2(double x, double2 y) {
  double2 r;

  r.x  = x + y.x;
  double v = r.x - x;
  r.y = (x - (r.x - v)) + (y.x - v) + y.y;

  return r;
}

static inline double2 ddadd_d2_d2_d2(double2 x, double2 y) {
  // |x| >= |y|

  double2 r;

#ifndef NDEBUG
  if (!(checkfp(x.x) || checkfp(y.x) || xfabs(x.x) >= xfabs(y.x))) fprintf(stderr, "[ddadd_d2_d2_d2 : %g %g]", x.x, y.x);
#endif

  r.x = x.x + y.x;
  r.y = x.x - r.x + y.x + x.y + y.y;

  return r;
}

static inline double2 ddadd2_d2_d2_d2(double2 x, double2 y) {
  double2 r;

  r.x  = x.x + y.x;
  double v = r.x - x.x;
  r.y = (x.x - (r.x - v)) + (y.x - v);
  r.y += x.y + y.y;

  return r;
}

static inline double2 ddsub_d2_d2_d2(double2 x, double2 y) {
  // |x| >= |y|

  double2 r;

#ifndef NDEBUG
  if (!(checkfp(x.x) || checkfp(y.x) || xfabs(x.x) >= xfabs(y.x))) fprintf(stderr, "[ddsub_d2_d2_d2 : %g %g]", x.x, y.x);
#endif

  r.x = x.x - y.x;
  r.y = x.x - r.x - y.x + x.y - y.y;

  return r;
}

static inline double2 dddiv_d2_d2_d2(double2 n, double2 d) {
  double t = 1.0 / d.x;
  double dh  = upper(d.x), dl  = d.x - dh;
  double th  = upper(t  ), tl  = t   - th;
  double nhh = upper(n.x), nhl = n.x - nhh;

  double2 q;

  q.x = n.x * t;

  double u = -q.x + nhh * th + nhh * tl + nhl * th + nhl * tl +
    q.x * (1 - dh * th - dh * tl - dl * th - dl * tl);

  q.y = t * (n.y - q.x * d.y) + u;

  return q;
}

static inline double2 ddmul_d2_d_d(double x, double y) {
  double xh = upper(x), xl = x - xh;
  double yh = upper(y), yl = y - yh;
  double2 r;

  r.x = x * y;
  r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl;

  return r;
}

static inline double2 ddmul_d2_d2_d(double2 x, double y) {
  double xh = upper(x.x), xl = x.x - xh;
  double yh = upper(y  ), yl = y   - yh;
  double2 r;

  r.x = x.x * y;
  r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.y * y;

  return r;
}

static inline double2 ddmul_d2_d2_d2(double2 x, double2 y) {
  double xh = upper(x.x), xl = x.x - xh;
  double yh = upper(y.x), yl = y.x - yh;
  double2 r;

  r.x = x.x * y.x;
  r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.x * y.y + x.y * y.x;

  return r;
}

static inline double2 ddsqu_d2_d2(double2 x) {
  double xh = upper(x.x), xl = x.x - xh;
  double2 r;

  r.x = x.x * x.x;
  r.y = xh * xh - r.x + (xh + xh) * xl + xl * xl + x.x * (x.y + x.y);

  return r;
}

static inline double2 ddrec_d2_d(double d) {
  double t = 1.0 / d;
  double dh = upper(d), dl = d - dh;
  double th = upper(t), tl = t - th;
  double2 q;

  q.x = t;
  q.y = t * (1 - dh * th - dh * tl - dl * th - dl * tl);

  return q;
}

static inline double2 ddrec_d2_d2(double2 d) {
  double t = 1.0 / d.x;
  double dh = upper(d.x), dl = d.x - dh;
  double th = upper(t  ), tl = t   - th;
  double2 q;

  q.x = t;
  q.y = t * (1 - dh * th - dh * tl - dl * th - dl * tl - d.y * t);

  return q;
}

static inline double2 ddsqrt_d2_d2(double2 d) {
  double t = sqrt(d.x + d.y);
  return ddscale_d2_d2_d(ddmul_d2_d2_d2(ddadd2_d2_d2_d2(d, ddmul_d2_d_d(t, t)), ddrec_d2_d(t)), 0.5);
}

//

static inline double atan2k(double y, double x) {
  double s, t, u;
  int q = 0;

  if (x < 0) { x = -x; q = -2; }
  if (y > x) { t = x; x = y; y = -t; q += 1; }

  s = y / x;
  t = s * s;

  u = -1.88796008463073496563746e-05;
  u = u * t + (0.000209850076645816976906797);
  u = u * t + (-0.00110611831486672482563471);
  u = u * t + (0.00370026744188713119232403);
  u = u * t + (-0.00889896195887655491740809);
  u = u * t + (0.016599329773529201970117);
  u = u * t + (-0.0254517624932312641616861);
  u = u * t + (0.0337852580001353069993897);
  u = u * t + (-0.0407629191276836500001934);
  u = u * t + (0.0466667150077840625632675);
  u = u * t + (-0.0523674852303482457616113);
  u = u * t + (0.0587666392926673580854313);
  u = u * t + (-0.0666573579361080525984562);
  u = u * t + (0.0769219538311769618355029);
  u = u * t + (-0.090908995008245008229153);
  u = u * t + (0.111111105648261418443745);
  u = u * t + (-0.14285714266771329383765);
  u = u * t + (0.199999999996591265594148);
  u = u * t + (-0.333333333333311110369124);

  t = u * t * s + s;
  t = q * (M_PI/2) + t;

  return t;
}

double xatan2(double y, double x) {
  double r = atan2k(xfabs(y), x);

  r = mulsign(r, x);
  if (xisinf(x) || x == 0) r = M_PI/2 - (xisinf(x) ? (sign(x) * (M_PI  /2)) : 0);
  if (xisinf(y)          ) r = M_PI/2 - (xisinf(x) ? (sign(x) * (M_PI*1/4)) : 0);
  if (             y == 0) r = (sign(x) == -1 ? M_PI : 0);

  return xisnan(x) || xisnan(y) ? NAN : mulsign(r, y);
}

static double2 atan2k_u1(double2 y, double2 x) {
  double u;
  double2 s, t;
  int q = 0;

  if (x.x < 0) { x.x = -x.x; x.y = -x.y; q = -2; }
  if (y.x > x.x) { t = x; x = y; y.x = -t.x; y.y = -t.y; q += 1; }

  s = dddiv_d2_d2_d2(y, x);
  t = ddsqu_d2_d2(s);
  t = ddnormalize_d2_d2(t);

  u = 1.06298484191448746607415e-05;
  u = mla(u, t.x, -0.000125620649967286867384336);
  u = mla(u, t.x, 0.00070557664296393412389774);
  u = mla(u, t.x, -0.00251865614498713360352999);
  u = mla(u, t.x, 0.00646262899036991172313504);
  u = mla(u, t.x, -0.0128281333663399031014274);
  u = mla(u, t.x, 0.0208024799924145797902497);
  u = mla(u, t.x, -0.0289002344784740315686289);
  u = mla(u, t.x, 0.0359785005035104590853656);
  u = mla(u, t.x, -0.041848579703592507506027);
  u = mla(u, t.x, 0.0470843011653283988193763);
  u = mla(u, t.x, -0.0524914210588448421068719);
  u = mla(u, t.x, 0.0587946590969581003860434);
  u = mla(u, t.x, -0.0666620884778795497194182);
  u = mla(u, t.x, 0.0769225330296203768654095);
  u = mla(u, t.x, -0.0909090442773387574781907);
  u = mla(u, t.x, 0.111111108376896236538123);
  u = mla(u, t.x, -0.142857142756268568062339);
  u = mla(u, t.x, 0.199999999997977351284817);
  u = mla(u, t.x, -0.333333333333317605173818);

  t = ddmul_d2_d2_d(t, u);
  t = ddmul_d2_d2_d2(s, ddadd_d2_d_d2(1, t));
  t = ddadd2_d2_d2_d2(ddmul_d2_d2_d(dd(1.570796326794896557998982, 6.12323399573676603586882e-17), q), t);

  return t;
}

double xatan2_u1(double y, double x) {
  double2 d = atan2k_u1(dd(xfabs(y), 0), dd(x, 0));
  double r = d.x + d.y;

  r = mulsign(r, x);
  if (xisinf(x) || x == 0) r = M_PI/2 - (xisinf(x) ? (sign(x) * (M_PI  /2)) : 0);
  if (xisinf(y)          ) r = M_PI/2 - (xisinf(x) ? (sign(x) * (M_PI*1/4)) : 0);
  if (             y == 0) r = (sign(x) == -1 ? M_PI : 0);

  return xisnan(x) || xisnan(y) ? NAN : mulsign(r, y);
}

