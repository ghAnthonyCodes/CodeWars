#include <cstdint>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdlib>

using namespace std;

using i64 = int64_t;
using u64 = uint64_t;

u64 M;

class Mu64_t {
  
  public:
    u64 n;
    Mu64_t(u64 n) {
      this->n = n;
    };
    friend Mu64_t operator+(const Mu64_t &a, const Mu64_t &b) {
      return Mu64_t(((a.n % M) + (b.n % M)) % M);
    }
    friend Mu64_t operator+(const Mu64_t &a, const u64 &b) {
      return Mu64_t(((a.n % M) + (b % M)) % M);
    }
    friend Mu64_t operator-(const Mu64_t &a, const u64 &b) {
      return Mu64_t(((a.n % M) - (b % M) + M) % M);
    }
    friend Mu64_t operator-(const Mu64_t &a, const Mu64_t &b) {
      return Mu64_t(((a.n % M) - (b.n % M) + M) % M);
    }
    friend Mu64_t operator*(const Mu64_t &a, const Mu64_t &b) {
      return Mu64_t(((a.n % M) * (b.n % M)) % M);
    }
};

/******************************************************************************/
/*                                                                            */
/* solve_square_mod                                                           */
/*                                                                            */
/******************************************************************************/
u64 solve_square_mod(u64 n, u64 l) {

  Mu64_t x = Mu64_t(n);
  Mu64_t y = Mu64_t(l);
  Mu64_t z = Mu64_t((l + 1) / 2);
  Mu64_t q = Mu64_t(l / 2);
  Mu64_t k = Mu64_t(n / 2);
  
  if (l >= n - 1) 
    return 0;

  if (l == 0)
    return (k * x * (x - 1)).n;

  if (l & 1)
    return (x * (k * (x - 1) - (x - 1 - y) * y - y * z)).n;
  
  return (x * (k * (x - 1) - (x - 1 - y) * y - q * (y + 1))).n;
}

/******************************************************************************/
/*                                                                            */
/* solve_side_lobe_mod                                                        */
/*                                                                            */
/******************************************************************************/
u64 solve_side_lobe_mod(u64 a, u64 n, u64 l) {

  Mu64_t x = Mu64_t(n);
  Mu64_t y = Mu64_t(l);
  Mu64_t z = Mu64_t(a);
  Mu64_t k = Mu64_t(n / 2);
  Mu64_t s = Mu64_t(2 * n - l);
  Mu64_t q = Mu64_t((2 * n - l - 1) / 2);
  Mu64_t r = Mu64_t((2 * n - l) / 2);
  Mu64_t g = Mu64_t(2 * n - l - 1);

  if (l >= 2 * n - 1)
    return 0;

  if (l <= n)
    return (z * (2 * k * (2 * x - 1) - k * (x - 1)) - y * z * x).n;

  // If l is odd, then 2n - l - 1 is even
  if (l & 1)
    return (z * s * q).n;
  
  // n < l < 2n - 1
  return (z * r * g).n;
}

/******************************************************************************/
/*                                                                            */
/* dive                                                                       */
/*                                                                            */
/******************************************************************************/
u64 dive(u64 r, u64 c, u64 l) {
  
  // Start building the answer
  u64 answer = 0;

  // Smallest case
  if (r == 1 && c == 1)
    return 0;
  
  // Swap if row > c, to reduce logic (the problem is symmetric)
  if (r > c) 
    swap(r, c);

  // If only columns remain
  if (r == 1) {
    if (l >= c - 1)
      return 0;  
    if ((c - l) % 2 == 0)
      return (Mu64_t((c - l) / 2) * Mu64_t(c - l - 1)).n;
    return (Mu64_t(c - l) * Mu64_t((c - l - 1) / 2)).n;
  }
    
  // Largest subsquare side length
  u64 n = pow(2, (u64)log2(r));

  // Solve main square
  answer = (answer + solve_square_mod(n, l)) % M;

  // Solve right lobe up to maximum of next power of 2 (2n) 
  if (c > n)
    answer = (answer + solve_side_lobe_mod(min(n, c - n), n, l)) % M;

  // Solve bottom lobe down to the next power of 2 (2n)
  if (r > n)
    answer = (answer + solve_side_lobe_mod(min(n, r - n), n, l)) % M;

  // Solve bottom right corner of main square
  if (c > n && r > n)
    answer = (answer + dive(min(n, r - n), min(n, c - n), l)) % M;

  // There can only be, at max, 1 side lobe remaining
  // Further, we no longer need to worry about the corner nonesense
  if (c <= r) 
    return answer;

  // Next power of 2
  n *= 2;
  while (c > 2 * n - 1) {
    answer = (answer + solve_side_lobe(r, n, l)) % M;
    n *= 2;
  }

  // Exit conditions
  if (c <= n || l >= 2 * n -1) 
    return answer;

  // Whatever is left from c, we will shift dive
  Mu64_t x = Mu64_t(n);
  Mu64_t y = Mu64_t(l);
  Mu64_t z = Mu64_t(r);
  Mu64_t q = Mu64_t(c);
  
  if (l <= n) {
    answer = (answer + ((x - y) * z * (q - x)).n) % M;
    return (answer + dive(r, c - n, 0)) % M;
  }
  
  // n < l < 2n-1
  return (answer + dive(r, c - n, l - n)) % M;
}

/******************************************************************************/
/*                                                                            */
/* elder_age                                                                  */
/*                                                                            */
/******************************************************************************/
u64 elder_age(u64 m, u64 n, u64 l, u64 t) {
  M = t;
  return dive(m, n, l);
}
