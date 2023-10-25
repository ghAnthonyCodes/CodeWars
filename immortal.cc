#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cassert>
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
      return Mu64_t(((a.n%M) + (b.n%M))%M);
    }
    friend Mu64_t operator+(const Mu64_t &a, const u64 &b) {
      return Mu64_t(((a.n%M) + (b%M))%M);
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
      return (Mu64_t((c - l) / 2) * Mu64_t(c - l -1)).n;
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
  u64 answer = dive(m, n, l) % t;
  return answer;
}

int main() {
  for (int i=0; i<10000000; i++) {
    u64 right_lobe;
    u64 right_lobe_mod;
    u64 bottom_lobe;
    u64 bottom_lobe_mod;

    M = 1 + rand() % 10000;
    u64 r = 1 + rand() % 10000;
    u64 c = 1 + rand() % 10000;
    u64 l = 1 + rand() % 1000;

    // Largest subsquare side length
    u64 n = pow(2, (u64)log2(min(r, c)));
    // printf("Testing: square(%llu, %llu) %% %llu \n", n, l, M);
    
    u64 square = solve_square(n, l) % M;
    u64 square_mod = solve_square_mod(n, l);

    // Solve right lobe up to maximum of next power of 2 (2n) 
    if (c > n) {
      right_lobe = solve_side_lobe(min(n, c - n), n, l) % M;
      right_lobe_mod = solve_side_lobe_mod(min(n, c - n), n, l);
    }
    // printf("%llu == %llu?\n", right_lobe, right_lobe_mod);

    // Solve bottom lobe down to the next power of 2 (2n)
    if (r > n) {
      bottom_lobe = solve_side_lobe(min(n, r - n), n, l) % M;
      bottom_lobe_mod = solve_side_lobe_mod(min(n, r - n), n, l);
    }
    // printf("L: %llu ?= R: %llu\n", left, right);

    assert(square == square_mod);
    assert(right_lobe == right_lobe_mod);
    assert(bottom_lobe == bottom_lobe_mod);
  }
  printf("Everything passed!\n");

  bool enable_4x4 = true;
  bool enable_8x8 = true;
  bool enable_32x32 = true;
  bool enable_odds = true;

  // Test cases 4x4
  if (enable_4x4) {
    assert(elder_age(4, 4, 0, 1000) == 24);
    assert(elder_age(4, 4, 1, 1000) == 12);
    assert(elder_age(4, 4, 2, 1000) == 4);
    assert(elder_age(4, 4, 3, 1000) == 0);
    assert(elder_age(4, 4, 4, 1000) == 0);
    assert(elder_age(4, 4, 123431, 1000) == 0);
  }

  // Test cases 8x8
  if (enable_8x8) {
    assert(elder_age(8, 8, 0, 1000) == 224);
    assert(elder_age(8, 8, 1, 1000) == 168);
    assert(elder_age(8, 8, 2, 1000) == 120);
    assert(elder_age(8, 8, 3, 1000) == 80);
    assert(elder_age(8, 8, 4, 1000) == 48);
    assert(elder_age(8, 8, 5, 1000) == 24);
    assert(elder_age(8, 8, 6, 1000) == 8);
    assert(elder_age(8, 8, 7, 1000) == 0);
    assert(elder_age(8, 8, 8, 1000) == 0);
    assert(elder_age(8, 8, 123431, 1000) == 0);
  }
    
  // Tough cases 32x32
  if (enable_32x32) {
    assert(elder_age(32, 32, 0, 123) == 5);
    assert(elder_age(32, 32, 13, 10000) == 5472);
  }

  // Odd shape cases
  if (enable_odds) {
    assert(elder_age(11, 10, 0, 10000) == 691);
    assert(elder_age(11, 10, 1, 10000) == 591);
    assert(elder_age(11, 10, 2, 10000) == 501);
    assert(elder_age(11, 10, 3, 10000) == 420);
    assert(elder_age(11, 10, 4, 10000) == 348);
    assert(elder_age(11, 10, 5, 10000) == 284);
    assert(elder_age(11, 10, 6, 10000) == 228);
    assert(elder_age(11, 10, 7, 10000) == 180);
    assert(elder_age(11, 10, 8, 10000) == 140);
    assert(elder_age(11, 10, 9, 10000) == 105);
    assert(elder_age(11, 10, 14, 10000) == 5);
    assert(elder_age(11, 10, 15, 10000) == 0);
    assert(elder_age(11, 10, 16, 10000) == 0);

    assert(elder_age(10, 11, 0, 10000) == 691);
    assert(elder_age(10, 11, 1, 10000) == 591);
    assert(elder_age(10, 11, 2, 10000) == 501);
    assert(elder_age(10, 11, 3, 10000) == 420);
    assert(elder_age(10, 11, 4, 10000) == 348);
    assert(elder_age(10, 11, 5, 10000) == 284);
    assert(elder_age(10, 11, 6, 10000) == 228);
    assert(elder_age(10, 11, 7, 10000) == 180);
    assert(elder_age(10, 11, 8, 10000) == 140);
    assert(elder_age(10, 11, 9, 10000) == 105);
    assert(elder_age(10, 11, 14, 10000) == 5);
    assert(elder_age(10, 11, 15, 10000) == 0);
    assert(elder_age(10, 11, 16, 10000) == 0);
  }

  // Wide cases
  assert(elder_age(3, 8, 0, 10000) == 84);
  assert(elder_age(3, 9, 0, 10000) == 111);
  assert(elder_age(3, 9, 1, 1000) == 87);
  assert(elder_age(3, 9, 2, 1000) == 66);
  assert(elder_age(3, 9, 3, 1000) == 48);
  assert(elder_age(3, 9, 4, 1000) == 33);
  assert(elder_age(3, 9, 5, 1000) == 21);
  assert(elder_age(3, 9, 6, 1000) == 12);
  assert(elder_age(3, 9, 7, 1000) == 6);
  assert(elder_age(3, 9, 8, 1000) == 3);
  assert(elder_age(3, 9, 9, 1000) == 1);
  assert(elder_age(3, 9, 10, 1000) == 0);
  assert(elder_age(3, 9, 15, 1000) == 0);
  assert(elder_age(3, 9, 16, 1000) == 0);
  assert(elder_age(3, 9, 17, 1000) == 0);

  assert(elder_age(3, 15, 0, 10000) == 318);
  assert(elder_age(3, 15, 1, 1000) == 276);
  assert(elder_age(3, 15, 2, 1000) == 237);
  assert(elder_age(3, 15, 3, 1000) == 201);
  assert(elder_age(3, 15, 4, 1000) == 168);
  assert(elder_age(3, 15, 5, 1000) == 138);
  assert(elder_age(3, 15, 6, 1000) == 111);
  assert(elder_age(3, 15, 7, 1000) == 87);
  assert(elder_age(3, 15, 8, 1000) == 66);
  assert(elder_age(3, 15, 9, 1000) == 48);
  assert(elder_age(3, 15, 10, 1000) == 33);
  assert(elder_age(3, 15, 11, 1000) == 21);
  assert(elder_age(3, 15, 12, 1000) == 12);
  assert(elder_age(3, 15, 13, 1000) == 6);
  assert(elder_age(3, 15, 14, 1000) == 2);
  assert(elder_age(3, 15, 15, 1000) == 0);
  assert(elder_age(3, 15, 16, 1000) == 0);
  assert(elder_age(3, 15, 17, 1000) == 0);

  assert(elder_age(6, 36, 5, 10000) == 2822);
  assert(elder_age(6, 36, 11, 10000) == 1832);
  assert(elder_age(6, 36, 25, 10000) == 362);
  assert(elder_age(6, 36, 32, 10000) == 68);

  // Tall cases
  assert(elder_age(8, 3, 0, 10000) == 84);
  assert(elder_age(9, 3, 0, 10000) == 111);
  assert(elder_age(9, 3, 1, 1000) == 87);
  assert(elder_age(9, 3, 2, 1000) == 66);
  assert(elder_age(9, 3, 3, 1000) == 48);
  assert(elder_age(9, 3, 4, 1000) == 33);
  assert(elder_age(9, 3, 5, 1000) == 21);
  assert(elder_age(9, 3, 6, 1000) == 12);
  assert(elder_age(9, 3, 7, 1000) == 6);
  assert(elder_age(9, 3, 8, 1000) == 3);
  assert(elder_age(9, 3, 9, 1000) == 1);
  assert(elder_age(9, 3, 10, 1000) == 0);
  assert(elder_age(9, 3, 15, 1000) == 0);
  assert(elder_age(9, 3, 16, 1000) == 0);
  assert(elder_age(9, 3, 17, 1000) == 0);

  assert(elder_age(15, 3, 0, 10000) == 318);
  assert(elder_age(15, 3, 1, 1000) == 276);
  assert(elder_age(15, 3, 2, 1000) == 237);
  assert(elder_age(15, 3, 3, 1000) == 201);
  assert(elder_age(15, 3, 4, 1000) == 168);
  assert(elder_age(15, 3, 5, 1000) == 138);
  assert(elder_age(15, 3, 6, 1000) == 111);
  assert(elder_age(15, 3, 7, 1000) == 87);
  assert(elder_age(15, 3, 8, 1000) == 66);
  assert(elder_age(15, 3, 9, 1000) == 48);
  assert(elder_age(15, 3, 10, 1000) == 33);
  assert(elder_age(15, 3, 11, 1000) == 21);
  assert(elder_age(15, 3, 12, 1000) == 12);
  assert(elder_age(15, 3, 13, 1000) == 6);
  assert(elder_age(15, 3, 14, 1000) == 2);
  assert(elder_age(15, 3, 15, 1000) == 0);
  assert(elder_age(15, 3, 16, 1000) == 0);
  assert(elder_age(15, 3, 17, 1000) == 0);

  assert(elder_age(36, 6, 5, 10000) == 2822);
  assert(elder_age(36, 6, 11, 10000) == 1832);
  assert(elder_age(36, 6, 25, 10000) == 362);
  assert(elder_age(36, 6, 32, 10000) == 68);
  
  assert(elder_age(28827050410, 35165045587, 7109602, 13719506) == 5456283);

  // Success message
  printf("All test cases pass.\n");
  return 0;
}