#include "mini-gmp.h"

#ifndef ECDH
#define ECDH
typedef unsigned char boolean;

typedef struct EPoint {
    mpz_t x;
    mpz_t y;
} EPoint;

void init_point(EPoint* P);
void clear_point(EPoint** P);
void copy_point(EPoint* P, EPoint* Q);
boolean point_equal(EPoint* P, EPoint* Q);
int mpz_legendre(const mpz_t n, const mpz_t p);
void fit(EPoint* R, mpz_t a, mpz_t b, mpz_t p);
void generate_secret_key(mpz_t a, const mpz_t p);
int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p);
void to_weierstrass(EPoint* G, mpz_t p, mpz_t a, mpz_t b);
void double_and_add(EPoint* R, EPoint* P, mpz_t n, mpz_t p, mpz_t a, mpz_t b);
void point_addition(EPoint* R, EPoint* P, EPoint* Q, mpz_t p, mpz_t a, mpz_t b);
#endif
