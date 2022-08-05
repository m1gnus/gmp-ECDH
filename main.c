#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mini-gmp.h"

typedef unsigned char boolean;

typedef struct EPoint {
    mpz_t x;
    mpz_t y;
} EPoint;

void init_point(EPoint* P) {
    mpz_init2(P->x, 512); mpz_init2(P->y, 512);
}

void clear_point(EPoint** P) {
    mpz_clear((*P)->x);
    mpz_clear((*P)->y);
    free(*P); (*P) = NULL; }

void copy_point(EPoint* P, EPoint* Q) {
    mpz_set(P->x, Q->x);
    mpz_set(P->y, Q->y);
}

boolean point_equal(EPoint* P, EPoint* Q) {
    return (mpz_cmp(P->x, Q->x) == 0) && (mpz_cmp(P->y, Q->y) == 0);
}

void to_weierstrass(EPoint* G, mpz_t p, mpz_t a, mpz_t b) {
    mpz_t tmp;
    mpz_init2(tmp, 512);

    /* x = x/b + A/(3*b) */
    mpz_invert(tmp, b, p);
    mpz_mul(G->x, G->x, tmp);
    mpz_mod(G->x, G->x, p);
    mpz_mul_ui(tmp, b, 3);
    mpz_invert(tmp, tmp, p);
    mpz_mul(tmp, tmp, a);
    mpz_mod(tmp, tmp, p);
    mpz_add(G->x, G->x, tmp);
    mpz_mod(G->x, G->x, p);

    /* y = y/b */
    mpz_invert(tmp, b, p);
    mpz_mul(G->y, G->y, tmp);
    mpz_mod(G->y, G->y, p);
}

void point_addition(EPoint* R, EPoint* P, EPoint* Q, mpz_t p, mpz_t a, mpz_t b) {
    mpz_t l1, l2, lambda, a4, B;
    mpz_init2(l1, 512);
    mpz_init2(l2, 512);
    mpz_init2(lambda, 512);
    mpz_init2(a4, 512);
    mpz_init2(B, 512);

    if (point_equal(P, Q)) {
        /* a4 = (3-A^2)/(3B^2) */
        mpz_mul(a4, a, a);
        mpz_mul_si(a4, a4, -1);
        mpz_add_ui(a4, a4, 3);
        mpz_mod(a4, a4, p);

        mpz_mul(B, b, b);
        mpz_mul_ui(B, B, 3);
        mpz_invert(B, B, p);

        mpz_mul(a4, a4, B);
        mpz_mod(a4, a4, p);

        /* ( 3*(x_p^2) + a ) * (2 * y_p)^-1 mod p */
        mpz_mul(l1, P->x, P->x);
        mpz_mul_ui(l1, l1, 3);
        mpz_add(l1, l1, a4);

        mpz_mul_ui(l2, P->y, 2);
        mpz_invert(l2, l2, p);
    }
    else {
        /* y_q - y_p */
        mpz_sub(l1, Q->y, P->y);

        /* (x_q - x_p)^(-1) */
        mpz_sub(l2, Q->x, P->x);
        mpz_invert(l2, l2, p);

    }
    /* l = (l1 * l2) mod p */
    mpz_mul(lambda, l1, l2);
    mpz_mod(lambda, lambda, p);

    /* l^2 */
    mpz_mul(l2, lambda, lambda);
    mpz_mod(l2, l2, p);

    /* x_r */
    mpz_sub(R->x, l2, P->x);
    mpz_sub(R->x, R->x, Q->x);
    mpz_mod(R->x, R->x, p);

    /* y_r */
    mpz_sub(R->y, P->x, R->x);
    mpz_mul(R->y, R->y, lambda);
    mpz_sub(R->y, R->y, P->y);
    mpz_mod(R->y, R->y, p);

    mpz_clear(l1);
    mpz_clear(l2);
    mpz_clear(a4);
    mpz_clear(B);
    mpz_clear(lambda);
}

void double_and_add(EPoint* R, EPoint* P, mpz_t n, mpz_t p, mpz_t a, mpz_t b) {
    mpz_t odd;
    mpz_init(odd);

    EPoint *appo = malloc(sizeof *appo), *appo2 = malloc(sizeof *appo2);
    init_point(appo);
    init_point(appo2);

    copy_point(appo2, P);
    mpz_set_si(R->x, -1);
    mpz_set_ui(R->y, 0);

    while (mpz_cmp_ui(n, 0) > 0) {
        mpz_mod_ui(odd, n, 2);

        if (mpz_cmp_ui(odd, 1) == 0) {
            if (mpz_cmp_si(R->x, 0) < 0) {
                copy_point(R, appo2);
            }
            else{
                copy_point(appo, R);
                point_addition(R, appo, appo2, p, a, b);
            }
            mpz_sub_ui(n, n, 1);
        }
        copy_point(appo, appo2);
        point_addition(appo2, appo, appo, p, a, b);
        mpz_divexact_ui(n, n, 2);
    }

    mpz_clear(odd);
    clear_point(&appo);
    clear_point(&appo2);
}

void generate_secret_key(mpz_t a, const mpz_t p) {
    /* Genera array di 32 byte e importa finchè non ottiene un numero minore di p */
    unsigned char random_array[32];
    mpz_set(a, p);

    while (mpz_cmp(a, p) >= 0) {
        for (unsigned char i=0; i<32; i++)
            random_array[i] = rand();
        mpz_import(a, 32, 1, sizeof random_array[0], 0, 0, random_array);
        printf("\r\nsecret: ");
        mpz_out_str(stdout, 10, a);
        printf("\r\n");
    }
}

int mpz_legendre(const mpz_t n, const mpz_t p) {
    mpz_t tmp;
    mpz_init2(tmp, 512);

    mpz_sub_ui(tmp, p, 1);
    mpz_divexact_ui(tmp, tmp, 2);
    mpz_powm(tmp, n, tmp, p);

    return (mpz_cmp_ui(tmp, 1) == 0) ? 1 : -1;
}

int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p) {
    unsigned int i, s;

      if (mpz_legendre(n, p) < 0) {
        return 0;
      }

      if(mpz_divisible_p(n, p)) {
          mpz_set_ui(q, 0);
          return 1;
      }

      if(mpz_tstbit(p, 1) == 1) {         /* p = 3 (mod 4) ?                  */
          mpz_set(q, p);
          mpz_add_ui(q, q, 1);
          mpz_fdiv_q_2exp(q, q, 2);
          mpz_powm(q, n, q, p);           /* q = n ^ ((p 1) / 4) (mod p)      */
          return 1;
      }

      unsigned int S=0;
      mpz_t R, M, c, z, t, Q, exp, tmp, b;
      mpz_init2(Q, 512);
      mpz_init2(R, 512);
      mpz_init2(M, 512);
      mpz_init2(c, 512);
      mpz_init2(z, 512);
      mpz_init2(t, 512);
      mpz_init2(exp, 512);
      mpz_init2(tmp, 512);
      mpz_init2(b, 512);

      mpz_set(Q, p);
      mpz_sub_ui(Q, Q, 1);
      while (mpz_tstbit(Q, S) == 0) {
          S++;
      }
      mpz_fdiv_q_2exp(Q, Q, S);
      mpz_set_ui(M, S);

      mpz_set_ui(z, 2);
      while(mpz_legendre(z, p) == 1) mpz_add_ui(z, z, 1);
      mpz_powm(c, z, Q, p);
      mpz_powm(t, n, Q, p);
      mpz_set(exp, Q);
      mpz_add_ui(exp, exp, 1);
      mpz_divexact_ui(exp, exp, 2);
      mpz_powm(R, n, exp, p);

      while (1) {
        if (!mpz_cmp_ui(t, 0)) mpz_set_ui(q, 0);
        else if (!mpz_cmp_ui(t, 1)) {
            mpz_set(q, R);
            return 1;
        }

        mpz_set(tmp, t);
        unsigned int i;
        for (i=0; i<mpz_get_ui(M); i++) {
            mpz_powm_ui(tmp, tmp, 2, p);
            if (!mpz_cmp_ui(tmp, 1)) break;
        }

        if (!mpz_cmp_ui(M, i)) return 1;

        mpz_sub_ui(exp, M, i);
        mpz_sub_ui(exp, exp, 1);
        mpz_set_ui(b, 2);
        mpz_set(tmp, p);
        printf("\r\nEXP (%d): ", i);
        mpz_out_str(stdout, 10, M);
        printf("\r\n");
        mpz_sub_ui(tmp, tmp, 1);
        printf("\r\nEXP: ");
        mpz_out_str(stdout, 10, exp);
        printf("\r\n");
        mpz_powm(exp, b, exp, tmp);
        mpz_powm(b, c, exp, p);
        mpz_set_ui(M, i);
        mpz_powm_ui(c, b, 2, p);
        mpz_mul(t, t, c);
        mpz_mod(t, t, p);
        mpz_mul(R, R, b);
        mpz_mod(R, R, p);
      }

      return 0;
}

void fit(EPoint* R, mpz_t a, mpz_t b, mpz_t p) {
    mpz_t a4, B, tmp, tmp2;
    mpz_init2(a4, 512);
    mpz_init2(B, 512);
    mpz_init2(tmp, 512);
    mpz_init2(tmp2, 512);


    /* a4 = (3-A^2)/(3B^2) */
    mpz_mul(a4, a, a);
    mpz_mul_si(a4, a4, -1);
    mpz_add_ui(a4, a4, 3);
    mpz_mod(a4, a4, p);

    mpz_mul(B, b, b);
    mpz_mul_ui(B, B, 3);
    mpz_invert(B, B, p);

    mpz_mul(a4, a4, B);
    mpz_mod(a4, a4, p);

    /* B = (2*A^3 - 9*A)/(27*B^3) */
    mpz_powm_ui(B, a, 3, p);
    mpz_mul_ui(B, B, 2);
    mpz_mod(B, B, p);
    mpz_mul_ui(tmp, a, 9);
    mpz_sub(B, B, tmp);
    mpz_mod(B, B, p);
    mpz_powm_ui(tmp, b, 3, p);
    mpz_mul_ui(tmp, tmp, 27);
    mpz_invert(tmp, tmp, p);
    mpz_mul(B, B, tmp);
    mpz_mod(B, B, p);

    mpz_powm_ui(tmp, R->x, 3, p);
    mpz_mul(tmp2, R->x, a4);
    mpz_mod(tmp2, tmp2, p);
    mpz_add(tmp, tmp, tmp2);
    mpz_add(tmp, tmp, B);
    mpz_mod(tmp, tmp, p);

    if (!mpz_sqrtm(R->y, tmp, p)) {
        printf("\r\n\r\nTMP: ");
        mpz_out_str(stdout, 10, tmp);
        printf("\r\n");
        mpz_out_str(stdout, 10, R->y);
        mpz_mul_si(tmp, tmp, -1);
        mpz_mod(tmp, tmp, p);
        if (!mpz_sqrtm(R->y, tmp, p))
            printf("\r\n\r\nTMP: ");
            mpz_out_str(stdout, 10, tmp);
            printf("\r\n");
            mpz_out_str(stdout, 10, R->y);
            perror("fit");
            exit(1);
    }
}

int main() {
    srand(time(0));

    /* public module p=2^255-19 */
    mpz_t p;
    mpz_init2(p, 512);
    mpz_set_ui(p, 2);
    mpz_pow_ui(p, p, 255);
    mpz_sub_ui(p, p, 19);

    mpz_t a,b;
    mpz_init2(a, 512);
    mpz_init2(b, 512);
    mpz_set_ui(a, 0x76d06);
    mpz_set_ui(b, 0x01);

    /* generator */
    EPoint* G = malloc(sizeof *G);
    init_point(G);
    mpz_set_ui(G->x, 0x09);
    mpz_set_str(G->y, "20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);

    EPoint* R = malloc(sizeof *R);
    init_point(R);

    to_weierstrass(G, p, a, b);
    printf("G.x: ");
    mpz_out_str(stdout, 10, G->x);
    printf("\r\nG.y: ");
    mpz_out_str(stdout, 10, G->y);

    printf("\r\n\r\nFitting the point:");
    fit(G, a, b, p);
    printf("\r\nG.x: ");
    mpz_out_str(stdout, 10, G->x);
    printf("\r\nG.y: ");
    mpz_out_str(stdout, 10, G->y);
    printf("\r\n\r\n");

    point_addition(R, G, G, p, a, b);
    printf("\r\nR.x: ");
    mpz_out_str(stdout, 10, R->x);
    printf("\r\nR.y: ");
    mpz_out_str(stdout, 10, R->y);

    EPoint* P = malloc(sizeof *G);
    init_point(P);
    mpz_set_str(P->x, "42413865998592550293572449399763517288723132413484395897391043645513011644945", 10);
    mpz_set_str(P->y, "21996169931478149623282952155306601567177220461061887349857934195604197012838", 10);

    EPoint* Q = malloc(sizeof *G);
    init_point(Q);
    mpz_set_str(Q->x, "38784272806674503386194598642122606159134019533217578547625273097149674728900", 10);
    mpz_set_str(Q->y, "6403266695435485875536658491776352978668566976402459773937556907604509847680", 10);

    point_addition(R, P, Q, p, a, b);
    printf("\r\nR.x: ");
    mpz_out_str(stdout, 10, R->x);
    printf("\r\nR.y: ");
    mpz_out_str(stdout, 10, R->y);

    mpz_t exp;
    mpz_init2(exp, 512);
    mpz_set_ui(exp, 1111);
    printf("\r\n\r\nTesting double and add:");
    double_and_add(R, G, exp, p, a, b);
    printf("\r\nR.x: ");
    mpz_out_str(stdout, 10, R->x);
    printf("\r\nR.y: ");
    mpz_out_str(stdout, 10, R->y);

    mpz_t alice, bob;
    mpz_init2(alice, 512);
    mpz_init2(bob, 512);

    generate_secret_key(alice, p);
    generate_secret_key(bob, p);

    EPoint *PA = malloc(sizeof *PA), *PB = malloc (sizeof *PB);
    init_point(PA);
    init_point(PB);

    double_and_add(PA, G, alice, p, a, b);
    double_and_add(PB, G, bob, p, a, b);

    EPoint *SA = malloc(sizeof *SA), *SB = malloc(sizeof *SB);
    init_point(SA);
    init_point(SB);

    double_and_add(SA, PB, alice, p, a, b);
    double_and_add(SB, PA, bob, p, a, b);

    printf("\r\nSA.x: ");
    mpz_out_str(stdout, 10, R->x);
    printf("\r\nSA.y: ");
    mpz_out_str(stdout, 10, R->y);

    printf("\r\n\r\nSB.x: ");
    mpz_out_str(stdout, 10, R->x);
    printf("\r\nSB.y: ");
    mpz_out_str(stdout, 10, R->y);

    mpz_clear(p);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(exp);
    clear_point(&G);
    clear_point(&R);
    clear_point(&P);
    clear_point(&Q);
    clear_point(&PA);
    clear_point(&PB);
    clear_point(&SA);
    clear_point(&SB);

    return 0;
}
