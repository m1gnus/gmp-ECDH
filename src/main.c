#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mini-gmp.h"
#include "ecdh.h"


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
