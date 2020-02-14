#ifndef _ytlu_randomgenerator_1_
#define _ytlu_randomgenerator_1_
#define RANDOM_RHO  0x7a35ef05UL  // rho = 8t-3, 8-byte unsigned long
#define RANDOM_MAX  ((double)0x40000000)  // maximan of congruence sequence
#define RANDOM_HALF ((double)0x20000000)
#define NRANDOM(nn) ((nn *= RANDOM_RHO) >> 2)
#define RRANDOM(nn) ( (((nn *= RANDOM_RHO) >> 2)-RANDOM_HALF) / RANDOM_HALF)
#endif
