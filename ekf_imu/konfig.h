#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>



/* State Space dimension */
#define SS_X_LEN    (4)
#define SS_Z_LEN    (6)
// #define SS_Z_LEN    (3)
#define SS_U_LEN    (3)
#define SS_DT_MILIS (10)                            /* 10 ms */
#define SS_DT       float_prec(SS_DT_MILIS/1000.)   /* Sampling time */

#define MATRIX_MAXIMUM_SIZE     (10)        /* Ukuran maksimum matrix yang diperbolehkan */

#define PAKAI_FLOAT         1
#define PAKAI_DOUBLE        2
#define PRESISI_FPU         (PAKAI_DOUBLE)

#if (PRESISI_FPU == PAKAI_FLOAT)
    #define float_prec      float
    #define float_prec_ZERO (1e-8)
#elif (PRESISI_FPU == PAKAI_DOUBLE)
    #define float_prec      double
    #define float_prec_ZERO (1e-15)
#else
    #error("Kepresisian FPU belum didefinisi!");
#endif

// -0.183381	-0.047919	-0.981873	
// -0.011493	-0.051165	0.031833	
//  0.980140	-0.180927	32.621292
#define ModelSistem_ACC_Z0      (-1)

#define ModelSistem_B0x_Awal    (1)//-0.081183)
#define ModelSistem_B0y_Awal    (1)//0.980140)
#define ModelSistem_B0z_Awal    (0.0)//(-0.180927)
// 42.441505	-48.060120	-2.635646	


/* Aktifkan definisi ini untuk implementasi kode di MCU (std lib tidak digunakan) */
#define SISTEM_PC                   1
#define SISTEM_EMBEDDED_NO_PRINT    2
#define SISTEM_EMBEDDED_ARDUINO     3

#define SISTEM_IMPLEMENTASI         SISTEM_EMBEDDED_ARDUINO

#endif // MAIN_H
