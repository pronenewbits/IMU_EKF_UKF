/*********************************************************************************************
 *  Fungsi Continuous-time Extended Kalman Filter (yang dihitung secara diskrit)
 *    (Karena umumnya dinamakan Continuous-discrete EKF, maka implementasi ini bernama CD_EKF)
 *
 *
 *         Formulasi plant yang diestimasi:
 *              x*(t) = f[x(t), u(t)] + w(t)          ; x=Nx1,    u=Mx1             ...{CD_EKF_1}
 *              z(t)  = h[x(t), u(t)] + v(t)          ; z=Zx1                       ...{CD_EKF_2}
 *
 *        dimana:
 *          x*(t): Perubahan x(t) terhadap waktu ==> dx/dt   : Nx1
 *          x(t) : State Variable sistem waktu ke-t          : Nx1
 *          z(t) : Keluaran sistem waktu ke-t                : Zx1
 *          u(t) : Masukan sistem waktu ke-t                 : Mx1
 *          w(t) : Noise proses, asumsi AWGN, kovarian Q     : Nx1
 *          v(t) : Noise pengukuran, asumsi AWGN, kovarian R : Nx1
 *
 *          f(..), h(..) adalah model (non-linear) dari sistem yang ingin diprediksi
 *
 *********************************************************************************************
 *         Algoritma Continuous-Discrete Extended-kalman-filter:
 *          Init:
 *                  P (k=0|k=0) = Identitas * Kovarian(P(k=0)), biasanya diisi dg nilai yg besar
 *                  X~(k=0|k=0) = Ekspektasi(X(k=0)), biasanya = 0
 *          Lakukan iterasi:
 *              Hitung paramater linearisasi A:
 *                  A = d(f(..))/dx |x~(k-1|k-1),u(k) = Jacobian(f,x~(k-1|k-1),u(k))...{CD_EKF_3}
 *
 *              Prediksi:
 *                  X*~(k|k-1) = f[X~(k-1|k-1), u(k)]                               ...{CD_EKF_4}
 *                  P*(k|k-1)  = A*P(k-1|k-1) + P(k-1|k-1)*A' + Q                   ...{CD_EKF_5}
 *                  X~(k|k-1)  = X~(k-1|k-1) + X*~(k|k-1)*dt    (integral orde 1)   ...{CD_EKF_6}
 *                  P(k|k-1)   = P(k-1|k-1) + P*(k|k-1)*dt      (integral orde 1)   ...{CD_EKF_7}
 *
 *              Hitung paramater linearisasi C:
 *                  C = d(f(..))/dx |x~(k|k-1),u(k) = Jacobian(h,x~(k|k-1),u(k))    ...{CD_EKF_8}
 *
 *              Update:
 *                  e(k)    = z(k) - h[X~(k|k-1), u(k)]                             ...{CD_EKF_9}
 *                  S       = C*P(k|k-1)*C' + R                                     ...{CD_EKF_10}
 *                  K       = P(k|k-1)*C'*(S^-1)                                    ...{CD_EKF_11}
 *                  X~(k|k) = X~(k|k-1) + K*e(k)                                    ...{CD_EKF_12}
 *                  P(k|k)  = (I - K*C)*P(k|k-1)                                    ...{CD_EKF_13}
 *
 *        *Catatan tambahan:
 *              - Perhatikan persamaan {CD_EKF_5} BERBEDA dengan versi diskritnya!!!!
 *              - Dengan asumsi masukan plant ZOH, u(k) = u(k|k-1),
 *                  Dengan asumsi tambahan observer dijalankan sebelum pengendali, u(k|k-1) = u(k-1),
 *                  sehingga u(k) [untuk perhitungan kalman] adalah nilai u(k-1) [dari pengendali].
 *              - Notasi yang benar adalah u(k|k-1), tapi disini menggunakan notasi u(k) untuk
 *                  menyederhanakan penulisan rumus.
 *              - Untuk memperjelas. Pada beberapa textbook, Matrix A == F, C == H, dan
 *                  fungsi non-linear f(..) sebagai transformasi X dan h(..) sebagai transformasi Z.
 *                  Perbedaan dengan algoritma pada implementasi ini hanya sebatas notasi saja.
 *              - Algoritma di atas mengasumsikan pada proses diskritisasi, lineasisasi
 *                  matrix A berdasarkan pendekatan:
 *                      X~(t=t_now) ~~ X~(k-1|k-1)      [linearisasi A pada X(t-1)]
 *                  Jika linearisasi menggunakan pendekatan:
 *                      X~(t=t_now) ~~ X~(k|k-1)        [linearisasi A pada X(t)]
 *                          maka urutan algoritma menjadi:
 *                              {CD_EKF_4} -> {CD_EKF_6} -> {CD_EKF_3} -> {CD_EKF_5} -> {CD_EKF_7}
 *              - Implementasi {CD_EKF_4) - {CD_EKF_7} di atas menggunakan implementasi naif yang mudah
 *                  untuk diimplementasikan,.
 *                  Sebagai contoh yang lebih robust & cepat, bisa dipakai pendekatan:
 *                    Jorgensen, J.B. A Computationally Efficient and Robust Implementation of
 *                      the Continuous-Discrete Extended Kalman Filter. 2007.
 *              - Pada contoh di atas X~(k=0|k=0) = [0]. Untuk mempercepat konvergensi bisa digunakan
 *                  informasi plant-spesific. Misal pada implementasi EKF untuk memfilter sensor
 *                  IMU (Inertial measurement unit) dengan X = [quaternion], dengan asumsi IMU
 *                  awalnya menghadap ke atas tanpa rotasi, X~(k=0|k=0) = [1, 0, 0, 0]'
 *
 *        Variabel:
 *          X_kira(k)    : X~(k) = x_estimasi(k) kalman filter   : Nx1
 *          X_dot_kira(k): X*~(k) = dX~(k)/dt                    : Nx1
 *          P(k)         : P(k) = matrix kovarian kalman filter  : NxN
 *          P_dot(k)     : P*(k) = dP(k)/dt                      : NxN
 *          A(k)         : Linearisasi dari fungsi non-linear f  : NxN
 *          C(k)         : Linearisasi dari fungsi non-linear h  : ZxN
 *          Q            : Matrix kovarian dari w(k)             : NxN
 *          R            : Matrix kovarian dari v(k)             : ZxZ
 *
 *
 ********************************************************************************************/

#include "cd_ekf.h"



CD_EKF::CD_EKF(const float_prec PInit, const float_prec QInit, const float_prec RInit)
{
    X_est.vIsiNol();
    X_est[0][0] = 1.0;       /* Quaternion(k = 0) = [1 0 0 0]' */

    P.vIsiDiagonal(PInit);
    Q.vIsiDiagonal(QInit);
    R.vIsiDiagonal(RInit);

}

void CD_EKF::vReset(const float_prec PInit, const float_prec QInit, const float_prec RInit)
{
    X_est.vIsiNol();
    X_est[0][0] = 1.0;       /* Quaternion(k = 0) = [1 0 0 0]' */
    X_dot_est.vIsiNol();

    P.vIsiDiagonal(PInit);
    Q.vIsiDiagonal(QInit);
    R.vIsiDiagonal(RInit);
}

void CD_EKF::vUpdate(Matrix &Z, Matrix &U)
{
    /* Dijalankan satu kali per waktu cuplik */


    /* ======================= Hitung paramater linearisasi A ======================= */
    /* A = d(F(..))/dx |x~(k-1|k-1),u(k)                                ...{CD_EKF_3} */
    vCalculateJacobianF(A, X_est, U);

    /* ================================== Prediksi ================================== */
    /* X*~(k|k-1) = F[X~(k-1|k-1), u(k)]                                ...{CD_EKF_4} */
    vUpdateNonlinearX(X_dot_est, X_est, U);

    /* P*(k|k-1)  = A*P(k-1|k-1) + P(k-1|k-1)*A' + Q                    ...{CD_EKF_5} */
    P_dot = (A*P) + (P*(A.Transpose())) + Q;

    /* X~(k|k-1)  = X~(k-1|k-1) + X*~(k|k-1)*dt    (Euler Method)       ...{CD_EKF_6} */
    X_est = X_est + (X_dot_est * SS_DT);

    /* P(k|k-1)   = P(k-1|k-1) + P*(k|k-1)*dt      (Euler Method)       ...{CD_EKF_7} */
    P = P + (P_dot * SS_DT);




    /* ======================= Hitung paramater linearisasi C ======================= */
    /* C = d(F(..))/dx |x~(k|k-1),u(k)                                  ...{CD_EKF_8} */
    vCalculateJacobianH(C, X_est, U);

    /* ================================== Koreksi =================================== */
    /* e(k)    = z(k) - H[X~(k|k-1), u(k)]                              ...{CD_EKF_9} */
    vUpdateNonlinearZ(Z_est, X_est, U);
    Err = Z - Z_est;

    /* S       = C*P(k|k-1)*C' + R                                      ...{CD_EKF_10} */
    S = (C*P*(C.Transpose())) + R;

    /* K       = P(k|k-1)*C'*(S^-1)                                     ...{CD_EKF_11} */
    Gain = P*(C.Transpose())*(S.Invers());

    /* X~(k|k) = X~(k|k-1) + K*e(k)                                     ...{CD_EKF_12} */
    X_est = X_est + (Gain*Err);

    /* P(k|k)  = (I - K*C)*P(k|k-1)                                     ...{CD_EKF_13} */
    Matrix Identitas = Matrix(SS_X_LEN, SS_X_LEN);
    Identitas.vSetIdentitas();
    P = (Identitas - (Gain*C))*P;




    /* ======= Tambahan untuk plant berupa sistem inersia berbasis quaternion ======= */

    /* # Lakukan normalisasi untuk memastikan batasan |q| = 1       --> Quaternion normalization
     * X = X/math.sqrt(X[0]**2 + X[1]**2 + X[2]**2 + X[3]**2)
     */
    if (!X_est.bNormVector()) {
        /* System error, reset EKF (atau sekalian reset MCU via watchdog?) */
        vReset(1000., Q[0][0], R[0][0]);
    }
}

void CD_EKF::vUpdateNonlinearX(Matrix &X_dot, Matrix &X, Matrix &U)
{
    float_prec q0, q1, q2, q3;
    float_prec p, q, r;

    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];

    p = U[0][0];
    q = U[1][0];
    r = U[2][0];

    /* Kode python box_quaternion_tanpa_magnetometer_EKF_vIMU6DOF+HMC.py:
     *  q0_dot = 1/2. * (  0   - p*q1 - q*q2 - r*q3)
     *  q1_dot = 1/2. * ( p*q0 +   0  + r*q2 - q*q3)
     *  q2_dot = 1/2. * ( q*q0 - r*q1 +  0   + p*q3)
     *  q3_dot = 1/2. * ( r*q0 + q*q1 - p*q2 +  0  )
     */

    X_dot[0][0] = 0.5 * (+0.00 -p*q1 -q*q2 -r*q3);
    X_dot[1][0] = 0.5 * (+p*q0 +0.00 +r*q2 -q*q3);
    X_dot[2][0] = 0.5 * (+q*q0 -r*q1 +0.00 +p*q3);
    X_dot[3][0] = 0.5 * (+r*q0 +q*q1 -p*q2 +0.00);
}

void CD_EKF::vUpdateNonlinearZ(Matrix &Z_est, Matrix &X, Matrix &U)
{
    float_prec q0, q1, q2, q3;
    float_prec q0_2, q1_2, q2_2, q3_2;

    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];

    q0_2 = q0 * q0;
    q1_2 = q1 * q1;
    q2_2 = q2 * q2;
    q3_2 = q3 * q3;

    /* Kode python box_quaternion_EKF_vIMU6DOF+HMC.py:
     *     DCM     = numpy.array([[(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                        2*(q1*q2+q0*q3),                        2*(q1*q3-q0*q2)],
     *                           [                   2*(q1*q2-q0*q3),     (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                        2*(q2*q3+q0*q1)],
     *                           [                   2*(q1*q3+q0*q2),                        2*(q2*q3-q0*q1),     (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]],
     *                           "float32")
     * 
     *  G_proj_sens = DCM * [0 0 -g]              --> Proyeksi gravitasi ke sensors
     *  M_proj_sens = DCM * [Mx My 0]             --> Proyeksi medan magnet ke sensors
     */
    Z_est[0][0] = (-2*q1*q3 +2*q0*q2) * ModelSistem_ACC_Z0;

    Z_est[1][0] = (-2*q2*q3 -2*q0*q1) * ModelSistem_ACC_Z0;

    Z_est[2][0] = (-(q0_2) +(q1_2) +(q2_2) -(q3_2)) * ModelSistem_ACC_Z0;

    Z_est[3][0] = (+(q0_2)+(q1_2)-(q2_2)-(q3_2)) * Mag_B0[0]
                   +2*(q1*q2+q0*q3) * Mag_B0[1]
                   +2*(q1*q3-q0*q2) * Mag_B0[2];

    Z_est[4][0] = 2*(q1*q2-q0*q3) * Mag_B0[0]
                   +(+(q0_2)-(q1_2)+(q2_2)-(q3_2)) * Mag_B0[1]
                   +2*(q2*q3+q0*q1) * Mag_B0[2];

    Z_est[5][0] = 2*(q1*q3+q0*q2) * Mag_B0[0]
                   +2*(q2*q3-q0*q1) * Mag_B0[1]
                   +(+(q0_2)-(q1_2)-(q2_2)+(q3_2)) * Mag_B0[2];
}

void CD_EKF::vCalculateJacobianF(Matrix &A, Matrix &X, Matrix &U)
{
    /* Kode python box_quaternion_tanpa_magnetometer_EKF_vIMU6DOF+HMC.py:
     *  q0_dot = 1/2. * (  0   - p*q1 - q*q2 - r*q3)
     *  q1_dot = 1/2. * ( p*q0 +   0  + r*q2 - q*q3)
     *  q2_dot = 1/2. * ( q*q0 - r*q1 +  0   + p*q3)
     *  q3_dot = 1/2. * ( r*q0 + q*q1 - p*q2 +  0  )
     */
    /* Kode python box_quaternion_tanpa_magnetometer_EKF_vIMU6DOF+HMC.py:
        F     = numpy.array([[    0, -1/2.*p, -1/2.*q, -1/2.*r],
                             [1/2.*p,      0,  1/2.*r, -1/2.*q],
                             [1/2.*q, -1/2.*r,      0,  1/2.*p],
                             [1/2.*r,  1/2.*q, -1/2.*p,      0]],
                             "float64")
    */
    float_prec p, q, r;

    p = U[0][0];
    q = U[1][0];
    r = U[2][0];

    A[0][0] =  0.000;
    A[1][0] =  0.5*p;
    A[2][0] =  0.5*q;
    A[3][0] =  0.5*r;

    A[0][1] = -0.5*p;
    A[1][1] =  0.000;
    A[2][1] = -0.5*r;
    A[3][1] =  0.5*q;

    A[0][2] = -0.5*q;
    A[1][2] =  0.5*r;
    A[2][2] =  0.000;
    A[3][2] = -0.5*p;

    A[0][3] = -0.5*r;
    A[1][3] = -0.5*q;
    A[2][3] =  0.5*p;
    A[3][3] =  0.000;
}

void CD_EKF::vCalculateJacobianH(Matrix &C, Matrix &X, Matrix &U)
{
    /* Kode python box_quaternion_EKF_vIMU6DOF+HMC.py:
     *     DCM     = numpy.array([[(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                        2*(q1*q2+q0*q3),                        2*(q1*q3-q0*q2)],
     *                           [                   2*(q1*q2-q0*q3),     (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                        2*(q2*q3+q0*q1)],
     *                           [                   2*(q1*q3+q0*q2),                        2*(q2*q3-q0*q1),     (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]],
     *                           "float32")
     * 
     *  G_proj_sens = DCM * [0 0 -g]              --> Proyeksi gravitasi ke sensors
     *  M_proj_sens = DCM * [Mx My 0]             --> Proyeksi medan magnet ke sensors
     */

    float_prec q0, q1, q2, q3;

    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];

    C[0][0] =  2*q2 * ModelSistem_ACC_Z0;
    C[1][0] = -2*q1 * ModelSistem_ACC_Z0;
    C[2][0] = -2*q0 * ModelSistem_ACC_Z0;
    C[3][0] =  2*q0*Mag_B0[0]+2*q3*Mag_B0[1]-2*q2*Mag_B0[2];
    C[4][0] = -2*q3*Mag_B0[0]+2*q0*Mag_B0[1]+2*q1*Mag_B0[2];
    C[5][0] =  2*q2*Mag_B0[0]-2*q1*Mag_B0[1]+2*q0*Mag_B0[2];

    C[0][1] = -2*q3 * ModelSistem_ACC_Z0;
    C[1][1] = -2*q0 * ModelSistem_ACC_Z0;
    C[2][1] =  2*q1 * ModelSistem_ACC_Z0;
    C[3][1] =  2*q1*Mag_B0[0]+2*q2*Mag_B0[1]+2*q3*Mag_B0[2];
    C[4][1] =  2*q2*Mag_B0[0]-2*q1*Mag_B0[1]+2*q0*Mag_B0[2];
    C[5][1] =  2*q3*Mag_B0[0]-2*q0*Mag_B0[1]-2*q1*Mag_B0[2];

    C[0][2] =  2*q0 * ModelSistem_ACC_Z0;
    C[1][2] = -2*q3 * ModelSistem_ACC_Z0;
    C[2][2] =  2*q2 * ModelSistem_ACC_Z0;
    C[3][2] = -2*q2*Mag_B0[0]+2*q1*Mag_B0[1]-2*q0*Mag_B0[2];
    C[4][2] =  2*q1*Mag_B0[0]+2*q2*Mag_B0[1]+2*q3*Mag_B0[2];
    C[5][2] =  2*q0*Mag_B0[0]+2*q3*Mag_B0[1]-2*q2*Mag_B0[2];

    C[0][3] = -2*q1 * ModelSistem_ACC_Z0;
    C[1][3] = -2*q2 * ModelSistem_ACC_Z0;
    C[2][3] = -2*q3 * ModelSistem_ACC_Z0;
    C[3][3] = -2*q3*Mag_B0[0]+2*q0*Mag_B0[1]+2*q1*Mag_B0[2];
    C[4][3] = -2*q0*Mag_B0[0]-2*q3*Mag_B0[1]+2*q2*Mag_B0[2];
    C[5][3] =  2*q1*Mag_B0[0]+2*q2*Mag_B0[1]+2*q3*Mag_B0[2];
}


