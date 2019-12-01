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

#ifndef CD_EKF_H
#define CD_EKF_H

#include "konfig.h"
#include "matrix.h"

class CD_EKF
{
public:
    CD_EKF(const float_prec PInit, const float_prec QInit, const float_prec RInit);
    void vReset(const float_prec PInit, const float_prec QInit, const float_prec RInit);
    void vUpdate(Matrix &Z, Matrix &U);
    Matrix BacaDataX() { return X_est; }
    Matrix BacaDataP() { return P; }
    Matrix BacaDataErr() { return Z_est; }

protected:
    void vUpdateNonlinearX(Matrix &X_dot, Matrix &X, Matrix &U);
    void vUpdateNonlinearZ(Matrix &Z_est, Matrix &X, Matrix &U);
    void vCalculateJacobianF(Matrix &A, Matrix &X, Matrix &U);
    void vCalculateJacobianH(Matrix &C, Matrix &X, Matrix &U);

private:
    Matrix X_est{SS_X_LEN, 1};
    Matrix X_dot_est{SS_X_LEN, 1};
    Matrix P{SS_X_LEN, SS_X_LEN};
    Matrix P_dot{SS_X_LEN, SS_X_LEN};
    Matrix A{SS_X_LEN, SS_X_LEN};
    Matrix C{SS_Z_LEN, SS_X_LEN};
    Matrix Z_est{SS_Z_LEN, 1};
    Matrix Err{SS_Z_LEN, 1};
    Matrix Q{SS_X_LEN, SS_X_LEN};
    Matrix R{SS_Z_LEN, SS_Z_LEN};
    Matrix S{SS_Z_LEN, SS_Z_LEN};
    Matrix Gain{SS_X_LEN, SS_Z_LEN};

    const float_prec Mag_B0[3] = {ModelSistem_B0x_Awal, ModelSistem_B0y_Awal, ModelSistem_B0z_Awal};
};

#endif // CD_EKF_H
