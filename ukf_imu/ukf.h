/**********************************************************************************************************************
 *  Fungsi Discrete Unscented-kalman-filter (dihitung secara diskrit), Ref:
 *          Van der. Merwe, .. (2004). Sigma-Point Kalman Filters for Probabilistic
 *          Inference in Dynamic State-Space Models (Ph.D. thesis). Oregon Health &
 *          Science University.
 *
 *         Formulasi plant yang diestimasi:
 *              x*(k) = f[x(k-1), u(k)] + w(k)          ; x=Nx1,    u=Mx1               ...{UKF_1}
 *              z(k)  = h[x(k), u(k)] + v(k)            ; z=Zx1                         ...{UKF_2}
 *
 *        dimana:
 *          x*(k)   : Update x(k-1) terhadap waktu cuplik       : Nx1
 *          x(k-1)  : State Variable sistem waktu ke-k          : Nx1
 *          z(k)    : Keluaran sistem waktu ke-k                : Zx1
 *          u(k)    : Masukan sistem waktu ke-k                 : Mx1
 *          w(k)    : Noise proses, asumsi AWGN, kovarian Q     : Nx1
 *          v(k)    : Noise pengukuran, asumsi AWGN, kovarian R : Nx1
 *
 *          f(..), h(..) adalah model (non-linear) dari sistem yang ingin diprediksi
 *
 **********************************************************************************************************************
 *         Algoritma Discrete Unscented-kalman-filter:
 *          Init:
 *                  P (k=0|k=0) = Identitas * Kovarian(P(k=0)), biasanya diisi dg nilai yg besar
 *                  X~(k=0|k=0) = Ekspektasi(X(k=0)), biasanya = 0
 *                  Wc, Wm      = First order & second order weight
 *          Lakukan iterasi:
 *              Kalkulasi sigma point:
 *                  XSigma = [X~ Y+sqrt(C*P) Y-sqrt(C*P)    ; Y = [X~ ... X~], Y=NxN    ...{UKF_3}
 * 
 * 
 *              Unscented Transform XSigma [f,XSigma,U,Wm,Wc,Q] -> [X~,XSigma,P,dSigX]:
 *                  XSigma(k|k+1) = f(XSigma(k-1|k-1), u(k))                            ...{UKF_4}
 *                  X~(k|k+1) = sum(Wm(i) * XSigma(k|k+1)(i))    ; i = 1 ... (2N+1)     ...{UKF_5}
 * 
 *                  dSigX = XSigma(k|k+1)(i) - Y(k)  ; Y(k) = [X~(k|k-1) .. X~(k|k-1)]
 *                                                   ; Y(k)=Nx(2N+1)                    ...{UKF_6}
 *                  P(k|k-1) = sum(Wc(i)*dSigX(i)*dSigX(i)') + Q   ; i=1...(2N+1)       ...{UKF_7}
 *
 * 
 *              Unscented Transform ZSigma [h,XSigma,U,Wm,Wc,R] -> [Z~,ZSigma,Pz,dSigZ]:
 *                  ZSigma(k|k+1) = h(XSigma(k|k-1), u(k))                              ...{UKF_4}
 *                  Z~(k|k+1) = sum(Wm(i) * ZSigma(k|k+1)(i))    ; i = 1 ... (2N+1)     ...{UKF_5}
 * 
 *                  dSigZ = ZSigma(k|k+1)(i) - Y(k)  ; Y(k) = [Z~(k|k-1) .. Z~(k|k-1)]
 *                                                   ; Y(k)=Zx(2N+1)                    ...{UKF_6}
 *                  Pz(k|k-1) = sum(Wc(i)*dSigZ(i)*dSigZ(i)') + R   ; i=1...(2N+1)      ...{UKF_7}
 * 
 * 
 *              Update the estimated system:
 *                  CrossCov(k) = sum(Wc(i)*dSigX(i)*dSigZ(i)')     ; i=1...(2N+1)      ...{UKF_8}
 *                  K           = CrossCov(k) * dSigZ^-1                                ...{UKF_9}
 *                  X~(k|k)     = X~(k|k-1) + K * (Z(k) - Z~(k|k+1))                    ...{UKF_10}
 *                  P(k|k)      = P(k|k-1) - K*dSigZ*K'                                 ...{UKF_11}
 *
 *        *Catatan tambahan:
 *              - Perhatikan persamaan f pada {UKF_4} adalah update versi diskrit!!!!   X~(k+1) = f(X~(k),u(k))
 *              - Dengan asumsi masukan plant ZOH, u(k) = u(k|k-1),
 *                  Dengan asumsi tambahan observer dijalankan sebelum pengendali, u(k|k-1) = u(k-1),
 *                  sehingga u(k) [untuk perhitungan kalman] adalah nilai u(k-1) [dari pengendali].
 *              - Notasi yang benar adalah u(k|k-1), tapi disini menggunakan notasi u(k) untuk
 *                  menyederhanakan penulisan rumus.
 *              - Pada contoh di atas X~(k=0|k=0) = [0]. Untuk mempercepat konvergensi bisa digunakan
 *                  informasi plant-spesific. Misal pada implementasi EKF untuk memfilter sensor
 *                  IMU (Inertial measurement unit) dengan X = [quaternion], dengan asumsi IMU
 *                  awalnya menghadap ke atas tanpa rotasi, X~(k=0|k=0) = [1, 0, 0, 0]'
 *
 *        Variabel:
 *          X_kira(k)    : X~(k) = X_Estimasi(k) kalman filter   : Nx1
 *          X_dot_kira(k): X*~(k) = dX~(k)/dt                    : Nx1
 *          P(k)         : P(k) = matrix kovarian kalman filter  : NxN
 *          P_dot(k)     : P*(k) = dP(k)/dt                      : NxN
 *          A(k)         : Linearisasi dari fungsi non-linear f  : NxN
 *          C(k)         : Linearisasi dari fungsi non-linear h  : ZxN
 *          Q            : Matrix kovarian dari w(k)             : NxN
 *          R            : Matrix kovarian dari v(k)             : ZxZ
 *
 *
 **********************************************************************************************************************/

#ifndef UKF_H
#define UKF_H

#include "konfig.h"
#include "matrix.h"

class UKF
{
public:
    UKF(const float_prec PInit, const float_prec QInit, const float_prec RInit);
    void vReset(const float_prec PInit, const float_prec QInit, const float_prec RInit);
    void vUpdate(Matrix &Z, Matrix &U);
    Matrix BacaDataX() { return X_Est; }
    Matrix BacaDataP() { return P; }
    Matrix BacaDataErr() { return Err; }

protected:
    typedef  void (UKF::*UpdateNonLinear)(Matrix &X_dot, Matrix &X, Matrix &U);
    bool bCalculateSigmaPoint();
    bool bUnscentedTransform(Matrix &Out, Matrix &OutSigma, Matrix &P, Matrix &DSig,
                             UpdateNonLinear _vFuncNonLinear,
                             Matrix &InpSigma, Matrix &InpVector,
                             Matrix &_Wm, Matrix &_Wc, Matrix &_CovNoise);
    void vUpdateNonlinearX(Matrix &X_Next, Matrix &X, Matrix &U);
    void vUpdateNonlinearZ(Matrix &Z_est, Matrix &X, Matrix &U);

private:
    Matrix X_Est{SS_X_LEN, 1};
    Matrix X_Sigma{SS_X_LEN, (2*SS_X_LEN + 1)};
    
    Matrix Z_Est{SS_Z_LEN, 1};
    Matrix Z_Sigma{SS_Z_LEN, (2*SS_X_LEN + 1)};
    
    Matrix P{SS_X_LEN, SS_X_LEN};
    Matrix P_Chol{SS_X_LEN, SS_X_LEN};
    
    Matrix DX{SS_X_LEN, (2*SS_X_LEN + 1)};
    Matrix DZ{SS_Z_LEN, (2*SS_X_LEN + 1)};
    
    Matrix PZ{SS_Z_LEN, SS_Z_LEN};
    Matrix CrossCov{SS_X_LEN, SS_Z_LEN};
    
    Matrix Wm{1, (2*SS_X_LEN + 1)};
    Matrix Wc{1, (2*SS_X_LEN + 1)};
    
    Matrix Q{SS_X_LEN, SS_X_LEN};
    Matrix R{SS_Z_LEN, SS_Z_LEN};

    Matrix Err{SS_Z_LEN, 1};
    Matrix Gain{SS_X_LEN, SS_Z_LEN};
    float_prec UKFconst;

    const float_prec Mag_B0[3] = {ModelSistem_B0x_Awal, ModelSistem_B0y_Awal, ModelSistem_B0z_Awal};
};

#endif // UKF_H
