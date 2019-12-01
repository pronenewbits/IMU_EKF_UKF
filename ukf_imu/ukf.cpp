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

#include "ukf.h"



UKF::UKF(const float_prec PInit, const float_prec QInit, const float_prec RInit)
{
    float_prec _alpha   = 1e-3;     // Default, tunable
    float_prec _ki      = 0.0;      // Default, tunable
    float_prec _beta    = 2.0;      // Default, tunable
    float_prec _lambda;

    X_Est.vIsiNol();
    X_Sigma.vIsiNol();
    Z_Est.vIsiNol();
    Z_Sigma.vIsiNol();
    P.vIsiNol();
    P_Chol.vIsiNol();
    PZ.vIsiNol();
    DX.vIsiNol();
    DZ.vIsiNol();
    CrossCov.vIsiNol();
    Wm.vIsiNol();
    Wc.vIsiNol();
    Q.vIsiNol();
    R.vIsiNol();
    Err.vIsiNol();
    Gain.vIsiNol();


    X_Est.vIsiNol();
    X_Est[0][0] = 1.0;       /* Quaternion(k = 0) = [1 0 0 0]' */

    P.vIsiDiagonal(PInit);
    Q.vIsiDiagonal(QInit);
    R.vIsiDiagonal(RInit);

    _lambda = (_alpha*_alpha)*(SS_X_LEN+_ki) - SS_X_LEN;
    UKFconst = SS_X_LEN + _lambda;

    Wm[0][0] = _lambda/UKFconst;
    for (int32_t _i = 1; _i < Wm.i32getKolom(); _i++) {
        Wm[0][_i] = 0.5/UKFconst;
    }
    Wc = Wm.Salin();
    Wc[0][0] = Wc[0][0] + (1.0-(_alpha*_alpha)+_beta);

    UKFconst = sqrt(UKFconst);
}

void UKF::vReset(const float_prec PInit, const float_prec QInit, const float_prec RInit)
{
    X_Est.vIsiNol();
    X_Est[0][0] = 1.0;       /* Quaternion(k = 0) = [1 0 0 0]' */

    P.vIsiDiagonal(PInit);
    Q.vIsiDiagonal(QInit);
    R.vIsiDiagonal(RInit);
}

void UKF::vUpdate(Matrix &Z, Matrix &U)
{
    /* Dijalankan satu kali per waktu cuplik */


    /*  XSigma = [X~ Y+sqrt(C*P) Y-sqrt(C*P)    ; Y = [X~ ... X~], Y=NxN    ...{UKF_3}  */
    bCalculateSigmaPoint();
    /* Unscented Transform XSigma [f,XSigma,U,Wm,Wc,Q] -> [X~,XSigma,P,dSigX]: ======== */
    bUnscentedTransform(X_Est, X_Sigma, P, DX, (&UKF::vUpdateNonlinearX), X_Sigma, U, Wm, Wc, Q);
    /* Unscented Transform ZSigma [h,XSigma,U,Wm,Wc,R] -> [Z~,ZSigma,Pz,dSigZ]: ======= */
    bUnscentedTransform(Z_Est, Z_Sigma, PZ, DZ, (&UKF::vUpdateNonlinearZ), X_Sigma, U, Wm, Wc, R);


    /* Update the estimated system: =================================================== */
    /*  CrossCov(k) = sum(Wc(i)*dSigX(i)*dSigZ(i)')     ; i=1...(2N+1)      ...{UKF_8}  */
    for (int32_t _i = 0; _i < DX.i32getBaris(); _i++) {
        for (int32_t _j = 0; _j < DX.i32getKolom(); _j++) {
            DX[_i][_j] *= Wc[0][_j];
        }
    }
    CrossCov = DX * (DZ.Transpose());

    /*  K           = CrossCov(k) * dSigZ^-1                                ...{UKF_9}  */
    Matrix PZ_Inv = PZ.Invers();
    if (!PZ_Inv.bCekMatrixValid()) {
        /* return false; */
        this->vReset(1000., Q[0][0], R[0][0]);
    }
    Gain = CrossCov * PZ_Inv;

    /*  X~(k|k)     = X~(k|k-1) + K * (Z(k) - Z~(k|k+1))                    ...{UKF_10}  */
    Err = Z - Z_Est;
    X_Est = X_Est + (Gain*Err);

    /*  P(k|k)      = P(k|k-1) - K*dSigZ*K'                                 ...{UKF_11}  */
    P = P - (Gain * PZ * Gain.Transpose());




    /* ======= Tambahan untuk plant berupa sistem inersia berbasis quaternion ======= */
    /* # Lakukan normalisasi untuk memastikan batasan |q| = 1       --> Quaternion normalization
     * X = X/math.sqrt(X[0]**2 + X[1]**2 + X[2]**2 + X[3]**2)
     */
    /* HATI-HATI!! Jika inisialisasi salah, (X(k=0) = [0]), hitungan akan NAN (pembagian dengan nol)!!  */
    if (!X_Est.bNormVector()) {
        /* System error, reset EKF (atau sekalian reset MCU via watchdog?) */
        vReset(1000., Q[0][0], R[0][0]);
    }
}

bool UKF::bCalculateSigmaPoint(void)
{
    /* XSigma = [X~ Y+sqrt(C*P) Y-sqrt(C*P)    ; Y = [X~ ... X~], Y=NxN    ...{UKF_3}  */
    
    
    P_Chol = P.CholeskyDec();
    if (!P_Chol.bCekMatrixValid()) {
        /* System Fail */
        return false;
    }
    P_Chol = P_Chol * UKFconst;

    /* _xSigma = [x Y+C*sqrt(P) Y-C*sqrt(P)], dimana
     *  Y = [x x ... x] dengan ukuran (len(x) x len(x))     (misal x = 3x1, Y = 3x3)
     *  Sehingga didapat len(_xSigma) = 2len(x) + 1
     */
    Matrix _Y(SS_X_LEN, SS_X_LEN);
    for (int32_t _i = 0; _i < SS_X_LEN; _i++) {
        _Y = _Y.InsertVector(X_Est, _i);
    }
    X_Sigma.vIsiNol();
    /* Set _xSigma = [x 0 0] */
    X_Sigma = X_Sigma.InsertVector(X_Est, 0);
    /* Set _xSigma = [x Y+C*sqrt(P) 0] */
    X_Sigma = X_Sigma.InsertSubMatrix((_Y + P_Chol), 0, 1);
    /* Set _xSigma = [x Y+C*sqrt(P) Y-C*sqrt(P)] */
    X_Sigma = X_Sigma.InsertSubMatrix((_Y - P_Chol), 0, (1+SS_X_LEN));

    return true;
}

bool UKF::bUnscentedTransform(Matrix &Out, Matrix &OutSigma, Matrix &P, Matrix &DSig,
                              UpdateNonLinear _vFuncNonLinear,
                              Matrix &InpSigma, Matrix &InpVector,
                              Matrix &_Wm, Matrix &_Wc, Matrix &_CovNoise)
{
    /* XSigma(k|k+1) = f(XSigma(k-1|k-1), u(k))                            ...{UKF_4}  */
    /* X~(k|k+1) = sum(Wm(i) * XSigma(k|k+1)(i))    ; i = 1 ... (2N+1)     ...{UKF_5}  */
    Out.vIsiNol();
    for (int32_t _j = 0; _j < InpSigma.i32getKolom(); _j++) {
        /* Transformasi non-linear per kolom */
        Matrix _AuxSigma1(InpSigma.i32getBaris(), 1);
        Matrix _AuxSigma2(OutSigma.i32getBaris(), 1);
        for (int32_t _i = 0; _i < InpSigma.i32getBaris(); _i++) {
            _AuxSigma1[_i][0] = InpSigma[_i][_j];
        }
        (this->*(_vFuncNonLinear))(_AuxSigma2, _AuxSigma1, InpVector);      /* Welp... */
        /* Gabungkan hasil transformasi non-linear ke matrix outSigma */
        OutSigma = OutSigma.InsertVector(_AuxSigma2, _j);

        _AuxSigma2 = _AuxSigma2 * _Wm[0][_j];
        Out = Out + _AuxSigma2;
    }

    /* dSigX = XSigma(k|k+1)(i) - Y(k)  ; Y(k) = [X~(k|k-1) .. X~(k|k-1)]
    /*                                  ; Y(k)=Nx(2N+1)                    ...{UKF_6}  */
    Matrix _AuxSigma1(OutSigma.i32getBaris(), OutSigma.i32getKolom());
    for (int32_t _j = 0; _j < OutSigma.i32getKolom(); _j++) {
        _AuxSigma1 = _AuxSigma1.InsertVector(Out, _j);
    }
    DSig = OutSigma - _AuxSigma1;

    /* Px(k|k-1) = sum(Wc(i)*dSigX(i)*dSigX(i)') + Q   ; i=1...(2N+1)      ...{UKF_7}  */
    _AuxSigma1 = DSig.Salin();
    for (int32_t _i = 0; _i < DSig.i32getBaris(); _i++) {
        for (int32_t _j = 0; _j < DSig.i32getKolom(); _j++) {
            _AuxSigma1[_i][_j] *= _Wc[0][_j];
        }
    }
    P = (_AuxSigma1 * (DSig.Transpose())) + _CovNoise;

    return true;
}

void UKF::vUpdateNonlinearX(Matrix &X_Next, Matrix &X, Matrix &U)
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

    X_Next[0][0] = (0.5 * (+0.00 -p*q1 -q*q2 -r*q3))*SS_DT + q0;
    X_Next[1][0] = (0.5 * (+p*q0 +0.00 +r*q2 -q*q3))*SS_DT + q1;
    X_Next[2][0] = (0.5 * (+q*q0 -r*q1 +0.00 +p*q3))*SS_DT + q2;
    X_Next[3][0] = (0.5 * (+r*q0 +q*q1 -p*q2 +0.00))*SS_DT + q3;
}

void UKF::vUpdateNonlinearZ(Matrix &Z_est, Matrix &X, Matrix &U)
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
    Z_est[0][0] = (2*q1*q3 -2*q0*q2) * ModelSistem_ACC_Z0;

    Z_est[1][0] = (2*q2*q3 +2*q0*q1) * ModelSistem_ACC_Z0;

    Z_est[2][0] = (+(q0_2) -(q1_2) -(q2_2) +(q3_2)) * ModelSistem_ACC_Z0;

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


