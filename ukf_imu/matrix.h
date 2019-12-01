/************************************************************************************
 *
 * Class Matrix
 *  Berisi kumpulan kode yang digunakan untuk melakukan operasi matrix.
 *
 *  Catatan:
 *    - Operasi indexing matrix dimulai dari 0, dengan format matrix[baris][kolom]
 *    - Data matrix disimpan dalam bentuk array 2 dimensi.
 *    - Setiap matrix yang menggunakan memori MATRIX_MAXIMUM_SIZE^2 (di variabel
 *       f32data), dan ukuran baris & kolom harus lebih kecil dari macro
 *       MATRIX_MAXIMUM_SIZE. Pendekatan ini digunakan untuk menghindari malloc
 *       (optimasi lebih lanjut bisa dilakukan untuk mengurangi penggunaan memori).
 * 
 * Class Matrix Versioning:
 *    v0.2 (2019-11-30), {PNb}: 
 *      -Fungsi yang disupport:
 *          - Operator ==
 *          - Normalisasi matrix
 *          - Cholesky Decomposition
 *          - InsertSubMatrix
 *          - InsertVector
 *
 *    v0.1 (2019-11-29), {PNb}: 
 *      -Fungsi yang disupport:
 *          - Operasi matrix dasar
 *          - Invers
 *          - Cetak
 * 
 *************************************************************************************/


#ifndef MATRIX_H
#define MATRIX_H

#include "konfig.h"

#if (SISTEM_IMPLEMENTASI == SISTEM_PC)
    #include <iostream>
    #include <iomanip>      // std::setprecision

    using namespace std;
#elif (SISTEM_IMPLEMENTASI == SISTEM_EMBEDDED_ARDUINO)
    #include <Wire.h>
#endif

#define MATRIX_PAKAI_BOUND_CHECKING

class Matrix
{
public:
    Matrix(const int32_t _i32baris, const int32_t _i32kolom)
    {
        this->i32baris = _i32baris;
        this->i32kolom = _i32kolom;
    }
    
    bool bCekMatrixValid() {
        /* Index terakhir untuk buffer jika ada kode yg buffer overflow 1 index */
        if ((this->i32baris > 0) && (this->i32baris < MATRIX_MAXIMUM_SIZE) && (this->i32kolom > 0) && (this->i32kolom < MATRIX_MAXIMUM_SIZE)) {
            return true;
        } else {
            return false;
        }
    }
    
    int32_t i32getBaris() { return this->i32baris; }
    int32_t i32getKolom() { return this->i32kolom; }
    
    /* Ref: https://stackoverflow.com/questions/6969881/operator-overload */
    class Proxy {
    public:
        Proxy(float_prec* _array, int32_t _maxKolom) : _array(_array) { this->_maxKolom = _maxKolom; }

        /* Modifikasi agar lvalue modifiable, ref:
         * https://stackoverflow.com/questions/6969881/operator-overload#comment30831582_6969904
         * (I know this is so dirty, but it makes the code so FABULOUS :D)
         */
        float_prec & operator[](int32_t _kolom) {
            #if defined(MATRIX_PAKAI_BOUND_CHECKING)
                if (_kolom >= this->_maxKolom) {
                    #if (SISTEM_IMPLEMENTASI == SISTEM_PC)
                        cout << "Matrix index out-of-bounds (kolom: " << _kolom << ")"<< endl;
                    #elif (SISTEM_IMPLEMENTASI == SISTEM_EMBEDDED_ARDUINO)
                        Serial.println("Matrix index out-of-bounds kolom");
                    #else
                        /* Silent function */
                    #endif
                    while(1);
                }
            #endif
            return _array[_kolom];
        }
    private:
        float_prec* _array;
        int32_t _maxKolom;
    };
    Proxy operator[](int32_t _baris) {
        #if defined(MATRIX_PAKAI_BOUND_CHECKING)
            if (_baris >= this->i32baris) {
                #if (SISTEM_IMPLEMENTASI == SISTEM_PC)
                    cout << "Matrix index out-of-bounds (baris: " << _baris << ")"<< endl;
                #elif (SISTEM_IMPLEMENTASI == SISTEM_EMBEDDED_ARDUINO)
                    Serial.println("Matrix index out-of-bounds baris");
                #else
                    /* Silent function */
                #endif
                while(1);
            }
        #endif
        return Proxy(f32data[_baris], this->i32kolom);      /* Parsing data kolom untuk bound checking */
    }

    bool operator == (Matrix _pembanding) {
        if ((this->i32baris != _pembanding.i32baris) || (this->i32kolom != _pembanding.i32getKolom())) {
            return false;
        }

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                if (fabs((*this)[_i][_j] - _pembanding[_i][_j]) > float_prec_ZERO) {
                    return false;
                }
            }
        }
        return true;
    }

    Matrix operator + (Matrix _penjumlah) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] + _penjumlah[_i][_j];
            }
        }
        return _outp;
    }

    Matrix operator - (Matrix _pengurang) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] - _pengurang[_i][_j];
            }
        }
        return _outp;
    }

    Matrix operator * (Matrix _pengali) {
        Matrix _outp(this->i32baris, _pengali.i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < _pengali.i32kolom; _j++) {
                _outp[_i][_j] = 0.0;
                for (int32_t _k = 0; _k < this->i32kolom; _k++) {
                    _outp[_i][_j] += ((*this)[_i][_k] * _pengali[_k][_j]);
                }
            }
        }
        return _outp;
    }

    Matrix operator * (float_prec _scalar) {
        Matrix _outp(this->i32baris, this->i32kolom);

        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j] * _scalar;
            }
        }
        return _outp;
    }

    void vRoundingElementToZero(const int32_t _i, const int32_t _j) {
        if (fabs((*this)[_i][_j]) < float_prec_ZERO) {
            (*this)[_i][_j] = 0.0;
        }
    }

    void vRoundingMatrixToZero() {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                if (fabs((*this)[_i][_j]) < float_prec_ZERO) {
                    (*this)[_i][_j] = 0.0;
                }
            }
        }
    }

    void vIsiHomogen(const float_prec _data) {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                (*this)[_i][_j] = _data;
            }
        }
    }

    void vIsiNol() {
        this->vIsiHomogen(0.0);
    }

    void vIsiRandom(const int32_t _batasAtas, const int32_t _batasBawah) {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                (*this)[_i][_j] = float_prec((rand() % (_batasAtas - _batasBawah + 1)) + _batasBawah);
            }
        }
    }

    void vIsiDiagonal(const float_prec _data) {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                if (_i == _j) {
                    (*this)[_i][_j] = _data;
                } else {
                    (*this)[_i][_j] = 0.0;
                }
            }
        }
    }

    void vSetIdentitas() {
        this->vIsiDiagonal(1.0);
    }

    /* Masukkan vektor ke matrix pada posisi kolom _posKolom
     * Contoh: A = Matrix 3x3, B = Vector 3x1
     *
     *  C = A.InsertVector(B, 1);
     *
     *  A = [A00  A01  A02]     B = [B00]
     *      [A10  A11  A12]         [B10]
     *      [A20  A21  A22]         [B20]
     *
     *  C = [A00  B00  A02]
     *      [A10  B10  A12]
     *      [A20  B20  A22]
     */
    Matrix InsertVector(Matrix _Vector, const int32_t _posKolom) {
        Matrix _outp(this->i32kolom, this->i32baris);
        if ((_Vector.i32baris != this->i32baris) || (_posKolom >= this->i32kolom)) {
            /* Return false */
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp = this->Salin();
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            _outp[_i][_posKolom] = _Vector[_i][0];
        }
        return _outp;
    }

    /* Masukkan submatrix ke matrix pada posisi baris _posBaris & posisi kolom _posKolom
     * Contoh: A = Matrix 4x4, B = Matrix 2x3
     *
     *  C = A.InsertSubMatrix(B, 1, 1);
     *
     *  A = [A00  A01  A02  A03]    B = [B00  B01  B02]
     *      [A10  A11  A12  A13]        [B10  B11  B12]
     *      [A20  A21  A22  A23]
     *      [A30  A31  A32  A33]
     *
     *
     *  C = [A00  A01  A02  A03]
     *      [A10  B00  B01  B02]
     *      [A20  B10  B11  B12]
     *      [A30  A31  A32  A33]
     */
    Matrix InsertSubMatrix(Matrix _subMatrix, const int32_t _posBaris, const int32_t _posKolom) {
        Matrix _outp(this->i32kolom, this->i32baris);
        if (((_subMatrix.i32baris+_posBaris) > this->i32baris) || ((_subMatrix.i32baris+_posKolom) > this->i32kolom)) {
            /* Return false */
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp = this->Salin();
        for (int32_t _i = 0; _i < _subMatrix.i32baris; _i++) {
            for (int32_t _j = 0; _j < _subMatrix.i32kolom; _j++) {
                _outp[_i + _posBaris][_j + _posKolom] = _subMatrix[_i][_j];
            }
        }
        return _outp;
    }

    Matrix Transpose() {
        Matrix _outp(this->i32kolom, this->i32baris);
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_j][_i] = (*this)[_i][_j];
            }
        }
        return _outp;
    }

    bool bNormVector() {
        float_prec _normM = 0.0;
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _normM = _normM + ((*this)[_i][_j] * (*this)[_i][_j]);
            }
        }
        if (_normM < float_prec_ZERO) {
            return false;
        }
        _normM = sqrt(_normM);
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                (*this)[_i][_j] /= _normM;
            }
        }
        return true;
    }
    
    Matrix Salin() {
        Matrix _outp(this->i32baris, this->i32kolom);
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                _outp[_i][_j] = (*this)[_i][_j];
            }
        }
        return _outp;
    }

    /* Melakukan operasi invers matrix dengan menggunakan algoritma Gauss-Jordan */
    Matrix Invers() {
        Matrix _outp(this->i32baris, this->i32kolom);
        Matrix _temp(this->i32baris, this->i32kolom);
        _outp.vSetIdentitas();
        _temp = this->Salin();


        /* Eliminasi Gauss... */
        for (int32_t _j = 0; _j < (_temp.i32baris)-1; _j++) {
            for (int32_t _i = _j+1; _i < _temp.i32baris; _i++) {
                if (fabs(_temp[_j][_j]) < float_prec_ZERO) {
                    // return false;    /* Matrix non-invertible */
                    _outp.i32baris = -1;
                    _outp.i32kolom = -1;
                    return _outp;
                }

                float_prec _tempfloat = _temp[_i][_j] / _temp[_j][_j];

                for (int32_t _k = 0; _k < _temp.i32kolom; _k++) {
                    _temp[_i][_k] -= (_temp[_j][_k] * _tempfloat);
                    _outp[_i][_k] -= (_outp[_j][_k] * _tempfloat);

                    _temp.vRoundingElementToZero(_i, _k);
                    _outp.vRoundingElementToZero(_i, _k);
                }

            }
        }

#if (1)
        /* Sampai sini seharusnya matrix _temp adalah matrix segitiga atas, tapi karena
         * keterbatasan kepresisian (rounding error), bisa jadi segitiga bawahnya
         * bukan 0 semua, jadikan 0 --> berguna untuk dekomposisi LU
         */
        for (int32_t _i = 1; _i < _temp.i32baris; _i++) {
            for (int32_t _j = 0; _j < _i; _j++) {
                _temp[_i][_j] = 0.0;
            }
        }
#endif


        /* Jordan... */
        for (int32_t _j = (_temp.i32baris)-1; _j > 0; _j--) {
            for (int32_t _i = _j-1; _i >= 0; _i--) {
                if (fabs(_temp[_j][_j]) < float_prec_ZERO) {
                    // return false;    /* Matrix non-invertible */
                    _outp.i32baris = -1;
                    _outp.i32kolom = -1;
                    return _outp;
                }

                float_prec _tempfloat = _temp[_i][_j] / _temp[_j][_j];
                _temp[_i][_j] -= (_temp[_j][_j] * _tempfloat);
                _temp.vRoundingElementToZero(_i, _j);

                for (int32_t _k = (_temp.i32baris - 1); _k >= 0; _k--) {
                    _outp[_i][_k] -= (_outp[_j][_k] * _tempfloat);
                    _outp.vRoundingElementToZero(_i, _k);
                }
            }
        }


        /* Normalisasi */
        for (int32_t _i = 0; _i < _temp.i32baris; _i++) {
            if (fabs(_temp[_i][_i]) < float_prec_ZERO) {
                // return false;    /* Matrix non-invertible */
                _outp.i32baris = -1;
                _outp.i32kolom = -1;
                return _outp;
            }

            float_prec _tempfloat = _temp[_i][_i];
            _temp[_i][_i] = 1.0;

            for (int32_t _j = 0; _j < _temp.i32baris; _j++) {
                _outp[_i][_j] /= _tempfloat;
            }
        }
        return _outp;
    }

    /* Melakukan operasi Dekomposisi Cholesky pada matrix dengan menggunakan algoritma Cholesky-Crout.
     *
     *      A = L*L'     ; A = matrix riil, positif definit, dan simetrik dengan ukuran MxM
     *
     *      L = A.CholeskyDec();
     *
     *      CATATAN! NOTE! Untuk menghemat komputansi, pengecekan matrix simetrik TIDAK dilakukan.
     *          Karena pemrosesan dilakukan pada segitiga kiri bawah dari matrix _A, maka
     *          diasumsikan segitiga atas dari _A juga merupakan simetrik dari segitiga bawah.
     *          (sebagai catatan, Scilab & MATLAB yang menggunakan Lapack routines DPOTRF
     *           memproses submatrix segitiga atas dari _A dan berperilaku kebalikan dengan
     *           fungsi ini, namun tetap valid secara matematis).
     *
     */
    Matrix CholeskyDec()
    {
        float_prec _tempFloat;

        Matrix _outp(this->i32baris, this->i32kolom);
        if (this->i32baris != this->i32kolom) {
            _outp.i32baris = -1;
            _outp.i32kolom = -1;
            return _outp;
        }
        _outp.vIsiHomogen(0.0);
        for (int32_t _j = 0; _j < this->i32kolom; _j++) {
            for (int32_t _i = _j; _i < this->i32baris; _i++) {
                _tempFloat = (*this)[_i][_j];
                if (_i == _j) {
                    for (int32_t _k = 0; _k < _j; _k++) {
                        _tempFloat = _tempFloat - (_outp[_i][_k] * _outp[_i][_k]);
                    }
                    if (_tempFloat < float_prec_ZERO) {
                        /* Matrix tidak positif definit */
                        _outp.i32baris = -1;
                        _outp.i32kolom = -1;
                        return _outp;
                    }
                    _outp[_i][_i] = sqrt(_tempFloat);
                } else {
                    for (int32_t _k = 0; _k < _j; _k++) {
                        _tempFloat = _tempFloat - (_outp[_i][_k] * _outp[_j][_k]);
                    }
                    if (fabs(_outp[_j][_j]) < float_prec_ZERO) {
                        /* Matrix tidak positif definit */
                        _outp.i32baris = -1;
                        _outp.i32kolom = -1;
                        return _outp;
                    }
                    _outp[_i][_j] = _tempFloat / _outp[_j][_j];
                }
            }
        }
        return _outp;
    }

#if (SISTEM_IMPLEMENTASI == SISTEM_PC)
    void vCetak() {
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            cout << "[ ";
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                cout << std::fixed << std::setprecision(3) << (*this)[_i][_j] << " ";
            }
            cout << "]";
            cout << endl;
        }
        cout << endl;
    }
#elif (SISTEM_IMPLEMENTASI == SISTEM_EMBEDDED_ARDUINO)
    void vCetak() {
        char _bufSer[10];
        for (int32_t _i = 0; _i < this->i32baris; _i++) {
            Serial.print("[ ");
            for (int32_t _j = 0; _j < this->i32kolom; _j++) {
                snprintf(_bufSer, sizeof(_bufSer)-1, "%2.2f ", (*this)[_i][_j]);
                Serial.print(_bufSer);
            }
            Serial.println("]");
        }
        Serial.println("");
    }
#else
    #warning("Fungsi Matrix.vCetak() tidak berfungsi");
    
    void vCetak() {}     /* Silent function */
#endif

private:
    int32_t i32baris;
    int32_t i32kolom;
    float_prec f32data[MATRIX_MAXIMUM_SIZE][MATRIX_MAXIMUM_SIZE] = {{0}};
};

#endif // MATRIX_H
