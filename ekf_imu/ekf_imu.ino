#include <Wire.h>
#include <elapsedMillis.h>
#include "konfig.h"
#include "matrix.h"
#include "cd_ekf.h"
#include "MPU9250.h"


#if (SISTEM_IMPLEMENTASI == SISTEM_PC)
    #include <iostream>
    #include <iomanip>      // std::setprecision
    using namespace std;
#endif


float_prec BIAS_MAG[3] = {8.9254665, 8.040476,  -25.126487};
// float_prec BIAS_MAG[3] = {0,0,0};


elapsedMillis timerLed, timerEKF;
uint64_t u64lamaEKF;

CD_EKF EKF_IMU(1000., 1e-6, 0.0015);
Matrix Z(SS_Z_LEN, 1);
Matrix U(SS_U_LEN, 1);
Matrix dataQuaternion(SS_X_LEN, 1);
/* The command from the PC */
char cmd;
char bufferTxSer[100];

/* An MPU9250 object with the MPU-9250 sensor on I2C bus 0 with address 0x68 */
MPU9250 IMU(Wire,0x68);

void setup() {
    /* serial to display data */
    Serial.begin(115200);
    while(!Serial) {}

    /* start communication with IMU */
    int status = IMU.begin();
    if (status < 0) {
        Serial.println("IMU initialization unsuccessful");
        Serial.println("Check IMU wiring or try cycling power");
        Serial.print("Status: ");
        Serial.println(status);
        while(1) {}
    }
    IMU.setAccelRange(MPU9250::ACCEL_RANGE_2G);
    IMU.setGyroRange(MPU9250::GYRO_RANGE_2000DPS);
    IMU.setDlpfBandwidth(MPU9250::DLPF_BANDWIDTH_184HZ);
    IMU.setSrd(19);
}

void serialFloatPrint(float f) {
    byte * b = (byte *) &f;
    for (int i = 0; i < 4; i++) {
        byte b1 = (b[i] >> 4) & 0x0f;
        byte b2 = (b[i] & 0x0f);

        char c1 = (b1 < 10) ? ('0' + b1) : 'A' + b1 - 10;
        char c2 = (b2 < 10) ? ('0' + b2) : 'A' + b2 - 10;

        Serial.print(c1);
        Serial.print(c2);
    }
}


void loop() {
    if (timerEKF > SS_DT_MILIS) {
        /* Freq EKF = 50 Hz */
        /* Update data mentah IMU */
        IMU.readSensor();
        Z[0][0] = IMU.getAccelX_mss();
        Z[1][0] = IMU.getAccelY_mss();
        Z[2][0] = IMU.getAccelZ_mss();
        Z[3][0] = IMU.getMagX_uT()-BIAS_MAG[0];
        Z[4][0] = IMU.getMagY_uT()-BIAS_MAG[1];
        Z[5][0] = IMU.getMagZ_uT()-BIAS_MAG[2];
        U[0][0] = IMU.getGyroX_rads();
        U[1][0] = IMU.getGyroY_rads();
        U[2][0] = IMU.getGyroZ_rads();
        /* Normalisasi */
        float_prec _normG = (Z[0][0] * Z[0][0]) + (Z[1][0] * Z[1][0]) + (Z[2][0] * Z[2][0]);
        _normG = sqrt(_normG);
        Z[0][0] = Z[0][0] / _normG;
        Z[1][0] = Z[1][0] / _normG;
        Z[2][0] = Z[2][0] / _normG;
        float_prec _normM = (Z[3][0] * Z[3][0]) + (Z[4][0] * Z[4][0]) + (Z[5][0] * Z[5][0]);
        _normM = sqrt(_normM);
        Z[3][0] = Z[3][0] / _normM;
        Z[4][0] = Z[4][0] / _normM;
        Z[5][0] = Z[5][0] / _normM;
        
        /* Update Kalman */
        u64lamaEKF = micros();
        EKF_IMU.vUpdate(Z, U);
        u64lamaEKF = (micros() - u64lamaEKF);
//         snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "Lama EKF = %lu us", (uint32_t)u64lamaEKF);
//         Serial.println(bufferTxSer); 
        timerEKF = 0;
    }
    
    if (Serial.available()) {
        cmd = Serial.read();
        if (cmd == 'v') {
            snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "Kode EKF di Teensy 4.0, diadopsi agar bisa terhubung ke demo Processing FreeIMU");
            Serial.print(bufferTxSer);
            Serial.print('\n');
        } else if (cmd == 'q') {
            dataQuaternion = EKF_IMU.BacaDataX();

            while (!Serial.available());
            uint8_t count = Serial.read();
            for (uint8_t i = 0; i < count; i++) {
                
                serialFloatPrint(dataQuaternion[0][0]);
                Serial.print(",");
                serialFloatPrint(dataQuaternion[1][0]);
                Serial.print(",");
                serialFloatPrint(dataQuaternion[2][0]);
                Serial.print(",");
                serialFloatPrint(dataQuaternion[3][0]);
                Serial.print(",");
                serialFloatPrint((float)u64lamaEKF);
                Serial.print(",");
                Serial.println("");
            }
        }
    }
}

