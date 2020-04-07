# Kode Teensy 4.0 untuk Extended Kalman Filter & Unscented Kalman Filter


# Hello, it seems you got here from ancient link, please go to [this repo for the latest updated code \(and more description\)!](https://github.com/pronenewbits/Arduino_AHRS_System).

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;


Koneksi Teensy 4.0 dengan MPU-9250:

![Alt text](MPU9250_Connection.png "MPU9250 Connection")





Jangan lupa kalibrasi magnetometer untuk mendapatkan hasil yang terbaik.

Contoh output magnetometer yang belum dikalibrasi:
![Alt text](2019-12-01_magneto_gabung_nonKalib.png "Uncalibrated Magnetometer")



Setelah dikalibrasi hard-iron:

![Alt text](2019-12-01_magneto_gabung_Kalib.png "Calibrated Magnetometer")

(TODO: implementasi kalibrasi soft-iron).



Untuk persamaan dinamik sistem IMU bisa dilihat di bawah ini:
![Alt text](Quaternion_IMU_Equation.png "Quaternion_IMU_Equation")
