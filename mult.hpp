#include <cmath>
void mult(double *a, double *b, double *c, int rowa, int cola, int rowb, int colb) {
    if(cola != rowb) {
        printf("error: Incorrect matrix sizes!");
        return;
    }
    for (int i = 0; i < cola*rowb; i++) {
        c[i] = 0;
        if (fabs(a[i]) < 1e-90) {
            a[i] = 0.;
        }
        if(fabs(b[i]) < 1e-90) {
            b[i] = 0. ;
        }
    }
    int t, q, r;
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    int v3 = rowa % 3;
    int h3 = colb % 3;
    for(r = 0; r < v3; r++) {
        for(t = 0; t < h3; t++) {
            double sum = 0;
            for(q = 0; q < cola; q++) {
                sum += a[r * cola + q] * b[q * colb + t];
            }
            c[r * colb + t] += sum;
        }
        for(; t < colb; t += 3) {
            s00 = 0;
            s01 = 0;
            s02 = 0;
            for(q = 0; q < cola; q++) {
                s00 += a[r * cola + q] * b[q * colb + t];
                s01 += a[r * cola + q] * b[q * colb + t + 1];
                s02 += a[r * cola + q] * b[q * colb + t + 2];
            }
            c[r * colb + t] += s00;
            c[r * colb + t + 1] += s01;
            c[r * colb + t + 2] += s02;
        }
    }
    for(;r < rowa; r += 3) {
        for(t = 0; t < h3; t++) {
            s00 = 0;
            s10 = 0;
            s20 = 0;
            for(q = 0; q < colb; q++) {
                s00 += a[r * cola + q] * b[q * colb + t];
                s10 += a[(r + 1) * cola + q] * b[q * colb + t];
                s20 += a[(r + 2) * cola + q] * b[q * colb + t];
            }
            c[r * colb + t] += s00;
            c[(r + 1) * colb + t] += s10;
            c[(r + 2) * colb + t] += s20;
        }
        for(; t < colb; t+=3) {
            s00 = 0;
            s01 = 0;
            s02 = 0;
            s10 = 0;
            s11 = 0;
            s12 = 0;
            s20 = 0;
            s21 = 0;
            s22 = 0;
            for(q = 0; q < cola; q++) {
                s00 += a[r * cola + q] * b[q * colb + t];
                s01 += a[r * cola + q] * b[q * colb + t + 1];
                s02 += a[r * cola + q] * b[q * colb + t + 2];
                s10 += a[(r + 1) * cola + q] * b[q * colb + t];
                s11 += a[(r + 1) * cola + q] * b[q * colb + t + 1];
                s12 += a[(r + 1) * cola + q] * b[q * colb + t + 2];
                s20 += a[(r + 2) * cola + q] * b[q * colb + t];
                s21 += a[(r + 2) * cola + q] * b[q * colb + t + 1];
                s22 += a[(r + 2) * cola + q] * b[q * colb + t + 2];
            }
            c[r * colb + t] += s00;
            c[r * colb + t + 1] += s01;
            c[r * cola + t + 2] += s02;
            c[(r + 1) * colb + t] += s10;
            c[(r + 1) * colb + t + 1] += s11;
            c[(r + 1) * colb + t + 2] += s12;
            c[(r + 2) * colb + t] += s20;
            c[(r + 2) * colb + t + 1] += s21;
            c[(r + 2) * colb + t + 2] += s22;
        }
    }
}
