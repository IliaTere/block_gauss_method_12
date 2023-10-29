#include <iostream>
#include <math.h>
#include <float.h>
#include "mult.hpp"
#include <iomanip>
//#include "reader.hpp"
#define UNUSED(x) (void)(x)
#define eps -DBL_MAX

using namespace std;

void get_block (double *matr, double *block , int n, int m, int i , int j )
{  
    int ii=0,jj=0;
    int k = n / m;
    int l = n - k * m;
    int v = (j < k ? m : l), h = (i < k ? m : l);
    //std::cout << "v= " << v << " h= " << h << std::endl;
    //block = new double[v * h];
    double *ma_bl = matr + i * n * m+ j * m;
    for (ii = 0; ii < h; ii++)
    {
        for (jj = 0; jj < v; jj++)
        {
            block[ii *m+jj]=ma_bl[ii *n+jj];
        }   
    }
}
void put_block (double *matr, double *block , int n, int m, int i , int j )
{
int ii=0,jj=0;
int k = n / m;
int l = n - k * m;
int v = (j < k ? m : l), h = (i < k ? m : l);
double *ma_bl = matr + i * n * m+ j * m;
for (ii = 0; ii < h; ii++)
{
    for (jj = 0; jj < v; jj++)
    {
        ma_bl[ii *n+jj]=block[ii *m+jj];
    }   
}
}

double norma(double* block, int m) {
    int i,j;
    double norm_m=0, temp;
    for(i=0;i<m;i++)  
    {
        temp=0;
        for(j=0;j<m;j++)
            temp+= fabs(block[j*m + i]);
        if(temp>norm_m)
            norm_m=temp;
    }
    return norm_m;
} 

    int treug(double* a,double* b,int n)
{
    for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                b[i*n+j] = 0;
        for(int i=0;i<n;i++)
            b[i*n+i] = 1;
    double *c;
    c = new double[n];

    int i,j,k,t;
    double p;

    //Проверяем на каждом шаге, что матрица сводится к треугольному виду.
    //Иначе ранк не максимален и нет обратной матрицы.
    for(i=0;i<n;i++)
    {
        t = -1;
        for(j=i;j<n;j++)
            if(a[j*n+i]>1e-18 || a[j*n+i]<-1e-18) t = 1;
        if(t==-1)
        {
            return -1;
        }

        p = a[i*n+i];
        t = i;

    //Метод Гаусса с выбором главного элемента по строке
        for(j=i;j<n;j++)
        {
            if(fabs(a[j*n+i])>fabs(p))
            {
                p = a[j*n+i];
                t = j;
            }
        }

        for(j=0;j<n;j++) c[j] = a[t*n+j];
        for(j=0;j<n;j++) a[t*n+j] = a[i*n+j];
        for(j=0;j<n;j++) a[i*n+j] = c[j];

        for(j=0;j<n;j++) c[j] = b[t*n+j];
        for(j=0;j<n;j++) b[t*n+j] = b[i*n+j];
        for(j=0;j<n;j++) b[i*n+j] = c[j];

        p = a[i*n+i];

        for(j=0;j<n;j++)
        {
            a[i*n+j] = a[i*n+j]/p;
            b[i*n+j] = b[i*n+j]/p;
        }

        for(j=i+1;j<n;j++)
        {
            p = a[j*n+i];
            for(k=0;k<n;k++)
            {
                a[j*n+k] = a[j*n+k] - a[i*n+k]*p;
                b[j*n+k] = b[j*n+k] - b[i*n+k]*p;
            }
        }

    }
    delete[]c;
    return 1;
}

void diag(double* a,double* b,int n)
{
    int i,j,k;
    double p;
    for(i=n-1;i>0;i--)
    {
        for(j=0;j<i;j++)
        {
            p = a[j*n+i];
            for(k=0;k<n;k++)
            {
                a[j*n+k] = a[j*n+k] - a[i*n+k]*p;
                b[j*n+k] = b[j*n+k] - b[i*n+k]*p;
            }
        }
    }
}

double ravno(double x)
{
    if(fabs(x)<1e-7){
        return 0.0;}
    else {return x;}
}

void ravnoBlock(double *block, int m)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            block[j*m + i] = ravno(block[j*m + i]);
        }
    }
}

int findmax(double* matr, double* block, int n, int m,int l, int j) { // TODO:  подавать на вход inverse
    int k = -1;
    int t;
    double tmp;
    double max_n = DBL_MAX;
    double* inverse = new double[m * m];
    for (int i = l; i < n/m; i++) {
        get_block(matr, block, n ,m , i, j);
        t = treug(block, inverse, m);
        if (t==-1) {
            continue;
        }
        diag(block,inverse,m);
        tmp = norma(inverse, m);
        if (tmp < max_n) {
            max_n = tmp;
            k = i;
        }
        
    }
    delete [] inverse;
    return k;
}

void swap_rows(double* matr, int k, int l, int n, int m) {
    double* tmp = new double[m*n];
    for (int i = 0; i < m*n; i++) {
        tmp[i] = matr[k*m * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*k * n + i] = matr[m*l * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*l * n + i] = tmp[i];
    }
    delete[] tmp;
}

void subtraction(double* a, double* b, int m) {
    for (int i = 0; i < m*m; i++) {
        a[i] -= b[i];
    }    
}

int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, double* block2) {
    for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                solution[i*n+j] = 0;
        for(int i=0;i<n;i++)
            solution[i*n+i] = 1;

    int t = -1; // ToDo: можем итерироваться до size который менятеся в зависимости от n, m
    
    for (int i = 0; i < n/m ; i++) {
        t = findmax(matr, block, n, m, i, i);
        if (t == -1) {
            throw("Exception");
            std::cout << "Необратимая матрица"  << std::endl;
        }
        if (t != i) {
        swap_rows(matr, i, t, n, m);
        swap_rows(solution, i, t, n, m);
        }
        get_block(matr, block, n, m, i, i);

        treug(block, inverse, m);
        diag(block,inverse,m);
        for (int j=i; j < n/m; j++) {
            get_block(matr, block, n, m, i, j);
            mult(inverse, block, tmp, m, m, m ,m);
            put_block(matr, tmp, n,m, i, j);
            get_block(solution, block, n, m, i, j);
            mult(inverse, block, tmp, m, m, m ,m);
            //ravnoBlock(tmp, m);
            put_block(solution, tmp, n,m, i, j);
        }                                       // Правильно +
        
        
        for (int f=i+1; f < n/m; f++) {

            
            get_block(matr, block, n, m, f, i); // То на что умножаем

            for(int j=i;j<n/m;j++) {
                get_block(matr, block1, n, m, i, j); // То что умножаем
                mult(block, block1, tmp, m , m ,m , m);
                get_block(matr, block2, n, m, f, j); // должно быть (f j)
                subtraction(block2, tmp, m);
                put_block(matr, block2, n,m, f, j);

                get_block(solution, block1, n, m, i, j);
                mult(block, block1, tmp, m , m ,m , m);
                get_block(solution, block2, n, m, f, j);
                subtraction(block2, tmp, m);
                put_block(solution, block2, n,m, f, j);
            }
        }
    }
    //Обратный ход
        for (int i = n/m-1; i > 0; i--) {
            for (int j = i-1; j > -1; j--) {
               
                get_block(matr,block1,n,m,j,i); //block1 - то на что умножаем строку
                
                for (int k = 0; k < n/m; k++) {
                    get_block(solution ,block, n ,m , i, k);
                    mult(block1, block, tmp, m, m, m, m);
                    
                    get_block(solution ,block2, n,m, i-1, k);
                
                    subtraction(block2, tmp, m);
                    //ravnoBlock(block2, m);
                    put_block(solution,block2, n,m, i-1, k);
                }  
            }
        }    
    
    printf("-----------------------ОТВЕТ\n");
    PrintDouble(matr, n, n);
    PrintDouble(solution, n , n);
    return 0;
}
