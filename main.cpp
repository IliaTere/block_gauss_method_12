#include <iostream>
#include <string>
#include <fstream>
#include "reader.hpp"
#include "solve.hpp"
#include "formula.hpp"
#include <cstring>
#include <time.h>
#include "residual.hpp"
// #include <fenv.h>

bool isNumber(std::string& str)
    {
    std::string::iterator it = std::begin(str);
    while (it != str.end() && std::isdigit(*it)) {
        it++;
    }
    return !str.empty() && it == str.end();
    }

int main(int argc, char **argv)
{
//     feenableexcept(FE_DIVBYZERO | FE_INVALID
//   | FE_OVERFLOW | FE_UNDERFLOW);
    if (argc < 5 || argc > 6) {
        std::cout << "error: To many(few) arguments";
        return -1;
    }
    for (int i = 1; i < 5; i++) {
        std::string str(argv[i]);
        if (isNumber(str) == false) {
            std::cout << "error: Invalid argument: " << str << std::endl;
            return -4;
        }
    }
    std::string s1(argv[1]);
    std::string s2(argv[2]);
    std::string s3(argv[3]);
    int s;
    int n = std::stoi(s1); // Размерность матрицы
    int m = std::stoi(s2); // Размерность блока
    int r = std::stoi(s3); // Кол-во выводимых значений в матрице
    clock_t start;
    clock_t end;
    if ( m == 0  || m > n) {
        printf("invalid block\n");
        return -8;
    }
	double* matr = new double[n*n];
    double* x = new double[n*n];
    if (strcmp(argv[4],"0") == 0) {
        if (argc != 6) {
            std::cout << "error: File not found" << std::endl;
            delete[] matr;
            delete[] x;
            return -2;
        }
        std::string name(argv[5]);
		int t = read_ff(name , matr , n*n);
        if(t != 0) {
            
            delete[] matr;
            delete[] x;
            return -5;
        }
        std::string tmp(argv[4]);
        s = stoi(tmp);
    }
    if ((strcmp(argv[4],"0") != 0)) 
    {
        if (argc > 5) {
            std::cout << "error: To many(few) arguments" << std::endl;
            delete[] matr;
            delete[] x;
            return -6;
        }
        std::string tmp(argv[4]);
        s = stoi(tmp);
        if (s < 1 || s > 4) {
            std::cerr << "error: Parametr s is not a valid number" << std::endl;
            delete[] matr;
            delete[] x;
            return -7;
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                matr[n*i + j] = formula(s,n,i+1,j+1);
            }
        }
    }
    
    PrintDouble(matr, n, r);
    double* block = new double[m*m];
    double* inverse = new double[m*m];
    double* temp = new double[m*m];
    double* temp1 = new double[m*m];
    int* block_index = new int[m];
    double* matrtmp = new double[n*n];
    if (n <= 11000) {
        for (int i = 0; i < n*n; i++)
            matrtmp[i] = matr[i];
    }
    double matrix_norm = norma(matr, n);
    initialize_matrix(x, n);
    start = clock();
    printf("\n-------------------\n");
    int sd = solve(n, m, matr, block, x, inverse, temp, temp1, matrix_norm, block_index);
    end = clock();
    printf("\n-------------------\n");
    double t1 = (double)(end - start) / (double)CLOCKS_PER_SEC;


    switch(sd)
    {
        case -1:
            printf(
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
                argv[0], 12, -1., -1. ,0. , 0., s, n, m);
                break;
        case 2:
            printf(
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
                argv[0], 12, -1., -1. ,0., 0. , s, n, m);
                break;
        case 0:
            double res1 = 0. ;
            double res2 = 0. ;
            PrintDouble(x, n, r);
            if (n <= 11000)
            {
                start = clock();
                res1 = find_residual(matrtmp, x, block, inverse, n, m);
                res2 = find_residual(x, matrtmp, block, inverse, n, m);
                end = clock();
                double t2 = (double)(end - start) / (double)CLOCKS_PER_SEC;
                printf(
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
                argv[0], 12, res1, res2, t1, t2, s, n, m);
            } else {
                printf(
                "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
                argv[0], 12, 0. , 0. , t1, 0. , s, n, m);
            }
            break;
    }
    delete[] block_index;
    delete[] matrtmp;

    delete[] block;
    delete[] inverse;
    delete[] temp;
    delete[] temp1; 
    delete[] matr;
    delete[] x;
    return 0;
}
