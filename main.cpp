#include <iostream>
#include <string>
#include <fstream>
#include "reader.hpp"
#include "solve.hpp"
#include "formula.hpp"
#include <cstring>
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
    int n = std::stoi(s1); // Размерность матрицы
    int m = std::stoi(s2); // Размерность блока
    int r = std::stoi(s3); // Кол-во выводимых значений в матрице
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
        std::cout << std::endl;
        PrintDouble(matr,n,r);
        double* block = new double[m*m];
        double* inverse = new double[m*m];
        double* temp = new double[m*m];
        double* temp1 = new double[m*m];
        double* temp2 = new double[m*m];
        int v = solve(n, m, matr, block, x, inverse, temp, temp1, temp2);
        if(v == 0) {}
        delete[] block;
        delete[] inverse;
        delete[] temp;
        delete[] temp2;
        delete[] temp1;
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
        int s = stoi(tmp);
        if (s < 1 || s > 4) {
            std::cerr << "error: Pasametr s is not a valid number" << std::endl;
            delete[] matr;
            delete[] x;
            return -7;
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                matr[n*i + j] = formula(s,n,i+1,j+1);
            }
        }
        PrintDouble(matr, n, r);
        std::cout << "--------------------------------" << std::endl;
        double* block = new double[m*m];
        double* inverse = new double[n*n];
        int i,j;
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                inverse[i*n+j] = 0;
        for(i=0;i<n;i++)
            inverse[i*n+i] = 1;
        double* temp = new double[m*m];
        double* temp1 = new double[m*m];
        double* temp2 = new double[m*m];

        int sd = solve(n, m, matr, block, x, inverse, temp, temp1, temp2);
        if(sd==-1) {}
        delete[] block;
        delete[] inverse;
        delete[] temp;
        delete[] temp2;
        delete[] temp1;
        //return 1;
    }
    delete[] matr;
    delete[] x;
    return 0;
}
