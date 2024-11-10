#include <iostream>
#include <stdlib.h>
#include <cmath>
using namespace std;

double formula(int s, int n, int i, int j) {
    if(s == 1)
        return (double)(n-max(i,j) + 1);
    if(s == 2)
        return (double)max(i,j);
    if(s == 3)
        return (double)abs(i-j);
    if(s == 4)
        return 1/(double)(i+j-1);
    return 0;   
}
