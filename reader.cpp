#include"functions.hpp"

bool is_double(const std::string& s) {
    if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+') && (s[0] != '.'))) return false;
    char* p;
    strtod(s.c_str(), &p);
    return (*p == 0);
}

int read_ff(const std::string& filename, double* result, int n) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error opening file: " << filename << std::endl;
        return -1;
    }
    std::string line;
    int i = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            if (is_double(token)) {
                if (i > n-1) {
                    printf("Incorrect size\n");
                    return -3;
                }
                result[i] = std::stod(token);
                i++;
            } else {
                std::cout << "Invalid double value: " << token << std::endl;
                return -2;
            }
        }
    }
    if (i != n) {
        std::cout << "Invalid length" << std::endl;
        return -3;
    }
    return 0;
}


void PrintDouble(double* matrix, int n, int r) {
	for(int i = 0; i < std::min(n, r); i++) {
        for(int j = 0; j < std::min(n, r); j++) {
            printf("%10.3e ", matrix[i*n+j]);
        }
        std::cout << std::endl;
    }
	printf("\n");
}
