#!/bin/bash

# Компилируем программу
g++ -I/usr/include/eigen3 test.cpp -o test

# Запускаем программу с разными аргументами
./test 1000 5
./test 1000 3
./test 1000 10
./test 1000 23
./test 1000 40
./test 1000 41
./test 1000 60
./test 2000 50
