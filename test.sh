#!/bin/bash
for input_formula in {1..4} # изменяем размер матрицы
do
    for block_size in {1..10} # изменяем размер блока
    do
        for matrix_size in {10..20} # изменяем размер вывода
        do
            ./a.out $matrix_size $block_size 0 $input_formula # запускаем программу с разными аргументами
        done
    done
done