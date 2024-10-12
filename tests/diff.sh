#!/bin/bash

# Проверяем, что переданы два аргумента
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file1> <file2>"
    exit 1
fi

file1="$1"
file2="$2"
output_file="differences.txt"

# Очищаем файл с отличиями перед началом
> "$output_file"

# Сравниваем строки в двух файлах
while IFS= read -r line1; do
    # Проверяем, есть ли строка в файле 2
    if ! grep -Fxq "$line1" "$file2"; then
        echo "$line1" >> "$output_file"
    fi
done < "$file1"

while IFS= read -r line2; do
    # Проверяем, есть ли строка в файле 1
    if ! grep -Fxq "$line2" "$file1"; then
        echo "$line2" >> "$output_file"
    fi
done < "$file2"

echo "Differences have been written to $output_file"