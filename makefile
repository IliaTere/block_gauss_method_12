# Переменные
CC = g++
CFLAGS = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
CFLAGS_LITE = -O3 -fstack-protector-all -g -W -Wall -Wextra -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-property-attribute-mismatch
# Список всех исходных файлов
SOURCES = main.cpp mult.cpp reader.cpp solve.cpp formula.cpp residual.cpp

# Список всех объектных файлов
OBJECTS = $(SOURCES:.cpp=.o)

# Целевой исполняемый файл
TARGET = a.out

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS)

# Правило для компиляции каждого исходного файла
%.o: %.cpp functions.hpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
# Правило для компиляции без флагов -Wunused
lite:
	$(CC) $(CFLAGS_LITE) -c main.cpp -o main.o
	$(CC) $(CFLAGS_LITE) -o $(TARGET) main.o
# Правило для запуска тестов
test: $(TARGET)
	rm tests/$(TARGET)
	cp $(TARGET) tests/$(TARGET)
	cd tests && ./test_small.sh
