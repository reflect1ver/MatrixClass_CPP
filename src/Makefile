CC=gcc
FLAGS=-std=c++17 -pedantic -Wall -Werror -Wextra
CHECKFLAGS=-lgtest -pedantic-errors -pthread #-lgmock
REPORTDIR=gcov_report
GCOV=--coverage
OS = $(shell uname)

all: test

s21_matrix_oop.a:
	$(CC) $(FLAGS) -c s21_matrix_oop.cpp -o matrix_oop.o
	ar rcs s21_matrix_oop.a matrix_oop.o

test: clean
	$(CC) $(FLAGS) $(GCOV) -c s21_matrix_oop.cpp
	$(CC) $(FLAGS) -c test.cpp $(CHECKFLAGS)
	$(CC) $(FLAGS) $(GCOV) -o matrix_test test.o s21_matrix_oop.o $(CHECKFLAGS) -lstdc++
	./matrix_test

gcov_report: test
	lcov -t "Unit-tests of matrix_oop" -o s21_matrix_oop.info -c -d .
	genhtml -o $(REPORTDIR) s21_matrix_oop.info
	open -a "Safari" ./$(REPORTDIR)/index.html

clang:
	clang-format -n *.h s21_matrix_oop.cpp test.cpp

check:
	cppcheck *.cpp && cppcheck --enable=all --suppress=missingIncludeSystem --language=c++ *.h

clean:
	@rm -rf ./*.o ./*.a ./a.out ./*.gcno ./*.gcda ./$(REPORTDIR) *.info ./*.info report matrix_test matrix_oop

valgrind: all
	valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all -v ./matrix_test