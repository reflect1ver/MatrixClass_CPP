#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <math.h>

#include <iostream>

#define E8 0.00000001

class S21Matrix {
 private:
  int _rows, _cols;
  double** _p;
  void Free();
  void mem();
  S21Matrix GetMinorMatrix(int rows, int cols);
  void Fill_by_number(double n);

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  bool EqMatrix(const S21Matrix& other) const noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(const double num, const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  S21Matrix operator=(const S21Matrix& other);
  S21Matrix operator+=(const S21Matrix& other);
  S21Matrix operator-=(const S21Matrix& other);
  S21Matrix operator*=(const double num);
  S21Matrix operator*=(const S21Matrix& other);
  double& operator()(int row, int col);

  int GetRows();
  int GetCols();
  void SetRows(int rows);
  void SetCols(int cols);
};

#endif
