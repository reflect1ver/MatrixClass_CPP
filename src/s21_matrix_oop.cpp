#include "s21_matrix_oop.h"

//////////////// КОНСТРУКТОРЫ И ДЕСТРУКТОР ////////////////
S21Matrix::S21Matrix() {
  _rows = 3;
  _cols = 3;
  mem();
}

S21Matrix::S21Matrix(int rows, int cols) : _rows(rows), _cols(cols) {
  if (_rows <= 0 || _cols <= 0) {
    throw std::out_of_range(
        "Incorrect input, matrices should have  cols and rows > 0");
  } else {
    mem();
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : _rows(other._rows), _cols(other._cols) {
  mem();
  for (int i = 0; i < _rows; i++) {
    for (int k = 0; k < _cols; k++) {
      _p[i][k] = other._p[i][k];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  _p = other._p;
  _rows = other._rows;
  _cols = other._cols;
  other._p = nullptr;
  other._rows = 0;
  other._cols = 0;
}

S21Matrix::~S21Matrix() { Free(); }

////////////////////// ТРЕБУЕМЫЕ ФУНКЦИИ //////////////////////

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _p[i][j] = _p[i][j] + other._p[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _p[i][j] = _p[i][j] - other._p[i][j];
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  bool res = true;
  if ((_rows != other._rows) || (_cols != other._cols)) {
    res = false;
  } else {
    for (int i = 0; i < _rows; i++) {
      for (int j = 0; j < _cols; j++) {
        if (fabs(_p[i][j] - other._p[i][j]) > 1e-7) {
          res = false;
          break;
        }
      }
      if (res == false) {
        break;
      }
    }
  }
  return (res);
}

void S21Matrix::MulNumber(const double num) {
  if (!_p) {
    throw std::out_of_range(
        "Incorrect input, matrices should have not null pointer");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _p[i][j] = _p[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (!_p && !other._p) {
    throw std::out_of_range(
        "Incorrect input, matrices should have not null pointer");
  } else if (_cols != other._rows) {
    throw std::out_of_range(
        "The number of columns of the first matrix is ​​not equal to the "
        "value of the second row of the matrix");
  } else {
    double res = 0;
    S21Matrix result(_rows, _cols);
    for (int i = 0; i < _rows; i++) {
      for (int j = 0; j < other._cols; j++) {
        for (int k = 0; k < _cols; k++) {
          res = res + _p[i][k] * other._p[k][j];
        }
        result._p[i][j] = res;
        res = 0;
      }
    }
    *this = result;
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(_cols, _rows);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      res._p[j][i] = _p[i][j];
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (_p == NULL) {
    throw std::out_of_range("Null pointer at source");
  } else if (_cols != _rows) {
    throw std::out_of_range("Wrong Matrix");
  } else {
    S21Matrix res(_rows, _cols);
    if (_rows == 1 && _cols == 1) {
      res(0, 0) = 1;
    } else {
      S21Matrix minor;
      double determ = 0;
      for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
          minor = GetMinorMatrix(i, j);
          determ = minor.Determinant();
          res(i, j) = determ * pow(-1, i + j);
        }
      }
    }
    return res;
  }
}

double S21Matrix::Determinant() {
  double res = 0;
  if (_p == NULL) {
    throw std::out_of_range("Null pointer at source");
  } else if (_cols != _rows || _cols <= 0) {
    throw std::out_of_range("It's not cube matrix");
  } else if (_cols == 1) {
    res = _p[0][0];
  } else if (_cols == 2) {
    res = _p[0][0] * _p[1][1] - _p[0][1] * _p[1][0];
  } else {
    int degree = 1;
    res = 0;
    for (int j = 0; j < _cols; j++) {
      double recursion_res = 0;
      S21Matrix recursion_matrix = GetMinorMatrix(0, j);
      recursion_res = recursion_matrix.Determinant();
      res += degree * _p[0][j] * recursion_res;
      degree = -degree;
    }
  }
  return res;
}

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix res;
  if (_p == NULL) {
    throw std::out_of_range("Null pointer at source");
  } else if (_rows != _cols) {
    throw std::out_of_range("It's not cube matrix");
  } else if (_rows == 1 && _cols == 1) {
    if (_p[0][0] == 0) {
      throw std::out_of_range("Determinant is 0");
    } else {
      S21Matrix buf_res(1, 1);
      buf_res._p[0][0] = 1 / _p[0][0];
      res = buf_res;
    }
  } else {
    double Deter = Determinant();
    if (Deter == 0) {
      throw std::out_of_range("Determinant is 0");
    } else {
      S21Matrix buf_calc_compl = CalcComplements();
      buf_calc_compl = buf_calc_compl.Transpose();
      buf_calc_compl.MulNumber(1 / Deter);
      res = buf_calc_compl;
    }
  }
  return res;
}

//////////////////// ОПЕРАТОРЫ ///////////////////

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  // creating result matrix
  S21Matrix res(_rows, _cols);
  res.SumMatrix(other);
  res.SumMatrix(*this);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  // creating result matrix
  S21Matrix res(_rows, _cols);
  res.SumMatrix(*this);
  res.SubMatrix(other);
  return res;
}

// ============================================================ //

S21Matrix S21Matrix::operator*(const double num) { return num * *this; }

S21Matrix operator*(const double num, const S21Matrix& other) {
  S21Matrix res(other._rows, other._cols);
  res = other;
  res.MulNumber(num);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res = *this;
  res.MulMatrix(other);
  return res;
}

// ============================================================ //

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

S21Matrix S21Matrix::operator=(const S21Matrix& other) {
  Free();
  _rows = other._rows;
  _cols = other._cols;
  mem();
  for (int i = 0; i < _rows; i++) {
    for (int k = 0; k < _cols; k++) {
      _p[i][k] = other._p[i][k];
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

double& S21Matrix::operator()(int row, int col) {
  if (row >= _rows || col >= _cols) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  return _p[row][col];
}

/////////////////// ACCESSTORS ///////////////////

int S21Matrix::GetRows() { return _rows; }

int S21Matrix::GetCols() { return _cols; }

/////////////////// MUTATORS ///////////////////

void S21Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::out_of_range("Rows cannot be <=0");
  }
  S21Matrix res(rows, _cols);
  for (int i = 0; i < (_rows < rows ? _rows : rows); i++) {
    for (int k = 0; k < _cols; k++) {
      res._p[i][k] = _p[i][k];
    }
  }
  *this = res;
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::out_of_range("Rows cannot be <=0");
  }
  S21Matrix res(_rows, cols);
  for (int i = 0; i < _rows; i++) {
    for (int k = 0; k < (_cols < cols ? _cols : cols); k++) {
      res._p[i][k] = _p[i][k];
    }
  }
  *this = res;
}

/////////////////// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ///////////////////

void S21Matrix::mem() {
  _p = new double*[_rows];
  for (int i = 0; i < _rows; i++) _p[i] = new double[_cols];
  Fill_by_number(0.0);
}

void S21Matrix::Free() {
  if (_p) {
    for (int i = 0; i < _rows; i++) delete[] _p[i];
    delete[] _p;
    _p = nullptr;
    _rows = 0;
    _cols = 0;
  }
}

S21Matrix S21Matrix::GetMinorMatrix(int rows, int cols) {
  if (_p == NULL || _rows <= 0 || _cols <= 0) {
    throw std::out_of_range(
        "Incorrect input, source matrices should have not null pointer and "
        "cols and rows > 0");
  } else {
    S21Matrix res(_rows - 1, _cols - 1);
    int counter_minor_rows = 0;
    int counter_minor_cols = 0;
    for (int i = 0; i < _rows; i++) {
      for (int j = 0; j < _cols; j++) {
        if (i != rows && j != cols) {
          res._p[counter_minor_rows][counter_minor_cols] = _p[i][j];
          ++counter_minor_cols;
          if (counter_minor_cols == res._cols) {
            counter_minor_cols = 0;
            ++counter_minor_rows;
          }
        }
      }
    }
    return res;
  }
}

void S21Matrix::Fill_by_number(double n) {
  for (int i = 0; i < _rows; i++) {
    for (int k = 0; k < _cols; k++) {
      _p[i][k] = n;
    }
  }
}
