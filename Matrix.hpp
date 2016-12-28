//
//  Matrix.hpp
//
//
//  Created by Ali Jbaily on 9/25/16.
//
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <string.h>
#include "ArrayX.cpp"

template <typename T>
class Matrix {
    ArrayX<T>* matrix;
    
public:
    static int defaultSize;
    int count = 0;
    void add(T, int);
    Matrix<T>();
    Matrix<T>(int);
    ~Matrix();
    //    Matrix(Matrix &obj);  // copy constructor
    ArrayX<T>& operator [](int);
    Matrix& operator=(const Matrix&);
    
private:
    void doubleSize(int);
};

#endif /* Matrix_hpp */


// 3 arrays (size, max, matrix)and 1 int (count)