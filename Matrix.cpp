//
//  Matrix.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 9/25/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "Matrix.hpp"

using namespace std;

template <typename T>
int Matrix<T>::defaultSize = 1;

template <typename T>
Matrix<T>::Matrix(){
    matrix = new ArrayX<T>[defaultSize];
    count = defaultSize;
}

template <typename T>
Matrix<T>::Matrix(int t) {
    matrix = new ArrayX<T>[t];
    count = t;
}

template <typename T>
Matrix<T>::~Matrix() {
    if (count > 0) {
        delete [] matrix;
        count = 0;
    }
}

template <typename T>
ArrayX<T>& Matrix<T>::operator [](int i) {
    return matrix[i];
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& obj) {
    if (this != &obj) {
        if (count > 0)
            delete [] matrix;
        count = obj.count;
        matrix = new ArrayX<T>[count];
        for (int i = 0; i < count; i++)
            matrix[i] = obj[i];
    }
    return *this;
}

template <typename T>
void Matrix<T>::add(T p, int index) {
    matrix[index].add(p);
}