#pragma once
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include "iostream"

using namespace Imagine;
using namespace std;



// returns a vector containing all pixels with distance less or equal to d away from p on image I
vector<IntPoint2> W (const Matrix<double>& I, const IntPoint2& p, int d);

// returns the partial x derivative of an image I
double dx(const Matrix<double>& I, const IntPoint2&  p);

// returns the partial y derivative of an image I
double dy(const Matrix<double>& I, const IntPoint2&  p);

// returns the partial t derivative of an image I
double dt(const Matrix<double>& I1, const Matrix<double>& I2, const IntPoint2&  p);

// returns the matrix A of a pixel p
Matrix<double> A(const Matrix<double>& I, const IntPoint2& p, int d);

// returns the vector b of a pixel p
Matrix<double> b(const Matrix<double>& I1, const Matrix<double>& I2, const IntPoint2&  p, int d);

Vector<double> V(const Matrix<double>& I1, const Matrix<double>& I2, const IntPoint2& p, int d);
