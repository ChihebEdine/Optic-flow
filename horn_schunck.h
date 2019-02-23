#pragma once
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include "iostream"

using namespace Imagine;
using namespace std;

//calculates the mean of V in p
double mean(const IntPoint2& p,  const Matrix<double>& V );

//calculates the fields Vx et Vy after n iterations
void HS(int n, double lambda, const Matrix<double>& I1, const Matrix<double>& I2, Matrix<double>& Vx, Matrix<double>& Vy);

void HS2(int n, double e , int d, double lambda, const Matrix<double>& I1, const Matrix<double>& I2, Matrix<double>& Vx, Matrix<double>& Vy);
