#include "lucas_kanade.h"



// returns a vector containing all pixels with distance less or equal to d away from p on image I
vector<IntPoint2> W (const Matrix<double>& I, const IntPoint2& p, int d)
{
    vector<IntPoint2> VP;
    for(int r = p.y()-d; r<=p.y()+d; r++)
        for(int c = p.x(); c<=p.x()+d; c++)
            if( r>=0 && r<I.nrow() && c>=0 && c<I.ncol())
            {
                IntPoint2 b (c,r);
                VP.push_back(b);
            }
    return VP;
}


// returns the partial x derivative of an image I
double dx(const Matrix<double>& I, const IntPoint2&  p)
{
    int left = p.x()-1;
    int right = p.x()+1;
    if(left<0)
        left +=1;
    if(right>= I.ncol())
        right-=1;

    return (I(p.y(), right)-I(p.y(),left))/(right - left);
}

// returns the partial y derivative of an image I
double dy(const Matrix<double>& I, const IntPoint2&  p)
{
    int up = p.y()-1;
    int down = p.y()+1;
    if(up<0)
        up+=1;
    if(down>=I.nrow())
        down-=1;
    return (I(down,p.x())-I(up,p.x()))/(down-up);
}

// returns the partial t derivative of an image I
double dt(const Matrix<double>& I1, const Matrix<double>& I2, const IntPoint2&  p)
{
    return I2(p.y(),p.x())-I1(p.y(),p.x());
}


// returns the matrix A of a pixel p
Matrix<double> A(const Matrix<double>& I, const IntPoint2&  p, int d)
{
    vector<IntPoint2> WP = W(I, p, d);
    int n = WP.size();
    Matrix<double> AP(n,2);
    for(int r=0; r<n; r++)
    {
        IntPoint2 p = WP[r];
        AP(r,0) = dx(I,p);
        AP(r,1) = dy(I,p);
    }
    return AP;
}


// returns the vector b of a pixel p
Matrix<double> b(const Matrix<double>& I1, const Matrix<double>& I2, const IntPoint2&  p, int d)
{
    vector<IntPoint2> WP = W(I1, p, d);
    int n = WP.size();
    Matrix<double> bP(n,1);
    for(int r=0; r<n; r++)
    {
        IntPoint2 p = WP[r];
        bP(r,0) = - dt(I1,I2,p);
    }
    return bP;
}

// returns the velocity field of two consecutive frames I1 and I2 using the Lucas-Kanade's method
Vector<double> V(const Matrix<double>& I1, const Matrix<double>& I2, const IntPoint2& p, int d)
{
    Matrix<double> A_P = A(I1,p,d);
    Matrix<double> b_P = b(I1,I2,p,d);
    Vector<double> V_P = linSolve(A_P,b_P);
    if(isnan(V_P[0]) ||  isnan(V_P[1]) || abs(V_P[0]) == INFINITY || abs(V_P[1]) == INFINITY || sqrt(V_P[0]*V_P[0]+V_P[1]*V_P[1])>10)
    {
        V_P[0] = 0;
        V_P[1] = 0;

    }
    return V_P;

}
