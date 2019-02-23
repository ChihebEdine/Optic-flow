#include "horn_schunck.h"
#include "lucas_kanade.h"

//calculates the mean of V in p
double norm2 (const Matrix<double>& Vx, const Matrix<double>& Vy)
{
    double N = 0;
    for(int x = 0; x<Vx.ncol() ; x++)
        for(int y = 0; y<Vx.nrow(); y++)
            N += Vx(y,x)*Vx(y,x) + Vy(y,x)*Vy(y,x);

    return sqrt(N);
}


double mean(const IntPoint2& p,  const Matrix<double>& V )
{
    int x1,x2,y1,y2;

    if(p.x() >= 2)
        x1 = p.x() - 2;
    else if (p.x() == 1)
        x1 = p.x() - 1;
    else
        x1 = p.x();

    if(p.x() <= V.ncol()- 3)
        x2 = p.x() + 2;
    else if (p.x() == V.ncol()- 2)
        x2 = p.x() + 1;
    else
        x2 = p.x();

    if(p.y() >= 2)
        y1 = p.y() - 2;
    else if (p.y() == 1)
        y1 = p.y() - 1;
    else
        y1 = p.y();

    if(p.y() <= V.nrow()- 3)
        y2 = p.y() + 2;
    else if (p.y() == V.nrow()- 2)
        y2 = p.y() + 1;
    else
        y2 = p.y();

    return (1/4)*(V(y1, p.x()) + V(y2,p.x())+ V(p.y(),x1) + V(p.y(),x2));
}


//calculates the fields Vx et Vy after n iterations
void HS(int n, double lambda, const Matrix<double>& I1, const Matrix<double>& I2, Matrix<double>& Vx, Matrix<double>& Vy)
{

    Matrix<double> Vx_old(I1.nrow(),I1.ncol());
    Vx_old.fill(0);
    Vx = Vx_old.clone();

    Matrix<double> Vy_old(I1.nrow(),I1.ncol());
    Vy_old.fill(0);
    Vy = Vy_old.clone();

    for(int i = 0; i<n; i++)
    {
        for(int x = 0; x<I1.ncol() ; x++)
            for(int y = 0; y<I1.nrow(); y++)
            {
                IntPoint2 p(x,y);
                double dxI = dx(I1,p);
                double dyI = dy(I1,p);
                double dtI = dt(I1,I2,p);
                double Vxm = mean(p,Vx_old);
                double Vym = mean(p,Vy_old);
                double C = (Vxm*dxI + Vym*dyI + dtI)/(dxI*dxI + dyI*dyI + lambda);
                Vx(y,x) = Vxm - C * dxI;
                Vy(y,x) = Vym - C * dyI;
            }
        Vx_old = Vx.clone();
        Vy_old = Vy.clone();
    }
}


void HS2(int n, double e, int d, double lambda, const Matrix<double>& I1, const Matrix<double>& I2, Matrix<double>& Vx, Matrix<double>& Vy){

    Matrix<double> Vx_old(I1.nrow(),I1.ncol());
    Matrix<double> Vy_old(I1.nrow(),I1.ncol());

    //LK initial speed
    for(int x = 0 ; x<I1.ncol() ; x++)
    {
        for(int y = 0 ; y<I1.nrow() ; y++)
        {
            IntPoint2 p(x,y);
            Vector<double> v = V(I1,I2, p, d);
            Vx_old(y,x) = v[0];
            Vy_old(y,x) = v[1];
            cout<<v[0]<<','<<v[1]<<endl;
        }
    }

    Vx = Vx_old.clone();
    Vy = Vy_old.clone();
    int i = 0;

    while(i<n)
    {
        for(int x = 0; x<I1.ncol() ; x++)
            for(int y = 0; y<I1.nrow(); y++)
            {
                IntPoint2 p(x,y);
                double dxI = dx(I1,p);
                double dyI = dy(I1,p);
                double dtI = dt(I1,I2,p);
                double Vxm = mean(p,Vx_old);
                double Vym = mean(p,Vy_old);
                double C = (Vxm*dxI + Vym*dyI + dtI)/(dxI*dxI + dyI*dyI + lambda);
                Vx(y,x) = Vxm - C * dxI;
                Vy(y,x) = Vym - C * dyI;
            }

        if(norm2(Vx-Vx_old,Vy-Vy_old)<e)
            break;

        Vx_old = Vx.clone();
        Vy_old = Vy.clone();
        i++;
    }
}








