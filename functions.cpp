#include "functions.h"

Image<byte> comb(const Image<byte>& I1, double alpha, const Image<byte>& I2){
    Image<byte> I(I1.width(), I1.height());
    for ( int x = 0 ; x<I.width() ; x++)
        for(int y = 0; y<I.height(); y++)
            I(x,y) = (byte) (int((1-alpha)*I1(x,y) + alpha*I2(x,y)));
    return I;
}

//transforms an image into a matrix
Matrix<double> transform (const Image<byte>& I)
{
    Matrix<double> M(I.height(), I.width());
    for ( int i = 0 ; i<I.height() ; i++)
        for(int j = 0; j<I.width(); j++)
            M(i,j) = (double)I(j,i)/255;
    return M;
}

//prints a 2D (scaled) vector on active window from pixel p
void show_vector(const Vector<double>& v, const IntPoint2& p, double scale, Color col )
{
    IntPoint2 p2(int(scale*v[0]), int(scale*v[1]));
    drawArrow(p,p+p2,col);
}

//cuts douwn a sub image
Image<byte> round_img( const Image<byte>& I, const IntPoint2& p, int width, int height)
{
    Image<byte> I_new(width, height);
    for(int x = p.x() ; x< p.x() + width ; x++)
    {
        for(int y = p.y() ; y<p.y() + height ; y++)
        {
            I_new(x-p.x(),y-p.y()) = I(x,y);
        }
    }
    return I_new;
}


//selects a sub image using the mouse
void select(Image<byte> I, IntPoint2& p, int& w, int& h)
{
    Window window = openWindow(I.width(),I.height(),"Select a sub image : 2 left clicks with the mouse");
    display(I);

    IntPoint2 p1,p2,Pmin,Pmax,p3;

    while(true)
    {
        if(getMouse(p3)==3)
            break;
        display(I);
        p1 = p3;
        getMouse(p2);
        Pmin.x() = min(p1.x(),p2.x());
        Pmin.y() = min(p1.y(),p2.y());
        Pmax.x() = max(p1.x(),p2.x());
        Pmax.y() = max(p1.y(),p2.y());
        w = Pmax.x() - Pmin.x();
        h = Pmax.y() - Pmin.y();
        drawRect(Pmin, w, h, GREEN);
    }
    closeWindow(window);
    p = Pmin;
}

//prints the volocity vector on the selected spot using Lucas-Kanade's method (with distance d to neighbours)
void show_velocity_LK(Image<byte> I1_rounded, Image<byte> I2_rounded, int d, int scale, Color col)
{
    Matrix<double> I1 = transform(I1_rounded);
    Matrix<double> I2 = transform(I2_rounded);

    Window w1 = openWindow(I1_rounded.width(), I1_rounded.height(), "t = t1, click here");
    setActiveWindow(w1);
    //display(I1_rounded);
    display(comb(I1_rounded,0.2,I2_rounded));

    Window w2 = openWindow(I2_rounded.width(), I2_rounded.height(), "t = t1 + dt");
    setActiveWindow(w2);
    display(I2_rounded);

    IntPoint2 p;

    while(true)
    {
        setActiveWindow(w1);
        if(getMouse(p) == 3)
            break;
        Vector<double> v = V(I1,I2, p, d);
        show_vector(v, p, scale, col);
        setActiveWindow(w2);
        show_vector(v, p, scale, col);
    }
    closeWindow(w1);
    closeWindow(w2);
}

//prints a diagram of the volocity field using Lucas-Kanade's method
void velocity_vector_field_LK(Image<byte> I1_rounded, Image<byte> I2_rounded, int d, int density, int scale, Color col)
{
    Matrix<double> I1 = transform(I1_rounded);
    Matrix<double> I2 = transform(I2_rounded);

    Window w1 = openWindow(I1_rounded.width(), I1_rounded.height(), "Image");
    setActiveWindow(w1);
    display(I1_rounded);

    Window w2 = openWindow(I1_rounded.width(), I1_rounded.height(), "Velocity Vector Field : Lucas-Kanade method");
    setActiveWindow(w2);

    for(int x = 0 ; x<I1_rounded.width() ; x+=density)
    {
        for(int y = 0 ; y<I1_rounded.height() ; y+=density)
        {
            IntPoint2 p(x,y);
            Vector<double> v = V(I1,I2, p, d);
            show_vector(v, p, scale, col);
        }
    }
    if(click()==3)
    {
        closeWindow(w1);
        closeWindow(w2);
    }
}

void LK_colorMap(Image<byte> I1_rounded, Image<byte> I2_rounded, int d)
{
    Matrix<double> I1 = transform(I1_rounded);
    Matrix<double> I2 = transform(I2_rounded);

    Window w1 = openWindow(I1_rounded.width(), I1_rounded.height(), "Image");
    setActiveWindow(w1);
    display(I1_rounded);

    Window w2 = openWindow(I1_rounded.width(), I1_rounded.height(), "Velocity Vector Field : Lucas-Kanade method");
    setActiveWindow(w2);
    Image<byte> CMap(I1_rounded.width(), I1_rounded.height());

    for(int x = 0 ; x<I1_rounded.width() ; x++)
    {
        for(int y = 0 ; y<I1_rounded.height() ; y++)
        {
            IntPoint2 p(x,y);
            Vector<double> v = V(I1,I2, p, d);
            CMap(x,y) = (byte)(int(255-255*sqrt(v[0]*v[0]+v[1]*v[1])/10));
        }
    }
    display(CMap);

    if(click()==3)
    {
        closeWindow(w1);
        closeWindow(w2);
    }
}

//prints a diagram of the volocity field using Horn-Schunck's method starting from a null vector field
void velocity_vector_field_HS(Image<byte> I1_rounded, Image<byte> I2_rounded,
                              double lambda, int n, int density, int scale, Color col)
{
    Matrix<double> I1 = transform(I1_rounded);
    Matrix<double> I2 = transform(I2_rounded);

    Window w1 = openWindow(I1_rounded.width(), I1_rounded.height(), "Image");
    setActiveWindow(w1);
    display(I1_rounded);

    Window w2 = openWindow(I1_rounded.width(), I1_rounded.height(), "Velocity Vector Field : Horn-Schunck method");
    setActiveWindow(w2);

    Matrix<double> Vx(I1.nrow(), I1.ncol());
    Vx.fill(0);
    Matrix<double> Vy(I1.nrow(), I1.ncol());
    Vy.fill(0);

    //HS(n, lambda, I1, I2, Vx, Vy);
    HS2(n, 0.0001 , 2, lambda, I1, I2, Vx, Vy);


    for(int x = 0; x<Vx.ncol() ; x+=density)
        for(int y = 0; y<Vx.nrow(); y+=density)
        {
            IntPoint2 p(x,y);
            IntPoint2 p2(int(scale*Vx(y,x)), int(scale*Vy(y,x)));
            drawArrow(p,p+p2,col);
        }

    if(click()==3)
    {
        closeWindow(w1);
        closeWindow(w2);
    }

}

void HS_colorMap(Image<byte> I1_rounded, Image<byte> I2_rounded,double lambda, int n)
{
    Matrix<double> I1 = transform(I1_rounded);
    Matrix<double> I2 = transform(I2_rounded);

    Window w1 = openWindow(I1_rounded.width(), I1_rounded.height(), "Image");
    setActiveWindow(w1);
    display(I1_rounded);

    Window w2 = openWindow(I1_rounded.width(), I1_rounded.height(), "Velocity Vector Field : Horn-Schunck method");
    setActiveWindow(w2);
    Image<byte> CMap(I1_rounded.width(), I1_rounded.height());


    Matrix<double> Vx(I1.nrow(), I1.ncol());
    Vx.fill(0);
    Matrix<double> Vy(I1.nrow(), I1.ncol());
    Vy.fill(0);

    //HS(n, lambda, I1, I2, Vx, Vy);
    HS2(n, 0.0001 , 2, lambda, I1, I2, Vx, Vy);

    for(int x = 0; x<Vx.ncol() ; x++)
        for(int y = 0; y<Vx.nrow(); y++)
        {
            CMap(x,y) = (byte)(int(255-255*10*sqrt(Vx(y,x)*Vx(y,x)+Vy(y,x)*Vy(y,x))/10));
        }
    display(CMap);

    if(click()==3)
    {
        closeWindow(w1);
        closeWindow(w2);
    }

}
