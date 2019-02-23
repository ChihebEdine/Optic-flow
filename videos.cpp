#include "videos.h"


vector< Matrix<double> > Load_data(int n) {
    vector< Matrix<double> >  data;
    for(int i = 0; i<n; i++) {
        Image<byte> I;
        load(I,srcPath("frame"+  to_string(i) +".png"));
        data.push_back(transform(I));
    }
    return data;

}

Image<byte> makeImage(byte* II, int w, int h){
    Image<byte> I(w,h);
    for(int x = 0 ; x<w ; x++)
        for(int y = 0 ; y<h ; y++){
            I(x,y) = II[y*w+x];
        }
    return I;
}


Image<byte> single_optic_flow(const Window& W,  const Matrix<double>& M1, const Matrix<double>& M2, int method){
    byte* I;
    int w = 0;
    int h = 0;

    int density = 8;
    int scale = 1;
    Color col = BLACK;

    int d = 3;

    int n = 10;
    double lambda = 0.01;

    setActiveWindow(W);

    if(method == 0) {
        // Lucas Kanade

        for(int x = 0 ; x<M1.ncol() ; x+=density)
        {
            for(int y = 0 ; y<M1.nrow() ; y+=density)
            {
                IntPoint2 p(x,y);
                Vector<double> v = V(M1,M2, p, d);
                show_vector(v, p, scale, col);
            }
        }
        captureWindow(I,w,h);
        clearWindow();
        return(makeImage(I,w,h));
    }

    else {
        // Horn Schunk
        Matrix<double> Vx(M1.nrow(), M1.ncol());
        Vx.fill(0);
        Matrix<double> Vy(M1.nrow(), M1.ncol());
        Vy.fill(0);

        HS(n, lambda, M1, M2, Vx, Vy);

        for(int x = 0; x<Vx.ncol() ; x+=density)
            for(int y = 0; y<Vx.nrow(); y+=density)
            {
                IntPoint2 p(x,y);
                IntPoint2 p2(int(scale*Vx(y,x)), int(scale*Vy(y,x)));
                drawArrow(p,p+p2,col);
            }
        captureWindow(I,w,h);
        clearWindow();
        return(makeImage(I,w,h));
    }

}



vector< Image<byte> > Optic_flow(const vector< Matrix<double> >& data, int method ) {
    vector< Image<byte> > OpticFlow;
    Window w = openWindow(data[0].ncol(), data[0].nrow(), "Optic Flow");

    for(int i = 0 ; i<data.size()-1; i++){
        Image<byte> I  = single_optic_flow(w,data[i],data[i+1], method);
        OpticFlow.push_back(I);
        save(I, "Oframe" + to_string(i) + ".png");
    }

    closeWindow(w);

    return OpticFlow;
}


void Save_data(const vector< Image<byte> >& images) {
    for(int i = 0 ; i<images.size() ; i++) {
        save(images[i], "Oframe" + to_string(i) + ".png");
    }
}
