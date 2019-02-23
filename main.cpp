#include "functions.h"
#include "videos.h"

int main()
{

    Image<byte> I_1, I_2;
    load(I_1,srcPath("frame2.png"));
    load(I_2,srcPath("frame3.png"));

    IntPoint2 p;
    int w,h;
    select(I_1,p,w,h);

    Image<byte> I1_rounded = round_img(I_1,p,w,h);
    Image<byte> I2_rounded = round_img(I_2,p,w,h);

    //Parameters
    int density = 8;
    int scale = 10;
    int n = 100; //horn-schunck iteration number
    double lambda = 0.1;
    int d = 5; //the distance used in lucas-kanade's method

    //Graphs
    //show_velocity_LK(I1_rounded, I2_rounded,d,scale,RED);
    //velocity_vector_field_LK(I1_rounded, I2_rounded,d,density,scale,RED);
    //LK_colorMap(I1_rounded, I2_rounded, d);


    //velocity_vector_field_HS(I1_rounded, I2_rounded,lambda,n,density,scale,BLUE);
    HS_colorMap(I1_rounded, I2_rounded,lambda, n);

    //rcon d
    //sortir une imge % à rcond pour voir les endroits ou la mat est inv
    // staitliserer

//================= Videos ===================================================
//    vector< Matrix<double> > data = Load_data(300);
//    vector< Image<byte> > images = Optic_flow(data, 1);
//    //Save_data(images);



    return 0;
}
