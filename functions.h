#pragma once

#include "lucas_kanade.h"
#include "horn_schunck.h"

//transforms an image into a matrix
Matrix<double> transform (const Image<byte> &I);

//prints a 2D (scaled) vector on active window from pixel p
void show_vector(const Vector<double>& v, const IntPoint2& p, double scale = 10, Color col = RED );

//cuts douwn a sub image
Image<byte> round_img( const Image<byte>& I, const IntPoint2& p, int width, int height);


//selects a sub image using the mouse
void select(Image<byte> I, IntPoint2& p, int& w, int& h);

//prints the volocity vector on the selected spot using Lucas-Kanade's method (with distance d to neighbours)
void show_velocity_LK(Image<byte> I1_rounded, Image<byte> I2_rounded, int d = 50, int scale = 10, Color col = RED);

//prints a diagram of the volocity field using Lucas-Kanade's method
void velocity_vector_field_LK(Image<byte> I1_rounded, Image<byte> I2_rounded, int d = 50, int density = 8, int scale = 10, Color col = RED);
void LK_colorMap(Image<byte> I1_rounded, Image<byte> I2_rounded, int d);

//prints a diagram of the volocity field using Horn-Schunck's method starting from a null vector field
void velocity_vector_field_HS(Image<byte> I1_rounded, Image<byte> I2_rounded,
                              double lambda = 0.1, int n = 20, int density = 8, int scale = 10, Color col = BLUE);
void HS_colorMap(Image<byte> I1_rounded, Image<byte> I2_rounded,double lambda, int n);
