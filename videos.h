#pragma once

#include "functions.h"
#include <vector>


vector< Matrix<double> > Load_data(int n);


vector< Image<byte> > Optic_flow(const vector< Matrix<double> >& data, int method );


void Save_data(const vector< Image<byte> >& images);
