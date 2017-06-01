//
//  l1_stat.hpp
//  osc
//
//  Created by Futao Zhang on 23/12/2016.
//  Copyright © 2016 Futao Zhang. All rights reserved.
//

#ifndef l1_stat_hpp
#define l1_stat_hpp

#include "l0_stat.h"
#endif /* l1_stat_hpp */

MatrixXd reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table = false);
void lin(VectorXd &y, MatrixXd &C, VectorXd &x, vector<double> &rst);
double var(VectorXd &x);
void cov(MatrixXd &mat, MatrixXd &covar);
