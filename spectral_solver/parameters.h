//
//  parameters.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef spectral_solver_parameters_h
#define spectral_solver_parameters_h

//#define PUNISHTERM //Add punish term can sometimes speed up
//#define MULTIPROCESS

const int Ntheta=64;
const int Nr=82;
const double StepT=1.0/16384/32;
const int NumField=1;
const double radius=1; //5;

#ifdef PUNISHTERM
const double punish=50;
#endif

//for solver
const double tolerance=1e-10;


const int logicNr=Nr-1;
const int Nrp=Nr/2;
const double PI=3.14159265358979323846264338328;
const int matrixH=Nrp*NumField;
const int matrixW=Ntheta*NumField;
const int NumPoints=Nrp*Ntheta;


#endif
