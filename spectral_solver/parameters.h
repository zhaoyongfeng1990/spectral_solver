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
#include <cmath>

const int Ntheta=64;
const int Nr=81*5+1;
const double StepT=1.0/16384/32;
const int NumField=3;
const double radius=3;

const double Alpha=2.08;
const double Beta=2.08;
const double Gamma=0.7;
const double Kn2=100;
const double Drho=0.0045*3.6;
const double Drho0=0.0001*3.6;
const double Dh=0.004*3.6;
const double Dn=0.008*3.6;
const double Kh=4;
const double M=20;

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
const int totalPoints=matrixH*Ntheta;

const int aliasingr=2*floor(Nr/3);
const int aliasingt=floor(Ntheta/3);

const double rp2=1.0/radius/radius;

#endif
