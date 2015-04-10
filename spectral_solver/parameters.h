//
//  parameters.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef spectral_solver_parameters_h
#define spectral_solver_parameters_h

 //Add punish term can sometimes speed up
//#define MULTIPROCESS
#define LINEAR_TEST_MODEL
//#define FU_MODEL
#include <cmath>


#ifdef FU_MODEL
#define PUNISHTERM
const int Ntheta=16;
const int Nr=81*5+1;
const long double iniStepT=1.0/16384/32;
const int NumField=3;
const long double radius=3;

const long double Alpha=2.08;
const long double Beta=2.08;
const long double Gamma=0.7;
const long double Kn2=100;
const long double Drho=0.0045*3.6;
const long double Drho0=0.0001*3.6;
const long double Dh=0.004*3.6;
const long double Dn=0.008*3.6;
const long double Kh=4;
const long double M=20;
#ifdef PUNISHTERM
const long double punish=50;
#endif
#endif


#ifdef LINEAR_TEST_MODEL
//#define PUNISHTERM
const int Ntheta=64;
const int Nr=82;
const long double iniStepT=1.0/16384/8;
const int NumField=3;
const long double radius=1;
#ifdef PUNISHTERM
const long double punish=50;
#endif
#endif



//for solver
const long double tolerance=1e-8;


const int logicNr=Nr-1;
const int Nrp=Nr/2;
const long double PI=3.1415926535897932384626433832795028841971693993751;
const int matrixH=Nrp*NumField;
const int matrixW=Ntheta*NumField;
const int NumPoints=Nrp*Ntheta;
const int totalPoints=matrixH*Ntheta;

const int aliasingr=2*floor(Nr/3);
const int aliasingt=floor(Ntheta/3);

const long double rp2=radius*radius;

const int IncreaseTimes=0;

#endif
