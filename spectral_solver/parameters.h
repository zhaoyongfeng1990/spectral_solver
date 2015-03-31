//
//  parameters.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef spectral_solver_parameters_h
#define spectral_solver_parameters_h

//#define LINEAR_TEST_MODEL
#define FU_MODEL
//#define INTERACTION_MODIFIED_FU

#include <cmath>

#ifdef FU_MODEL
//#define PUNISHTERM
const int Ntheta=16;
const int Nr=81*9+1;
const double StepT=1.0/16384/16;
const int NumField=1; //3;
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
const double punish=30;
#endif

#endif

#ifdef LINEAR_TEST_MODEL

//#define PUNISHTERM
const int Ntheta=64;
const int Nr=81+1;
const double StepT=1.0/16384/8;
const int NumField=3;
const double radius=1;

#ifdef PUNISHTERM
const double punish=100;
#endif

#endif

#ifdef INTERACTION_MODIFIED_FU
const int Ntheta=64;
const int Nr=81*5+1;
const double StepT=1.0/16384/16;
const int NumField=5;
const double radius=3;

const double Dh1=0.004*3.6;
const double Dh2=0.004*3.6;
const double Dn=0.008*3.6;
const double Drho=0.0001*3.6;
const double Kh1=4;
const double Kh2=4;
const double M=20;
const double Gamma1=0.9;
const double Gamma2=0.7;
const double Kn1=100;
const double Kn2=100;
const double Alpha1=2.08;
const double Alpha2=2.08;
const double Beta1=2.08;
const double Beta2=2.08;
const double Drho1=0;
const double Drho2=0;
const double KC1=0.0089*3.6*3;
const double KC2=1.5;
const double KD1=0.0089*3.6*3;
const double KD2=1.5;
#endif

//for solver
const double tolerance=1e-8;


const int logicNr=Nr-1;
const int Nrp=Nr/2;
const double PI=3.14159265358979323846264338328;
const int matrixH=Nrp*NumField;
const int matrixW=Ntheta*NumField;
const int NumPoints=Nrp*Ntheta;
const int totalPoints=matrixH*Ntheta;

const int aliasingr=Nr-4; //2*floor(Nr/3);
const int aliasingt=floor(Ntheta/3);

const double rp2=1.0/radius/radius;
#endif
