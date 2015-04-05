//
//  matrix.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/4/2.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "matrix.h"

template <typename DType>
matrix<DType>::matrix()
{
    sizex=0;
    sizey=0;
    NumE=0;
    data=0;
}

template <typename DType>
matrix<DType>::matrix(const int size1, const int size2)
{
    sizex=size2;
    sizey=size1;
    NumE=size1*size2;
    data=new DType[NumE];
}

template <typename DType>
matrix<DType>::matrix(const matrix<DType>& cp)
{
    sizex=cp.sizex;
    sizey=cp.sizey;
    NumE=sizex*sizey;
    data=new DType[NumE];
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=cp.data[iter];
    }
}

template <typename DType>
matrix<DType>::~matrix()
{
    if (data!=0)
    {
        delete [] data;
    }
}