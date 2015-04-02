//
//  matrix.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/4/2.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef __spectral_solver__matrix__
#define __spectral_solver__matrix__

//#include <stdio.h>

template <typename DType>
class matrix
{
public:
    matrix();
    matrix(const int size1, const int size2);
    matrix(const matrix<DType>& cp);
    
    ~matrix();
    
    DType& operator[](const int idx);
    matrix& operator=(const matrix<DType>& cp);
    void operator*=(const double factor);
    void operator/=(const double factor);
    
    void eleMultiply(const matrix<DType>& matrix);
    void eleDivide(const matrix<DType>& matrix);
    DType& ele(const int x, const int y);
    
    DType* data;
    int sizex;
    int sizey;
    int NumE;
};

#endif /* defined(__spectral_solver__matrix__) */
