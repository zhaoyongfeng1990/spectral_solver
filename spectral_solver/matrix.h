//
//  matrix.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/4/2.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef __spectral_solver__matrix__
#define __spectral_solver__matrix__

#define CHECK_FOR_ROBUSTNESS

template <typename DType>
struct complex
{
    DType data[2];
};

template <typename DType>
class math_vector
{
public:
    math_vector();
    math_vector(int size1);
    math_vector(math_vector<DType> const& cp);
    
    ~math_vector();
    
    /*inline*/ void alloc(int size1);
    
    /*inline*/ DType& operator[](int idx);
    /*inline*/ math_vector<DType>& operator=(math_vector<DType> const& cp);
    
    /*inline*/ void operator*=(DType factor);
    /*inline*/ void operator/=(DType factor);
    /*inline*/ void operator*=(int factor);
    /*inline*/ void operator/=(int factor);
    /*inline*/ void operator*=(double factor);
    /*inline*/ void operator/=(double factor);
    
    /*inline*/ void operator+=(math_vector<DType> const& matrixB);
    /*inline*/ void operator-=(math_vector<DType> const& matrixB);
    
    /*inline*/ void eleMultiply(math_vector<DType> const& matrixB);
    /*inline*/ void eleDivide(math_vector<DType> const& matrixB);
    
    /*inline*/ void add(math_vector<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void minus(math_vector<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void eleMultiply(math_vector<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void eleDivide(math_vector<DType> const& mA, math_vector<DType> const& mB);
    
    DType* data;
    int size;
};

template <typename DType>
class math_vector_complex
{
public:
    math_vector_complex();
    math_vector_complex(int size1);
    math_vector_complex(math_vector<DType> const& cp);
    math_vector_complex(math_vector_complex<DType> const& cp);
    math_vector_complex(math_vector<DType> const& Re, math_vector<DType> const& Im);
    
    ~math_vector_complex();
    
    /*inline*/ void alloc(int size1);
    
    /*inline*/ DType& operator[](int idx);
    /*inline*/ math_vector_complex<DType>& operator=(math_vector<DType> const& cp);
    /*inline*/ math_vector_complex<DType>& operator=(math_vector_complex<DType> const& cp);
    
    /*inline*/ void operator*=(DType factor);
    /*inline*/ void operator/=(DType factor);
    /*inline*/ void operator*=(int factor);
    /*inline*/ void operator/=(int factor);
    /*inline*/ void operator*=(double factor);
    /*inline*/ void operator/=(double factor);
    
    /*inline*/ void operator*=(const complex<DType>& factor);
    /*inline*/ void operator/=(const complex<DType>& factor);
    
    /*inline*/ void operator+=(math_vector_complex<DType> const& matrixB);
    /*inline*/ void operator-=(math_vector_complex<DType> const& matrixB);
    /*inline*/ void operator+=(math_vector<DType> const& matrixB);
    /*inline*/ void operator-=(math_vector<DType> const& matrixB);
    
    /*inline*/ void eleMultiply(math_vector<DType> const& matrixB);
    /*inline*/ void eleDivide(math_vector<DType> const& matrixB);
    /*inline*/ void eleMultiply(math_vector_complex<DType> const& matrixB);
    /*inline*/ void eleDivide(math_vector_complex<DType> const& matrixB);
    
    /*inline*/ void add(math_vector<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void minus(math_vector<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void eleMultiply(math_vector<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void eleDivide(math_vector<DType> const& mA, math_vector<DType> const& mB);
    
    /*inline*/ void add(math_vector<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void minus(math_vector<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void minus(math_vector_complex<DType> const& mA, math_vector<DType> const& mB);
    /*inline*/ void eleMultiply(math_vector<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void eleDivide(math_vector<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void eleDivide(math_vector_complex<DType> const& mA, math_vector<DType> const& mB);
    
    /*inline*/ void add(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void minus(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void eleMultiply(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB);
    /*inline*/ void eleDivide(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB);
    
    /*inline*/ complex<DType> ele(int idx);
    /*inline*/ DType& ReEle(int idx);
    /*inline*/ DType& ImEle(int idx);
    
    DType* data;
    int size;
    int NumD;
};

template <typename DType>
class matrix
{
public:
    matrix();
    matrix(int size1, int size2);
    matrix(const matrix<DType>& cp);
    
    ~matrix();
    
    /*inline*/ void alloc(int size1, int size2);
    
    /*inline*/ DType& operator[](int idx);
    /*inline*/ matrix<DType>& operator=(const matrix<DType>& cp);
    
    /*inline*/ void operator*=(DType factor);
    /*inline*/ void operator/=(DType factor);
    /*inline*/ void operator*=(int factor);
    /*inline*/ void operator/=(int factor);
    /*inline*/ void operator*=(double factor);
    /*inline*/ void operator/=(double factor);
    
    /*inline*/ void operator+=(const matrix<DType>& matrixB);
    /*inline*/ void operator-=(const matrix<DType>& matrixB);
    
    /*inline*/ void eleMultiply(const matrix<DType>& matrixB);
    /*inline*/ void eleDivide(const matrix<DType>& matrixB);
    
    /*inline*/ void add(const matrix<DType>& mA, const matrix<DType>& mB);
    /*inline*/ void minus(const matrix<DType>& mA, const matrix<DType>& mB);
    /*inline*/ void eleMultiply(const matrix<DType>& mA, const matrix<DType>& mB);
    /*inline*/ void eleDivide(const matrix<DType>& mA, const matrix<DType>& mB);
    
    /*inline*/ DType& ele(int y, int x);
    
    /*inline*/ void setRow(int Idx, DType number);
    /*inline*/ void setCol(int Idx, DType number);
    /*inline*/ void scaleRow(int Idx, DType number);
    /*inline*/ void scaleCol(int Idx, DType number);
    /*inline*/ void allSet(DType number);
    
    DType* data;
    int sizex;
    int sizey;
    int NumE;
};

template <typename DType>
class matrix_complex
{
public:
    matrix_complex();
    matrix_complex(int size1, int size2);
    matrix_complex(const matrix<DType>& cp);
    matrix_complex(matrix_complex<DType> const& cp);
    matrix_complex(const matrix<DType>& Re, const matrix<DType>& Im);
    
    ~matrix_complex();
    
    /*inline*/ void alloc(int size1, int size2);
    
    /*inline*/ DType& operator[](int idx);
    /*inline*/ matrix_complex<DType>& operator=(const matrix<DType>& cp);
    /*inline*/ matrix_complex<DType>& operator=(matrix_complex<DType> const& cp);
    
    /*inline*/ void operator*=(DType factor);
    /*inline*/ void operator/=(DType factor);
    /*inline*/ void operator*=(int factor);
    /*inline*/ void operator/=(int factor);
    /*inline*/ void operator*=(double factor);
    /*inline*/ void operator/=(double factor);
    
    /*inline*/ void operator*=(const complex<DType>& factor);
    /*inline*/ void operator/=(const complex<DType>& factor);
    
    /*inline*/ void operator+=(matrix_complex<DType> const& matrixB);
    /*inline*/ void operator-=(matrix_complex<DType> const& matrixB);
    /*inline*/ void operator+=(const matrix<DType>& matrixB);
    /*inline*/ void operator-=(const matrix<DType>& matrixB);
    
    /*inline*/ void eleMultiply(const matrix<DType>& matrixB);
    /*inline*/ void eleDivide(const matrix<DType>& matrixB);
    /*inline*/ void eleMultiply(matrix_complex<DType> const& matrixB);
    /*inline*/ void eleDivide(matrix_complex<DType> const& matrixB);
    
    /*inline*/ void add(const matrix<DType>& mA, const matrix<DType>& mB);
    /*inline*/ void minus(const matrix<DType>& mA, const matrix<DType>& mB);
    /*inline*/ void eleMultiply(const matrix<DType>& mA, const matrix<DType>& mB);
    /*inline*/ void eleDivide(const matrix<DType>& mA, const matrix<DType>& mB);
    
    /*inline*/ void add(const matrix<DType>& mA, matrix_complex<DType> const& mB);
    /*inline*/ void minus(const matrix<DType>& mA, matrix_complex<DType> const& mB);
    /*inline*/ void minus(matrix_complex<DType> const& mA, const matrix<DType>& mB);
    /*inline*/ void eleMultiply(const matrix<DType>& mA, matrix_complex<DType> const& mB);
    /*inline*/ void eleDivide(const matrix<DType>& mA, matrix_complex<DType> const& mB);
    /*inline*/ void eleDivide(matrix_complex<DType> const& mA, const matrix<DType>& mB);
    
    /*inline*/ void add(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB);
    /*inline*/ void minus(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB);
    /*inline*/ void eleMultiply(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB);
    /*inline*/ void eleDivide(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB);
    
    /*inline*/ complex<DType> ele(int y, int x);
    /*inline*/ DType& ReEle(int y, int x);
    /*inline*/ DType& ImEle(int y, int x);
    
    /*inline*/ void setRow(int Idx, DType Re, DType Im);
    /*inline*/ void setCol(int Idx, DType Re, DType Im);
    /*inline*/ void scaleRow(int Idx, DType Re, DType Im);
    /*inline*/ void scaleCol(int Idx, DType Re, DType Im);
    
    DType* data;
    int sizex;
    int sizey;
    int NumE;
    
    int width;
    int NumD;
};

template <typename DType>
matrix<DType>::matrix()
{
    sizex=0;
    sizey=0;
    NumE=0;
    data=0;
}

template <typename DType>
matrix<DType>::matrix(int size1, int size2)
{
    sizex=size2;
    sizey=size1;
    NumE=size1*size2;
    data=new DType[NumE];
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::alloc(int size1, int size2)
{
    sizex=size2;
    sizey=size1;
    NumE=size1*size2;
    if (data!=0)
    {
        delete [] data;
    }
    data=new DType[NumE];
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
matrix<DType>::matrix(matrix<DType> const& cp)
{
    sizex=cp.sizex;
    sizey=cp.sizey;
    NumE=cp.NumE;
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

template <typename DType>
/*inline*/ DType& matrix<DType>::operator[](int idx)
{
    return data[idx];
}

template <typename DType>
/*inline*/ matrix<DType>& matrix<DType>::operator=(matrix<DType> const& cp)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=cp.sizex;
        sizey=cp.sizey;
        NumE=cp.NumE;
        data=new DType[NumE];
    }
    else if (sizex!=cp.sizex || sizey!=cp.sizey)
    {
        delete [] data;
        sizex=cp.sizex;
        sizey=cp.sizey;
        NumE=cp.NumE;
        data=new DType[NumE];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=cp.data[iter];
    }
    return *this;
}

template <typename DType>
/*inline*/ void matrix<DType>::operator*=(DType factor)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator/=(DType factor)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator*=(int factor)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator/=(int factor)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator*=(double factor)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator/=(double factor)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator+=(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]+=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::operator-=(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]-=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::eleMultiply(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]*=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::eleDivide(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]/=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::add(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=mA.data[iter]+mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::minus(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=mA.data[iter]-mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::eleMultiply(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=mA.data[iter]*mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::eleDivide(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        data=new DType[NumE];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=mA.data[iter]/mB.data[iter];
    }
}

template <typename DType>
/*inline*/ DType& matrix<DType>::ele(int y, int x)
{
    return data[x+y*sizex];
}

template <typename DType>
/*inline*/ void matrix<DType>::setRow(int Idx,DType number)
{
    for (int iter=Idx*sizex; iter<Idx*sizex+sizex; ++iter)
    {
        data[iter]=number;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::setCol(int Idx,DType number)
{
    for (int iter=Idx; iter<NumE+Idx; iter+=sizex)
    {
        data[iter]=number;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::scaleRow(int Idx,DType number)
{
    for (int iter=Idx*sizex; iter<Idx*sizex+sizex; ++iter)
    {
        data[iter]*=number;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::scaleCol(int Idx,DType number)
{
    for (int iter=Idx; iter<NumE+Idx; iter+=sizex)
    {
        data[iter]*=number;
    }
}

template <typename DType>
/*inline*/ void matrix<DType>::allSet(DType number)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter]=number;
    }
}

template <typename DType>
matrix_complex<DType>::matrix_complex()
{
    sizex=0;
    sizey=0;
    NumE=0;
    width=0;
    NumD=0;
    data=0;
}

template <typename DType>
matrix_complex<DType>::matrix_complex(int size1, int size2)
{
    sizex=size2;
    sizey=size1;
    NumE=size1*size2;
    width=sizex*2;
    NumD=2*NumE;
    data=new DType[NumD];
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::alloc(int size1, int size2)
{
    sizex=size2;
    sizey=size1;
    NumE=size1*size2;
    width=sizex*2;
    NumD=2*NumE;
    if (data!=0)
    {
        delete [] data;
    }
    data=new DType[NumD];
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
matrix_complex<DType>::matrix_complex(matrix<DType> const& cp)
{
    sizex=cp.sizex;
    sizey=cp.sizey;
    NumE=cp.NumE;
    width=sizex*2;
    NumD=2*NumE;
    data=new DType[NumD];
    for (int iter=0; iter<NumE; ++iter)
    {
        data[2*iter]=cp.data[iter];
    }
}

template <typename DType>
matrix_complex<DType>::matrix_complex(matrix<DType> const& Re, matrix<DType> const& Im)
{
    sizex=Re.sizex;
    sizey=Re.sizey;
    NumE=Re.NumE;
    width=sizex*2;
    NumD=2*NumE;
    data=new DType[NumD];
    for (int iter=0; iter<NumE; ++iter)
    {
        data[2*iter]=Re.data[iter];
        data[2*iter+1]=Im.data[iter];
    }
}

template <typename DType>
matrix_complex<DType>::matrix_complex(matrix_complex<DType> const& cp)
{
    sizex=cp.sizex;
    sizey=cp.sizey;
    NumE=cp.NumE;
    width=cp.width;
    NumD=cp.NumD;
    data=new DType[NumD];
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=cp.data[iter];
    }
}

template <typename DType>
matrix_complex<DType>::~matrix_complex()
{
    if (data!=0)
    {
        delete [] data;
    }
}

template <typename DType>
/*inline*/ DType& matrix_complex<DType>::operator[](int idx)
{
    return data[idx];
}

template <typename DType>
/*inline*/ matrix_complex<DType>& matrix_complex<DType>::operator=(matrix<DType> const& cp)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=cp.sizex;
        sizey=cp.sizey;
        NumE=cp.NumE;
        width=cp.width;
        NumD=cp.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=cp.sizex || sizey!=cp.sizey)
    {
        delete [] data;
        sizex=cp.sizex;
        sizey=cp.sizey;
        NumE=cp.NumE;
        width=cp.width;
        NumD=cp.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=cp.data[iter];
    }
    return *this;
}

template <typename DType>
/*inline*/ matrix_complex<DType>& matrix_complex<DType>::operator=(matrix_complex<DType> const& cp)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=cp.sizex;
        sizey=cp.sizey;
        NumE=cp.NumE;
        width=cp.width;
        NumD=cp.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=cp.sizex || sizey!=cp.sizey)
    {
        delete [] data;
        sizex=cp.sizex;
        sizey=cp.sizey;
        NumE=cp.NumE;
        width=cp.width;
        NumD=cp.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=cp.data[iter];
    }
    return *this;
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator*=(DType factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator/=(DType factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator*=(int factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator/=(int factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator*=(double factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator/=(double factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator*=(const complex<DType>& factor)
{
    DType& fRe=factor.data[0];
    DType& fIm=factor.data[1];
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType NewRe=Re*fRe-Im*fIm;
        DType NewIm=Re*fIm+Im*fRe;
        Re=NewRe;
        Im=NewIm;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator/=(const complex<DType>& factor)
{
    DType& fRe=factor.data[0];
    DType& fIm=factor.data[1];
    DType d=fRe*fRe+fIm*fIm;
    fRe=fRe/d;
    fIm=-fIm/d;
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType NewRe=Re*fRe-Im*fIm;
        DType NewIm=Re*fIm+Im*fRe;
        Re=NewRe;
        Im=NewIm;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator+=(matrix_complex<DType> const& matrixB)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]+=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator-=(matrix_complex<DType> const& matrixB)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]-=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator+=(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]+=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::operator-=(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]-=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleMultiply(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]*=matrixB.data[iter];
        data[iter*2+1]*=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleDivide(matrix<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]/=matrixB.data[iter];
        data[iter*2+1]/=matrixB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleMultiply(matrix_complex<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType& fRe=matrixB.data[2*iter];
        DType& fIm=matrixB.data[2*iter+1];
        DType NewRe=Re*fRe-Im*fIm;
        DType NewIm=Re*fIm+Im*fRe;
        Re=NewRe;
        Im=NewIm;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleDivide(matrix_complex<DType> const& matrixB)
{
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType& fRe=matrixB.data[2*iter];
        DType& fIm=matrixB.data[2*iter+1];
        DType d=fRe*fRe+fIm*fIm;
        DType NewRe=Re*fRe+Im*fIm;
        DType NewIm=Im*fRe-Re*fIm;
        Re=NewRe/d;
        Im=NewIm/d;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::add(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]+mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::minus(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]-mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleMultiply(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]*mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleDivide(matrix<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=sizex*sizey;
        width=sizex*2;
        NumD=2*NumE;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]/mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::add(matrix<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mB.sizex || sizey!=mB.sizey)
    {
        delete [] data;
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]+mB.data[iter*2];
        data[iter*2+1]=mB.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::minus(matrix<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mB.sizex || sizey!=mB.sizey)
    {
        delete [] data;
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]-mB.data[iter*2];
        data[iter*2+1]=-mB.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::minus(matrix_complex<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter*2]-mB.data[iter];
        data[iter*2+1]=mA.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleMultiply(matrix<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mB.sizex || sizey!=mB.sizey)
    {
        delete [] data;
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter]*mB.data[iter*2];
        data[iter*2+1]=mA.data[iter]*mB.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleDivide(matrix_complex<DType> const& mA, matrix<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        data[iter*2]=mA.data[iter*2]/mB.data[iter];
        data[iter*2+1]=mA.data[iter*2+1]/mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleDivide(matrix<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mB.sizex || sizey!=mB.sizey)
    {
        delete [] data;
        sizex=mB.sizex;
        sizey=mB.sizey;
        NumE=mB.NumE;
        width=mB.width;
        NumD=mB.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& BRe=mB.data[iter*2];
        DType& BIm=mB.data[iter*2+1];
        DType d=BRe*BRe+BIm*BIm;
        data[iter*2]=mA.data[iter]*BRe/d;
        data[iter*2+1]=-mA.data[iter]*BIm/d;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::add(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=mA.data[iter]+mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::minus(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=mA.data[iter]-mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleMultiply(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& Re=mA.data[2*iter];
        DType& Im=mA.data[2*iter+1];
        DType& fRe=mB.data[2*iter];
        DType& fIm=mB.data[2*iter+1];
        data[2*iter]=Re*fRe-Im*fIm;
        data[2*iter+1]=Re*fIm+Im*fRe;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::eleDivide(matrix_complex<DType> const& mA, matrix_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
    else if (sizex!=mA.sizex || sizey!=mA.sizey)
    {
        delete [] data;
        sizex=mA.sizex;
        sizey=mA.sizey;
        NumE=mA.NumE;
        width=mA.width;
        NumD=mA.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumE; ++iter)
    {
        DType& Re=mA.data[2*iter];
        DType& Im=mA.data[2*iter+1];
        DType& fRe=mB.data[2*iter];
        DType& fIm=mB.data[2*iter+1];
        DType d=fRe*fRe+fIm*fIm;
        data[2*iter]=(Re*fRe+Im*fIm)/d;
        data[2*iter+1]=(Im*fRe-Re*fIm)/d;
    }
}

template <typename DType>
/*inline*/ complex<DType> matrix_complex<DType>::ele(int y, int x)
{
    complex<DType> result;
    result.data[0]=data[2*x+y*width];
    result.data[1]=data[2*x+1+y*width];
    return result;
}

template <typename DType>
/*inline*/ DType& matrix_complex<DType>::ReEle(int y, int x)
{
    return data[2*x+y*width];
}

template <typename DType>
/*inline*/ DType& matrix_complex<DType>::ImEle(int y, int x)
{
    return data[2*x+1+y*width];
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::setRow(int Idx, DType Re, DType Im)
{
    for (int iter=Idx*width; iter<Idx*width+width; iter+=2)
    {
        data[iter]=Re;
        data[iter+1]=Im;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::setCol(int Idx, DType Re, DType Im)
{
    for (int iter=Idx*2; iter<NumD; iter+=width)
    {
        data[iter]=Re;
        data[iter+1]=Im;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::scaleRow(int Idx, DType Re, DType Im)
{
    for (int iter=Idx*width; iter<Idx*width+width; iter+=2)
    {
        DType& dRe=data[iter];
        DType& dIm=data[iter+1];
        DType NewRe=dRe*Re-dIm*Im;
        DType NewIm=dRe*Im+dIm*Re;
        dRe=NewRe;
        dIm=NewIm;
    }
}

template <typename DType>
/*inline*/ void matrix_complex<DType>::scaleCol(int Idx, DType Re, DType Im)
{
    for (int iter=Idx*2; iter<NumD; iter+=width)
    {
        DType& dRe=data[iter];
        DType& dIm=data[iter+1];
        DType NewRe=dRe*Re-dIm*Im;
        DType NewIm=dRe*Im+dIm*Re;
        dRe=NewRe;
        dIm=NewIm;
    }
}

template <typename DType>
math_vector<DType>::math_vector()
{
    size=0;
    data=0;
}

template <typename DType>
math_vector<DType>::math_vector(int size1)
{
    size=size1;
    data=new DType[size];
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::alloc(int size1)
{
    size=size1;
    if (data!=0)
    {
        delete [] data;
    }
    data=new DType[size];
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
math_vector<DType>::math_vector(math_vector<DType> const& cp)
{
    size=cp.size;
    data=new DType[size];
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=cp.data[iter];
    }
}

template <typename DType>
math_vector<DType>::~math_vector()
{
    if (data!=0)
    {
        delete [] data;
    }
}

template <typename DType>
/*inline*/ DType& math_vector<DType>::operator[](int idx)
{
    return data[idx];
}

template <typename DType>
/*inline*/ math_vector<DType>& math_vector<DType>::operator=(math_vector<DType> const& cp)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=cp.size;
        data=new DType[size];
    }
    else if (size!=cp.size)
    {
        delete [] data;
        size=cp.size;
        data=new DType[size];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=cp.data[iter];
    }
    return *this;
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator*=(DType factor)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator/=(DType factor)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator*=(int factor)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator/=(int factor)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator*=(double factor)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator/=(double factor)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator+=(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]+=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::operator-=(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]-=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::eleMultiply(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]*=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::eleDivide(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]/=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::add(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        data=new DType[size];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        data=new DType[size];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=mA.data[iter]+mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::minus(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        data=new DType[size];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        data=new DType[size];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=mA.data[iter]-mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::eleMultiply(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        data=new DType[size];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        data=new DType[size];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=mA.data[iter]*mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector<DType>::eleDivide(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        data=new DType[size];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        data=new DType[size];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter]=mA.data[iter]/mB.data[iter];
    }
}

template <typename DType>
math_vector_complex<DType>::math_vector_complex()
{
    size=0;
    NumD=0;
    data=0;
}

template <typename DType>
math_vector_complex<DType>::math_vector_complex(int size1)
{
    size=size1;
    NumD=2*size;
    data=new DType[NumD];
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::alloc(int size1)
{
    size=size1;
    NumD=2*size;
    if (data!=0)
    {
        delete [] data;
    }
    data=new DType[NumD];
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=0;
    }
}

template <typename DType>
math_vector_complex<DType>::math_vector_complex(math_vector<DType> const& cp)
{
    size=cp.size;
    NumD=2*size;
    data=new DType[NumD];
    for (int iter=0; iter<size; ++iter)
    {
        data[2*iter]=cp.data[iter];
    }
}

template <typename DType>
math_vector_complex<DType>::math_vector_complex(math_vector<DType> const& Re, math_vector<DType> const& Im)
{
    size=Re.size;
    NumD=2*size;
    data=new DType[NumD];
    for (int iter=0; iter<size; ++iter)
    {
        data[2*iter]=Re.data[iter];
        data[2*iter+1]=Im.data[iter];
    }
}

template <typename DType>
math_vector_complex<DType>::math_vector_complex(math_vector_complex<DType> const& cp)
{
    size=cp.size;
    NumD=2*size;
    data=new DType[NumD];
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=cp.data[iter];
    }
}

template <typename DType>
math_vector_complex<DType>::~math_vector_complex()
{
    if (data!=0)
    {
        delete [] data;
    }
}

template <typename DType>
/*inline*/ DType& math_vector_complex<DType>::operator[](int idx)
{
    return data[idx];
}

template <typename DType>
/*inline*/ math_vector_complex<DType>& math_vector_complex<DType>::operator=(math_vector<DType> const& cp)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=cp.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=cp.size)
    {
        delete [] data;
        size=cp.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=cp.data[iter];
    }
    return *this;
}

template <typename DType>
/*inline*/ math_vector_complex<DType>& math_vector_complex<DType>::operator=(math_vector_complex<DType> const& cp)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=cp.size;
        NumD=cp.NumD;
        data=new DType[NumD];
    }
    else if (size!=cp.size)
    {
        delete [] data;
        size=cp.size;
        NumD=cp.NumD;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=cp.data[iter];
    }
    return *this;
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator*=(DType factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator/=(DType factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator*=(int factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator/=(int factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator*=(double factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]*=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator/=(double factor)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]/=factor;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator*=(const complex<DType>& factor)
{
    DType& fRe=factor.data[0];
    DType& fIm=factor.data[1];
    for (int iter=0; iter<size; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType NewRe=Re*fRe-Im*fIm;
        DType NewIm=Re*fIm+Im*fRe;
        Re=NewRe;
        Im=NewIm;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator/=(const complex<DType>& factor)
{
    DType& fRe=factor.data[0];
    DType& fIm=factor.data[1];
    DType d=fRe*fRe+fIm*fIm;
    fRe=fRe/d;
    fIm=-fIm/d;
    for (int iter=0; iter<size; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType NewRe=Re*fRe-Im*fIm;
        DType NewIm=Re*fIm+Im*fRe;
        Re=NewRe;
        Im=NewIm;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator+=(math_vector_complex<DType> const& math_vectorB)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]+=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator-=(math_vector_complex<DType> const& math_vectorB)
{
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]-=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator+=(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]+=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::operator-=(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]-=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleMultiply(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]*=math_vectorB.data[iter];
        data[iter*2+1]*=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleDivide(math_vector<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]/=math_vectorB.data[iter];
        data[iter*2+1]/=math_vectorB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleMultiply(math_vector_complex<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType& fRe=math_vectorB.data[2*iter];
        DType& fIm=math_vectorB.data[2*iter+1];
        DType NewRe=Re*fRe-Im*fIm;
        DType NewIm=Re*fIm+Im*fRe;
        Re=NewRe;
        Im=NewIm;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleDivide(math_vector_complex<DType> const& math_vectorB)
{
    for (int iter=0; iter<size; ++iter)
    {
        DType& Re=data[2*iter];
        DType& Im=data[2*iter+1];
        DType& fRe=math_vectorB.data[2*iter];
        DType& fIm=math_vectorB.data[2*iter+1];
        DType d=fRe*fRe+fIm*fIm;
        DType NewRe=Re*fRe+Im*fIm;
        DType NewIm=Im*fRe-Re*fIm;
        Re=NewRe/d;
        Im=NewIm/d;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::add(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]+mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::minus(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]-mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleMultiply(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]*mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleDivide(math_vector<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]/mB.data[iter];
        data[iter*2+1]=0;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::add(math_vector<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]+mB.data[iter*2];
        data[iter*2+1]=mB.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::minus(math_vector<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]-mB.data[iter*2];
        data[iter*2+1]=-mB.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::minus(math_vector_complex<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter*2]-mB.data[iter];
        data[iter*2+1]=mA.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleMultiply(math_vector<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter]*mB.data[iter*2];
        data[iter*2+1]=mA.data[iter]*mB.data[iter*2+1];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleDivide(math_vector_complex<DType> const& mA, math_vector<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        data[iter*2]=mA.data[iter*2]/mB.data[iter];
        data[iter*2+1]=mA.data[iter*2+1]/mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleDivide(math_vector<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        DType& BRe=mB.data[iter*2];
        DType& BIm=mB.data[iter*2+1];
        DType d=BRe*BRe+BIm*BIm;
        data[iter*2]=mA.data[iter]*BRe/d;
        data[iter*2+1]=-mA.data[iter]*BIm/d;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::add(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=mA.data[iter]+mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::minus(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<NumD; ++iter)
    {
        data[iter]=mA.data[iter]-mB.data[iter];
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleMultiply(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        DType& Re=mA.data[2*iter];
        DType& Im=mA.data[2*iter+1];
        DType& fRe=mB.data[2*iter];
        DType& fIm=mB.data[2*iter+1];
        data[2*iter]=Re*fRe-Im*fIm;
        data[2*iter+1]=Re*fIm+Im*fRe;
    }
}

template <typename DType>
/*inline*/ void math_vector_complex<DType>::eleDivide(math_vector_complex<DType> const& mA, math_vector_complex<DType> const& mB)
{
#ifdef CHECK_FOR_ROBUSTNESS
    if (data==0)
    {
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
    else if (size!=mA.size)
    {
        delete [] data;
        size=mA.size;
        NumD=2*size;
        data=new DType[NumD];
    }
#endif
    for (int iter=0; iter<size; ++iter)
    {
        DType& Re=mA.data[2*iter];
        DType& Im=mA.data[2*iter+1];
        DType& fRe=mB.data[2*iter];
        DType& fIm=mB.data[2*iter+1];
        DType d=fRe*fRe+fIm*fIm;
        data[2*iter]=(Re*fRe+Im*fIm)/d;
        data[2*iter+1]=(Im*fRe-Re*fIm)/d;
    }
}

template <typename DType>
/*inline*/ complex<DType> math_vector_complex<DType>::ele(int idx)
{
    complex<DType> result;
    result.data[0]=data[2*idx];
    result.data[1]=data[2*idx+1];
    return result;
}

template <typename DType>
/*inline*/ DType& math_vector_complex<DType>::ReEle(int idx)
{
    return data[2*idx];
}

template <typename DType>
/*inline*/ DType& math_vector_complex<DType>::ImEle(int idx)
{
    return data[2*idx+1];
}



#endif /* defined(__spectral_solver__matrix__) */
