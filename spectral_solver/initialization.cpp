//
//  initialization.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/11.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"

void solver::initialization()
{
    if (0==cRank)
    {
        for (int itert=0; itert<Ntheta; ++itert)
        {
            for (int iterr=0; iterr<Nrp; ++iterr)
            {
#ifdef FU_MODEL
                gsl_matrix_set(Fields, iterr, itert, 2*exp(-r->data[iterr]*r->data[iterr]/(0.04/radius/radius)));
                //gsl_matrix_set(Fields, 2*Nrp+iterr, itert, 15);
#endif
#ifdef LINEAR_TEST_MODEL
                gsl_matrix_set(Fields, iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
                gsl_matrix_set(Fields, Nrp+iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
                gsl_matrix_set(Fields, 2*Nrp+iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
#endif
#ifdef INTERACTION_MODIFIED_FU
                gsl_matrix_set(Fields, iterr, itert, 2*exp(-r->data[iterr]*r->data[iterr]/(0.04/radius/radius)));
                gsl_matrix_set(Fields, Nrp+iterr, itert, 2*exp(-r->data[iterr]*r->data[iterr]/(0.04/radius/radius)));
                gsl_matrix_set(Fields, 4*Nrp+iterr, itert, 15);
#endif
            }
        }
    }
}
