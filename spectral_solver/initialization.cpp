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
                Fields.ele(iterr, itert)=2*exp(-r[iterr]*r[iterr]/(0.04/radius/radius));
                Fields.ele(2*Nrp+iterr, itert)=15;
#endif
#ifdef LINEAR_TEST_MODEL
                //Fields.ele(iterr, itert)=sin(theta[itert]);
                Fields.ele(iterr, itert)=exp(-10*r[iterr]*r[iterr]);
                Fields.ele(Nrp+iterr, itert)=exp(-10*r[iterr]*r[iterr]);
                Fields.ele(2*Nrp+iterr, itert)=exp(-10*r[iterr]*r[iterr]);
#endif
#ifdef INTERACTION_MODIFIED_FU
                Fields.ele(iterr, itert)=2*exp(-r[iterr]*r[iterr]/(0.04/radius/radius));
                Fields.ele(Nrp+iterr, itert)=2*exp(-r[iterr]*r[iterr]/(0.04/radius/radius));
                Fields.ele(4*Nrp+iterr, itert)=15;
#endif
            }
        }
    }
}
