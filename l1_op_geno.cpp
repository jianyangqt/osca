//
//  l1_op_geno.cpp
//  osc
//
//  Created by Futao Zhang on 31/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "l1_op_geno.h"

void _calc_ld(double* ld,float* snpData,int rsize,int csize) {
    
    vector<float> tmpX(rsize);
    vector<float> tmpX2(rsize);
    
#ifndef __APPLE__
    setNbThreads(thread_num);
    //int Num_CPU_Cores;
#if defined _WIN64 || defined _WIN32
    /*
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    Num_CPU_Cores=si.dwNumberOfProcessors;
     */
#else
   // Num_CPU_Cores=sysconf(_SC_NPROCESSORS_CONF);
#endif
   // omp_set_num_threads(Num_CPU_Cores);
    omp_set_num_threads(thread_num);
#endif
    
    #pragma omp parallel for
    for(int i=0;i<rsize;i++){
        for(int j=0;j<csize;j++){
            tmpX[i] += snpData[i*csize+j];
            tmpX2[i] += snpData[i*csize+j]*snpData[i*csize+j];
        }
    }
    float tmp0;
    #pragma omp parallel for private(tmp0)
    for(int i=0;i<rsize;i++){
        for(int j=0;j<i;j++){
            tmp0=0;
            for(int k=0;k<csize;k++) tmp0 += snpData[i*csize+k]*snpData[j*csize+k];
            ld[i*(i-1)/2+j] =
            (tmp0-tmpX[i]*tmpX[j]/csize)/sqrtf((tmpX2[i]-tmpX[i]*tmpX[i]/csize)*(tmpX2[j]-tmpX[j]*tmpX[j]/csize));
            
        }
    }
    
}
