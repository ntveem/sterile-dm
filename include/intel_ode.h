/******************************************************************************* 
!                              INTEL CONFIDENTIAL 
!   Copyright(C) 2007-2008 Intel Corporation. All Rights Reserved. 
!   The source code contained  or  described herein and all documents related to 
!   the source code ("Material") are owned by Intel Corporation or its suppliers 
!   or licensors.  Title to the  Material remains with  Intel Corporation or its 
!   suppliers and licensors. The Material contains trade secrets and proprietary 
!   and  confidential  information of  Intel or its suppliers and licensors. The 
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and 
!   treaty  provisions. No part of the Material may be used, copied, reproduced, 
!   modified, published, uploaded, posted, transmitted, distributed or disclosed 
!   in any way without Intel's prior express written permission. 
!   No license  under any  patent, copyright, trade secret or other intellectual 
!   property right is granted to or conferred upon you by disclosure or delivery 
!   of the Materials,  either expressly, by implication, inducement, estoppel or 
!   otherwise.  Any  license  under  such  intellectual property  rights must be 
!   express and approved by Intel in writing.
! 
!******************************************************************************
!   
!  Header file for Intel(R) ODE Solvers
! 
!*******************************************************************************/

#ifndef _INTEL_ODE_H_
#define _INTEL_ODE_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void dodesol(int*,int*,double*,double*,double*,void*,void*,\
			 double*,double*,double*,double*,double*,int*,int*);
void dodesol_rkm9st(int*,int*,double*,double*,double*,void*,\
					double*,double*,double*,double*,double*,int*);
void dodesol_mk52lfn(int*,int*,double*,double*,double*,void*,\
					 double*,double*,double*,double*,double*,int*,int*);
void dodesol_mk52lfa(int*,int*,double*,double*,double*,void*,void*,\
					 double*,double*,double*,double*,double*,int*,int*);
void dodesol_rkm9mkn(int*,int*,double*,double*,double*,void*,\
					 double*,double*,double*,double*,double*,int*,int*);
void dodesol_rkm9mka(int*,int*,double*,double*,double*,void*,void*,\
					 double*,double*,double*,double*,double*,int*,int*);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _INTEL_ODE_H_ */
