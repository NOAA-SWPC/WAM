#if 0
>>
>> Make this header file available as ESMFVersionDefine.h in order to build
>> NEMS against an ESMF 5_2_0r_beta_snapshot installation.
>>
#endif

#undef ESMF_3

#ifndef ESMF_MAJOR_VERSION
#define ESMF_MAJOR_VERSION 5
#define ESMF_MINOR_VERSION 2
#define ESMF_REVISION 0
#define ESMF_PATCHLEVEL 0
#endif

#include "./ESMFVersionLogic.h"
