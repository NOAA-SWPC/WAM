#if 0
>>
>> Make this header file available as ESMFVersionDefine.h in order to build
>> NEMS against an ESMF 4_0_0rp2 installation.
>>
#endif

#undef ESMF_3

#ifndef ESMF_MAJOR_VERSION
#define ESMF_MAJOR_VERSION 4
#define ESMF_MINOR_VERSION 0
#define ESMF_REVISION 0
#define ESMF_PATCHLEVEL 2
#endif

#include "./ESMFVersionLogic.h"
