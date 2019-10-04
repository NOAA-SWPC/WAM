#if 0
>>
>> Make this header file available as ESMFVersionDefine.h in order to build
>> NEMS against an ESMF 3 installation.
>>
#endif

#define ESMF_3

#ifndef ESMF_MAJOR_VERSION
#define ESMF_MAJOR_VERSION 3
#define ESMF_MINOR_VERSION 1
#define ESMF_REVISION 0
#define ESMF_PATCHLEVEL 4
#endif

#include "./ESMFVersionLogic.h"
