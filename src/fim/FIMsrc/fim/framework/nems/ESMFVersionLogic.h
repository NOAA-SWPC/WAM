#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)

#define esmf_logfounderror esmf_logmsgfounderror
#define ESMF_LogFoundError ESMF_LogMsgFoundError
#define ESMF_LogFoundAllocError ESMF_LogMsgFoundAllocError

#define STATENAME stateName
#define FIELDNAME name
#define FIELDNAMELIST nameList
#define FIELDCOUNT nameCount
#define DSTPET dst
#define SRCPET src
#define FARRAYPTR fptr

#else

#define STATENAME name
#define FIELDNAME fieldName
#define FIELDNAMELIST fieldNameList
#define FIELDCOUNT fieldCount
#define DSTPET dstPet
#define SRCPET srcPet
#define FARRAYPTR farrayPtr

#endif
