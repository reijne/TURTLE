#if !defined(_WINF2C_H_)
#define _WINF2C_H_

/*
 * $Id: winf2c.h,v 1.1.1.1 2000-10-26 16:29:23 psh Exp $
 */

typedef struct{
        char *string;
        int  len;
}_fcd;

#define _fcdlen(x)   (x).len
#define _fcdtocp(x)  (x).string
#define _cptofcd(str,len) str,len

#define FATR __stdcall

#endif /* _WINF2C_H_ */
