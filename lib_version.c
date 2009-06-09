#include <stdio.h>
#include <math.h>
#include <string.h>
#include "idl_export.h"

IDL_VPTR lib_version(int argc, IDL_VPTR argv[]) {
/* Declare variables from IDL */
IDL_VPTR vers;				// Strings

vers = IDL_StrToSTRING("clean tools 20090609\0");

return(vers);

}
