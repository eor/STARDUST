#ifndef __log_h__
#define __log_h__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#define clean_errno() (errno == 0 ? "None" : strerror(errno))

#define log_debug(M, ...) fprintf(fpDebugLog, "DEBUG: " M "\n", ##__VA_ARGS__); fflush(fpDebugLog)

#define log_main(M, ...)  fprintf(fpMainLog,  "MAIN:  " M "\n", ##__VA_ARGS__); fflush(fpMainLog)

#define log_error(M, ...) fprintf(fpMainLog,  "ERROR: (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__); fflush(fpMainLog) 

#define check(A, M, ...) if(!(A)) { log_error(M, ##__VA_ARGS__); errno=0; goto error; }

#endif
