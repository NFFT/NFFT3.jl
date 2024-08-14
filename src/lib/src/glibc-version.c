#include <gnu/libc-version.h>
const char* glibc_version(){
  return(gnu_get_libc_version());
}
