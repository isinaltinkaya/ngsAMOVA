#include "shared.h"
#include <time.h>
#include <ctype.h>

int strIsNumeric(const char* val) {
    for (size_t i = 0; i < strlen(val); i++) {
        if (!isdigit(val[i])) {
            return 0;
        }
    }
    return 1;
}

using size_t = decltype(sizeof(int));

/// @brief get current time
/// @return time as char*
char* get_time() {
    time_t current_time;
    struct tm* local_time;
    current_time = time(NULL);
    local_time = localtime(&current_time);
    return (asctime(local_time));
}

