@ECHO   OFF
GCC %1.c -o %1.exe -Wall -Werror -Wextra -ggdb -lm