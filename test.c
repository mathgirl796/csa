#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

int main(){

    char* a = malloc(2);
    char* b = malloc(3);
    free(a);
    a = b;
    free(b);
    printf("%ld\n", LONG_MAX);
}