#include <stdio.h>
#include <limits.h>


int main(){

    int a = 1;
    int b = 2;
    const int *p = &a;
    p = &b;
    return 0;
}