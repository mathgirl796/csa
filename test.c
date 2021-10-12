#include <stdio.h>

int main(){

    unsigned char a = 124;
    printf("%d\n", a << 4 >> 6);
    printf("%-10s\n", "123");
    printf("%-.3s\n", "qwertyuioplkjhgfdsaaaaaaaaaaaaaaa");
    return 0;
}