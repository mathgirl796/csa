#ifndef DR_HON_H
#define DR_HON_H


/*
    build SA and psi for a compressedString

    input:
        compressedString: 2bit encode byte array, see more in dr_tools.h
        length: the number of base in your compressedString, not size of byte array
    output: two long* pointers RET
        RET[0] point to a LONG type array of size LENGTH + 1, storing SA of the compressedString, 
            whose 0-th element is SA[0], correspond to the '$' symbol
        RET[1] point to a LONG type array of size LENGTH + 1, storing psi of the compressedString, 
            whose 0-th element is psi[0], correspond to the '$' symbol
*/
long** HonSaPsi(const char* compressedString, long length);








#endif