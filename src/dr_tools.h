#ifndef DR_TOOLS_H
#define DR_TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "dr_sort.h"
#include "utlib/uthash.h"
/*
    simply print "Hello world!" to default output stream.
*/
void Hello(void);

/*
    compress an uncased acgt string to a 2bit per base format.
    *length will be change to number of real compressed letters
    
    input: 
        acgtString: string of 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'.
        *length: length of the string, compressing will stop when it reaches the length or the end of the string.
    output: 2bit format string of acgt
*/
char* CompressBase(const char* acgtString, long* length);

/*
    decompress a 2bit format string to an UPPER case ACGT string

    input:
        compressedString: 2bit per base format string
        start: 0-base start pos of the dna sequance
        length: length of the sequence you want to decompress from the start
    output: 1byte encoded string of ACGT
    warning: SEGMENT FAULT! no check in this function about start and length
*/
char* DeCompressBase(const char* compressedString, long start, long length);

/*
    retrive a charater in dna string

    input:
        compressedString: 2bit per base format string
        pos: 0-base base position you want to retrieve in the compressedString
    output: one of 'A', 'C', 'G', 'T'
    warning: SEGMENT FAULT! no check in this function about pos
*/
char RetrieveCompressBase(const char* compressedString, long pos);

/*
    function: return (long)cell(byteStringLength / 4.0)
    implemented by logic, without float type
*/
long Compute2BitStrLength(long byteStringLength);

/*
    quick sort a compressed acgt string to get its SA.

    input: 
        cprString: compressed acgt string
        start: 0-base start pos of the dna sequance
        length: length of the sequence you want to decompress from the start
    output: array of type LONG, represent a suffix array
    warning: SEGMENT FAULT! no check in this function about start and length
*/
// long* QuickSortSaFromCs(const char* cprString, long start, long length);


/*
    return the number of reads and lengthes of each fasta string 
    in the input file in order without read the file into memory

    input:
        filePath: path of a fasta format file
    output: array of type LONG, first element represents the number of reads
        others represent lengthes of each fasta string
        in the input file in order
*/
long* CountFasta(const char* filePath);

/*
    build SA for compressedString

    input:
        compressedString: see function CompressBase
        start_pos: 0-base start position you want to begin your sort in original string
        length: the number of base start from start_pos you want to sort
    output: LONG array SA whose size is length + 1, with SA[0] = length correspond to '$'
        it is a permutation of [0, length]
*/
long* BuildSA_QuickSort(const char* compressedString, long start_pos, long length);
long* BuildSA_QuickSort_CompareToEnd(const char* compressedString, long start_pos, long length, long compressedStringLength);
/*
    build Psi array from an SA

    input: 
        compressedString: source string of SA, used to do binary search
        start_pos: start position in the compressedString where your SA comes from
        SA: output of function BuildSA_QuickSort
        length: length of origin string where your SA comes from

    output: LONG array Psi whose size is the same as SA, Psi[0] = SA^(-1)[0]
        it is a permutation of [0, length]
*/
long* BuildPsi_BinarySearch(const char* compressedString, long start_pos, const long* SA, long length);
long* BuildPsi_BinarySearch_CompareToEnd(const char* compressedString, long start_pos, const long* SA, long length, long compressedStringLength);

/*
    find RANK[target] in a psi array

    input:
        psi: a permutation of [0, length]
        length: length of string where your psi comes from
        target: target in [0, length], represents start_pos in a acgt string 
    output:
        RANK[target] of the SA correspond to your psi
*/
long RetriveRankFromPsi(const long* psi, long target);

/*
    Build acgtString from a psi array

    input:
        psi: an array whose length is (length + 1)
        length: the length of acgtString correspond to psi, often len(psi) - 1
        characterSet: substring of 'ACGT'
        l_x: 0-base start pos of ACGT region in your psi
    output:
        char array of 'A','C','G','T' whose length is length

    instructions:
        if your characterSet is shorter than the detected characterSet of psi, output would be strange
        please make sure your characterSet is right especially when its length is less than 4
*/
char* BuildStrFromPsi(const long* psi, long length, const char* characterSet, const long l_x[4]); 


struct read {
    char* name;
    char* seq;
    long length;
    UT_hash_handle hh;
};

void add_read(struct read** reads, char* name, char* seq, long length);

/*
    read a fastq file into two separate hash tables

    input:
        filePath:   path to a fastq format file
    output:
        reads:      hashtable, key is seq name, value is read (if you don't need it , just convey a NULL pointer)
        qualities:  hashtable, key is seq name, value is quality string (if you don't need it , just convey a NULL pointer)
*/
void ReadFastq(char* filePath, struct read** reads, struct read** qualities);


/*
    return the number of reads and lengthes of each fasta string 
    in the input file in order without read the file into memory

    input:
        filePath: path of a fasta format file
    output: array of type LONG, first element represents the number of reads
        others represent lengthes of each fasta string
        in the input file in order
*/
long* CountFastq(const char* filePath);


/*
    read a fasta file into a hash table

    input:
        filePath:   path to a fastq format file
    output:
        reads:      hashtable, key is seq name, value is read
*/
void ReadFasta(char* filePath, struct read** reads);

/*
    read a line of string from file

    input:
        f: a fp opened with "r"
    output:
        length: length of readed line (if you don't need it , just convey a NULL pointer)
        string: a pointer to memory where the string stores (we don't need a buffer)
*/
char* readline(FILE* f, long* plength);

#endif