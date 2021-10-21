#include "dr_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "dr_sort.h"

void Hello(void) {
    printf("Hello, world!\n");
}

long* CountFasta(const char* filePath) {
    // 打开文件
    FILE* f = fopen(filePath, "r");
    if (f == NULL) {
        perror("CountFasta");
        exit(1);
    }
    // 逐个读取文件字符，根据情况执行相应行为
    long maxSize = 1;
    long pos = 0;
    long* ret = malloc(maxSize * sizeof(long));
    while (1) {
        char c = fgetc(f);
        // putchar(c);
        /* 读到文件结尾：结束循环 */
        if (c == EOF) break; 
        /* 读到">"：读掉这一行, 游标加1, 游标处初始化为0 */ 
        else if (c == '>') { 
            while ((c = fgetc(f)) != '\n' && c != EOF);
            pos ++;
            // printf("maxSize: %ld\npos: %ld\n", maxSize, pos);
            if (pos == maxSize && c != EOF) {
                maxSize *= 2;
                ret = realloc(ret, maxSize * sizeof(long)); // 这里千万记得乘sizeof(long)
            }
            ret[pos] = 0;
        }
        /* 遇到空白符：跳过 */
        else if (isspace(c)) { 
            continue;
        }
        /* 其他情况：正常计数 */
        else if (pos >= 1) { 
            ret[pos] ++;
            // printf("%ld\n", ret[pos]);
        }
    }
    ret[0] = pos;
    ret = realloc(ret, (pos + 1) * sizeof(long));  // 回收未使用的临时堆内存
    fclose(f);
    return ret;
}


long Compute2BitStrLength(long byteStringLength) {
    return (byteStringLength % 4 == 0) ? byteStringLength / 4 : byteStringLength / 4 + 1;
}

// rule: ['a': b'00', 'c': b'01', 'g': b'10', 't': b'11']
// acgtString "123456789" will be converted to "'1234','5678','9xxx'"
char* CompressBase(const char* acgtString, long* length) {
    long strLength = *length;
    long realStringLength = strLength > *length ? strLength : *length;
    long retStringLength = Compute2BitStrLength(realStringLength);
    char* ret = malloc(retStringLength * sizeof(char));

    long sPos = -1; // pos of acgtString
    long rPos = -1; // pos of ret

    while(1) {
        if (++rPos > retStringLength) break;
        char buffer = 0;
        for (int i = 0; i < 4; i++) {
            if (++sPos >= realStringLength) break;
            char select = 0;
            switch(acgtString[sPos]){
                case 'A':
                case 'a': select = 0; break;
                case 'C':
                case 'c': select = 1; break;
                case 'G': 
                case 'g': select = 2; break;
                case 'T': 
                case 't': select = 3; break;
                default:
                    fprintf(stderr, "CompressBase: unexpected character %d at string pos %ld\n", acgtString[sPos], sPos);
                    exit(1);
            }
            buffer += select << (3-i) * 2;
        }
        ret[rPos] = buffer;
    }

    *length = realStringLength;
    return ret;
}

char* DeCompressBase(const char* compressedString, long start, long length) {
    long cPos = start / 4; // compressedString pos
    char* ret = malloc(length * sizeof(char) + 1); // +1 for '\0'

    long rPos = -1; // ret pos
    int i = start % 4;
    while(1) {
        if (++rPos >= length) break;

        char code = (unsigned char) // 每次移位后都要强转为无符号数，因为移位表达式的值是有符号数
            ((unsigned char)compressedString[cPos] << 2 * i) // 保证进行的是逻辑右移
             >> 6;
        if (++i == 4) {
            i = 0;
            cPos++;
        }

        switch (code)
        {
            case 0: ret[rPos] = 'A'; break;
            case 1: ret[rPos] = 'C'; break;
            case 2: ret[rPos] = 'G'; break;
            case 3: ret[rPos] = 'T'; break;
            default: fprintf(stderr, "DeCompressBase: decode fault: %d!\n", code); exit(1);
        }
    }
    ret[rPos] = '\0';
    return ret;
}

char RetrieveCompressBase(const char* compressedString, long pos) {
    long cPos = pos / 4; // compressedString pos
    int i = pos % 4;
    char ret;

    char code = (unsigned char) // 每次移位后都要强转为无符号数，因为移位表达式的值是有符号数
        ((unsigned char)compressedString[cPos] << 2 * i) // 保证进行的是逻辑右移
        >> 6;

    switch (code)
    {
        case 0: ret = 'A'; break;
        case 1: ret = 'C'; break;
        case 2: ret = 'G'; break;
        case 3: ret = 'T'; break;
        default: fprintf(stderr, "RetrieveCompressBase: decode fault: %d!\n", code); exit(1);
    }
    return ret;
}

const char* BuildSA_QuickSort_CompressedArray = NULL;
long BuildSA_QuickSort_MaxCmpPos = 0;
int BuildSA_QuickSort_CmpFunc(long l, long r) {
    while (l < BuildSA_QuickSort_MaxCmpPos && r < BuildSA_QuickSort_MaxCmpPos) {
        char lchar = RetrieveCompressBase(BuildSA_QuickSort_CompressedArray, l);
        char rchar = RetrieveCompressBase(BuildSA_QuickSort_CompressedArray, r);
        if (lchar < rchar) return -1;
        else if (lchar > rchar) return 1;
        else {l++; r++;}
    }
    if (l > r) return -1;
    else if (l < r) return 1;
    else {
        fprintf(stderr, "BuildSA_QuickSort_CmpFunc: equal happens when compare suffixes!\n");
        exit(1);
    }
}

long* BuildSA_QuickSort(const char* compressedString, long start_pos, long length) {
    BuildSA_QuickSort_CompressedArray = compressedString;
    BuildSA_QuickSort_MaxCmpPos = start_pos + length;
    long* quickSortResult = QuickSort(start_pos, start_pos+length, BuildSA_QuickSort_CmpFunc, SimpleQuickSortPartition);
    for (long i = 0; i < length; ++i) quickSortResult[i] -= start_pos;
    long* SA = malloc(sizeof(long) * (length + 1));
    SA[0] = length;
    memcpy(SA + 1, quickSortResult, sizeof(long) * length);    
    free(quickSortResult);
    return SA;
}

long* BuildSA_QuickSort_CompareToEnd(const char* compressedString, long start_pos, long length, long compressedStringLength) {
    BuildSA_QuickSort_CompressedArray = compressedString;
    BuildSA_QuickSort_MaxCmpPos = compressedStringLength;
    long* quickSortResult = QuickSort(start_pos, start_pos+length, BuildSA_QuickSort_CmpFunc, SimpleQuickSortPartition);
    for (long i = 0; i < length; ++i) quickSortResult[i] -= start_pos;
    long* SA = malloc(sizeof(long) * (length + 1));
    SA[0] = length;
    memcpy(SA + 1, quickSortResult, sizeof(long) * length);    
    free(quickSortResult);
    return SA;
}

const char* BuildPsi_BinarySearch_CompressedArray = NULL;
long BuildPsi_BinarySearch_CompressedString_MaxCmpPos = 0;
const long* BuildPsi_BinarySearch_SA = NULL;
long BuildPsi_BinarySearch_CompressedString_StartPos = 0;
// target指向字符串的某个下标， pos是SA数组的某个下标
int BuildPsi_BinarySearch_CmpFunc(const void* target, long pos) {
    long l = *(long*)target;
    long r = BuildPsi_BinarySearch_SA[pos] + BuildPsi_BinarySearch_CompressedString_StartPos;
    if (l == r) return 0;
    // printf("%ld %ld\n", l, r); // 打印调试信息
    while (l < BuildPsi_BinarySearch_CompressedString_MaxCmpPos && r < BuildPsi_BinarySearch_CompressedString_MaxCmpPos) {
        char lchar = RetrieveCompressBase(BuildPsi_BinarySearch_CompressedArray, l);
        char rchar = RetrieveCompressBase(BuildPsi_BinarySearch_CompressedArray, r);
        if (lchar < rchar) return -1;
        else if (lchar > rchar) return 1;
        else {l++; r++;}
    }
    if (l > r) return -1;
    else if (l < r) return 1;
    else return 0;
}
long* BuildPsi_BinarySearch(const char* compressedString, long start_pos, const long* SA, long length) {
    BuildPsi_BinarySearch_CompressedArray = compressedString;
    BuildPsi_BinarySearch_CompressedString_StartPos = start_pos;
    BuildPsi_BinarySearch_CompressedString_MaxCmpPos = start_pos + length;
    BuildPsi_BinarySearch_SA = SA;
    long SA_size = length + 1;
    long* psi = malloc(sizeof(long) * (length + 1));
    long target = start_pos;
    /* psi[0]是特殊情况 */
    psi[0] = BinarySearch(0, SA_size, &target, BuildPsi_BinarySearch_CmpFunc);
    for (long i = 1; i < SA_size; ++i) {
        target = SA[i] + 1 + start_pos;
        psi[i] = BinarySearch(1, SA_size, &target, BuildPsi_BinarySearch_CmpFunc);
        if (psi[i] == BINARY_SEARCH_NOT_FOUND) psi[i] = start_pos;
    }
    return psi;
}

long* BuildPsi_BinarySearch_CompareToEnd(const char* compressedString, long start_pos, const long* SA, long length, long compressedStringLength) {
    BuildPsi_BinarySearch_CompressedArray = compressedString;
    BuildPsi_BinarySearch_CompressedString_StartPos = start_pos;
    BuildPsi_BinarySearch_CompressedString_MaxCmpPos = compressedStringLength;
    BuildPsi_BinarySearch_SA = SA;
    long SA_size = length + 1;
    long* psi = malloc(sizeof(long) * (length + 1));
    long target = start_pos;
    /* psi[0]是特殊情况 */
    psi[0] = BinarySearch(1, SA_size, &target, BuildPsi_BinarySearch_CmpFunc);
    for (long i = 1; i < SA_size; ++i) {
        target = SA[i] + 1 + start_pos;
        psi[i] = BinarySearch(1, SA_size, &target, BuildPsi_BinarySearch_CmpFunc);
        if (psi[i] == BINARY_SEARCH_NOT_FOUND) psi[i] = start_pos;
    }
    return psi;
}

long RetriveRankFromPsi(const long* psi, long target) {
    long rank = 0;
    for (int i = 0; i < target + 1; ++i) {
        rank = psi[rank];
    }
    return rank;
}


char* BuildStrFromPsi(const long* psi, long length, const char* characterSet, const long l_x[4])
{
    long pos = psi[0];
    char* retString = malloc(sizeof(char) * length);

    for (long i = 0; i < length; ++i) {
        if (l_x[0] <= pos && pos < l_x[1]) retString[i] = characterSet[0];
        else if (l_x[1] <= pos && pos < l_x[2]) retString[i] = characterSet[1];
        else if (l_x[2] <= pos && pos < l_x[3]) retString[i] = characterSet[2];
        else if (l_x[3] <= pos && pos <= length) retString[i] = characterSet[3];
        else retString[i] = 'N';
        pos = psi[pos];
    }

    return retString;
}