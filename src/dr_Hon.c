#include "dr_Hon.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dr_sort.h"
#include "dr_tools.h"


const char* HonSaPsi_BinarySearch_X_Region_CompressedArray = NULL;
const long* HonSaPsi_BinarySearch_X_Region_SA = NULL;
long HonSaPsi_BinarySearch_X_Region_StartPos = 0;
// target should be either of 'A', 'C', 'G' and 'T'
int HonSaPsi_BinarySearch_X_Region_CmpFunc(const void* target, long pos) {
    char t = *(char*)target;
    char chAtPos = RetrieveCompressBase(
        HonSaPsi_BinarySearch_X_Region_CompressedArray, 
        HonSaPsi_BinarySearch_X_Region_SA[pos] + HonSaPsi_BinarySearch_X_Region_StartPos
    );
    if (t > chAtPos) return 1;
    else if (t < chAtPos) return -1;
    else if (t == chAtPos) return 0;
}
long** HonSaPsi(const char* compressedString, long length) {
    /* 拒绝处理空串 */
    if (length <= 0) {fprintf(stderr, "HonSaPsi: null string was inputed!\n"); exit(1);}

    /* 计算一些元信息 */
    long lNormalSegmentLength = (long)log2((double)length);
    long lNumNormalSegment = length / lNormalSegmentLength;
    long lNormalSegmentTotalLength = lNumNormalSegment * lNormalSegmentLength;
    long lLastSegmentLength = length == lNormalSegmentTotalLength? lNormalSegmentLength : length - lNormalSegmentTotalLength;
    printf("length:\t\t\t%ld\nnormal seg length:\t%ld\nnum normal seg:\t\t%ld\nlast seg length:\t%ld\nnormal seg total len:\t%ld\n", 
        length, lNormalSegmentLength, lNumNormalSegment, lLastSegmentLength, lNormalSegmentTotalLength);

    /* 定义用于增量构建的变量 */
    long* SA_Tprime = NULL;
    long* psi_Tprime = NULL;
    long* SA_Ti = NULL;
    long* psi_Ti = NULL;
    long* order = malloc(sizeof(long) * (lNormalSegmentLength + 1));
    long new_block_start_pos;
    long TprimeLength;

    // 处理第一块
    SA_Tprime = BuildSA_QuickSort(compressedString, length - lLastSegmentLength, lLastSegmentLength);
    psi_Tprime = BuildPsi_BinarySearch(compressedString, length - lLastSegmentLength, SA_Tprime, lLastSegmentLength + 1);
    /*调试信息*/for (long i = 0; i < lLastSegmentLength + 1; ++i) {printf("%ld ", SA_Tprime[i]);} printf("\n");
    new_block_start_pos = length - lLastSegmentLength - lNormalSegmentLength;
    TprimeLength = lLastSegmentLength + 1;
    // 寻找初始块最长后缀在Tprime的SA中的rank，时间复杂度为O(logn) (线性搜索，搜索域长度为logn)
    long lastOrder;
    for (long i = 1; i < lLastSegmentLength + 1; ++i) { 
        if (SA_Tprime[i] == 0) {lastOrder = i; break;}
    }
    
    for (; new_block_start_pos >= 0; new_block_start_pos -= lNormalSegmentLength) {
        /* sort the suffixes suf_1, suf_2, ..., suf_l */
        SA_Ti = BuildSA_QuickSort(compressedString, new_block_start_pos, lNormalSegmentLength);

        /* for every suf_i, calculate order(suf_i, T') */
        /////* binary search to find l_x and r_x */
        order[0] = 0;
        HonSaPsi_BinarySearch_X_Region_CompressedArray = compressedString;
        HonSaPsi_BinarySearch_X_Region_SA = SA_Tprime;
        HonSaPsi_BinarySearch_X_Region_StartPos = new_block_start_pos + lNormalSegmentLength;
        char target[4] = {'A', 'C', 'G', 'T'};
        long tag[8]; // correspond to l_a, r_a, l_c, r_c, l_g, r_g, l_t, r_t
        for (int i = 0; i < 4; ++i) {
            tag[i*2] = BinarySearchLeftBound(1, TprimeLength, &target[i], HonSaPsi_BinarySearch_X_Region_CmpFunc);
            tag[i*2 + 1] = BinarySearchRightBound(1, TprimeLength, &target[i], HonSaPsi_BinarySearch_X_Region_CmpFunc);
            // /*打印调试信息*/printf("%ld %ld ", tag[i*2], tag[i*2 + 1]);
        }
        /////* caculate order */
        for (long i = lNormalSegmentLength; i >= 1; --i) {
            char c = RetrieveCompressBase(compressedString, new_block_start_pos + i - 1);
            int tag_pos;
            for (tag_pos = 0; tag_pos < 4; ++tag_pos){
                if (target[tag_pos] == c) break;
            }
            BinarySeachRightBoundLessEqualTarget(tag[tag_pos], tag[tag_pos+1], &lastOrder, TODOCmpFunc);
        }

        /* computer the pai function for T^iT' */


        TprimeLength += lNormalSegmentLength;

        break;
    }
}