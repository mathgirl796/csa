#include "dr_Hon.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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

const long* HonSaPsi_Caculate_Order_BinarySearch_Psi_PsiArray = NULL;
int HonSaPsi_Caculate_Order_BinarySearch_Psi_CmpFunc(const void* target, long pos) {
    long t = *(long*)target;
    long eleAtPos = HonSaPsi_Caculate_Order_BinarySearch_Psi_PsiArray[pos];
    if (t > eleAtPos) return 1;
    else if (t < eleAtPos) return -1;
    else return 0;
}
long** HonSaPsi(const char* compressedString, long length) {
    /* 拒绝处理空串 */
    if (length <= 0) {fprintf(stderr, "HonSaPsi: null string was inputed!\n"); exit(1);}

    /* 计算一些元信息 */
    long lNormalSegmentLength = (long)log2((double)length);
    // long lNormalSegmentLength = 5000;
    long lNumNormalSegment = length / lNormalSegmentLength;
    long lNormalSegmentTotalLength = lNumNormalSegment * lNormalSegmentLength;
    long lLastSegmentLength = length == lNormalSegmentTotalLength? lNormalSegmentLength : length - lNormalSegmentTotalLength;
    printf("length:\t\t\t%ld\nnormal seg length:\t%ld\nnum normal seg:\t\t%ld\nlast seg length:\t%ld\nnormal seg total len:\t%ld\n", 
        length, lNormalSegmentLength, lNumNormalSegment, lLastSegmentLength, lNormalSegmentTotalLength);

    /* 定义用于增量构建的变量 */
    long* SA_Tprime = NULL; // 处理第一块时初始化
    long* psi_Tprime = NULL; // 处理第一块时初始化
    long* SA_Ti = NULL;  // 每次迭代会被赋值
    long* psi_Ti = NULL;  // 每次迭代会被赋值
    long* order = malloc(sizeof(long) * (lNormalSegmentLength + 1));
    long new_block_start_pos;
    long TprimeLength;
    long* tempPsi_TiTprime = NULL;

    // 处理第一块
    SA_Tprime = BuildSA_QuickSort(compressedString, length - lLastSegmentLength, lLastSegmentLength);
    psi_Tprime = BuildPsi_BinarySearch(compressedString, length - lLastSegmentLength, SA_Tprime, lLastSegmentLength);
    // /*调试信息*/printf("first psi_Tprime");for (long i = 0; i < lLastSegmentLength + 1; ++i) {printf("%ld ", psi_Tprime[i]);} printf("\n");
    new_block_start_pos = length - lLastSegmentLength - lNormalSegmentLength;
    TprimeLength = lLastSegmentLength + 1;
    // 寻找初始块最长后缀在Tprime的SA中的rank，时间复杂度为O(logn) (线性搜索，搜索域长度为logn)
    long lastOrder;
    for (long i = 1; i < lLastSegmentLength + 1; ++i) { 
        if (SA_Tprime[i] == 0) {lastOrder = i; break;}
    }
    
    int debug_flag = 1;
    long print_seg = 10000;
    long print_count = new_block_start_pos / print_seg;
    int log_flag = 0; // 每隔print_seg长度，打印一轮log。log_flag会在对应轮开头设为1，结尾设为0
    time_t time_mid, time_lr=0, time_f=0, time_g=0, time_order=0, time_merge=0;
    if (debug_flag == 1) {
        printf("start_pos\tlf\t\torder\t\tf\t\tg\t\tmerge\n");
    }
    for (; new_block_start_pos >= 0; new_block_start_pos -= lNormalSegmentLength) {
        /* sort the suffixes suf_1, suf_2, ..., suf_l */
        SA_Ti = BuildSA_QuickSort_CompareToEnd(compressedString, new_block_start_pos, lNormalSegmentLength, length);
        psi_Ti = BuildPsi_BinarySearch_CompareToEnd(compressedString, new_block_start_pos, SA_Ti, lNormalSegmentLength, length);

        long print_count_a = new_block_start_pos / print_seg;
        if (print_count_a < print_count) {
            print_count = print_count_a;
            log_flag = 1;
            
            printf("%ld\t\t", new_block_start_pos);
        }
        // /*调试信息*/ printf("************************** 新块起始位置:%ld ************************\n", new_block_start_pos);
        //             printf("新块内容: ");for(int i = new_block_start_pos;  i < new_block_start_pos + lNormalSegmentLength; ++i){printf("%c", RetrieveCompressBase(compressedString, i));}printf("\n");
        //             printf("即将处理: ");for(int i = new_block_start_pos;  i < length ; ++i){printf("%c", RetrieveCompressBase(compressedString, i));}printf("\n");        
        // /* 调试信息：psi_Ti的内容 */printf("psi_Ti\t");for(int i=0;i<lNormalSegmentLength+1;++i){printf("%ld ",psi_Ti[i]);}printf("\n");
        // /* 调试信息：psi_Ti的内容 */printf("SA_Ti\t");for(int i=0;i<lNormalSegmentLength+1;++i){printf("%ld ",SA_Ti[i]);}printf("\n");

        time_mid = clock();
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
        }
        if (log_flag == 1 && debug_flag) {
            printf("%lfs\t", (double)(time_lr) / CLOCKS_PER_SEC);
            time_lr = 0;
        }
        else {
            time_lr += clock()-time_mid;
            time_mid = clock();
        }
        // /*调试信息: l_c,r_c的内容*/printf("lc,rc: ");for (int i = 0; i < 4; ++i) {printf("%ld %ld ", tag[i*2], tag[i*2 + 1]);}printf("\n");
        /////* caculate order */
        for (long i = lNormalSegmentLength; i >= 1; --i) {
            char c = RetrieveCompressBase(compressedString, new_block_start_pos + i - 1);
            int tag_pos;
            switch (c)
            {
            case 'A': tag_pos=0; break;
            case 'C': tag_pos=1; break;
            case 'G': tag_pos=2; break;
            case 'T': tag_pos=3; break;
            }
            // 在c区搜索最大的b
            HonSaPsi_Caculate_Order_BinarySearch_Psi_PsiArray = psi_Tprime;
            long max_b;
            if (tag[2*tag_pos] == BINARY_SEARCH_NOT_FOUND) {
                while(tag[2*tag_pos] == BINARY_SEARCH_NOT_FOUND && tag[2*tag_pos] >= 0){
                    tag_pos--;
                }
                if (tag_pos < 0) max_b = 0;
                else max_b = tag[2*tag_pos + 1];
            }
            else {
                max_b = BinarySeachRightBoundLessEqualTarget(tag[2*tag_pos], tag[2*tag_pos+1]+1, &lastOrder, HonSaPsi_Caculate_Order_BinarySearch_Psi_CmpFunc);
                if (max_b == BINARY_SEARCH_NOT_FOUND) {
                    /* c区所有psi值都比lastOrder大，则max_b = l_c - 1 */
                    max_b = tag[2*tag_pos] - 1; 
                }
            }
            // /*调试信息*/printf("suf%ld\tc:%c\tl:%ld\tr:%ld\tlastOrder:%ld\torder:%ld\t\n", i, c, tag[2*tag_pos], tag[2*tag_pos+1], lastOrder, max_b);
            order[i] = max_b;
            lastOrder = order[i];
        }
        if (log_flag == 1 && debug_flag) {
            printf("%lfs\t", (double)(time_order) / CLOCKS_PER_SEC);
            time_order = 0;
        }
        else {
            time_order += clock()-time_mid;
            time_mid = clock();
        }
        /* computer the pai function for T^iT' */
        //// compute f and g
        long* f = malloc(sizeof(long) * (TprimeLength + 1));
        long* g = malloc(sizeof(long) * (lNormalSegmentLength + 1));
        
        ////// compute f
        for (long j = 0; j <= TprimeLength; ++j) {
            long count = 0;
            for (long k = 1; k <= lNormalSegmentLength; ++k) {
                if (order[k] < j) ++count;
            }
            f[j] = j + count;
        }
        if (log_flag == 1 && debug_flag) {
            printf("%lfs\t", (double)(time_f) / CLOCKS_PER_SEC);
            time_f = 0;
        }
        else {
            time_f += clock()-time_mid;
            time_mid = clock();
        }
        ////// compute g
        long rank = 0;
        for (long j = 1; j <= lNormalSegmentLength; ++j) {
            rank = psi_Ti[rank];
            g[j] = order[j] + rank;
        }
        if (log_flag == 1 && debug_flag) {
            printf("%lfs\t", (double)(time_g) / CLOCKS_PER_SEC);
            time_g = 0;
        }
        else {
            time_g += clock()-time_mid;
            time_mid = clock();
        }
        // /*调试信息：order*/for (long j = 1; j <= lNormalSegmentLength; ++j){printf("order[%ld]:%ld ", j, order[j]);}printf("\n");
        // /*调试信息：f*/for (long j = 1; j < TprimeLength; ++j){printf("f[%ld]:%ld ", j, f[j]);}printf("\n");
        // /*调试信息：g*/for (long j = 1; j <= lNormalSegmentLength; ++j){printf("g[%ld]:%ld ", j, g[j]);}printf("\n");

        //// merge and update lastOrder
        tempPsi_TiTprime = malloc(sizeof(long) * (TprimeLength + lNormalSegmentLength + 1)); // 这个指针不要释放，它会被赋给一个上层变量，最后会被返回
        tempPsi_TiTprime[0] = g[1];
        tempPsi_TiTprime[g[lNormalSegmentLength]] = f[psi_Tprime[0]];
        for (long j = 1; j < TprimeLength; ++j) {
            tempPsi_TiTprime[f[j]] = f[psi_Tprime[j]];
        }
        for (long j = 1; j < lNormalSegmentLength; ++j) {
            tempPsi_TiTprime[g[j]] = g[j+1];
        }
        // /*调试信息*/printf("psi_Tprime\t");for(int i=0;i<TprimeLength;i++){printf("%ld ",psi_Tprime[i]);}printf("\n");
        free(psi_Tprime);
        psi_Tprime = tempPsi_TiTprime;
        lastOrder = psi_Tprime[0];
        TprimeLength += lNormalSegmentLength;
        /* 释放临时内存 */
        free(f);
        free(g);
        if (log_flag == 1 && debug_flag) {
            printf("%lfs\n", (double)(time_merge) / CLOCKS_PER_SEC);
            time_merge = 0;
        }
        else {
            time_merge += clock()-time_mid;
            time_mid = clock();
        }

        /* 根据merge后的psi_Tprime构造新的SA_Tprime */
        SA_Tprime = realloc(SA_Tprime, sizeof(long) * TprimeLength);
        long mergeSA_pos = 0;
        for (long i = 0; i < TprimeLength; ++i) {
            mergeSA_pos = psi_Tprime[mergeSA_pos];
            SA_Tprime[mergeSA_pos] = i;
        }

        // /*调试信息*/printf("merge后的psi_Tprime:\t");for (long i = 0; i < TprimeLength; ++i) {printf("%ld ", psi_Tprime[i]);}printf("\n");
        // /*调试信息*/printf("merge后的SA_Tprime:\t");for (long i = 0; i < TprimeLength; ++i) {printf("%ld ", SA_Tprime[i]);}printf("\n");
        // /*调试信息*/break;
        log_flag = 0;
        // getchar();
    }

    /* 释放临时变量 */
    free(SA_Ti);
    free(psi_Ti);
    free(order);
    
    /* 返回结果 */
    long** result = malloc(sizeof(long*) * 2);
    result[0] = SA_Tprime;
    result[1] = psi_Tprime;
    return result;
}