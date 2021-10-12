#include "dr_sort.h"
#include <stdlib.h>
#include <stdio.h>

/*
    递归地执行快速排序
*/
void QuickSortRecursive(long* array, long l, long r, int (*CmpFunc)(long, long), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long, long)));


long* QuickSort(long num_elements, int (*CmpFunc)(long, long), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long, long))) {
    /* 为将要返回的数组申请内存、初始化 */
    long* ret = malloc(num_elements * sizeof(long));
    for (long i = 0; i < num_elements; ++i) 
        ret[i] = i;

    QuickSortRecursive(ret, 0, num_elements-1, CmpFunc, QuickSortPartition);
    return ret;
}

void QuickSortRecursive(long* array, long l, long r, int (*CmpFunc)(long, long), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long, long))) {
    if (l < r) {
        long pivot = QuickSortPartition(array, l, r, CmpFunc);
        QuickSortRecursive(array, l, pivot-1, CmpFunc, QuickSortPartition);
        QuickSortRecursive(array, pivot+1, r, CmpFunc, QuickSortPartition);
    }
}

long SimpleQuickSortPartition(long* array, long l, long r, int (*CmpFunc)(long, long)) {
    long pivot = array[l];
    while (l < r) {
        while (l < r && CmpFunc(array[r], pivot) >= 0) {
            --r;
        }
        array[l] = array[r];
        while (l < r && CmpFunc(array[l], pivot) <= 0) {
            ++l;
        }
        array[r] = array[l];
    }
    array[l] = pivot;
    return l;
}
