#include "dr_sort.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>


/* .c文件内部函数声明 */
void QuickSortRecursive(long* array, long l, long r, int (*CmpFunc)(long, long), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long, long)));


long* QuickSort(long start_pos, long end_pos, int (*CmpFunc)(long, long), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long, long))) {
    /* 为将要返回的数组申请内存、初始化 */
    long num_elements = end_pos - start_pos;
    long* ret = malloc(num_elements * sizeof(long));
    for (long i = start_pos; i < end_pos; ++i) {
        ret[i - start_pos] = i;
    }
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

long BinarySearch(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos)) {
    return BinarySearchLeftBound(start_pos, end_pos, target, CmpFunc);
}

long BinarySearchLeftBound(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos)) {
    if (end_pos <= start_pos) return BINARY_SEARCH_NOT_FOUND;

    long left = start_pos;
    long right = end_pos;
    while (left < right) {
        long mid = left + (right - left) / 2;
        if (CmpFunc(target, mid) == 0) right = mid;
        else if (CmpFunc(target, mid) > 0) left = mid + 1;
        else if (CmpFunc(target, mid) < 0) right = mid;
    }
    
    if (left == end_pos || CmpFunc(target, left) !=0) return BINARY_SEARCH_NOT_FOUND;
    return left;
}

long BinarySearchRightBound(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos)) {
    if (end_pos <= start_pos) return BINARY_SEARCH_NOT_FOUND;

    long left = start_pos;
    long right = end_pos;
    while (left < right) {
        long mid = left + (right - left) / 2;
        if (CmpFunc(target, mid) == 0) left = mid + 1;
        else if (CmpFunc(target, mid) > 0) left = mid + 1;
        else if (CmpFunc(target, mid) < 0) right = mid;
    }
    
    if (left == start_pos || CmpFunc(target, left - 1) !=0) return BINARY_SEARCH_NOT_FOUND;
    return left - 1;
}


long BinarySeachRightBoundLessEqualTarget(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos)) {
    if (CmpFunc(target, start_pos) < 0) return BINARY_SEARCH_NOT_FOUND;

    long left = start_pos;
    long right = end_pos;
    while (left < right) {
        long mid = left + (right - left) / 2;
        if (CmpFunc(target, mid) == 0) left = mid + 1;
        else if (CmpFunc(target, mid) > 0) left = mid + 1;
        else if (CmpFunc(target, mid) < 0) right = mid;
    }

    return left - 1;
}
