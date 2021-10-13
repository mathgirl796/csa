#ifndef DR_SORT_H
#define DR_SORT_H
#include <limits.h>

#define BINARY_SEARCH_NOT_FOUND LONG_MAX

/*
    sort an index array [0, 1, 2, ..., num_elements-1] by CmpFun

    input:
        start_pos: start index of elements you want to sort in your CmpFunc
        end_pos: similar to start_pos, elements whose index in [start_pos, end_pos) will be sort
        CmpFun: input left_pos l and right_pos r, return 0 if elements[l] == elements[r], -1 if <, 1 if >
            l and r should be meaningful for all intergers between 0 and num_elements-1
        QuickSortPartition: how to do detailed partition
            some functions are provide for this arg in this lib such as SimpleQuickSortPartion
    output:
        array of type LONG of size num_elements, represent sort result

    warning: please insure your start_pos and end_pos legal

    example:
        input 5, CmpFun is string compare of "shit$"'s suffixes
        output: [4, 2, 1, 0, 3], represents: ["$", "it$", "hit$", "shit$", "t$"]
*/
long* QuickSort(long start_pos, long end_pos, int (*CmpFunc)(long pos_1, long pos_2), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long pos_1, long pos_2)));

/*
    Provided for default occation of QuickSort
    在[l, r]正数区间内选出一个数作为下标去索引pivot，以pivot为中间值划分后，返回pivot所在的下标
*/
long SimpleQuickSortPartition(long* array, long l, long r, int (*CmpFunc)(long pos_1, long pos_2));

/*
    binary search target in your continuous poses [start_pos, end_pos) 
        using your compare function CmpFunc
    the poses must correspond to an "ascending" list of elements

    hint: you can trick in your CmpFunc to make this function do binary 
        search in a list of "descending" elements 

    input:
        start_pos: start index of elements you want to search in your CmpFunc
        end_pos: similar to start_pos, elements whose index in [start_pos, end_pos) will be search
        target: the element to be search
        CmpFunc: function to compare target and element at pos
            it need to return -1 if target < element at pos, 0 if = and 1 if >
    output:
        if found, return the found element's pos which is within [start_pos, end_pos)
        if not found, return a constant long int BINARY_SEARCH_NOT_FOUND, which is LONG_MAX in limits.h
        if multiple target exists, return the LEFT bound of the target pos region
*/
long BinarySearch(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos));

/*
    exactly the same as BinarySearch
*/
long BinarySearchLeftBound(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos));

/*
    slightly different with BinarySearch
    THAT IS: if multiple target exists, return the RIGHT bound of the target pos region
*/
long BinarySearchRightBound(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos));


/*
    slightly different with BinarySearch
    THAT IS: return the RIGHT bound of the pos that CmpFunc(target, pos) >= 0
*/
long BinarySeachRightBoundLessEqualTarget(long start_pos, long end_pos, const void* target, int (*CmpFunc)(const void* target, long pos));


#endif