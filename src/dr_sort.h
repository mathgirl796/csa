#ifndef DR_SORT_H
#define DR_SORT_H


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
long* QuickSort(long start_pos, long end_pos, int (*CmpFunc)(long, long), long (*QuickSortPartition)(long* array, long l, long r, int (*CmpFunc)(long, long)));

/*
    Provided for default occation of QuickSort
    在[l, r]正数区间内选出一个数作为下标去索引pivot，以pivot为中间值划分后，返回pivot所在的下标
*/
long SimpleQuickSortPartition(long* array, long l, long r, int (*CmpFunc)(long, long));






#endif