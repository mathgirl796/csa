#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "utlib/uthash.h"
#include "dr_tools.h"
#include "dr_sort.h"
#include "dr_Hon.h"
#include "test.h"
#include "utlib/utstring.h"


int main() {

    char* filePath = "data/Ecoli_4x.fq1";
    char* testFilePath = "data/test.fq";

    struct read* reads;
    struct read* qualities;
    ReadFastq(testFilePath, &reads, &qualities);

    struct read *s, *tmp;
    HASH_ITER(hh, reads, s, tmp) {
        printf("%d\t%s\t%s\n", s->length, s->name, s->seq);
    }

    return 0;
}
