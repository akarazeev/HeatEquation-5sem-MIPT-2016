//
//  main.c
//  HeatEquation
//
//  Created by <anton.karazeev@gmail.com> on 18/10/16. MIPT.
//  Copyright © 2016 Anton Karazeev. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include "heat_equation.h"

#define DEBUG 0

int main(int argc, char** argv) {
    struct timeval t1, t2;
	gettimeofday(&t1, NULL);

    int i;
    clean_dumps();

    init_from_config();
    for (i = 1; i < iterations + 1; ++i) {
        if (i % period == 0) {
            printf("> cur_iter: %d\n", i);
            dump_to_file(i);
        }
        next_age();
    }

    gettimeofday(&t2, NULL);
    unsigned long int tv1 = t1.tv_sec*1e6 + t1.tv_usec;
	unsigned long int tv2 = t2.tv_sec*1e6 + t2.tv_usec;
	printf("TIME: %f\n", (tv2 - tv1)/(float)1e6);

    free_all();

    return 0;
}
