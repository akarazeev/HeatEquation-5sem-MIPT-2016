/**
 * main.c
 * Heat Equation
 *
 *          2016. MIPT
 * Created by <anton.karazeev@gmail.com> on 18/10/16.
 *
 *
 * Square Matrix:
 *   *-----------*
 *   |           |
 *   |           |
 *   |           |
 *   |           |
 *   |           |
 *   *-----------*
 *   <-- size  -->
 *
 * Passing arguments:
 * - dx, dy
 * - dt
 * - period
 * - file with initial matrix
 * - <prefix> of dump files: "<prefix>_iteration.txt"
 * - number of _border_ function
 * - number of _heat_ function
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "heat_equation.h"

int main(int argc, char** argv) {
    int i;
    clean_dumps();

    init_from_file();
    // init_from_func();
    for (i = 1; i < iterations + 1; ++i) {
        if (i % period == 0) {
            dump_to_file(i);
            mpi_dump(i);
        }
        next_age();
        // printf(">%d\n", i);
    }

    free_all();

    return 0;
}
