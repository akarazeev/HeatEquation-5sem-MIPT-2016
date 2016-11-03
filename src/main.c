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
 * - total_iterations
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

#define MASTER  0
#define NONE   -1
#define BEGIN  41
#define END    42
#define NEI1    1
#define NEI2    2

int      size;
int      period;
int      iterations;
int      taskid;
double   dt;
double   cur_iter = 0.0;

int   numworkers;
int   numtasks;   /* number of tasks */
int   min_number_rows;
int   extra_rows;
int   number_rows;
int   start_row;
int   worker_number;
int   source;
int   nei1;
int   nei2;

char   prefix[10];

double** grid[2];

#include "heat_equation.h"

static inline void dump_to_file(int iter) {
    MPI_File file;
    MPI_Status status;
    MPI_Datatype grid_type;

    MPI_Type_contiguous(size * number_rows, MPI_DOUBLE, &grid_type);
    MPI_Type_commit(&grid_type);

    char cur_file_name[4096];
    char file_path[4096];

    strcpy(file_path, "/Users/AntonKarazeev/Desktop/Paral/3_heat_equation/heat_equation/dump_files/");
    strcat(file_path, prefix);

    assert(period != 0);
    sprintf(cur_file_name, "_%03d", (int) iter / period);
    strcat(cur_file_name, ".txt");

    strcat(file_path, cur_file_name);

    puts(file_path);
    MPI_File_open(MPI_COMM_SELF, file_path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    if (taskid == MASTER + 1) {
        int tmp[2];
        tmp[0] = size;
        tmp[1] = iter;
        MPI_File_write_at(file, 0, &tmp, 2, MPI_INT, &status);
    }

    MPI_File_set_view(file, (2 * sizeof(int)) + (start_row * size * sizeof(double)), MPI_DOUBLE, grid_type, "native", MPI_INFO_NULL);
    int i;
    for (i = 0; i < number_rows; ++i) {
        MPI_File_write(file, &grid[0][start_row + i], size, MPI_DOUBLE, &status);
    }
    MPI_File_close(&file);
}

int main(int argc, char** argv) {
    int   i;          /* loop variable */

    MPI_Status status;

    int** test[2];
    for (i = 0; i < 2; ++i) {
        test[i] = malloc(2 * sizeof(int*));
    }
    for (i = 0; i < 2; ++i) {
        test[0][i] = malloc(sizeof(int));
        test[1][i] = malloc(sizeof(int));
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    assert(numtasks >= 2);

    numworkers = numtasks - 1;

    printf("%d %d\n", numtasks, taskid);

    if (taskid == MASTER) {

        /* Master Code */

        int tmp = 1;
        clean_dumps();      /* remove all dump files */
        init_from_config(); /* read all parameters */

        MPI_Bcast(&tmp, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

        min_number_rows = size / numworkers;
        extra_rows = size % min_number_rows;
        start_row = 0;

        for (i = 1; i <= numworkers; ++i) {
            number_rows = (i <= extra_rows) ? min_number_rows + 1 : min_number_rows;

            if (i == 1) {
                nei1 = NONE;
            } else  {
                nei1 = i - 1;
            }

            if (i == numworkers) {
                nei2 = NONE;
            } else {
                nei2 = i + 1;
            }

            // /* send info */
            worker_number = i;
            MPI_Send(&worker_number, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&start_row,     1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&number_rows,   1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&nei1,          1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&nei2,          1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&iterations,    1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&period,        1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&size,          1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&prefix,       10, MPI_CHAR, i, BEGIN, MPI_COMM_WORLD);

            MPI_Send(&grid[0][start_row][0], size, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);

            start_row += number_rows;
        }

        for (i = 1; i <= numworkers; ++i) {
            MPI_Recv(&start_row, 1, MPI_INT, i, END, MPI_COMM_WORLD, &status);
        }

        /* end of master code */

    } else {

        /* Worker Code */

        int tmp = 0;
        printf("tmp1: %d\n", tmp);
        MPI_Bcast(&tmp, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        printf("tmp2: %d\n", tmp);

        /* receive all info */
        MPI_Recv(&worker_number, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&start_row, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&number_rows, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&nei1, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&nei2, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&iterations, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&period, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&size, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&prefix, 10, MPI_CHAR, MASTER, BEGIN, MPI_COMM_WORLD, &status);

        printf("> start_row: %d\n", start_row);
        printf("size: %d\n", size);

        int k;
        grid[0] = malloc(size * sizeof(double*));
        grid[1] = malloc(size * sizeof(double*));

        for (k = start_row-1; k <= start_row + number_rows; ++k) {
            grid[0][k] = malloc(size * sizeof(double));
            grid[1][k] = malloc(size * sizeof(double));
        }

        printf(":( %d\n", taskid);
        MPI_Recv(&grid[0][start_row][0], size, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        printf(":) %d\n", taskid);

        int cur_iter;
        for (cur_iter = 1; cur_iter < iterations + 1; ++cur_iter) {
            printf("> cur_iter: %d\n", cur_iter);
            if (nei1 != NONE) {
                printf("> start_row1: %d\n", start_row);
                MPI_Send(&grid[0][start_row][0], size, MPI_DOUBLE, nei1, NEI2, MPI_COMM_WORLD);
                MPI_Recv(&grid[0][start_row-1][0], size, MPI_DOUBLE, nei1, NEI1, MPI_COMM_WORLD, &status);
            }
            if (nei2 != NONE) {
                printf("> start_row2: %d\n", start_row);
                MPI_Send(&grid[0][start_row + number_rows - 1][0], size, MPI_DOUBLE, nei2, NEI1, MPI_COMM_WORLD);
                MPI_Recv(&grid[0][start_row + number_rows][0], size, MPI_DOUBLE, nei2, NEI2, MPI_COMM_WORLD, &status);
            }

            if (cur_iter % period == 0) {
                dump_to_file(cur_iter);
            }

            //TODO: add update

        }
        /* end of worker code */

        MPI_Send(&start_row, 1, MPI_INT, MASTER, END, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    free_all();

    return 0;
}
