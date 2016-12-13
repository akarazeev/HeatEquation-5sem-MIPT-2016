/**
 * main.c
 * Heat Equation
 *
 *          2016. MIPT
 * Created by <anton.karazeev@gmail.com> on 18/10/16.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <dirent.h>
#include "heat_equation.h"

static inline void alloc_worker() {
    if (nei1 != NONE && nei2 != NONE) {
        start_r = start_row-1;
        end_r = start_row + number_rows;
    } else if (nei1 != NONE) {
        start_r = start_row-1;
        end_r = start_row + number_rows - 1;
    } else {
        start_r = start_row;
        end_r = start_row + number_rows;
    }

    p = (double*) malloc(size * (end_r - start_r + 1) * sizeof(double));
    grid[0] = (double**) malloc(size * sizeof(double*));
    grid[1] = (double**) malloc(size * sizeof(double*));

    // printf("start_r %d end_r %d _ ", start_r, end_r);
    // printf("start_row %d number_rows %d\n", start_row, number_rows);

    int i;
    for (i = start_r; i <= end_r; ++i) {
        // printf("i %d  ", i);
        grid[0][i] = &(p[(i - start_r) * size]);
        grid[1][i] = &(p[(i - start_r) * size]);
    }
    // printf("taskid=%d", taskid);
}

static inline void free_worker() {
    free(p);
    free(grid[0]);
    free(grid[1]);
}

static inline void free_master() {
    free(p);
    free(grid[0]);
    free(grid[1]);
}

static inline double sec_deriv_x(int x, int y, double t) {

    //FIXME: DO ot right!!

    // if (x + dx >= size) {
    //     assert(x - dx >= 0);
    //     return (border(x + dx, y, t) - (2.0 * grid[0][x][y]) + grid[0][x - dx][y]) / (2.0 * dx * dx);
    // } else if (x - dx < 0) {
    //     assert(x + dx < size);
    //     return (grid[0][x + dx][y] - (2.0 * grid[0][x][y]) + border(x - dx, y, t)) / (2.0 * dx * dx);
    // } else {
    if (x + dx >= size || x - dx < 0) {
        return 0;
    } else {
        if (cur_iter == 1) {
            printf("-at x and y-> %f\n", grid[0][x][y]);
            printf("-----> %f\n", (grid[0][x + dx][y] - (2.0 * grid[0][x][y]) + grid[0][x - dx][y]) / (2.0 * dx * dx));
        }
        return (grid[0][x + dx][y] - (2.0 * grid[0][x][y]) + grid[0][x - dx][y]) / (2.0 * dx * dx);
    }
    // }
    // return 1;
}

static inline double sec_deriv_y(int x, int y, double t) {
    // if (y + dy >= size) {
    //     assert(y - dy >= 0);
    //     // printf("b1_ y: %d\n", y);
    //     return (border(x, y + dy, t) - (2.0 * grid[0][x][y]) + grid[0][x][y - dy]) / (2.0 * dy * dy);
    // } else if (y - dy < 0) {
    //     assert(y + dy < size);
    //     // printf("b2_ y: %d\n", y);
    //     return (grid[0][x][y + dy] - (2.0 * grid[0][x][y]) + border(x, y - dy, t)) / (2.0 * dy * dy);
    // } else {
    //     // printf("y: %d\n", y);
    //     return (grid[0][x][y + dy] - (2.0 * grid[0][x][y]) + grid[0][x][y - dy]) / (2.0 * dy * dy);
    // }
    return 0;
    // return 0;
}

static inline void update_grid() {
    // int i;
    // int cur_task_size = number_rows * size;
    // for (i = 0; i < cur_task_size; ++i) {
    //     int x = start_row + (i / size);
    //     int y = i % size;

    //     assert(y < size);
    //     assert(x < size);

    //     grid[0][x][y] = grid[1][x][y];
    // }
    int x;
    int y;
    // puts("NEXT");
    for (x = start_row; x < start_row + number_rows; ++x) {
        for (y = 0; y < size; ++y) {
            grid[0][x][y] = grid[1][x][y];
        }
    }
}

static inline void next_age() {
    int i;
    double cur_time = (float) cur_iter * dt;
    int cur_task_size = number_rows * size;
    // printf("> cur_task_size: %d\n", cur_task_size);
    // printf("> start_row end: %d %d\n", start_row, start_row + number_rows);
    // for (i = 0; i < cur_task_size; ++i) {
    //     int x = start_row + (i / size);
    //     int y = i % size;

    //     assert(y < size);
    //     assert(x < size);

    //     // printf("> position: %d %d\n", x, y);
    //     grid[1-iz][x][y] += dt * (heat(x, y, cur_time) +
    //                 (a_squared * (sec_deriv_y(x, y, cur_time) +
    //                 sec_deriv_x(x, y, cur_time))));
    // }
    int x;
    int y;
    for (x = start_row; x < start_row + number_rows; ++x) {
        for (y = 0; y < size; ++y) {
            grid[1][x][y] += dt * (heat(x, y, cur_time) +
                        (a_squared * (sec_deriv_y(x, y, cur_time) +
                        sec_deriv_x(x, y, cur_time))));
            if (cur_iter == 1) {
                printf(">  sec_deriv:%f\n", sec_deriv_x(x, y, cur_time));
                printf(">> %f x:%d y:%d %d\n", grid[1][x][y], x, y, worker_number);
            }
        }
    }
    // update_grid();
}

int main(int argc, char** argv) {
    const double a_squared = 0.1;

    int      size;
    int      period;
    int      iterations;
    int      cur_iter;
    int      taskid;
    int      dx;
    int      dy;
    double   dt;

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
    int   side = 0;

    int start_r;
    int end_r;
    int permit = 0;

    char   prefix[PREFIX_LEN];

    double** grid[2];

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    assert(numtasks >= 2);

    numworkers = numtasks - 1;

    printf("%d %d\n", numtasks, taskid);

    if (taskid == MASTER) {

        /* Master Code */

        clean_dumps();      /* remove all dump files */
        init_from_config(); /* read all parameters */

        int tmp = 1;
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

            /* send info */

            worker_number = i;
            MPI_Send(&worker_number, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&start_row,     1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&number_rows,   1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&nei1,          1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&nei2,          1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&iterations,    1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&period,        1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&size,          1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&dx,            1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&dy,            1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&dt,            1, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&prefix, PREFIX_LEN, MPI_CHAR, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&grid[0][start_row][0], size * number_rows, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);

            start_row += number_rows;
        }

        /* this two fors are totally NEW */

        puts("LLLLLOOOOOOLLLLL1");
        for (i = 1; i <= numworkers; ++i) {
            puts("LLLLLOOOOOOLLLLL2");
            // MPI_Recv(&worker_number, 1, MPI_INT, i, PERMIT1, MPI_COMM_WORLD, &status);
        }
        // for (cur_iter = 1; cur_iter < iterations + 1; ++cur_iter) {
        //     for (i = 1; i <= numworkers; ++i) {
        //         puts("LLLLLOOOOOOLLLLL2");
        //         MPI_Recv(&worker_number, 1, MPI_INT, i, PERMIT1, MPI_COMM_WORLD, &status);
        //     }
        //     puts("LLLLLOOOOOOLLLLL3");
        //     permit = 1;
        //     for (i = 1; i <= numworkers; ++i) {
        //         MPI_Send(&permit,        1, MPI_INT, i, PERMIT1, MPI_COMM_WORLD);
        //     }
        //     for (i = 1; i <= numworkers; ++i) {
        //         MPI_Recv(&worker_number, 1, MPI_INT, i, PERMIT2, MPI_COMM_WORLD, &status);
        //     }
        //     puts("LLLLLOOOOOOLLLLL");
        //     for (i = 1; i <= numworkers; ++i) {
        //         MPI_Send(&permit,        1, MPI_INT, i, PERMIT2, MPI_COMM_WORLD);
        //     }
        // }

        for (i = 1; i <= numworkers; ++i) {
            MPI_Recv(&start_row, 1, MPI_INT,                                 i, END, MPI_COMM_WORLD, &status);
            // MPI_Recv(&number_rows, 1, MPI_INT,                               i, END, MPI_COMM_WORLD, &status);
            // MPI_Recv(&grid[1][start_row][0], size * number_rows, MPI_DOUBLE, i, END, MPI_COMM_WORLD, &status);
        }

        /* it's time to free master's grid */

        puts("lol");
        free_master();

        puts("> MASTER is DONE");

        /* decode binary files into txt files */

        DIR *dp;
        struct dirent *ep;
        dp = opendir ("./dump_files/");
        if (dp != NULL) {
            while ((ep = readdir(dp))) {
                if (!(strcmp(ep->d_name, ".") == 0 || strcmp(ep->d_name, "..") == 0)) {
                    char file_path[4096];
                    strcpy(file_path, "./dump_files/");
                    strcat(file_path, ep->d_name);

                    puts(file_path);

                    FILE* f = fopen(file_path, "rb");
                    int* fir = malloc(2 * sizeof(int));
                    fread(fir, sizeof(int), 2, f);
                    // printf("%d %d\n", fir[0], fir[1]);
                    double* x = malloc(size * size * sizeof(double));
                    fread(x, sizeof(double), size * size, f);
                    fclose(f);

                    strcat(file_path, ".txt");
                    f = fopen(file_path, "w");
                    // fprintf(f, "%d\n%d\n", fir[0], fir[1]);
                    for (i = 0; i < size * size; ++i) {
                        fprintf(f, "%f ", x[i]);
                        // if ((i+1) % size == 0) {
                        //     printf("\n");
                        // }
                    }
                    free(x);
                    free(fir);
                    fclose(f);
                }
            }
            closedir (dp);
        } else {
            perror ("Couldn't open the directory");
        }
        /* end of master code */

    } else {

        /* Worker Code */

        int tmp = 0;
        MPI_Bcast(&tmp, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

        /* receive all info */

        MPI_Recv(&worker_number, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&start_row,     1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&number_rows,   1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&nei1,          1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&nei2,          1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&iterations,    1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&period,        1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&size,          1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&dx,            1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&dy,            1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&dt,            1, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&prefix, PREFIX_LEN, MPI_CHAR, MASTER, BEGIN, MPI_COMM_WORLD, &status);

        printf("> start_row: %d\n", start_row);

        MPI_Datatype grid_type;
        MPI_Type_contiguous(size * number_rows * sizeof(double), MPI_BYTE, &grid_type);
        MPI_Type_commit(&grid_type);

        alloc_worker();

        MPI_Recv(&grid[0][start_row][0], size * number_rows, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);

        for (cur_iter = 1; cur_iter < iterations + 1; ++cur_iter) {
            // printf("> cur_iter: %d\n", cur_iter);
            // printf(" 1 %d taskid=%d", cur_iter, taskid);
            puts("keks1");
            if (nei1 != NONE) {
                MPI_Send(&grid[0][start_row][0],   size, MPI_DOUBLE, nei1, NEI2, MPI_COMM_WORLD);
                MPI_Recv(&grid[0][start_row-1][0], size, MPI_DOUBLE, nei1, NEI1, MPI_COMM_WORLD, &status);
            }
            puts("keks2");
            if (nei2 != NONE) {
                MPI_Send(&grid[0][start_row + number_rows - 1][0], size, MPI_DOUBLE, nei2, NEI1, MPI_COMM_WORLD);
                MPI_Recv(&grid[0][start_row + number_rows][0],     size, MPI_DOUBLE, nei2, NEI2, MPI_COMM_WORLD, &status);
            }
            // printf(" 3 %d taskid=%d", cur_iter, taskid);
            puts("keks3");
            if (cur_iter % period == 0) {
                // printf(">> %d A\n", worker_number);
                dump_to_file(cur_iter, grid_type);
            }

            /* it's time for next age */

            next_age();

            // printf(">> %d B\n", worker_number);
            puts("keks4");
            printf("taskid: %d\n", taskid);
            MPI_Send(&worker_number, 1, MPI_INT, MASTER, PERMIT1, MPI_COMM_WORLD);
            // MPI_Recv(&permit,        1, MPI_INT, MASTER, PERMIT1, MPI_COMM_WORLD, &status);
            // MPI_Barrier(MPI_COMM_WORLD);
            update_grid();
            // MPI_Send(&worker_number, 1, MPI_INT, MASTER, PERMIT2, MPI_COMM_WORLD);
            // MPI_Recv(&permit,        1, MPI_INT, MASTER, PERMIT2, MPI_COMM_WORLD, &status);

        }

        /* send back to master */

        MPI_Send(&start_row, 1, MPI_INT, MASTER, END, MPI_COMM_WORLD);
        // MPI_Send(&number_rows, 1, MPI_INT, MASTER, END, MPI_COMM_WORLD);
        // MPI_Send(&grid[1][start_row][0], size * number_rows, MPI_DOUBLE, MASTER, END, MPI_COMM_WORLD);

        /* it's time to free things up */

        free_worker();

        printf("> WORKER #%d is DONE\n", worker_number);

        /* end of worker code */

    }

    if (taskid == MASTER) {
        puts("111111111");
    }

    MPI_Finalize();

    return 0;
}
