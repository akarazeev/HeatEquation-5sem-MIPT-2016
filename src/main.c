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

#define MASTER  0
#define NONE   -1
#define BEGIN  41
#define END    42
#define NEI1    1
#define NEI2    2
#define PREFIX_LEN 10

#define DEBUG 0

int    numworkers;
int    numtasks;   /* number of tasks */
int    min_number_rows;
int    extra_rows;
int    number_rows;
int    start_row;
int    worker_number;
int    source;
int    nei1;
int    nei2;
int    size;
int    period;
int    iterations;
int    cur_iter;
int    taskid;
int    dx;
int    dy;
int    start_r;
int    end_r;
int    side = 0;
double dt;

char     prefix[PREFIX_LEN];
double*  p1;
double*  p2;
double** grid[2];

const double a_squared = 0.1;

#include "heat_equation.h"

static inline void dump_to_file(int iter, MPI_Datatype grid_type) {
    MPI_File file;
    MPI_Status status;

    char cur_file_name[4096];
    char file_path[4096];

    strcpy(file_path, "dump_files/");
    strcat(file_path, prefix);

    assert(period != 0);
    sprintf(cur_file_name, "_%03d", (int) iter / period);
    strcat(file_path, cur_file_name);

    MPI_File_open(MPI_COMM_SELF, file_path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    // if (worker_number == MASTER + 1) {
    //     int tmp[2];
    //     tmp[0] = size;
    //     tmp[1] = iter;
    //     MPI_File_write_at(file, 0, &tmp, 2, MPI_INT, &status);
    // }
    // printf("> taskid=%d DUMPING2\n", taskid);
    // MPI_File_set_view(file, (2 * sizeof(int)) + (start_row * size * sizeof(double)), MPI_DOUBLE, grid_type, "native", MPI_INFO_NULL);

    MPI_File_set_view(file, sizeof(double) + start_row * size * sizeof(double), MPI_BYTE, grid_type, "native", MPI_INFO_NULL);

    MPI_File_write(file, &grid[side][start_row][0], size * number_rows, MPI_DOUBLE, &status);
    MPI_File_close(&file);
}

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

    p1 = (double*) malloc(size * (end_r - start_r + 1) * sizeof(double));
    p2 = (double*) malloc(size * (end_r - start_r + 1) * sizeof(double));
    grid[0] = (double**) malloc(size * sizeof(double*));
    grid[1] = (double**) malloc(size * sizeof(double*));

    int i;
    for (i = start_r; i <= end_r; ++i) {
        grid[0][i] = &(p1[(i - start_r) * size]);
        grid[1][i] = &(p2[(i - start_r) * size]);
    }
}

static inline void free_worker() {
    free(p1);
    free(p2);
    free(grid[0]);
    free(grid[1]);
}

static inline void free_master() {
    free(p1);
    free(p2);
    free(grid[0]);
    free(grid[1]);
}

static inline double heat(int x, int y, double t) {
    if (x == size / 2 && y == size / 2) {
        return -1;
    } else {
        return 1;
    }
}

static inline double border(int x, int y, double t) {
    return 2;
}

static inline double sec_deriv_x(int x, int y, double t) {
    if (x + dx >= size) {
        assert(x - dx >= 0);
        return (border(x + dx, y, t) - (2.0 * grid[side][x][y]) + grid[side][x - dx][y]) / (2.0 * dx * dx);
    } else if (x - dx < 0) {
        assert(x + dx < size);
        return (grid[side][x + dx][y] - (2.0 * grid[side][x][y]) + border(x - dx, y, t)) / (2.0 * dx * dx);
    } else {
        return (grid[side][x + dx][y] - (2.0 * grid[side][x][y]) + grid[side][x - dx][y]) / (2.0 * dx * dx);
    }
}

static inline double sec_deriv_y(int x, int y, double t) {
    if (y + dy >= size) {
        assert(y - dy >= 0);
        return (border(x, y + dy, t) - (2.0 * grid[side][x][y]) + grid[side][x][y - dy]) / (2.0 * dy * dy);
    } else if (y - dy < 0) {
        assert(y + dy < size);
        return (grid[side][x][y + dy] - (2.0 * grid[side][x][y]) + border(x, y - dy, t)) / (2.0 * dy * dy);
    } else {
        return (grid[side][x][y + dy] - (2.0 * grid[side][x][y]) + grid[side][x][y - dy]) / (2.0 * dy * dy);
    }
}

static inline void next_age() {
    double cur_time = (float) cur_iter * dt;
    int cur_task_size = number_rows * size;
#if DEBUG
    printf("> cur_task_size: %d\n", cur_task_size);
    printf("> start_row end: %d %d\n", start_row, start_row + number_rows);
#endif

    int i;
    #pragma omp parallel for
    for (i = 0; i < cur_task_size; ++i) {
        int x = start_row + (i / size);
        int y = i % size;

        assert(y < size);
        assert(x < size);

        grid[1-side][x][y] = grid[side][x][y] + (dt * (heat(x, y, cur_time) +
                            (a_squared * (sec_deriv_y(x, y, cur_time) +
                            sec_deriv_x(x, y, cur_time)))));
    }
}

int main(int argc, char** argv) {
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    double begin = MPI_Wtime();

    assert(numtasks >= 2);

    numworkers = numtasks - 1;

#if DEBUG
    printf("%d %d\n", numtasks, taskid);
#endif

    if (taskid == MASTER) {

        /* Master Code */

        clean_dumps();      /* remove all dump files */
        init_from_config(); /* read all parameters */

#if DEBUG
        printf("%d %d\n", numworkers, size);
#endif
        assert(numworkers <= size);

        int tmp = 1;
        MPI_Bcast(&tmp, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

        min_number_rows = size / numworkers;
        extra_rows = size % numworkers;
        printf("min_number_rows:%d extra_rows:%d\n", min_number_rows, extra_rows);
        start_row = 0;

        int i;
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
            MPI_Send(&grid[side][start_row][0], size * number_rows, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);

            start_row += number_rows;
        }

        for (i = 1; i <= numworkers; ++i) {
            MPI_Recv(&start_row, 1, MPI_INT, i, END, MPI_COMM_WORLD, &status);
        }

        /* it's time to free master's grid */

        free_master();

#if DEBUG
        puts("> MASTER is DONE");
#endif

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

#if DEBUG
                    puts(file_path);
#endif

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
                    int i;
                    for (i = 0; i < size * size; ++i) {
                        fprintf(f, "%f ", x[i]);
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

#if DEBUG
        printf("> start_row: %d number_rows: %d\n", start_row, number_rows);
#endif

        MPI_Datatype grid_type;
        MPI_Type_contiguous(size * number_rows * sizeof(double), MPI_BYTE, &grid_type);
        MPI_Type_commit(&grid_type);

        alloc_worker();

        MPI_Recv(&grid[side][start_row][0], size * number_rows, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        for (cur_iter = 1; cur_iter < iterations + 1; ++cur_iter) {
            if (nei1 != NONE) {
                MPI_Send(&grid[side][start_row][0],   size, MPI_DOUBLE, nei1, NEI2, MPI_COMM_WORLD);
                MPI_Recv(&grid[side][start_row-1][0], size, MPI_DOUBLE, nei1, NEI1, MPI_COMM_WORLD, &status);
            }
            if (nei2 != NONE) {
                MPI_Recv(&grid[side][start_row + number_rows][0],     size, MPI_DOUBLE, nei2, NEI2, MPI_COMM_WORLD, &status);
                MPI_Send(&grid[side][start_row + number_rows - 1][0], size, MPI_DOUBLE, nei2, NEI1, MPI_COMM_WORLD);
            }
            if (cur_iter % period == 0) {
#if DEBUG
                printf("> cur_iter: %d\n", cur_iter);
#endif
                dump_to_file(cur_iter, grid_type);
            }

            /* it's time for next age */

            next_age();

            side = 1 - side;
        }

        /* send back to master */

        MPI_Send(&start_row, 1, MPI_INT, MASTER, END, MPI_COMM_WORLD);

        /* it's time to free things up */

        free_worker();

#if DEBUG
        printf("> WORKER #%d is DONE\n", worker_number);
#endif

        /* end of worker code */

    }

    double end = MPI_Wtime();

    if (taskid == MASTER) {
        printf("TIME: %lf\n", end - begin);
    }

    MPI_Finalize();

    return 0;
}
