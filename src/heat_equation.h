/**
 * heat_equation.h
 * Heat Equation
 *
 *          2016. MIPT
 * Created by <anton.karazeev@gmail.com> on 20/10/16.
 *
 */

#ifndef MODULE_H_
#define MODULE_H_

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <ctype.h>
// #include <omp.h>
#include "/usr/local/opt/llvm/lib/clang/3.8.1/include/omp.h"
#include <math.h>
#include <dirent.h>
#include "border_functions.h"
#include "heat_functions.h"

#define MASTER  0
#define NONE   -1
#define BEGIN  41
#define END    42
#define NEI1    1
#define NEI2    2
#define PERMIT1 -2
#define PERMIT2 -3

#define PREFIX_LEN 10

double   (*border_functions[5]) (int x, int y, double t);
double   (*heat_functions[4]) (int x, int y, double t);

int    border_number;
int    heat_number;

// static inline double heat(int x, int y, double t);
// static inline double border(int x, int y, double t);
// static inline double init(int x, int y);
// static inline void   init_from_config();
// static inline double sec_deriv_x(int x, int y, double t);
// static inline double sec_deriv_y(int x, int y, double t);
// static inline void   next_age();
// static inline void   update_grid();
// static inline void   clean_dumps();
// static inline void   dump_to_file(double*** grid, int cur_iter, MPI_Datatype grid_type);
// static inline void   free_all();
// static inline void   free_master();
// static inline int    is_number(char* str);

/* -- Implementation -- */

static inline double heat(int x, int y, double t, int size) {
    // return (*heat_functions[heat_number]) (x, y, t);
    if (x == size / 2 && y == size / 2) {
        // return -1;
    }
    return 0;
}

static inline double border(int x, int y, double t, int size) {
    // return (*border_functions[border_number]) (x, y, t);
    return 0;
}

static inline void free_master(double** grid) {
    free(grid[0]);
    free(grid[1]);
}

static inline void dump_to_file(double*** grid, int start_row, int size, int number_rows, int cur_iter, MPI_Datatype grid_type) {
    MPI_File file;
    MPI_Status status;

    char cur_file_name[4096];
    char file_path[4096];

    strcpy(file_path, "/Users/AntonKarazeev/Desktop/Paral/3_heat_equation/heat_equation/dump_files/");
    strcat(file_path, prefix);

    sprintf(cur_file_name, "_%03d", (int) iter / period);
    strcat(file_path, cur_file_name);

    MPI_File_open(MPI_COMM_SELF, file_path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    // if (worker_number == MASTER + 1) {
    //     int tmp[2];
    //     tmp[0] = size;
    //     tmp[1] = iter;
    //     MPI_File_write_at(file, 0, &tmp, 2, MPI_INT, &status);
    // }
    // MPI_File_set_view(file, (2 * sizeof(int)) + (start_row * size * sizeof(double)), MPI_DOUBLE, grid_type, "native", MPI_INFO_NULL);

    MPI_File_set_view(file, sizeof(double) + start_row * size * sizeof(double), MPI_BYTE, grid_type, "native", MPI_INFO_NULL);

    MPI_File_write(file, &grid[0][start_row][0], size * number_rows, MPI_DOUBLE, &status);
    MPI_File_close(&file);
}

static inline double init(int x, int y) {
    // if (pow(x - (size/2.0), 2) + pow(y - (size/2.0), 2) < pow(size/5.0, 2)) {
    //     return 0;
    // } else {
    //     return 1;
    // }
    if (x == 3 || x == 2) {
        return 20;
    } else {
        return 10;
    }
}

// static inline void fill_array_of_functions() {
//     border_functions[0] = border1;
//     border_functions[1] = border2;
//     border_functions[2] = border3;
//     border_functions[3] = border4;
//     border_functions[4] = border5;

//     heat_functions[0] = heat1;
//     heat_functions[1] = heat2;
//     heat_functions[2] = heat3;
//     heat_functions[3] = heat4;
// }

static inline void init_from_config(double*** grid, int* size, int* dx, int* dy, double* dt, int* iterations, int* period, ) {
    puts("> init_from_config");

    // FIXME make this work
    // fill_array_of_functions();

    FILE* f_config = fopen("res/config", "r");
    assert(f_config != NULL);

    fscanf(f_config, "%d %d", &dx, &dy);
    fscanf(f_config, "%lf", &dt);
    fscanf(f_config, "%d", &iterations);
    fscanf(f_config, "%d", &period);

    assert(period < iterations);

    char init_file_name[20];

    fscanf(f_config, "%s", init_file_name);
    fscanf(f_config, "%s", prefix);
    fscanf(f_config, "%d", &border_number);
    fscanf(f_config, "%d", &heat_number);
    heat_number--;
    border_number--;

    double* p;
    if (is_number(init_file_name) == 0) {
        puts("> fill_from_file");

        char init_file_path[30];
        strcpy(init_file_path, "res/");
        strcat(init_file_path, init_file_name);

        FILE* f_init = fopen(init_file_path, "r");
        assert(f_init != NULL);
        fscanf(f_init, "%d", &size);

        p = (double*) malloc(size * size * sizeof(double));
        grid[0] = (double**) malloc(size * sizeof(double*));
        grid[1] = (double**) malloc(size * sizeof(double*));

        int i;
        int j;
        for (i = 0; i < size; ++i) {
            grid[0][i] = &(p[i * size]);
            grid[1][i] = &(p[i * size]);
            for (j = 0; j < size; ++j) {
                double tmp;
                fscanf(f_init, "%lf", &tmp);
                grid[0][i][j] = tmp;
                grid[1][i][j] = tmp;
            }
        }
    } else {
        puts("> fill_with_func");

        size = atoi(init_file_name);

        assert(size > 0);

        p = (double*) malloc(size * size * sizeof(double));
        grid[0] = (double**) malloc(size * sizeof(double*));
        grid[1] = (double**) malloc(size * sizeof(double*));

        int k;
        int j;
        for (k = 0; k < size; ++k) {
            grid[0][k] = &(p[k * size]);
            grid[1][k] = &(p[k * size]);
            for (j = 0; j < size; ++j) {
                grid[0][k][j] = init(k, j);
                grid[1][k][j] = grid[0][k][j];
                // printf("pos: %d %d _ %f %f\n", k, j, grid[0][k][j], grid[1][k][j]);
            }
        }
    }

    puts  ("+------*------+");
    printf("| Iterations: %d\n", iterations);
    printf("| Period: %d\n", period);
    printf("| Size: %d\n", size);
    puts  ("+------*------+");
}

// static inline double sec_deriv_x(int x, int y, double t) {
//     if (x + dx >= size) {
//         assert(x - dx >= 0);
//         return (border(x + dx, y, t) - (2.0 * grid[0][x][y]) + grid[0][x - dx][y]) / (2.0 * dx * dx);
//     } else if (x - dx <= 0) {
//         assert(x + dx < size);
//         return (grid[0][x + dx][y] - (2.0 * grid[0][x][y]) + border(x - dx, y, t)) / (2.0 * dx * dx);
//     } else {
//         return (grid[0][x + dx][y] - (2.0 * grid[0][x][y]) + grid[0][x - dx][y]) / (2.0 * dx * dx);
//     }
// }

// static inline double sec_deriv_y(int x, int y, double t) {
//     if (y + dy >= size) {
//         assert(y - dy >= 0);
//         return (border(x, y + dy, t) - (2.0 * grid[0][x][y]) + grid[0][x][y - dy]) / (2.0 * dy * dy);
//     } else if (y - dy <= 0) {
//         assert(y + dy < size);
//         return (grid[0][x][y + dy] - (2.0 * grid[0][x][y]) + border(x, y - dy, t)) / (2.0 * dy * dy);
//     } else {
//         return (grid[0][x][y + dy] - (2.0 * grid[0][x][y]) + grid[0][x][y - dy]) / (2.0 * dy * dy);
//     }
// }

// static inline void next_age() {
//     int i;
//     #pragma omp parallel for
//     for (i = 0; i < size * size; ++i) {
//         int x = i / size;
//         int y = i % size;
//         grid[1][x][y] += dt * (heat(x, y, cur_time) +
//                     (a_squared * (sec_deriv_y(x, y, cur_time) +
//                      sec_deriv_x(x, y, cur_time))));
//     }
//     update_grid();
//     cur_time += dt;
// }

// static inline void update_grid() {
//     int i;
//     #pragma omp parallel for
//     for (i = 0; i < size * size; ++i) {
//         int x = i / size;
//         int y = i % size;
//         grid[0][x][y] = grid[1][x][y];
//     }
// }

static inline void clean_dumps() {
    puts("> clean_dumps");
    DIR* dp;
    struct dirent* ep;
    char* path = "dump_files/";

    dp = opendir (path);
    if (dp != NULL) {
        while ((ep = readdir(dp))) {
            if (strcmp(ep->d_name, ".") != 0 && strcmp(ep->d_name, "..") != 0) {
                char cur_path[30];
                strcpy(cur_path, path);
                strcat(cur_path, ep->d_name);
                remove(cur_path);
            }
        }
        closedir (dp);
    }
    else {
        perror ("Couldn't open the directory");
    }
}

static inline int is_number(char* str) {
    int i;
    for (i = 0; i < strlen(str); ++i) {
        if (isnumber(str[i]) == 0) {
            return 0;
        }
    }
    return 1;
}

#endif
