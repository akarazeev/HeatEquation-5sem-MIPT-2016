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

double   (*border_functions[5]) (int x, int y, double t);
double   (*heat_functions[4]) (int x, int y, double t);

int    border_number;
int    heat_number;

static inline double init(int x, int y);
static inline void   init_from_config();
static inline double sec_deriv_x(int x, int y, double t);
static inline double sec_deriv_y(int x, int y, double t);
static inline void   clean_dumps();
static inline void   free_all();
static inline int    is_number(char* str);

#include "border_functions.h"
#include "heat_functions.h"

/* -- Implementation -- */

// static inline double heat(int x, int y, double t) {
//     return (*heat_functions[heat_number]) (x, y, t);
// }
//
// static inline double border(int x, int y, double t) {
//     return (*border_functions[border_number]) (x, y, t);
// }

static inline double init(int x, int y) {
    if (pow(x - (size/2.0), 2) + pow(y - (size/2.0), 2) < pow(size/5.0, 2)) {
        return 10;
    } else {
        return -1;
    }
}

static inline void fill_array_of_functions() {
    border_functions[0] = border1;
    border_functions[1] = border2;
    border_functions[2] = border3;
    border_functions[3] = border4;
    border_functions[4] = border5;

    heat_functions[0] = heat1;
    heat_functions[1] = heat2;
    heat_functions[2] = heat3;
    heat_functions[3] = heat4;
}

static inline void init_from_config() {
    puts("> init_from_config");

    // FIXME make it work!
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

    if (is_number(init_file_name) == 0) {
        puts("> fill_from_file");

        char init_file_path[30];
        strcpy(init_file_path, "res/");
        strcat(init_file_path, init_file_name);

        FILE* f_init = fopen(init_file_path, "r");
        assert(f_init != NULL);
        fscanf(f_init, "%d", &size);

        p1 = (double*) malloc(size * size * sizeof(double));
        p2 = (double*) malloc(size * size * sizeof(double));
        grid[0] = (double**) malloc(size * sizeof(double*));
        grid[1] = (double**) malloc(size * sizeof(double*));

        int i;
        for (i = 0; i < size; ++i) {
            grid[0][i] = &(p1[i * size]);
            grid[1][i] = &(p2[i * size]);
            int j;
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

        p1 = (double*) malloc(size * size * sizeof(double));
        p2 = (double*) malloc(size * size * sizeof(double));
        grid[0] = (double**) malloc(size * sizeof(double*));
        grid[1] = (double**) malloc(size * sizeof(double*));

        int k;
        #pragma omp parallel for
        for (k = 0; k < size; ++k) {
            grid[0][k] = &(p1[k * size]);
            grid[1][k] = &(p2[k * size]);
            int j;
            for (j = 0; j < size; ++j) {
                grid[0][k][j] = init(k, j);
                grid[1][k][j] = grid[0][k][j];
            }
        }
    }

    puts  ("+------*------+");
    printf("| Iterations: %d\n", iterations);
    printf("| Period: %d\n", period);
    printf("| Size: %d\n", size);
    puts  ("+------*------+");
}

static inline void clean_dumps() {
    puts("> clean_dumps");
    DIR* dp;
    struct dirent* ep;
    char* path = "dump_files/";

    dp = opendir (path);
    if (dp != NULL) {
        while ((ep = readdir (dp))) {
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
