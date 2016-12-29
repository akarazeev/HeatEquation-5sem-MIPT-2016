//
//  heat_equation.h
//  HeatEquation
//
//  Created by <anton.karazeev@gmail.com> on 20/10/16. MIPT.
//  Copyright © 2016 Anton Karazeev. All rights reserved.
//

#ifndef MODULE_H_
#define MODULE_H_

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <ctype.h>
#include <omp.h>
#include <math.h>
#include <dirent.h>

#define trace(x) printf("%d\n", x);

double** grid;
double** grid_old;
double (*border_functions[5]) (int x, int y, double t);
double (*heat_functions[4]) (int x, int y, double t);

int dx;
int dy;
int size = -1;
int period;
int iterations;
int border_number;
int heat_number;
double dt;
char prefix[10];
double cur_time = 0.0;

const double a_squared = 0.1;

static inline double heat(int x, int y, double t);
static inline double border(int x, int y, double t);
static inline double init(int x, int y);
static inline void init_from_config();
static inline double sec_deriv_x(int x, int y, double t);
static inline double sec_deriv_y(int x, int y, double t);
static inline void next_age();
static inline void update_grid();
static inline void clean_dumps();
static inline void free_all();
static inline void dump_to_file(int cur_iteration);
static inline int is_number(char* str);

#include "border_functions.h"
#include "heat_functions.h"

/* -- Implementation -- */

static inline double heat(int x, int y, double t) {
    return (*heat_functions[heat_number]) (x, y, t);
}

static inline double border(int x, int y, double t) {
    return (*border_functions[border_number]) (x, y, t);
}

static inline double init(int x, int y) {
    // if (pow(x - (size/2.0), 2) + pow(y - (size/2.0), 2) < pow(size/5.0, 2)) {
    //     return 0;
    // } else {
    //     return 1;
    // }
    return 1;
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

    fill_array_of_functions();

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

        int i;
        grid = malloc(size * sizeof(double*));
        grid_old = malloc(size * sizeof(double*));
        for (i = 0; i < size; ++i) {
            grid[i] = malloc(size * sizeof(double));
            grid_old[i] = malloc(size * sizeof(double));
            int j;
            for (j = 0; j < size; ++j) {
                double tmp;
                fscanf(f_init, "%lf", &tmp);
                grid[i][j] = tmp;
                grid_old[i][j] = tmp;
            }
        }
    } else {
        puts("> fill_with_func");

        size = atoi(init_file_name);
        assert(size > 0);

        int i;
        grid = malloc(size * sizeof(double*));
        grid_old = malloc(size * sizeof(double*));
        for (i = 0; i < size; ++i) {
            grid[i] = malloc(size * sizeof(double));
            grid_old[i] = malloc(size * sizeof(double));
            int j;
            for (j = 0; j < size; ++j) {
                grid[i][j] = init(i, j);
                grid_old[i][j] = grid[i][j];
            }
        }
    }

    puts("+------*------+");
    printf("| Iterations: %d\n", iterations);
    printf("| Period: %d\n", period);
    printf("| Size: %d\n", size);
    puts("+------*------+");
}

static inline double sec_deriv_x(int x, int y, double t) {
    if (x + dx >= size) {
        assert(x - dx >= 0);
        return (border(x + dx, y, t) - (2.0 * grid_old[x][y]) + grid_old[x - dx][y]) / (2.0 * dx * dx);
    } else if (x - dx <= 0) {
        assert(x + dx < size);
        return (grid_old[x + dx][y] - (2.0 * grid_old[x][y]) + border(x - dx, y, t)) / (2.0 * dx * dx);
    } else {
        return (grid_old[x + dx][y] - (2.0 * grid_old[x][y]) + grid_old[x - dx][y]) / (2.0 * dx * dx);
    }
}

static inline double sec_deriv_y(int x, int y, double t) {
    if (y + dy >= size) {
        assert(y - dy >= 0);
        return (border(x, y + dy, t) - (2.0 * grid_old[x][y]) + grid_old[x][y - dy]) / (2.0 * dy * dy);
    } else if (y - dy <= 0) {
        assert(y + dy < size);
        return (grid_old[x][y + dy] - (2.0 * grid_old[x][y]) + border(x, y - dy, t)) / (2.0 * dy * dy);
    } else {
        return (grid_old[x][y + dy] - (2.0 * grid_old[x][y]) + grid_old[x][y - dy]) / (2.0 * dy * dy);
    }
}

static inline void next_age() {
    int i;
    #pragma omp parallel for
    for (i = 0; i < size * size; ++i) {
        int x = i / size;
        int y = i % size;
        grid[x][y] += dt * (heat(x, y, cur_time) +
                    (a_squared * (sec_deriv_y(x, y, cur_time) +
                    sec_deriv_x(x, y, cur_time))));
    }
    update_grid();
    cur_time += dt;
}

static inline void update_grid() {
    int i;
    #pragma omp parallel for
    for (i = 0; i < size * size; ++i) {
        int x = i / size;
        int y = i % size;
        grid_old[x][y] = grid[x][y];
    }
}

static inline void free_all() {
    int i;
    #pragma omp parallel for
    for (i = 0; i < size; ++i) {
        free(grid[i]);
        free(grid_old[i]);
    }
    free(grid);
    free(grid_old);
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

static inline void dump_to_file(int cur_iteration) {
#if DEBUG
    printf("> dump_to_file #%d/%d\n", cur_iteration, iterations);
#endif
    assert(ceil(log10(cur_iteration)) < 10);

    char cur_file_name[30];
    char file_path[50];

    strcpy(file_path, "dump_files/");
    strcat(file_path, prefix);

    sprintf(cur_file_name, "_%03d", (int) cur_iteration / period);
    strcat(cur_file_name, ".txt");

    strcat(file_path, cur_file_name);

    FILE* f = fopen(file_path, "w");
    assert(f != NULL);

    fprintf(f, "%d\n", size);
    fprintf(f, "%d\n", cur_iteration);
    int i;

    fprintf(f, "%f", grid[0][0]);
    for (i = 1; i < size * size; ++i) {
        fprintf(f, " %f", grid[i / size][i % size]);
    }

    fclose(f);
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
