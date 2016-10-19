/**
 * main.c
 * heat_equation
 *
 *          2016. MIPT
 * Created by <anton.karazeev@gmail.com> on 18/10/16.
 *
 * History:
 * 18.10 - created
 *
 *   *-----------*
 *   |           |
 *   |           |
 *   |           |
 *   |           |
 *   |           |
 *   *-----------*
 *   <-- size  -->
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <dirent.h>

#define trace(x) printf("%d\n", x);

double** grid;
int size = 80;
double cur_time = 0.0;
const double dt = 0.1;
const double temp = 100.0;
const double a_squared = 10.0;
const int dx = 1;
const int dy = 2;
const int period = 10;
const int iterations = 200;

static inline double heat(int x, int y, double t) {
    return 0;
}

static inline double border(int x, int y, double t) {
    return ((float)x/(float)size) * 100;
}

static inline double init(int x, int y) {
    if (pow(x - (size/2.0), 2) + pow(y - (size/2.0), 2) < pow(size/5.0, 2)) {
        return 100;
    } else {
        return 1;
    }
}

static inline void init_from_func() {
    puts("> init_from_func");

    int i;
    grid = malloc(size * sizeof(double*));
    for (i = 0; i < size; ++i) {
        grid[i] = malloc(size * sizeof(double));
        int j;
        for (int j = 0; j < size; ++j) {
            grid[i][j] = init(i, j);
        }
    }
}

static inline void init_from_file() {
    puts("> init_from_file");

    FILE* f = fopen("hnu.txt", "r");
    assert(f != NULL);
    fscanf(f, "%d", &size);
    printf("size: %d\n", size);

    int i;
    grid = malloc(size * sizeof(double*));
    for (i = 0; i < size; ++i) {
        grid[i] = malloc(size * sizeof(double));
        int j;
        for (int j = 0; j < size; ++j) {
            double tmp;
            fscanf(f, "%lf", &tmp);
            grid[i][j] = tmp;
        }
    }
}

static inline double sec_deriv_x(int x, int y, double t) {
    if (x + dx >= size) {
        assert(x - dx >= 0);
        return (border(x + dx, y, t) - (2.0 * grid[x][y]) + grid[x - dx][y]) / (2.0 * dx * dx);
    } else if (x - dx <= 0) {
        assert(x + dx < size);
        return (grid[x + dx][y] - (2.0 * grid[x][y]) + border(x - dx, y, t)) / (2.0 * dx * dx);
    } else {
        return (grid[x + dx][y] - (2.0 * grid[x][y]) + grid[x - dx][y]) / (2.0 * dx * dx);
    }
}

static inline double sec_deriv_y(int x, int y, double t) {
    if (y + dy >= size) {
        assert(y - dy >= 0);
        return (border(x, y + dy, t) - (2.0 * grid[x][y]) + grid[x][y - dy]) / (2.0 * dy * dy);
    } else if (y - dy <= 0) {
        assert(y + dy < size);
        return (grid[x][y + dy] - (2.0 * grid[x][y]) + border(x, y - dy, t)) / (2.0 * dy * dy);
    } else {
        return (grid[x][y + dy] - (2.0 * grid[x][y]) + grid[x][y - dy]) / (2.0 * dy * dy);
    }
}

static inline void next_age() {
    int i;
    for (i = 0; i < size * size; ++i) {
        int x = i / size;
        int y = i % size;
        grid[x][y] += dt * (heat(x, y, cur_time) + 
                    (a_squared * (sec_deriv_y(x, y, cur_time) + 
                    sec_deriv_x(x, y, cur_time))));
    }
    cur_time += dt;
}

static inline void clean_dumps() {
    puts("> cleand_dumps");
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
    printf("> dump_to_file #%d/%d\n", cur_iteration, iterations);
    assert(ceil(log10(cur_iteration)) < 10);

    char cur_file[30];
    sprintf(cur_file, "dump_files/dump_%03d.txt", cur_iteration);

    FILE* f = fopen(cur_file, "w");
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

int main(int argc, char** argv) {
    int i;

    clean_dumps();
    init_from_file();
    for (i = 1; i < iterations + 1; ++i) {
        if (i % period == 0) {
            dump_to_file(i);
        }
        next_age();
        // printf(">%d\n", i);
    }
    return 0;
}