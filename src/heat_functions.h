/**
 * heat_functions.h
 * Heat Equation
 *
 *          2016. MIPT
 * Created by <anton.karazeev@gmail.com> on 21/10/16.
 *
 */

static inline double heat1(int x, int y, double t) {
    if ((x < 0.3*size && x > 0.1*size) || (x < 0.9*size && x > 0.7*size)) {
        return -1;
    } else {
        return 0;
    }
}

static inline double heat2(int x, int y, double t) {
    return 0;
}