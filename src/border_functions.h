/**
 * border_functions.h
 * Heat Equation
 *
 *          2016. MIPT
 * Created by <anton.karazeev@gmail.com> on 21/10/16.
 *
 */

static inline double border1(int x, int y, double t) {
    if (y <= 0) {
        return 50.0 * exp(-pow(1.0 - (2.0 * (float)x/(float)size),2)/0.5);
    } else {
        return 0;
    }
}

static inline double border2(int x, int y, double t) {
    return 100;
}

static inline double border3(int x, int y, double t) {
    if (y <= 0) {
        return 3;
    } else {
        return 0;
    }
}