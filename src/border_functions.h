//
//  border_functions.h
//  HeatEquation
//
//  Created by <anton.karazeev@gmail.com> on 21/10/16. MIPT.
//  Copyright © 2016 Anton Karazeev. All rights reserved.
//

static inline double border1(int x, int y, double t) {
    if (y <= 0) {
        return 50.0 * exp(-pow(1.0 - (2.0 * (float)x/(float)size_x),2)/0.5);
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

static inline double border4(int x, int y, double t) {
    return 0;
}

static inline double border5(int x, int y, double t) {
    if (x >= 0 && x < size_x) {
        if (y > size_y) {
            return grid[x][size_x-1];
        } else if (y < 0) {
            return grid[x][0];
        } else {
            return grid[x][y];
        }
    } else if (x < 0) {
        assert(y >= 0 && y < size_y);
        return grid[0][y];
    } else if (x >= size_x) {
        assert(y >= 0 && y < size_y);
        return grid[size_x-1][y];
    } else {
        assert(0);
    }
}
