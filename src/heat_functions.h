//
//  heat_functions.h
//  HeatEquation
//
//  Created by <anton.karazeev@gmail.com> on 21/10/16. MIPT.
//  Copyright © 2016 Anton Karazeev. All rights reserved.
//

static inline double heat1(int x, int y, double t) {
    if (((x < 0.3*size_x && x > 0.1*size_x) || (x < 0.9*size_x && x > 0.7*size_x)) && (y <= 0.7*size_y && y >= 0.3*size_y)) {
        return -1;
    } else {
        return 0;
    }
}

static inline double heat2(int x, int y, double t) {
    return 0;
}

static inline double heat3(int x, int y, double t) {
    if (pow(x,2) + pow(y,2) <= pow(0.3*size_x,2)) {
    	return -0.5;
    } else {
    	return 0;
    }
}

static inline double heat4(int x, int y, double t) {
    if (pow(x-(0.5*size_x),2) + pow(y-(0.5*size_y),2) <= pow(0.3*size_x,2)) {
    	return -1;
    } else {
    	return 0;
    }
}
