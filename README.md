# MIPT. Modeling of Heat equation

![Demo](https://cloud.githubusercontent.com/assets/11920213/19579585/2648d77e-972a-11e6-90d7-b91923fd859d.gif)

***

#TODO:
 - _heat_ and _border_ functions choose by name
 - use iniparser to parse _config_ file
 - use PATH_MAX to limit _prefix_ size
 - try to READ NetCDF

***

There will be instructions soon...

***

### Structure of config file:
1. dx, dy
2. dt
3. total_iterations
4. period
5. file with initial matrix OR size of matrix (one number)
6. <prefix> of dump files: "<prefix>_iteration.txt"
7. number of _border_ function
8. number of _heat_ function
9. TODO: number of threads

### Border functions:
1. Left side heating
2. Unifrom heating around the perimeter
3. Uniform heating of the left side
4. Just zero
5. Continuous ...

### Heat functions:
1. Two plates
2. Just zero
3. Negative semicircle
4. Negative circle at the center