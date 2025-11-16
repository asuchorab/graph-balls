# graph-balls 
Computing graph dataset metrics.

A ball refers to a neighbourhood of a node or an edge, containing everything up to a given radius, analogous to ball in geometry.

## Instructions

Required C++17 and CMake, tested with MSVC, uses only portable functions, I haven't tried on GNU for some time, but works on MinGW.

All libraries are compiled in, except tbbmalloc, a memory allocation optimisation library, which is disabled in cmake by default.

To see an example usage, run the application in `runtime` directory, using arguments
```
graphballs -d 1 -R -t -o output/[name] ./datasets/Countries.csv
```

To see more info on arguments, use `graphballs --help`. 
