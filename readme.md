# graph-balls 
Computing graph dataset metrics.

A ball refers to a neighbourhood of a node or an edge, containing everything up to a given radius, analogous to ball in geometry.

## Instructions

Required C++17 and CMake, tested with MSVC, uses only portable functions, I haven't tried on GNU for some time, but works on MinGW.

All libraries are compiled in, except tbbmalloc, a memory allocation optimisation library, which is disabled in cmake by default.

To see an example usage, run the application in `runtime` directory, using arguments
```
graphballs -d 4 -R -t -o output/[name] ./datasets/Countries.csv
```

To see more info on arguments, use `graphballs --help`. 

The feature number 1 has been done, distinguishability of edges can be checked with with -e flag:
```
graphballs -e -d 4 -o output/[name] ./datasets/Countries.csv
```

The features number 2 and 3 are partially done, handling of multiedges and edge labels is controlled with -R and -l. Not working for edge-based tasks yet.
```
graphballs -d 4 -R -t -o output/[name] ./datasets/AristoV4.csv
graphballs -d 4 -t -o output/[name] ./datasets/AristoV4.csv
graphballs -d 4 -l -t -o output/[name] ./datasets/AristoV4.csv
```
