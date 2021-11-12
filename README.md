# OSCA  

## Requirement

    Eigen and MKL(distribute within intel oneAPI)

## compile  
Samplely, just run following command, the __osca__ would be compiled.
```
make  
```

This would search Eigen and MKL head file under `/usr/include` and search library 
under `/usr/lib64`.  
If you don't install these in these directory, you may need edit Makefile or 
specific following variable when run `make` command.
```
    make EIGEN_PATH="eigen head file path" MKL_INCLUDE="mkl head file path" \
    MKL_LIB="mkl library file path"
```

By default only the dynamic version would compiled, if you want compile static 
version, run:
```
make osca_static
```

After compile, a file named __osca__ would appared under same directory of Makefile, and
if you compile static version, the execuable file, named __osca_static__, would show 
up too. 
