# OSCA  

OSCA (OmicS-data-based Complex trait Analysis) is a software tool for the analysis of complex traits using multi-omics data and genetic analysis of molecular phenotypes.

## Requirement

- Eigen

- MKL(Math Kernel Library)
    MKL is distributed within intel oneAPI or you can get it at https://github.com/oneapi-src/oneMKL.

- Rmath library

## Instatll requirements

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
Eigen just contain head files, and do not need compile to install it. Here, we use version 3.3.7 when we compile osca.

- MKL
    You can read document and download MKL at this webpage: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit.  

- Rmath library
    I using version 3.6 as example to install this library.

    * download R  
    `wget https://cloud.r-project.org/src/base/R-3/R-3.6.3.tar.gz`  
    `tar -zxf R-3.6.3.tar.gz`  
    `cd R-3.6.3`  
    
    * configure  
    `./configure --prefix="PATH where you what put this lib"`  
    
    * compile and install  
    `cd src/nmath/standalone`  
    `make`  
    `make install`  
    
    Then the head files and library files would show up under prefix path.

## Compile  

Samplely, if all requirement is install under **/usr** just run following command, the **osca** would be compiled.
```
make  
```
This would search head file under `/usr/include` and search library 
under `/usr/lib64`.  

If you don't install requirements in other directory, you may need edit Makefile or 
specific following variable when run `make` command.

```
make EIGEN_PATH="eigen head file path" MKL_INCLUDE="mkl head file path" \
MKL_LIB="mkl library file path" RMath_INCLUDE="path where install Rmath head files" RMath_LIB="path where install Rmath library"
```

By default only the dynamic version would compiled, if you want compile static 
version, run:

```
make VIRABLE="VALUE" osca_static
```
Using `DEBUG=ON` to switch on debug mode.

After compile, a file named **osca** would appared under same directory of Makefile, and
if you compile static version, the execuable file, named **osca_static**, would show 
up too. 

## USAGE

See https://yanglab.westlake.edu.cn/software/osca/#Overview for its usage and data resources.