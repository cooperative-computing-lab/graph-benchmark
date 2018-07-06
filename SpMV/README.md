Hybrid SpMV using MPI and OpenMP Benchmark

A version of the Message Passing Interface (MPI) is required
in order to use this software. To build the SpMV benchmark:

    make

To run the benchmark:

    ./spmv -lm <input_matrix.mtx> --omp-threads <INT> 
        --use-barriers <true|false> --distribution-method 
            <split|balanced>
            
--use-barriers is a flag that facilitates the usage of MPI_Barrier
in order to obtain greater program phase time measurements. 
--distribution-method determines the work load distribution
pattern. "split" uses a disjoint sub matrix of the input matrix
and assigns it to an MPI process, whereas "balanced" uniformly 
distributes non-zeros to each MPI process. 



To obtain the original benhmark matrix suite we have provided 
bash script:

    sh spmv_benchmark_suite_prep.sh
    
Running this script will automatically download the original 
25 matrices tested (approximately 25GB) and prepare them for
use. The prepared matrices can be found:

    cd matrices/
    
A tool is provided that allows the use to convert a symmetric
matrix in Matrix Market Format (MMF) to its expanded form. To
build the matrix expander:

    cd utilities
    make

To utilize the Matrix_Expander:
    ./Matrix_Expander -i <input_file> -o <output_file>
