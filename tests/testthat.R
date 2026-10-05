# Limit BLAS/OpenMP thread usage before loading testthat/irtQ, so that the
# environment variables take effect before the underlying BLAS/OpenMP
# libraries are loaded into the session.
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OMP_THREAD_LIMIT = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")

library(testthat)
library(irtQ)

# Directly call the BLAS/OpenMP C-level thread control APIs at runtime, so the
# thread count is forced to 1 regardless of when the libraries were loaded
# (the Sys.setenv() calls above only take effect at load time).
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::omp_set_num_threads(1)   # force one OpenMP thread at runtime
  RhpcBLASctl::blas_set_num_threads(1)  # force one BLAS thread at runtime
}

test_check("irtQ")
