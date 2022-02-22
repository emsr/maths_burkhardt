void dgemm (char transa, char transb, int m, int n, int k, 
            double alpha, double a[], int lda, double b[], int ldb, double beta, 
            double c[], int ldc );
void dtrmm (char side, char uplo, char transa, char diag, int m, int n, 
            double alpha, double a[], int lda, double b[], int ldb );
void dtrsm (char side, char uplo, char transa, char diag, int m, int n, 
            double alpha, double a[], int lda, double b[], int ldb );

