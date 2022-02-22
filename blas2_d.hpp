void dgemv ( char trans, int m, int n, double alpha, double a[], int lda, 
  double x[], int incx, double beta, double y[], int incy );
void dger ( int m, int n, double alpha, double x[], int incx, double y[], 
  int incy, double a[], int lda );
void dtrmv ( char uplo, char trans, char diag, int n, double a[], int lda, 
  double x[], int incx );
