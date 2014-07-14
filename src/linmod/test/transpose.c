#include <stdio.h>
#include <stdlib.h>


// row-major printing
void show_matrix(int m, int n, double *x)
{
  int i, j;
  
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
      printf("%2g ", x[i*n + j]);
    
    putchar('\n');
  }
}



void rect_test()
{
  int i;
  double x[15];
  for (i = 0; i < 15; i++) x[i] = i + 1;
  
  puts("Row major:");
  show_matrix(5, 3, x);
  
  c2r(5, 3, x);
  
  puts("\nCol major:");
  show_matrix(3, 5, x);
  
  r2c(3, 5, x);
  
  puts("\nRow major:");
  show_matrix(5, 3, x);
}



void square_test()
{
  int i;
  double x[16];
  for (i = 0; i < 16; i++) x[i] = i + 1;
  
  puts("Row major:");
  show_matrix(4, 4, x);
  
  c2r(4, 4, x);
  
  puts("\nCol major:");
  show_matrix(4, 4, x);
}



int main(void)
{
  rect_test();
  
/*  square_test();*/
  
  return 0;
}

