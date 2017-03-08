/*created by Ahsen Dinc ad39724
SCC374C/394C: Parallel Computing for Science and Engineering
Assignment #0 Serial (CPU) code*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#define array(i,j) array[(i) + (j)*n]
#define x(i,j) x[(i) + (j)*n]
#define y(i,j) y[(i) + (j)*n]
double gtod_timer(void);

void initialize(float *array,  int n){
    srand(8);
    int i, j;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            array(i, j) = (rand()/(float)RAND_MAX);
        }
    }
}
//yani y (n-2)*(n-2)lik mi?
void smooth(float *x, float *y, int n, float a, float b, float c) {
    int i, j;
#pragma omp for nowait schedule(static) private(i,j)
   for (i = 1; i < n - 1; i++) {
        for (j = 1; j < n - 1; j++) {
            y(i, j) = (a * (x(i - 1, j - 1) + x(i - 1, j + 1) + x(i + 1, j - 1) + x(i + 1, j + 1))) +
                      (b * (x(i - 1, j + 0) + x(i + 1, j + 0) + x(i + 0, j - 1) + x(i + 0, j + 1))) +
                      (c * x(i + 0, j + 0));

        }
    }
}
void count (float *array, int n, float t, int *num){
    *num=0;
    int tmp_num=0;
    int i,j;
#pragma omp parallel num_threads(16)
{
printf("%d\n", omp_get_thread_num());
#pragma omp for schedule(dynamic) private(i,j) reduction(+:tmp_num)
    for(i=1; i<n-1; i++){
        for(j=1; j<n-1; j++){
            if (array(i,j) < t){
                //*num = *num + 1;
                tmp_num = tmp_num+1;
            }
        }
    }
	*num = tmp_num;
}
}
int main() {
    const float a=0.05, b=0.1 , c=0.4;
    float t=0.1;
    float *x, *y;
    int num_x, num_y;
    int n=32768+2;
    //int n=100;//to make the test faster
    float start, stop, t_allocx, t_allocy, t_initx;
    double  t_smooth, t_countx, t_county ;
    double time, t0, t1;
    //struct timespec tp;
    //clockid_t clk_id;
    //clk_id = CLOCK_REALTIME;

    start=(float)clock()/CLOCKS_PER_SEC;
    x= (float *)malloc(sizeof(float)*n*n);
    stop = (float)clock()/CLOCKS_PER_SEC;
    t_allocx =stop-start;

    start=(float)clock()/CLOCKS_PER_SEC;
    y= (float *)malloc(sizeof(float)*n*n);
    stop = (float)clock()/CLOCKS_PER_SEC;
    t_allocy =stop-start;

    start=(float)clock()/CLOCKS_PER_SEC;
    initialize(x,n);
    stop = (float)clock()/CLOCKS_PER_SEC;
    t_initx =stop-start;
	
    int nt;
    
#ifdef _OPENMP
#pragma omp parallel private(nt)
{ nt = omp_get_num_threads(); if(nt<1) printf("NO print, OMP warmup.\n"); }
#endif

     double start_time = omp_get_wtime();
#pragma omp parallel num_threads(16)
{
	printf("%d\n", omp_get_thread_num());
	smooth(x,y,n,a,b,c);    
}
    t_smooth = omp_get_wtime() - start_time;

    //start=(float)clock()/CLOCKS_PER_SEC;
    //#pragma omp parallel num_threads(16)
   //{
//    t0 = gtod_timer();
    start_time = omp_get_wtime();    
    count(x, n, t, &num_x);
   t_countx = omp_get_wtime() - start_time;
    //stop = (float)clock()/CLOCKS_PER_SEC;
    //t_countx =stop-start;
//	t_countx = t1-t0;
    //start=(float)clock()/CLOCKS_PER_SEC;
    start_time = omp_get_wtime();
    count(y, n, t, &num_y);
    //stop = (float)clock()/CLOCKS_PER_SEC;
    //t_county =stop-start;
    t_county = omp_get_wtime() - start_time;

    FILE* fp;
    fp = fopen("hw1_output", "w+");
    fprintf(fp, "created by Ahsen Dinc ad38724\nSCC374C/394C: Parallel Computing for Science and Engineering\nAssignment #0 Serial (CPU) code\n\n");
    fprintf(fp, "Summary\n");
    fprintf(fp, "-------\n\n");
    fprintf(fp, "Number of elements in a row/column       ::  %d\n\n", n);
    fprintf(fp, "Number of inner elements in a row/column ::  %d\n\n", n-2);
    fprintf(fp, "Total number of elements                 ::  %d\n\n", n*n);
    fprintf(fp, "Total number of inner elements           ::  %d\n\n", (n-2)*(n-2));
    fprintf(fp, "Memory (GB) used per array               ::  %.5f\n\n", 4*(n*n)/(float)(1024*1024*1024));
    fprintf(fp, "Threshold                                ::  %.2f\n\n", t);
    fprintf(fp, "Smoothing constants (a, b, c)            ::  %.2f %.2f %.2f\n\n", a, b, c);
    fprintf(fp, "Number of elements below threshold (X)   ::  %d\n\n", num_x);
    fprintf(fp, "Fraction of elements below threshold     ::  %9.5e \n\n", num_x/(float)((n-2)*(n-2)));
    fprintf(fp, "Number of elements below threshold (Y)   ::  %d\n\n", num_y);
    fprintf(fp, "Fraction of elements below threshold     ::  %9.5e\n\n\n", num_y/(float)((n-2)*(n-2)));

    //fprintf(fp, "Action       ::     time/s     Time resolution = %9.2e\n",clock_getres( clk_id, &tp));
    fprintf(fp, "Action       ::     time/s     Time resolution = %9.4e\n",1/(float)CLOCKS_PER_SEC);
 
    fprintf(fp, "------\n\n");
    fprintf(fp, "CPU: Alloc-X ::     %.4f\n\n", t_allocx);
    fprintf(fp, "CPU: Alloc-Y ::     %.4f\n\n", t_allocy);
    fprintf(fp, "CPU: Init-X  ::     %.4f\n\n", t_initx);
    fprintf(fp, "CPU: Smooth  ::     %.4f\n\n", t_smooth);
    fprintf(fp, "CPU: Count-X ::     %.4f\n\n", t_countx);
    fprintf(fp, "CPU: Count-Y ::     %.4f\n\n", t_county);

    fclose(fp);

    return 0;
}
