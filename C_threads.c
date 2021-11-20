#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"
#include <time.h>
#include <pthread.h>

/*code for timing is from 
stackoverflow.com/questions/13950290/clock-gettime-nanoseconds-calculation*/
struct timespec diff(struct timespec start, struct timespec end){
    struct timespec temp;
    if ((end.tv_nsec - start.tv_nsec) < 0) {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } 
    else{
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}

/*matrix in coo format is stores as array of coo structs
Useful for qsort() */
typedef struct coo{
    int i;
    int j;
    int value;
}coo;

/*Data to be passed to each thread*/
typedef struct Data{
    int *row, *col, *val;
    int start_row, end_row;
}Data;


/*compare to be used with qsort*/
int compareCoo(const void *a, const void *b);

/*Given the .mtx file this function extracts coordinates of non zeroes*/
int* getMatrixData(char* mat, int** I, int** J);

void cooToCsr(int dim,int nz,coo ** pairs,int ** csr_row,int ** csr_col );

void partitionArr(int size, int parts, int** indeces);

int roundN(float num){
    return (int)(num < 0 ? (num - 0.5) : (num + 0.5));
}

void* worker(void* allD);


int main(int argc, char *argv[]){

    int i, *data, *I, *J;
    data = getMatrixData(argv[1], &I, &J);
    //int max_threads = atoi(argv[2]);

    int dim = data[0];
    int nz = data[1];

    /* non zero elements are only given once in .mtx file */
    printf("\nNumber of rows is %d and number of non zeroes is: %d \n", dim, 2*nz);

    /*After getting data create coo form of array. Value is 1 and will be ignored for now.*/
    coo* pairs ;
    pairs = malloc(2*nz*sizeof(coo));

    /*Since matrix is symmetrical we need to account for twice as many elements*/
    for (int k=0; k<nz; k++){
        pairs[k].i = J[k];
        pairs[k].j = I[k];
        pairs[2*nz-k-1].i = I[k];
        pairs[2*nz-k-1].j = J[k];
    }

    /*Always free malloced memory*/
    free(I);
    free(J);

    /*using a custom comparator function coo is (i, j) sorted*/
    qsort((void*)pairs, 2*nz, sizeof(pairs[0]), compareCoo);

    printf("\n-------------------------------\n");
    
    /*COO to CSR conversion*/
    int *A_row, *A_col;
    A_col = malloc(2*nz*sizeof(int));
    A_row = malloc((dim+1)*sizeof(int));

    for(i=0; i<dim+1; i++){
        A_row[i] = 0;
    }

    cooToCsr(dim,nz,&pairs,&A_row,&A_col);

    free(pairs);


    /*-------------------------------------------*/
    char answer;
    printf("Do you want to print CSR? [Y/N]: ");
    

    if(scanf("%c", &answer) && answer == 'Y'){
        printf("Col_index is: ");
        for(i=0;i<2*nz;i++){
            printf(" %d ", A_col[i]);
        }
        printf("\nRow_index is: ");
        for(i=0;i<dim+1;i++)
            printf(" %d ", A_row[i]);
    }

    printf("\n-------------------------------\n");

    FILE *fp;
    fp = fopen("threadS.txt", "w");
    struct timespec tStart,tEnd;
    double minTime = 100;
    
    /*Outside for is for benchmarking*/
    for(int max_threads = 2; max_threads < 17; max_threads++){
        double mean_time = 0;
        for(int measurements = 0 ; measurements<10 ; measurements++){
            int* C_val;
            C_val = malloc(2*nz*sizeof(int));
            for (i = 0; i < 2*nz; ++i)  C_val[i] = 0 ;

            int dif,starti, startj, endi, endj ;

            /*----------------------------------------------*/
            int* partition; 
            partitionArr(dim, max_threads, &partition);
            clock_gettime(CLOCK_MONOTONIC, &tStart);

            /*Divide array in aprox. equal parts for better use of memory */
            
            pthread_t tid[max_threads];

            /*Fine point:
                at first glance we could use the same struct
                with different values for start and end but 
                there is no guarantee the thread will be created 
                before the for loop gets to the next i changing
                start and end.
            */
            Data allData[max_threads];

            for(i=0;i<max_threads;i++){
                allData[i].start_row = partition[i];
                allData[i].end_row = partition[i+1];
                allData[i].row = A_row;
                allData[i].col = A_col;
                allData[i].val = C_val;
                
                if( pthread_create(&tid[i], NULL, worker, (void*)&allData[i]) ){
                    printf("Error when creating thread: %d \n", i);
                    exit(1);
                }
            }

            /*Wait for all threads to finish*/
            for(i=0;i<max_threads;i++) pthread_join(tid[i], NULL);


            clock_gettime(CLOCK_MONOTONIC, &tEnd);
            /*--------------------------------------------------*/

   
            struct timespec tResult = diff(tStart,tEnd);
            double timeR = (double)(tResult.tv_sec+(double)tResult.tv_nsec/1000000000);

            mean_time += timeR;
            //fprintf(fp, "%f\n", timeR);
            //printf( "Benchmark %d done\n", measurements);
        }
        fprintf(fp, "%f\n", mean_time/10);
        printf( "Benchmark %d done\n", max_threads);
    }
    
    fclose(fp);
    return 0;
}



int* getMatrixData(char* mat, int** I, int** J){
    FILE *f;
    int ret_code, M, N, nz, i;   
    int* ret = malloc(2*sizeof(int));


    if ((f = fopen(mat, "r")) == NULL) {
            printf("Unexpected argument\n");
            exit(1);
    }
    
    /* find out size of sparse matrix */
    ret_code = mm_read_mtx_crd_size(f, &ret[0], &N, &ret[1]);

    /* reseve memory for matrices */
    *I = malloc(ret[1] * sizeof(int));
    *J = malloc(ret[1] * sizeof(int));

    for (i=0; i<ret[1]; i++){
        fscanf(f, "%d %d \n", &(*I)[i], &(*J)[i]);
        (*I)[i]--;  /* adjust from 1-based to 0-based */
        (*J)[i]--;    
    }
    return ret;
}

int compareCoo(const void *c1, const void *c2){
    const struct coo *a = c1;
    const struct coo *b = c2;

    if(a->i < b->i)
        return -1;
    else if(a->i > b->i)
        return 1;
    else if(a->i == b->i){
        if(a->j < b->j)
            return -1;
        else if(a->j > b->j)
            return 1;
        else 
            return 0;
    } 
    return 2;
}

void cooToCsr(int dim,int nz,coo ** pairs,int ** csr_row,int ** csr_col ){
    int temp = 0;
    (*csr_row)[dim] = 2*nz;
    for(int k=0; k<2*nz; k++){
        (*csr_col)[k] = (*pairs)[k].j;

        if((*pairs)[k].i != temp){
            temp+=(*pairs)[k].i - temp;
            (*csr_row)[temp] = k;
        }
    }

    for(int k = 0; k<dim+1 ; k++){
        if((*csr_row)[k] == 0) (*csr_row)[k] = (*csr_row)[k-1];
    }
}

/*partition the row collumn so that each thread has aprox. the same amount of work*/
void partitionArr(int size, int parts, int** indeces){
    float temp;
    int forparts = parts;
    *indeces = malloc( (parts+1) * sizeof(int));
    (*indeces)[0] = 0;
    (*indeces)[parts] = size;

    for(int i = 1; i< forparts+1; i++){
        
        temp = roundN( (float)size/(float)parts );
        size-=temp;
        parts-=1;

        (*indeces)[i] = temp+(*indeces)[i-1];
    }
}


void* worker(void* allD){
    Data *allData = (Data*) allD; 
    int dif,starti, startj, endi, endj ;
    int* A_row = allData->row;
    int* A_col = allData->col;
    int* C_val = allData->val;

    /*very similar to linear programm except rows are divided amongst the threads*/
    for(int row = allData->start_row ; row < allData->end_row; row++){
        for(int col = A_row[row] ; col < A_row[row+1]; col++){

            /*current index to compare*/
            starti = A_row[row];
            startj = A_row[ A_col[col] ];
            
            /*end of rows*/
            endi = A_row[row+1];
            endj = A_row[ A_col[col]+1 ];

            /*find common elements*/
            while(starti < endi && startj < endj){
                dif = A_col[starti] - A_col[startj]; 

                /*start with comparisons most likely to occur*/
                if(dif<0){   
                     starti++;  
                }
                else if(dif>0){
                    startj++;
                }
                else{
                    C_val[col]++;
                    starti++;
                    startj++;
                }     
            }
        }
    }
    return NULL;
}
