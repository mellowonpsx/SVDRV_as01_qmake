// need for print working directory
#define _GNU_SOURCE
#include <unistd.h>

// include
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <pthread.h>

// N is dimension of matrix, images and W ([N]x[N])
#define N 512
// TH_NUM is the number of pthread created
// my laptop has only 2 core so:
// 1 thread vs 2 has big improvement (~100% faster)
// 2 thread vs 8 thread has no improvement
#define TH_NUM 8

// define the struct for complex number storage
typedef struct
{
    double real;
    double img;
} complex;

// define the struct for passing data to pthread called function, will be explained better later
typedef struct
{
    complex** result;
    complex** matrix;
    complex** inputMatrix;
    int dim;
    int start_row;
    int end_row;
    int start_cols;
    int end_cols;
} thread_data;

// sum 2 complex number and return the result
complex complex_sum (complex a, complex b)
{
    complex res;
    res.real = a.real + b.real;
    res.img = a.img + b.img;
    return res;
}

// subtract 2 complex number and return the result
complex complex_diff (complex a, complex b)
{
    complex res;
    res.real = a.real - b.real;
    res.img = a.img - b.img;
    return res;
}

// multiply 2 complex number and return the result
// (x+yi)*(u+vi) = (x*u)+(x*v*i)+(y*i*u)+(y*v*-1)
// (a.real+a.img)*(b.real+b.img) = (a.real*b.real)+(a.real*b.img*i)+(a.img*i*b.real)+(a.img*b.img*-1)
// (a.real+a.img)*(b.real+b.img) = (a.real*b.real)-(a.img*b.img)+((a.real*b.img)+(a.img*b.real))*i
// checked, seems correct
complex complex_mul (complex a, complex b)
{
    complex res;
    res.real = (a.real*b.real)-(a.img*b.img);
    res.img = (a.real*b.img)+(a.img*b.real);
    return res;
}

// divide 2 complex number and retur result
// checked, seems correct
complex complex_div (complex a, complex b)
{
    complex res;
    complex complexConiugate;
    complex quad;
    complexConiugate.real = b.real;
    complexConiugate.img = -b.img;
    //res = a*complexConiugate;
    res = complex_mul(a, complexConiugate);
    //quad = b*complexConiugate;
    quad = complex_mul(b, complexConiugate);
    // prevent div by zero
    if(quad.real != 0.0)
    {
        res.real = res.real/quad.real;
        res.img = res.img/quad.real;
        return res;
    }
    else
    {
        // we decide to exit in case of div by zero
        exit(-1);
        res.img = 0.0;
        res.real = 0.0;
        return res;
    }
}

// printa a complex number in cool style
// examples:
// does not write 0.0000+3i but 3i
// does not write 0.0000-3i but -3i
// does not write 3+0.00000i but 3
// does not write 3 2i but 3+2i
void printComplex(const complex c)
{
    if(c.real==0)
    {
        if(c.img!=0)
        {
            printf("%fi", c.img);
        }
        else
        {
            printf("%f", c.real);
        }
    }
    else
    {
        printf("%f", c.real);
        if(c.img!=0)
        {
            if(c.img>0)
            {
                printf("+%fi", c.img);
            }
            else
            {
                printf("%fi", c.img);
            }
        }
    }
}

// print an entire matrix made of complex number, tabulated
void printComplexMatrix(complex** matrix,int dim)
{
    int i,j;
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            printComplex(matrix[i][j]);
            printf("\t");
        }
        printf("\n");
    }
}

// print an entire matrix made of double number, tabulated
void printMatrix(double** matrix,int dim)
{
    int i,j;
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            printf("%g\t",matrix[i][j]);
        }
        printf("\n");
    }
}

// version used without thread
// will explaind the multithread version
int matrix_multiply_matrix(complex** result, complex** matrix, complex** inputMatrix, int dim)
{
    int i,j,k;
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            result[i][j].real = 0;
            result[i][j].img = 0;
            for(k=0; k<dim; k++)
            {
                //result[i][j] = result[i][j]+(matrix[i][k]*inputMatrix[k][j]);
                result[i][j] = complex_sum(result[i][j], complex_mul(matrix[i][k], inputMatrix[k][j]));
            }
        }
    }
    return 0;
}

// calculate the euclidian square distance of two matrix made of double element
double matrix_euclidian_square_distance_matrix(double** matrix1, double** matrix2, int dim)
{
    int i,j;
    double distance = 0.0;
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            double diff = matrix1[i][j]-matrix2[i][j];
            distance += pow(diff,2.0);
        }
    }
    return distance;
}

// multiply 2 matrix row by column, used with pthread
// pthread accept only function that have ONE (void*) input and return (void*)
void* matrix_multiply_matrix_thread(void *threaddata)
{
    // thread_data is a struct that contain function parameters
    thread_data *data;
    // explicit casting needed
    data = (thread_data *) threaddata;
    // result is a pointer to a complex matrix that will contain the result of R = AxB moltiplication (row by column)
    complex** result = data->result;
    // matrix is a pointer to the first matrix (A)
    complex** matrix = data->matrix;
    // matrix is a pointer to the second matrix (B)
    complex** inputMatrix = data->inputMatrix;
    // dim = matrix dimention, we assume squared matrix
    int dim = data->dim;
    // start_row and end_row indicate the first and the last row to multiplicate
    // we use this to limitate the moltiplication on a subset of row
    // and we assign each subset of row to a different pthread
    int start_row = data->start_row;
    int end_row = data->end_row;
    // same ad start_row and end_row but not used
    int start_cols = data->start_cols;
    int end_cols = data->end_cols;
    int i,j,k;
    // matrix multiplication row by column
    for(i=start_row; i<end_row; i++)
    {
        for(j=start_cols; j<end_cols; j++)
        {
            result[i][j].real = 0;
            result[i][j].img = 0;
            for(k=0; k<dim; k++)
            {
                //result[i][j] = result[i][j]+(matrix[i][k]*inputMatrix[k][j]);
                result[i][j] = complex_sum(result[i][j], complex_mul(matrix[i][k], inputMatrix[k][j]));
            }
        }
    }
    // nothing to return
    return NULL;
}

// perform matrix visualization ordering on complex element matrix
int dft_matrix_visual_ordering(complex **matrix, int dim)
{
    int i,j;
    complex temp;
    for(i=0; i<(dim/2); i++)
    {
        for(j=0; j<(dim/2); j++)
        {
            //we switch first quadrant element with fourth quadrant element
            temp = matrix[i][j];
            matrix[i][j] = matrix[i+(dim/2)][j+(dim/2)];
            matrix[i+(dim/2)][j+(dim/2)] = temp;

            //we switch second quadrant element with third quadrant element
            temp = matrix[i][j+(dim/2)];
            matrix[i][j+(dim/2)] = matrix[i+(dim/2)][j];
            matrix[i+(dim/2)][j] = temp;
        }
    }
    return 0;
}

// absolute (module) of a complex number
// |c| = (real(c)^2+img(c)^2)^(1/2)
double complex_abs(complex num)
{
    return sqrt(pow(num.real,2.0)+pow(num.img,2.0));
}

// calcolate absolute value of each cell of a complex matrix,
// put the result in a result matrix.
// will memorize and return the max value of the result matrix (used for compact code)
double abs_matrix(double **result, complex **matrix, int dim)
{
    int i,j;
    double max = 0.0;
    for(i=0; i<dim; i++)
    {
        for(j=0;j<dim; j++)
        {
            result[i][j] = complex_abs(matrix[i][j]);
            if(result[i][j]> max)
                max = result[i][j];
        }
    }
    return max;
}

// perform the scale operation on a matrix, put the result in a new matrix
int scale_matrix(double **result, double **matrix, double C, int dim)
{
    int i,j;
    for(i=0; i<dim; i++)
    {
        for(j=0;j<dim; j++)
        {
            result[i][j] = C*log(1+matrix[i][j]);
        }
    }
    return 0;
}

// copy the content of a file to a complex matrix,
// the file we use will contain only double value
// but we store in complex matrix for elaboration
int file_to_complex_matrix(FILE *file, complex **matrix, int dim)
{
    int i,j;
    if(file == NULL)
    {
        return 1;
    }

    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            complex a;
            a.img = 0.0;
            double value_read = 0.0;
            int read = fscanf(file,"%lf",&value_read);
            if(!read)
            {
                return 1;
            }
            a.real = value_read;
            matrix[i][j] = a;
            //must be the same of (cause W is symetric)
            //inputMatrix[j][i] = a;
        }
    }
    return 0;
}

// copy the content of a file to a matrix,
// the file we use will contain only double value
int file_to_matrix(FILE *file, double **matrix, int dim)
{
    int i,j;
    if(file == NULL)
    {
        return 1;
    }

    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            double value_read = 0.0;
            int read = fscanf(file,"%lf",&value_read);
            if(!read)
            {
                return 1;
            }
            matrix[i][j] = value_read;
        }
    }
    return 0;
}

// perform the multithread elaboration of dft on a complex inputMatrix and return a pointer to the result matrix
double** dft(complex** inputMatrix, int dim)
{
    int i,j;
    pthread_t thread[TH_NUM];
    complex **W, **DFT_2d, **DFT_1d;
    double **DFT_2d_abs, **DFT_2d_scaled;

    // Allocate memory
    W = (complex**) malloc(sizeof(complex*)*dim);
    DFT_2d = (complex**) malloc(sizeof(complex*)*dim);
    DFT_1d = (complex**) malloc(sizeof(complex*)*dim);
    DFT_2d_abs = (double**) malloc(sizeof(double*)*dim);
    DFT_2d_scaled = (double**) malloc(sizeof(double*)*dim);
    for (i = 0; i < N; i++)
    {
      W[i] = (complex*) malloc(sizeof(complex)*dim);
      DFT_2d[i] = (complex*) malloc(sizeof(complex)*dim);
      DFT_1d[i] = (complex*) malloc(sizeof(complex)*dim);
      DFT_2d_abs[i] = (double*) malloc(sizeof(double)*dim);
      DFT_2d_scaled[i] = (double*) malloc(sizeof(double)*dim);
    }

    // generate W[dim][dim]
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            // W is symmetric, [i][j] or [j][i] not does not change anything
            double w_i_j = -(2.0*M_PI/N)*(double)i*(double)j;
            W[i][j].real = cos(w_i_j);
            W[i][j].img = sin(w_i_j);
        }
    }

    // we must create thread_data for passing data to the function called by pthread.
    // we must create a thread_data for each pthread cause otherwise the thread_data content
    // can be modified before the function called by pthread access to thread_data
    thread_data passing_thread_data[TH_NUM];
    int rc;

    // we calculate dft in this way:
    // DFT_1d = W x inputMatrix;
    // DFT_2d = DFT_1d x W;
    // DFT_2d = ordering(DFT_2d);
    // DFT_2d_abs = |DFT_2d|;
    // DFT_2d_scaled = scale(DFT_2d_abs);


    // we must finish all first scan (computation of DFT_1d) before startin second scan (computation of DFT_2d)
    // for each pthread (in this case 8)
    for(i=0; i<TH_NUM; i++)
    {
        passing_thread_data[i].result = DFT_1d;
        passing_thread_data[i].matrix = W;
        passing_thread_data[i].inputMatrix = inputMatrix;
        passing_thread_data[i].dim = dim;
        // we perform multiplication on a limited set of row
        // in this case 512/8 = 64
        passing_thread_data[i].start_row = (dim*i)/TH_NUM;
        passing_thread_data[i].end_row = (dim*(i+1))/TH_NUM;
        // we perform multiplication on all comumns
        passing_thread_data[i].start_cols = 0;
        passing_thread_data[i].end_cols = dim;

        // we must create a passing_thread_data for each pthread cause
        // if the value change before the pthread use the data
        // that pthread may use a wrong data

        // calling pthread function, now that function runs on different pthread and this execution continues
        rc = pthread_create( &thread[i], NULL, matrix_multiply_matrix_thread, (void *) &passing_thread_data[i]);
        if (rc)
        {
            //printf("ERROR; return code from pthread_create() is %d\n", rc);
            return NULL;
        }
    }

    for(i=0; i<TH_NUM; i++)
    {
        // we must wait alla the matrix is completed before starting new scan
        rc = pthread_join(thread[i],NULL);
        if (rc)
        {
            return NULL;
        }
    }

    // now i have DFT_1d = W x inputMatrix;

    // repeat each operation for the second scan (computation of DFT_2d)
    for(i=0; i<TH_NUM; i++)
    {
        passing_thread_data[i].result = DFT_2d;
        // Y = X*W performs 1d_dft on each cols
        passing_thread_data[i].matrix = DFT_1d;
        passing_thread_data[i].inputMatrix = W;
        passing_thread_data[i].dim = dim;
        passing_thread_data[i].start_row = (dim*i)/TH_NUM;
        passing_thread_data[i].end_row = (dim*(i+1))/TH_NUM;
        passing_thread_data[i].start_cols = 0;
        passing_thread_data[i].end_cols = dim;

        rc = pthread_create( &thread[i], NULL, matrix_multiply_matrix_thread, (void *) &passing_thread_data[i]);
        if (rc)
        {
            return NULL;
        }
    }

    for(i=0; i<TH_NUM; i++)
    {
        pthread_join(thread[i],NULL);
        if (rc)
        {
            return NULL;
        }
    }

    // now i have DFT_2d = DFT_1d x W;

    // order the dft matrix
    if(dft_matrix_visual_ordering(DFT_2d,dim))
    {
        printf("Error during matrix_visula_ordering");
        return NULL;
    }

    double max;

    // perform DFT_2d_abs = |DFT_2d|;
    // remember abd_matrix return max value od DFT_2d_abs matrix
    if((max = abs_matrix(DFT_2d_abs,DFT_2d,dim))<0)
    {
        return NULL;
    }

    // calculate scale value
    double C = 255.0/log(max);

    // perform DFT_2d_scaled = scale(DFT_2d_abs);
    if(scale_matrix(DFT_2d_scaled,DFT_2d_abs,C,dim))
    {
        return NULL;
    }

    //free memory!!
    for(i=0; i<dim; i++)
    {
        free(W[i]);
        free(DFT_2d[i]);
        free(DFT_1d[i]);
        free(DFT_2d_abs[i]);
    }
    free(W);
    free(DFT_2d);
    free(DFT_1d);
    free(DFT_2d_abs);

    // return DFT_value
    return DFT_2d_scaled;
}


int main()
{
    int i;
    // the minimum distance is the max bound of double number
    double min_distance = DBL_MAX;

    complex **img_matrix;
    double** img_dft_scaled;
    double** file_dft_scaled;

    // memory allocation
    img_matrix = (complex**) malloc(sizeof(complex*)*N);
    file_dft_scaled = (double**) malloc(sizeof(double*)*N);
    for (i = 0; i < N; i++)
    {
        img_matrix[i] = (complex*) malloc(sizeof(complex)*N);
        file_dft_scaled[i] = (double*) malloc(sizeof(double)*N);
    }

    // load from file
    FILE *imgFile = NULL;

    // get current directory, may not work on certain system
    char *current_directory = (char *)get_current_dir_name();
    printf("Current directory is: %s\n", current_directory);
    printf("Please put file image file in \"../As_01_material_04-03-14/image.txt\"\n");
    printf("Please put all database file in \"../As_01_material_04-03-14/db/\"\n");
    printf("Each file will be readed and compared\n");

    // open image file to elaborate
    imgFile = fopen("../As_01_material_04-03-14/image.txt","r");

    if(imgFile == NULL)
    {
        printf("Error during opening image file \n");
        return 1;
    }

    // transfer image data to complex matrix
    if(file_to_complex_matrix(imgFile,img_matrix,N))
    {
        printf("Error during file to matrix operation \n");
        return 1;
    }
    fclose(imgFile);

    //calculate dft of img
    img_dft_scaled = dft(img_matrix, N);

    // open database directory
    DIR *db_dir;
    const char db_dir_str[]="../As_01_material_04-03-14/db/";

    db_dir = opendir(db_dir_str);

    if(db_dir == NULL)
    {
        printf("Error during opening db directory %s: %s\n",
               db_dir_str,
               strerror(errno));
        return 1;
    }

    // a defined struct to read directory
    struct dirent* in_file;
    FILE *file;
    char *nearest_filename = NULL;

    // while i have file in directory
    while ((in_file = readdir(db_dir)))
    {
        // not working with all filesystem
        // if the file is a file (and not a directory)
        if(in_file->d_type != DT_REG)
            continue;

        // construct complete filename
        char *complete_filename;
        complete_filename = (char*)malloc(sizeof(char)*(strlen(db_dir_str)+strlen(in_file->d_name)+1));
        strcpy(complete_filename, db_dir_str);
        strcat(complete_filename,in_file->d_name);

        // open the file
        file = fopen(complete_filename, "r");
        if (file == NULL)
        {
            printf("Error during opening file %s: %s\n", complete_filename, strerror(errno));
            return 1;
        }

        // copy file content to a double matrix
        if(file_to_matrix(file, file_dft_scaled,N))
        {
            printf("Error during file to matrix operation\n");
            return 1;
        }

        // euclidian distance from this file and the calculated dft of img
        double distance = matrix_euclidian_square_distance_matrix(img_dft_scaled, file_dft_scaled, N);

        printf("File \"%s\" has distance: %f\n",in_file->d_name,distance);

        // remember closest (file with less distances)
        if(distance < min_distance)
        {
            free(nearest_filename);
            nearest_filename = (char*)malloc(sizeof(char)*(strlen(in_file->d_name)+1));
            min_distance = distance;
            strcpy(nearest_filename, in_file->d_name);
        }

        // close file
        fclose(file);
        free(complete_filename);
    }
    if(nearest_filename)
    {
        printf("the nearest file is: %s (distance %g)\n",nearest_filename, min_distance);
    }
    // when we exit the main, the program and memory will be freed, so we can skip memory deallocation
    return 0;
}
// add this to make another commit and another push
