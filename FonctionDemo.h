/*---------------------------------------------------*/
/* module  : FonctionDemo3.h                         */
/* auteur  : Max Mignotte                            */
/* revision: Francois Destrempes                     */
/* date    : 21/09/99--08/10/04                      */              
/* langage : C                                       */
/* labo    : DIRO                                    */
/*---------------------------------------------------*/

#ifndef FONCTIONDEMO_H
#define FONCTIONDEMO_H

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define SQUARE(X) ((X)*(X))
#define MAX(i,j)  ((i)>(j)?(i):(j))
#define MIN(i,j)  ((i)<(j)?(i):(j))

#define NBCHAR 200

#define FFT   1
#define IFFT -1
#define FFT2D 2

#define GREY_LEVEL 255
#define PI 3.141592654

#define WHITE 255
#define BLACK 0

#define TRUE 1
#define FALSE 0


typedef struct {
    int row;
    int col;
    float **data;  
} mat_t;

typedef struct {
    int size;
    float *data;  
} vect_t;



/*------------------------------------------------*/
/* PROTOTYPES DES FONCTIONS  ---------------------*/
/*------------------------------------------------*/
float*  fmatrix_allocate_1d(int);
float** fmatrix_allocate_2d(int,int);
void    free_fmatrix_1d(float*);
void    free_fmatrix_2d(float**);
float** LoadImagePgm(char*,int*,int*);
void    SaveImagePgm(char*,float**,int,int);
void    fourn(float*,unsigned long*,int,int);
void    FFTDD(float**,float**,int,int);
void    IFFTDD(float**,float**,int,int);
void    Mod(float**,float**,float**,int,int);
void    Mult(float**,float,int,int);
void    Recal(float**,int,int);
void    MultMatrix(float**,float**,float**,float**,float**,float**,int,int);
void    SquareMatrix(float**,float**,float**,float**,int,int);
float   funcgauss2D(float x,float y,float var);
void    compute_histo(float** mat,int lgth,int wdth,float* hist);
float   gaussian_noise(float var,float mean);
void    add_gaussian_noise(float** mat,int lgth,int wdth,float var);
void    add(float** matr, float** mat1, float** mat2,int lgth, int wdth);
void    substract(float** matr, float** mat1, float** mat2,int lgth, int wdth);
void    copy(float** mat_in, float** mat_out,int lgth, int wdth);
void    Recal_haar_step(float** image,int x0,int y0,int lgth,int wdth);
void    Recal_haar(float** image,int M,float** mat_tmp,int lgth,int wdth);

void haar1D(float* signal,float* work,int lgth);
void haar2D(float** image,int lgth,int wdth);
void haar2D_complete(float** image,float** haar,int M,int lgth,int wdth);
void ihaar1D(float* signal,float* work,int lgth);
void ihaar2D(float** image,int lgth,int wdth);
void ihaar2D_complete(float** haar,float** haar_inverse,int M,
		      int lgth,int wdth);


void save_and_show(char *name, float **mat, int lgth, int wdth);
float isnr(float **f, float **g, float **f_hat, int lgth, int wdth);
void copy_val(float **dst, int lgth, int wdth, float val);
void copy_zero(float **dst, int lgth, int wdth);
void mult(float **dst, int lgth, int wdth, float val);
void clean_image(float **mat, int lgth, int wdth);
void print_image(float **mat, int lgth, int wdth);
float thresholded_value(float val, float threshold);
void denoising_haar(float **f, float **g, float **image, int lgth, int wdth, float threshold, int nbLevels);
mat_t* allocate_mat(int row, int col);
vect_t* allocate_vect(int size);
void free_mat(mat_t* mat);
void free_vect(vect_t* vect);
mat_t* load_image(char* name);
int is_valid_index_mat(mat_t* mat, int i, int j);


#endif
