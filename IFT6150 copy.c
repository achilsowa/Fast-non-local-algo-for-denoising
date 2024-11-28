/*------------------------------------------------------*/
/* Prog    : TpIFT6150-3-A.c                            */
/* Auteur  :                                            */
/* Date    :                                            */
/* version :                                            */
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define NAME_IMG_IN "photograph"

#define NAME_IMG_OUT1 "photograph_bruite_B"
#define NAME_IMG_OUT2 "photograph_debruite_B"

#define TRUE 1
#define FALSE 0

void SaveAndShow(char *name, float **mat, int lgth, int wdth)
{
  char buff[100 + NBCHAR];
  SaveImagePgm(name, mat, lgth, wdth);

  sprintf(buff, "display %s.pgm&", name);
  // printf("cmd: %s\n", buff);
  system(buff);
}

/**
 * Calcule l'ISNR
 */
float ISNR(float **f, float **g, float **f_hat, int lgth, int wdth)
{
  int i, j;
  float f_g = 0.0, f_f_hat = 0.0;
  for (i = 0; i < lgth; ++i)
    for (j = 0; j < lgth; ++j)
    {
      f_g += SQUARE(f[i][j] - g[i][j]);
      f_f_hat += SQUARE(f[i][j] - f_hat[i][j]);
    }
  return 10 * log10f(f_g / f_f_hat);
}

/**
 * Initialise l'image dst a zero
 */
void CopyZero(float **dst, int lgth, int wdth)
{
  int i, j;
  for (i = 0; i < lgth; ++i)
    for (j = 0; j < wdth; ++j)
    {
      dst[i][j] = 0.0;
    }
}

/**
 * Nettoie l'image en mettant a zero les pixels negatifs et a 255 les pixels > 255
*/
void CleanImage(float **mat, int lgth, int wdth) {
  int i, j;
  for (i = 0; i < lgth; i++)
    for (j = 0; j < wdth; j++){
      if (mat[i][j] < 0) mat[i][j] = 0;
      if (mat[i][j] > 255) mat[i][j] = 255;
    }
}

/**
 * Valeur du seuillage
*/
float ThresholdedValue(float val, float threshold) {
  float fval = fabsf(val);
  if (fval < threshold) {
    return 0;
  }else {
    return 255 * (val/fval) * (fval-threshold)/(255-threshold);
  }
}

/**
 * Effectue une filtrage dans le domaine des ondelettes
 */
void DenoisingHaar(float **f, float **g, float **image, int lgth, int wdth, float threshold, int nbLevels){
  float **haar = fmatrix_allocate_2d(lgth, wdth);
  
  CopyZero(haar, lgth, wdth);
  haar2D_complete(g, haar, nbLevels, lgth, wdth);

  int i, j, mlgth = lgth/pow(2, nbLevels),  mwdth = wdth/pow(2, nbLevels);

  for(i=0; i<lgth; ++i)
    for(j=0; j<wdth; ++j) 
      if (i>= mlgth || j>= mwdth) {
        haar[i][j] = ThresholdedValue(haar[i][j], threshold);
      }
  
  ihaar2D_complete(haar, f, nbLevels, lgth, wdth);
  
  CleanImage(f, lgth, wdth);

  free_fmatrix_2d(haar);
}


int main(int argc, char **argv)
{
  int i, j, k, l;
  int length, width;
  int nbLevels;

  float threshold;
  float var;

  float **image; /* image d'entree */
  float **g;     /* image degradee */
  float **f;     /* image restoree */

  printf("Entrez la variance du bruit: ");
  scanf("%f", &var);
  // var = 900;

  printf("Entrez le nombre de niveaux a traiter : ");
  scanf("%d", &nbLevels);
  // nbLevels = 4;

  printf("Entrez le seuil : ");
  scanf("%f", &threshold);
  // threshold = 47;

  /* ouvrir l'image d'entree */
  image = LoadImagePgm(NAME_IMG_IN, &length, &width);

  /* ajouter du bruit a l'image d'entree (add_gaussian_noise) */
  g = fmatrix_allocate_2d(length, width);
  copy(g, image, length, width);
  add_gaussian_noise(g, length, width, var);

  /* debruiter l'image en seuillant les coefficients de Haar */
  f = fmatrix_allocate_2d(length, width);
  DenoisingHaar(f, g, image, length, width, threshold, nbLevels);
  printf("ISNR:\t %f\n", ISNR(image, g, f, length, width));
  /* afficher l'ISBN */

  /* sauvegarder les images */
  SaveAndShow(NAME_IMG_OUT1, g, length, width);

  SaveAndShow(NAME_IMG_OUT2, f, length, width);

  /* liberer la memoire */
  free_fmatrix_2d(image);
  free_fmatrix_2d(g);
  free_fmatrix_2d(f);


  /*retour sans probleme*/
  printf("\n C'est fini ... \n\n\n");
  return 0;
}
