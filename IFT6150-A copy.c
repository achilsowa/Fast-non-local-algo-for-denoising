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
#define NAME_IMG_IN "lenna"

#define NAME_IMG_NOISED "lenna_bruite_A"
#define NAME_IMG_DENOISED "lenna_debruite_A"


#define NEIGHBORHOOD_SIZE 2 //7
#define SEARCHING_SIZE 3 //21
#define MIN_SIMILARITY -1e6

void softmax(const float *input, float *output, int size, int hide) {
  float max_val = input[0];
  for (int i = 1; i < size; ++i) {
      if (input[i] > max_val) {
          max_val = input[i];
      }
  }

  float sum_exp = 0.0; 
  for (int i = 0; i < size; ++i) {
    sum_exp += exp(input[i] - max_val);
  }

  for (int i = 0; i < size; ++i) {
    output[i] = exp(input[i] - max_val) / sum_exp;
  }

  if (hide % 500 == 0) {
  printf("input: ");
  for(int i=0; i<size; ++i) printf("%.1f ", input[i]);
  printf("\noutput: ");
  for(int i=0; i<size; ++i) printf("%.2f ", output[i]);
  sum_exp = 0;
  for (int i = 0; i < size; ++i) sum_exp += output[i];
  printf("\t %f\n\n", sum_exp);
  }
}

/**
 * save the neihborhood of size (M, M), with 0 where the neihborhood might not be defined (padding)
*/
void neighborhood(mat_t* img, mat_t* n, int *pixel) {
  int M = n->row, i = pixel[0], j= pixel[1];
  int start_i = MAX(0, i-M/2), start_j = MAX(0, j-M/2);
  int end_i = MIN(img->row, i+M/2), end_j = MIN(img->col, j+M/2);
  
  copy_zero(n->data, n->row, n->col);
  
  for(int x=start_i; x<end_i; ++x)
    for(int y=start_j; y<end_j; ++y) 
      n->data[x-start_i][y-start_j] = img->data[x][y];
}

float neighborhood_distance(mat_t* n_i, mat_t* n_j) {
  float ss = 0;
  int M = n_i->row;
  for(int l=0; l<M; ++l)
    for(int k=0; k<M; ++k)
      ss += SQUARE(n_i->data[l][k] - n_j->data[l][k]);
  return ss;
}

/**
 * calcule la matrice -S(i,j)
*/
void similarity(mat_t* img, mat_t* s, int neighbor_size) {
  copy_val(s->data, s->row, s->col, MIN_SIMILARITY);
  mat_t* n_i = allocate_mat(neighbor_size, neighbor_size), *n_j = allocate_mat(neighbor_size, neighbor_size);
  int similarity_window_size = (int)sqrtf(s->col);
  int pos_i[2], pos_j[2], xi, yi, xj, yj;

  
  /**
   * Pour chaque pixel i, on evalue son voisinage n_i, et les n_j correspondant 0 la fenetre de similarite (soit 21x21) dans le papier.
   * Pour les pixels i de position < 21x21, comme (0, 0) on ne prend en compte que les pixels voisins valides, i.e. ceux de droite.
   * Pareillement pour les pixels du bords droit de l'image, on ne prend en compte que les pixels voisins valides, i.e. ceux de gauche.
  */
  for(int i=0; i<s->row; ++i) {
    xi = (int)i/img->row, yi = (int)i%img->col;
    pos_i[0] = xi, pos_i[1] = yi;
    neighborhood(img, n_i, pos_i);
    
    for(int j=0; j<s->col; ++j) {
      xj = xi + (int)j/similarity_window_size - similarity_window_size/2;
      yj = yi + (int)j%similarity_window_size - similarity_window_size/2;

      // Si le point (xj, yj) est dans les limites de l'image, on me a jour S(i, j), si non, on laisse la valeur par defaut
      if (xj < 0 || xj >= img->row || yj < 0 || yj >= img->col) {
        s->data[i][j] = MIN_SIMILARITY;
      } else {
        pos_j[0] = xj, pos_j[1] = yj;
        neighborhood(img, n_j, pos_j);
        s->data[i][j] = MAX(-neighborhood_distance(n_i, n_j), MIN_SIMILARITY);
      }
    }

    // int start_x = MAX(0, x-similar_window_size/2), start_y = MAX(0, y-similar_window_size/2);
    // int end_x = MIN(img->row, x+similar_window_size/2), end_y = MIN(img->col, y+similar_window_size/2);



    // // return;
    // // if (i%500 == 0) printf("pixel %d / %d\n", i, s->row);
    // for(pos_j[0]=start_x; pos_j[0]<end_x; ++pos_j[0])
    //   for(pos_j[1]=start_y; pos_j[1]<end_y; ++pos_j[1]) {
    //     neighborhood(img, n_j, pos_j);
    //     // printf("s = MAX(%f, %f) = %f\n", -neighborhood_distance(n_i, n_j), MIN_SIMILARITY, MAX(-neighborhood_distance(n_i, n_j), MIN_SIMILARITY));
    //     s->data[i][pos_j[0]*similar_window_size + pos_j[1]] = MAX(-neighborhood_distance(n_i, n_j), MIN_SIMILARITY);
    //     // s->data[i][(pos_j[0]-start_x)*similar_window_size + pos_j[1]-start_y] = -neighborhood_distance(n_i, n_j);
    //   }
  }
  printf("end");
  free_mat(n_i);
  free_mat(n_j);
}

/**
 * Evalue la matrice w(i,j), de dimension n^2*similar_window_size*2 a partir de S(i,j) de meme dimension 
 * et de l'equation $w(i,j) = \frac{1}{Z(i)}e^{-\frac{S(i,j)}{h^2}}
*/
void weights(mat_t* s, mat_t* w,  float h) {
  mat_t* s_h2 = allocate_mat(s->row, s->col);
  copy(s_h2->data, s->data, s->row, s->col);
  printf("1/h^2 = %f, %d, %d\n", 1/SQUARE(h), s_h2->row, s_h2->col);
  mult(s_h2->data, s_h2->row, s_h2->col, 1/SQUARE(h));
  for(int i=0; i<s->row; ++i) {
    softmax(s_h2->data[i], w->data[i], w->col, i);
  }
  // for(int i=0; i<)
  free_mat(s_h2);
}

void weighted_average(mat_t* img, mat_t* w, mat_t* output) {
  printf("a   ");
  int n_pixels = img->row * img->col, similar_window_size = (int)sqrtf(w->col);
  int pos_i[2], pos_j[2];
  copy_zero(output->data, output->row, output->col);
  printf("b   ");
  
  /**
   * Pour chaque pixel i, on evalue son voisinage n_i, et les n_j correspondant 0 la fenetre de similarite (soit 21x21) dans le papier.
   * Pour les pixels i de position < 21x21, comme (0, 0) on ne prend en compte que les pixels voisins valides, i.e. ceux de droite.
   * Pareillement pour les pixels du bords droit de l'image, on ne prend en compte que les pixels voisins valides, i.e. ceux de gauche.
  */
  int x, y, start_x, start_y, end_x, end_y;
  for(int i=0; i<n_pixels; ++i) {
    // printf("x ");
    x = pos_i[0]= (int)i/img->row, y = pos_i[1] = (int)i%img->col;
    start_x = MAX(0, x-similar_window_size/2), start_y = MAX(0, y-similar_window_size/2);
    end_x = MIN(img->row, x+similar_window_size/2), end_y = MIN(img->col, y+similar_window_size/2);
    // printf("%d, %d, %d, %d\n", x-similar_window_size/2, y-similar_window_size/2, x+similar_window_size/2, x+similar_window_size/2);
    if (i%500 == 0) printf("pixel %d / %d \t %d", i, n_pixels, similar_window_size);
    int totalcount = 0;
    float zi = 0;
    start_x = start_y = 0; end_x = w->row, end_y = w->col;
    for(pos_j[0]=start_x; pos_j[0]<end_x; ++pos_j[0])
      for(pos_j[1]=start_y; pos_j[1]<end_y; ++pos_j[1]) {
        // if (i%500 == 0) printf("%f, %f\n", zi, w->data[i][pos_j[0]*similar_window_size + pos_j[1]]);
        zi += w->data[i][pos_j[0]*similar_window_size + pos_j[1]];

        output->data[x][y] += w->data[i][pos_j[0]*similar_window_size + pos_j[1]] * img->data[pos_j[0]][pos_j[1]];
        totalcount++;
      }

    if (i%500 == 0) printf("\t (%d, %d) => (%d, %d) => (%d, %d) . (%d, %f)\n", start_x, start_y, x, y, end_x, end_y, totalcount ,zi);
  }

}

void denoise_image(mat_t* img, mat_t* output, float neighbor_size, float similarity_window_size, float h) {
  int n_pixels = img->row * img->col;
  mat_t* s = allocate_mat(n_pixels, SQUARE(similarity_window_size));
  mat_t* w = allocate_mat(n_pixels, SQUARE(similarity_window_size));
  printf("#1\n");
  similarity(img, s, neighbor_size);
  printf("#2\n");
  weights(s, w, h);
  printf("#3\n");
  weighted_average(img, w, output);
  printf("#4\n");
  
  clean_image(output->data, output->row, output->col);
  free_mat(s);
  free_mat(w);
}


int main(int argc, char **argv)
{
  int i, j, k, l;
  int length, width;
  int nbLevels;

  float threshold;
  float var=400;

  float K = atof(argv[1]);
  float h = K*sqrtf(var);
  
  mat_t *image = load_image(NAME_IMG_IN); /* image d'entree */
  mat_t *noised_image = load_image(NAME_IMG_NOISED);     /* image degradee */
  mat_t *denoised_image = allocate_mat(image->row, image->col);    /* image restoree */

  denoise_image(noised_image, denoised_image, NEIGHBORHOOD_SIZE, SEARCHING_SIZE, h);

  // Recal(denoised_image->data, denoised_image->row, denoised_image->col);
  save_and_show(NAME_IMG_DENOISED, denoised_image->data, denoised_image->row, denoised_image->row);
  save_and_show(NAME_IMG_NOISED, noised_image->data, noised_image->row, noised_image->row);
  
  printf("Runing with h = %f\n", h);
  
  /* liberer la memoire */
  free_mat(image);
  free_mat(noised_image);
  free_mat(denoised_image);


  /*retour sans probleme*/
  printf("\n C'est fini ... \n\n\n");
  return 0;
}
