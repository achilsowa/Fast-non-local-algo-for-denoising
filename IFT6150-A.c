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
#include <time.h>


#include "FonctionDemo.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define NAME_IMG_IN "photograph" //"lenna"

#define NAME_IMG_NOISED "photograph_bruite_B" //"lenna_bruite_A"
#define NAME_IMG_DENOISED "photograph_debruite_B" // "lenna_debruite_A"


#define NEIGHBORHOOD_SIZE 7
#define SEARCHING_SIZE 21
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

  if (hide % 500 == 0 && FALSE) {
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
 * return the position (xi, yi) of left pixel and (xj, yj) of rightmost pixel. possibly negative
 */ 
void neighborhood_positions(int* pixel, int *rect_pos, int M) {
  rect_pos[0] = pixel[0]-M/2;
  rect_pos[1] = pixel[1]-M/2;
  rect_pos[2] = pixel[0]+M/2;
  rect_pos[3] = pixel[1]+M/2;
}

/**
 * save the neihborhood of size (M, M), with 0 where the neihborhood might not be defined (padding)
*/
void neighborhood(mat_t* img, mat_t* n, int *pixel) {
  int M = n->row, x, y;
  copy_val(n->data, n->row, n->col, GREY_LEVEL/2);
  
  for(int i=0; i<M; ++i)
    for(int j=0; j<M; ++j) {
      x = pixel[0] + i - M/2, y = pixel[1] + j - M/2;
      // Si le point (xj, yj) est dans les limites de l'image, on met a jour n, si non, on laisse la valeur par defaut
      if (x < 0 || x >= img->row || y < 0 || y >= img->col) 
        n->data[i][j] = 0*GREY_LEVEL/2;
      else
        n->data[i][j] = img->data[x][y];
  }
}

float neighborhood_distance(mat_t* n_i, mat_t* n_j) {
  float ss = 0;
  int M = n_i->row;
  for(int l=0; l<M; ++l)
    for(int k=0; k<M; ++k)
      ss += SQUARE(n_i->data[l][k] - n_j->data[l][k]);
  return ss;
}

float neighborhood_distance_from_positions(mat_t* img, int* n_i_pos, int* n_j_pos, int M) {
  float ss = 0, vilk, vjlk;
  int nil, nik, njl, njk;
  for(int l=0; l<M; ++l)
    for(int k=0; k<M; ++k) {
      nil = n_i_pos[0]+l, nik = n_i_pos[1]+k;
      njl = n_j_pos[0]+l, njk = n_j_pos[1]+k;
      vilk = is_valid_index_mat(img, nil, nik) ? img->data[nil][nik] : 0;
      vjlk = is_valid_index_mat(img, njl, njk) ? img->data[njl][njk] : 0;
      ss += SQUARE(vilk-vjlk);
    }
  return ss;
}

/**
 * calcule la matrice -S(i,j)
*/
void similarity(mat_t* img, mat_t* s, int neighbor_size) {
  copy_val(s->data, s->row, s->col, MIN_SIMILARITY);
  mat_t* n_i = allocate_mat(neighbor_size, neighbor_size), *n_j = allocate_mat(neighbor_size, neighbor_size);
  int similarity_window_size = (int)sqrtf(s->col);
  int pos_i[2], pos_j[2], r_pos_i[4], r_pos_j[4], xi, yi, xj, yj;

  
  /**
   * Pour chaque pixel i, on evalue son voisinage n_i, et les n_j correspondant 0 la fenetre de similarite (soit 21x21) dans le papier.
   * Pour les pixels i de position < 21x21, comme (0, 0) on ne prend en compte que les pixels voisins valides, i.e. ceux de droite.
   * Pareillement pour les pixels du bords droit de l'image, on ne prend en compte que les pixels voisins valides, i.e. ceux de gauche.
  */
  for(int i=0; i<s->row; ++i) {
    xi = (int)i/img->row, yi = (int)i%img->col;
    pos_i[0] = xi, pos_i[1] = yi;
    // neighborhood(img, n_i, pos_i);
    neighborhood_positions(pos_i, r_pos_i, neighbor_size);
    
    for(int j=0; j<s->col; ++j) {
      xj = xi + (int)j/similarity_window_size - similarity_window_size/2;
      yj = yi + (int)j%similarity_window_size - similarity_window_size/2;

      // Si le point (xj, yj) est dans les limites de l'image, on met a jour S(i, j), si non, on laisse la valeur par defaut
      if (! is_valid_index_mat(img, xj, yj)) {
        s->data[i][j] = MIN_SIMILARITY;
      } else {
        pos_j[0] = xj, pos_j[1] = yj;
        // neighborhood(img, n_j, pos_j);
        neighborhood_positions(pos_j, r_pos_j, neighbor_size);
        // s->data[i][j] = MAX(-neighborhood_distance(n_i, n_j), MIN_SIMILARITY);
        s->data[i][j] = MAX(-neighborhood_distance_from_positions(img, r_pos_i, r_pos_j, neighbor_size), MIN_SIMILARITY);
      }
    }

  }
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
  // printf("1/h^2 = %f, %d, %d\n", 1/SQUARE(h), s_h2->row, s_h2->col);
  mult(s_h2->data, s_h2->row, s_h2->col, 1/SQUARE(h));
  for(int i=0; i<s->row; ++i) {
    softmax(s_h2->data[i], w->data[i], w->col, i);
  }
  // for(int i=0; i<)
  free_mat(s_h2);
}

void weighted_average(mat_t* img, mat_t* w, mat_t* output) {
  int n_pixels = img->row * img->col, similarity_window_size = (int)sqrtf(w->col);
  copy_zero(output->data, output->row, output->col);
  int print_every = 5000000;
  
  /**
   * Pour chaque pixel i, on evalue son voisinage n_i, et les n_j correspondant 0 la fenetre de similarite (soit 21x21) dans le papier.
   * Pour les pixels i de position < 21x21, comme (0, 0) on ne prend en compte que les pixels voisins valides, i.e. ceux de droite.
   * Pareillement pour les pixels du bords droit de l'image, on ne prend en compte que les pixels voisins valides, i.e. ceux de gauche.
  */
  int xi, yi, xj, yj;
  for(int i=0; i<w->row; ++i) {
    xi = (int)i/img->row, yi = (int)i%img->col;
    // if (i%print_every == 0) printf("pixel %d / %d \t %d", i, w->row, similarity_window_size);
    int totalcount = 0;
    float zi = 0;
    
    for(int j=0; j<w->col; ++j) {
      xj = (int)j/similarity_window_size, yj = (int)j%similarity_window_size;
      xj = xi + (int)j/similarity_window_size - similarity_window_size/2;
      yj = yi + (int)j%similarity_window_size - similarity_window_size/2;

      // if (i%print_every == 0) printf("\n%f, %f, (%d, %d), %d", zi, w->data[i][j], xj, yj, is_valid_index_mat(img, xj, yj));
      zi += w->data[i][j];

      if (is_valid_index_mat(img, xj, yj)) output->data[xi][yi] += w->data[i][j] * img->data[xj][yj]; 
      if (w->data[i][j] != 0) totalcount++;
    }

    // if (i%print_every == 0) printf("\t (%d, %d) . (%d, %f)\n", xi, yi, totalcount ,zi);
  }

}

void denoise_image(mat_t* img, mat_t* output, float neighbor_size, float similarity_window_size, float h) {
  time_t t1, t2, t3, t4, t5;
  int n_pixels = img->row * img->col;
  mat_t* s = allocate_mat(n_pixels, SQUARE(similarity_window_size));
  mat_t* w = allocate_mat(n_pixels, SQUARE(similarity_window_size));
  time(&t1);
  
  printf("\n#1 Calcul de la matrice S(i,j):\n");
  similarity(img, s, neighbor_size);
  time(&t2);
  printf("\t %.1fs\n", (float)(t2-t1));

  printf("#2 Calcul de la matrice w(i,j):\n");
  weights(s, w, h);
  time(&t3);
  printf("\t %.1fs\n", (float)(t3-t2));

  printf("#3 Calcul du moyennage par w(i,j):\n");
  weighted_average(img, w, output);
  time(&t4);
  printf("\t %.1fs\n", (float)(t4-t3));
  
  // print_image(output->data, output->row, output->col);
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
  float h = var/K; //best for 3.5 or 4
  
  mat_t *image = load_image(NAME_IMG_IN); /* image d'entree */
  mat_t *noised_image = load_image(NAME_IMG_NOISED);     /* image degradee */
  mat_t *denoised_image = allocate_mat(image->row, image->col);    /* image restoree */

  
  // print_image(noised_image->data, 3, 3); return 0;
  denoise_image(noised_image, denoised_image, NEIGHBORHOOD_SIZE, SEARCHING_SIZE, h);

  printf("ISNR:\t %f\n", isnr(image->data, noised_image->data, denoised_image->data, image->row, image->col));
  
  
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
