#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

#define PI 3.14159265

void PrintVect( double* vect, int taille){
  for(int i = 0; i < taille; i++){
    printf("histogramme[%d] = %.0f.\n", i, vect[i]);
  }
}

double* histogrammeCreat(int* taille, double* min, double** ecartsLocaux, int dimx, int dimy, int cumul){
  *min = 3000;
  double maxEcr = 0;
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      if(ecartsLocaux[i][j]<*min)
        *min = ecartsLocaux[i][j];
      if(ecartsLocaux[i][j]>maxEcr)
        maxEcr = ecartsLocaux[i][j];
    }
  }
  *taille = (int)(maxEcr/0.01) - (int)(*min/0.01) + 2;
  double *histogramme = NULL;
  histogramme = malloc(*taille*sizeof(double));
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      int indice =(int)((ecartsLocaux[i][j] - *min)/0.01);
      if(cumul == 1){
        for(int i = *taille; i >= indice; i--)
          histogramme[i]++;
      }
      else{
        histogramme[indice]++;
      }
    }
  }
  PrintVect(histogramme, *taille);
  return histogramme;
}

double estimBruit(double** ims, int dimx, int dimy, int tailleBlocs, double pourcentile) {
  double** masque = alloue_image_double(3, 3);
  masque[0][0] = 0; masque[0][1] = -((double)1/5); masque[0][2] = 0;
  masque[1][0] = -((double)1/5); masque[1][1] = 1; masque[1][2] = -((double)1/5);
  masque[2][0] = 0; masque[2][1] = -((double)1/5); masque[2][2] = 0;
  double** imsFiltree = alloue_image_double(dimx, dimy);
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      imsFiltree[i][j] = gradientPixel(ims, i, j, masque, dimx, dimy);
    }
  }
  // PrintMat(imsFiltree, dimx, dimy);
  double nbElemsBloc = pow(2*tailleBlocs+1, 2);
  double** ecartsLocaux = alloue_image_double(dimx, dimy);
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      double ecLocal = 0;
      double somme = 0;
      double sommeCarre = 0;
      for(int k = -tailleBlocs; k <= tailleBlocs; k++) {
        for(int l = -tailleBlocs; l <= tailleBlocs; l++) {
          somme += imsFiltree[(i+k+dimx)%dimx][(j+l+dimy)%dimy];
          sommeCarre += pow(imsFiltree[(i+k+dimx)%dimx][(j+l+dimy)%dimy], 2);
        }
      }
      double carreMoyenne = pow(somme / nbElemsBloc, 2);
      //printf("car = %f\n", carreMoyenne);
      double moyDesCarres = sommeCarre / nbElemsBloc;
      //printf("moy = %f\n", moyDesCarres);
      ecartsLocaux[i][j] = sqrt(moyDesCarres - carreMoyenne);
    }
  }
  PrintMat(ecartsLocaux, dimx, dimy);
  int taille;
  double min;
  double *hist = histogrammeCreat(&taille, &min, ecartsLocaux, dimx, dimy, 1);
  double seuil = dimx * dimy * pourcentile * 0.01;
  printf("seuil = %f\n", seuil);
  printf("taille = %d\n", taille);
  int i = 0;
  double cur = hist[i];
  while(i < taille && cur < seuil) {
    i++;
    cur = hist[i];
  }
  printf("seuil atteint à i = %d\n", i);
  double ecartFinal = 1.13 * (min + i * 0.01);
  return ecartFinal;
}
