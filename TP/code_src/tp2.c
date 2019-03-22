#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

#define PI 3.14159265

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
  return 1.0;
}
