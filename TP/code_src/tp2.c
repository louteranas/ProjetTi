#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

#define PI 3.14159265
#define max(X, Y) (((X) < (Y)) ? (Y) : (X))

void PrintVect( double* vect, int taille){
  for(int i = 0; i < taille; i++){
    printf("histogramme[%d] = %.0f.\n", i, vect[i]);
  }
}

// ----------------------------- ESTIMATION DU BRUIT ---------------------------------------------------

double* histogrammeCreat(int* taille, double* min, double** ecartsLocaux, int dimx, int dimy, int cumul, double pas){
  *min = 30000000;
  double maxEcr = 0;
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      if(ecartsLocaux[i][j]<*min)
        *min = ecartsLocaux[i][j];
      if(ecartsLocaux[i][j]>maxEcr)
        maxEcr = ecartsLocaux[i][j];
    }
  }
  *taille = (int)(maxEcr/pas) - (int)(*min/pas) + 2;
  double *histogramme = NULL;
  histogramme = malloc(*taille*sizeof(double));
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      int indice =(int)((ecartsLocaux[i][j] - *min)/pas);
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
  double *hist = histogrammeCreat(&taille, &min, ecartsLocaux, dimx, dimy, 1, 0.01);
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

// ----------------------------- FILTRE NLMEANS ---------------------------------------------------


double **nlmeans(double** ims, int dimx, int dimy, int sizeRegion, int sizePatch, double sigma, double h) {
  //Déclaration des variables temporaires
  double sommeNum = 0; double sommeDenom = 0;
  double** wRegion = NULL;
  double d; double w;
  double** output = alloue_image_double(dimx, dimy);
  //Pour tous les pixels de l'image
  printf("\n");
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      //Pour tous les pixels de la région
      wRegion = alloue_image_double(2*sizeRegion + 1, 2*sizeRegion + 1);
      for(int x = i - 2*sizeRegion; x <= i + 2*sizeRegion; x++) {
        for(int y = j - 2*sizeRegion; y <= j + 2*sizeRegion; y++) {
          d = calculD(ims, i, j, x, y, dimx, dimy, sizePatch);
          w = max(d-2*sigma*sigma, 0);
          w = exp(-(w/h*h));
          sommeNum += w * ims[(x+dimx)%dimx][(y+dimy)%dimy];
          sommeDenom += w;
        }
      }
      //Calcul final du pixel
      output[i][j] = sommeNum /sommeDenom;
      //On réinitialise les variable temporaires
      free(*wRegion); free(wRegion);
      sommeNum = 0; sommeDenom = 0;
    }
    printf("%d \n", i);
  }
  printf("\n");
  return output;
}

double calculD(double **ims, int X, int Y, int x, int y, int dimx, int dimy, int sizePatch) {
  double sommeTemp = 0; double sommeTot = 0;
  for(int i = -sizePatch; i <= sizePatch; i++) {
    for(int j = -sizePatch; j <= sizePatch; j++) {
      sommeTemp += pow(ims[(x+i+dimx)%dimx][(y+j+dimy)%dimy] - ims[(X+i+dimx)%dimx][(Y+j+dimy)%dimy], 2);
    }
    sommeTot += sommeTemp;
  }
  sommeTot = sommeTot / pow(2*sizePatch + 1, 2);
  return sommeTot;
}

// ----------------------------- FILTRE ADAPTATIF RECURSIF---------------------------------------------------

double** adaptRecursif(double** ims, int dimx, int dimy) {
  
  return ims;
}

// ----------------------------- FILTRE BILATERAL ---------------------------------------------------



// ----------------------------- FILTRE MEDIAN ---------------------------------------------------
