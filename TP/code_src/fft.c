#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"
#include "divers.c"

#define PI 3.14159265

#define SWAP(a,b) { tempr=(a);(a)=(b);(b)=tempr; }

static int fournFFT( double *data, long nn[], int ndim, int direct) {  // Direct =1 pour la transformee, -11 poru la transformee inverse
   int idim;
   unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
   unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
   double tempi,tempr;
   double theta,wi,wpi,wpr,wr,wtemp;

  if (direct!=1 && direct!=-1) return -1;

   // Total number of complex values.
   for (ntot=1,idim=0;idim<ndim;idim++)
      ntot *= nn[idim];
   --data;// Due to Numerical Recipies algorithm implementation :: array[1...n];

   nprev=1;
   // main loop over dimensions.

   for (idim=0;idim<ndim;idim++) {
      n=nn[idim];
      nrem=ntot/(n*nprev);
      ip1=nprev << 1; // multiplied by 2.
      ip2=ip1*n;
      ip3=ip2*nrem;
      i2rev=1;

      for (i2=1;i2<=ip2;i2+=ip1) { // Bit reversal section of the routine.
	 if (i2 < i2rev) {
	    for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	       for (i3=i1;i3<=ip3;i3+=ip2) {
		  i3rev=i2rev+i3-i2;
		  SWAP(data[i3],data[i3rev]);
		  SWAP(data[i3+1],data[i3rev+1]);
	       }
	    }
	 }
	 ibit=ip2 >> 1; // Divided by 2.
	 while (ibit >= ip1 && i2rev > ibit) {
	    i2rev -= ibit;
	    ibit >>= 1;
	 }
	 i2rev += ibit;
      }
      ifp1=ip1;

      while (ifp1 < ip2) {
	 ifp2=ifp1 << 1;
	 theta=direct*(M_PI*2)/(ifp2/ip1); // isign = +1 for transform
	 wtemp=sin(0.5*theta);
	 wpr=-2.0*wtemp*wtemp;
	 wpi=sin(theta);
	 wr=1.0;
	 wi=0.0;
	 for (i3=1;i3<=ifp1;i3+=ip1) {
	    for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	       for (i2=i1;i2<=ip3;i2+=ifp2) {
		  k1=i2;
		  k2=k1+ifp1;
		  tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
		  tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
		  data[k2]=data[k1]-tempr;
		  data[k2+1]=data[k1+1]-tempi;
		  data[k1] += tempr;
		  data[k1+1] += tempi;
	       }
	    }
	    wr=(wtemp=wr)*wpr-wi*wpi+wr;
	    wi=wi*wpr+wtemp*wpi+wi;
	 }
	 ifp1=ifp2;
      }
      nprev *= n;
   }

   return 1;
}

#undef SWAP

/*
 * Determines the next power of 2a from num.
 */
int nextpow2( double num ) { double f;
   int exp;
   f=frexp(fabs(num),&exp);
   if (f && f == 0.5) exp = exp-1;
   //return((int)pow(2.,exp));
   return((int)scalbn(1.,exp));
}

/**
 * Gray image.
 * Images MUST have their dimensions all powers of 2 -> reshape.
 */
static int allfft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int direct, int dimx, int dimy) {
   long nn[2];
   unsigned long n;
   double *data, *pf;
   int x,y;
   // Build a vector with couples (ims_reel[i][j],im_imag[i][j]).
   // return dimensions ->all power of 2.
   // New values are set to 0.
   nn[0] = dimx;
   nn[1] = dimy;
   n= 2;
   data=pf=(double*)calloc(2*dimx*dimy,sizeof(double));
   if (!data) { printf("Erreru allocation\n"); return -1; }
   for (y=0; y<dimy; y++) {
      for (x=0;x<dimx;x++) {
         *(pf++) = (double)ims_reel[y][x];
	 *(pf++) = (double)ims_imag[y][x];
      }
      pf+=2*(nn[0]-x);
    }
   pf+=2*((nn[1]-y)*nn[0]);

   fournFFT(data,nn,n,direct);
   pf=data;

   // Build outputs images.
   for (y=0; y<dimy; y++) {
     for (x=0;x<dimx; x++)
        if (direct==1) { imd_reel[y][x]= *(pf++); imd_imag[y][x]= *(pf++); }
        else { imd_reel[y][x]= *(pf++)/(nn[0]*nn[1]); imd_imag[y][x]= *(pf++)/(nn[0]*nn[1]); }
	 // Discard those points.
     pf+=2*(nn[0]-x);
   }
   pf+=2*((nn[1]-y)*nn[0]);
   free(data);
   return 1;
}

int ispowerof2(double num) { return (pow(2.,log2(num))==num); }

int fft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int dimx, int dimy) {
  if (pow(2.,log2(dimx))!=dimx || pow(2.,log2(dimy))!=dimy) { printf("Pas une puissance de 2.\n"); return -1;  }
  return allfft(ims_reel,ims_imag,imd_reel,imd_imag ,1,dimx,dimy);
}

int ifft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int dimx, int dimy) {
  if (pow(2.,log2(dimx))!=dimx || pow(2.,log2(dimy))!=dimy) { printf("Pas une puissance de 2.\n"); return -1;  }
  return allfft(ims_reel,ims_imag,imd_reel,imd_imag ,-1,dimx,dimy);
}


/* Shift les cadres de ims (reeel et imag) dans imd reeel et imaghinaire) */
int fftshift( double** imsr, double** imsi, double** imdr, double** imdi, int nl, int nc ) {
   int midx =nc >> 1;
   int midy = nl >> 1;
   int finx, finy;
   int x,y;

   // Cas impair -> +1
   finx= (nc % 2 == 1)?  midx+1 : midx;
   finy= (nl % 2 == 1)?  midy+1 : midy;

   // Cadre 1 -> Cadre 4.
   for (y=0; y < finy; y++)
      for (x=0; x < finx; x++) {
	 imdr[y][x] = imsr[y+midy][x+midx];
	 imdi[y][x] = imsi[y+midy][x+midx];
      }

   // Cadre 2 -> Cadre 4.
   for (y=0; y < finy; y++)
      for ( x=finx; x < nc; x++) {
	 imdr[y][x] = imsr[y+midy][x-finx];
	 imdi[y][x] = imsi[y+midy][x-finx];
      }


   // Cadre 3 -> cadre 2.
   for (y=finy ; y< nl; y++)
      for ( x=0; x < finx; x++) {
	 imdr[y][x] = imsr[y-finy][x+midx];
	 imdi[y][x] = imsi[y-finy][x+midx];
      }

   // Cadre 4 -> Cadre 1.
   for (y=finy; y < nl; y++)
      for ( x=finx; x < nc; x++) {
	 imdr[y][x] = imsr[y-finy][x-finx];
	 imdi[y][x] = imsi[y-finy][x-finx];
      }
    return 1;
 }

double** padimdforfft(double** im, int* pnl, int* pnc) {
  if (ispowerof2(*pnl) && ispowerof2(*pnc)) return im;
  else { double** res=NULL; int i,j,anl=*pnl,anc=*pnc;
    *pnl=nextpow2(*pnl); *pnc=nextpow2(*pnc);
    if( (res=alloue_image_double(*pnl,*pnc))==NULL) return NULL;
    for (i=0; i<anl; i++) for(j=0;j<anc; j++) res[i][j]=im[i][j];
    return res;
  }
}

double** padimucforfft(unsigned char** im, int* pnl, int* pnc) {
  double** res=NULL; int i,j,anl=*pnl,anc=*pnc;
  *pnl=nextpow2(*pnl); *pnc=nextpow2(*pnc);
  if( (res=alloue_image_double(*pnl,*pnc))==NULL) return NULL;
  for (i=0; i<anl; i++) for(j=0;j<anc; j++) res[i][j]=im[i][j];
  return res;
}

double** padimd(double** im, int nl, int nc, int anl, int anc) { int i,j;
  double **res;
  if( (res=alloue_image_double(nl,nc))==NULL) return NULL;
  for (i=0; i<anl; i++) for(j=0;j<anc; j++) res[i][j]=im[i][j];
     return res;
}

double** norme(double** real, double** imag, int nl, int nc) { int i;
  double** res;
  if( (res=alloue_image_double(nl,nc))==NULL) return NULL;
  for (i=0; i<nl*nc; i++) (*res)[i]=hypot((*real)[i],(*imag)[i]);
  return res;
}

double** phase(double** real, double** imag, int nl, int nc) { int i;
  double** res;
  if( (res=alloue_image_double(nl,nc))==NULL) return NULL;
  for (i=0; i<nl*nc; i++) (*res)[i]=atan2((*imag)[i],(*real)[i]);
  return res;
}

#define HAMMING(t,n) (0.56+0.44*cos(2*M_PI*(t)/(double)(n)))
double** hamming_double(double** im, double** res,int  nl, int nc){ int i,j;
   if( res==NULL && (res=alloue_image_double(nl,nc))==NULL) return NULL;
   for (i=0; i<nl; i++) for(j=0;j<nc; j++) res[i][j]=im[i][j]*HAMMING(i-nl/2,nl)*HAMMING(j-nc/2,nc);
   return res;
}

double** hamming_uc(unsigned char** im, int  nl, int nc){ int i,j;
   double **res;
   if( (res=alloue_image_double(nl,nc))==NULL) return NULL;
   for (i=0; i<nl; i++) for(j=0;j<nc; j++) res[i][j]=im[i][j]*HAMMING(i,nl)*HAMMING(j,nc);
   return res;
}

double** fftFiltreGaussien(double variance, int nl, int nc){
    double** fftFiltre = alloue_image_double(nl,nc);
    for(int u = 0; u<nl; u++){
        for(int v = 0; v<nc; v++){
            double norm = pow(((u - nl/2)/(double)nl), 2) + pow(((v - nc/2)/(double)nc), 2);
            fftFiltre[u][v] = exp(-2*pow(PI, 2)*pow(variance, 2)*norm);
        }
    }

    return fftFiltre;
}


void filtrageGaussienImg(double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , double variance, int dimx, int dimy){
    double ** tfImage_reel = alloue_image_double(dimx, dimy);
    double ** tfImage_imag = alloue_image_double(dimx, dimy);

    double ** shiftTfImage_reel = alloue_image_double(dimx, dimy);
    double ** shiftTfImage_imag = alloue_image_double(dimx, dimy);

    double ** produit_reel = alloue_image_double(dimx, dimy);
    double ** produit_imag = alloue_image_double(dimx, dimy);

    double ** shiftProduit_reel = alloue_image_double(dimx, dimy);
    double ** shiftProduit_imag = alloue_image_double(dimx, dimy);

    if(!fft(ims_reel, ims_imag, tfImage_reel, tfImage_imag , dimx, dimy)){
      printf("erreur de FFT");
    };
    if(!fftshift(tfImage_reel, tfImage_imag, shiftTfImage_reel, shiftTfImage_imag, dimx, dimy)){
      printf("erreur du premier shift");
    };

    double ** filtre = fftFiltreGaussien(variance, dimx, dimy);

    for(int u = 0; u<dimx; u++){
        for(int v = 0; v<dimy; v++){
            produit_reel[u][v] = filtre[u][v]*shiftTfImage_reel[u][v];
            produit_imag[u][v] = filtre[u][v]*shiftTfImage_imag[u][v];
        }
    }

    if(!fftshift(produit_reel, produit_imag, shiftProduit_reel, shiftProduit_imag, dimx, dimy)){
      printf("erreur du deuxieme shift");
    };

    if(!ifft(shiftProduit_reel, shiftProduit_imag, imd_reel, imd_imag, dimx, dimy)){
      printf("erreur de IFFT");
    };
    free(tfImage_reel);
    free(tfImage_imag);
    free(shiftTfImage_reel);
    free(shiftTfImage_imag);
    free(produit_reel);
    free(produit_imag);
    free(shiftProduit_reel);
    free(shiftProduit_imag);
}

double *PsnrFiltreGaussien(double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag, int dimx, int dimy){
  double* output = NULL;
  output = malloc(40*sizeof(double));
  int comp = 0;
  for(double variance = 1; variance < 20; variance+=0.5){
    filtrageGaussienImg(ims_reel, ims_imag, imd_reel, imd_imag , variance, dimx, dimy);
    output[comp] = psnr_double(ims_reel, imd_reel, dimx, dimy);
    printf("%f\n", output[comp]);
    comp+=1;
  }
  return output;
}

double cgimsLocal(double** ims, int posx, int posy, double resLissage, int dimMask, int dimx, int dimy){
  double frac = (1/(2*PI*pow(resLissage, 2)));
  double somme1 = 0;
  for(int i = -dimMask; i < dimMask; i++){
    double somme0 = 0;
    for(int j = -dimMask; j < dimMask; j++){
      somme0 += exp(-pow(j/resLissage,2) / 2) * ims[(posx + i + dimx) % dimx][(posy +j+ dimy) % dimy];
    }
    somme1 += somme0 * exp(-pow(i/resLissage,2) / 2);
  }
  return(somme1*frac);
}

double** cgims(double** ims, int dimMask, double resLissage, int dimx, int dimy){
  double** output = NULL;
  output=alloue_image_double(dimx,dimy);
  for(int x = 0; x <dimx; x++){
    for(int y = 0; y <dimy; y++){
      output[x][y] = cgimsLocal(ims, x, y, resLissage, dimMask, dimx, dimy);
    }
  }
  return output;
}

double* eqmConv(double variance, double** ims, int dimx, int dimy) {
  double dimMask;
  double* eqmData = NULL;
  eqmData = malloc(6 * sizeof(double));
  double** imsi = alloue_image_double(dimx, dimy);
  double** im1 = alloue_image_double(dimx, dimy);
  double** im1i = alloue_image_double(dimx, dimy);

  filtrageGaussienImg(ims, imsi, im1, im1i , variance, dimx, dimy);
  double** im2 = NULL;
  for(int i = 1; i<=6; i++) {
    dimMask = i * variance;
    im2 = cgims(ims, dimMask, variance, dimx, dimy);
    eqmData[i-1] = eqm_double(im1, im2, dimx, dimy);
  }
  free(imsi);
  free(im1);
  free(im1i);
  return eqmData;
}

void eqmTotal(double** ims1, double** ims2, int dimx, int dimy, int dimx2, int dimy2) {
  double* eqmData = NULL;
  eqmData = eqmConv(1, ims1, dimx, dimy);
  printf("Im1 => var = 1 : \n");
  printEqm(eqmData);
  eqmData = eqmConv(5, ims1, dimx, dimy);
  printf("Im1 => var = 5 : \n");
  printEqm(eqmData);
  eqmData = eqmConv(10, ims1, dimx, dimy);
  printf("Im1 => var = 10 : \n");
  printEqm(eqmData);
  eqmData = eqmConv(20, ims1, dimx, dimy);
  printf("Im1 => var = 20 : \n");
  printEqm(eqmData);

  printf("\n");
  eqmData = eqmConv(1, ims2, dimx2, dimy2);
  printf("Im2 => var = 1 : \n");
  printEqm(eqmData);
  eqmData = eqmConv(5, ims2, dimx2, dimy2);
  printf("Im2 => var = 5 : \n");
  printEqm(eqmData);
  eqmData = eqmConv(10, ims2, dimx2, dimy2);
  printf("Im2 => var = 10 : \n");
  printEqm(eqmData);
  eqmData = eqmConv(20, ims2, dimx2, dimy2);
  printf("Im2 => var = 20 : \n");
  printEqm(eqmData);
}

void printEqm(double* eqmData) {
  for(int i = 0; i < 6; i++) {
    printf("%f\n", eqmData[i]);
  }
  printf("\n");
}

double** lissageGaussien(double variance){
  double ** filtre = alloue_image_double(3, 3);
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      filtre[i][j] = (1/(PI*pow(variance, 4))) * (((pow(i, 2) + pow(j, 2))/(2*pow(variance, 2)))-1) * exp(-((pow(i, 2) + pow(j, 2))/(2*pow(variance, 2))));
    }
  }
  return filtre;
}

double ** logMasque(double variance){
  double ** lissage = lissageGaussien(variance);
  double ** masque = alloue_image_double(3, 3);
  masque[0][0] = lissage[0][1];
  masque[0][1] = lissage[0][0] - 4*lissage[0][1] + lissage[0][2];
  masque[0][2] = lissage[0][1];
  masque[1][0] = lissage[1][1];
  masque[1][1] = lissage[1][0] - 4*lissage[1][1] + lissage[1][2];
  masque[1][2] = lissage[1][1];
  masque[2][0] = lissage[2][1];
  masque[2][1] = lissage[2][0] - 4*lissage[2][1] + lissage[2][2];
  masque[2][2] = lissage[2][1];
  return masque;
}

double ** laplacienFiltre(double **ims, int dimx, int dimy){
  double ** masque = logMasque(4);
  double** output = alloue_image_double(dimx, dimy);
  printf("%d %d\n", dimx, dimy);
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      output[i][j] = gradientPixel(ims, i, j, masque, dimx, dimy);
    }
  }
  return output;
}

void compareEtSupprLapl(double ** ims, int dimx, int dimy, double** lapls) {
  for(int i = 1; i < dimx-1; i++) {
    for(int j = 1; j < dimy-1; j++) {
      if(lapls[i][j] > 0) {
        if(lapls[i+1][j] >= 0 && lapls[i][j+1] >= 0 && lapls[i][j-1] >= 0 && lapls[i-1][j] >= 0) {
          ims[i][j] = 0;
        }
      }
      else if(lapls[i][j] < 0) {
        if(lapls[i+1][j] <= 0 && lapls[i][j+1] <= 0 && lapls[i][j-1] <= 0 && lapls[i-1][j] <= 0) {
          ims[i][j] = 0;
        }
      }
      else {
        ims[i][j] = 0;
      }
    }
  }
}

void compareEtSuppr(double ** ims, int dimx, int dimy, double** grads, double** p1, double** p2) {
  for(int i = 1; i < dimx-1; i++) {
    for(int j = 1; j < dimy-1; j++) {
      if(grads[i][j] <= p1[i][j] || grads[i][j] <= p2[i][j]) {
        ims[i][j] = 0.;
      }
      else if(grads[i][j] < 10) {
        ims[i][j] = 0.;
      }
      else if(grads[i][j] < 80 && (grads[i-1][j] < 80 && grads[i+1][j] < 80 && grads[i][j-1] < 80 && grads[i-1][j+1] < 80 && grads[i][j+1] < 80 && grads[i-1][j-1] < 80 && grads[i+1][j+1] < 80 && grads[i+1][j-1] < 80)) {
        ims[i][j] = 0.;
      }
    }
  }
}

void calculGradiantP1P2(double ** ims, int dimx, int dimy, double** grads, double** p1, double** p2, double** masque, double** masque2){
  for(int i = 1; i < dimx-1; i++) {
    for(int j = 1; j < dimy-1; j++) {
      double gradx = gradientPixel(ims, i, j, masque, dimx, dimy);
      double grady = gradientPixel(ims, i, j, masque2, dimx, dimy);
      if(gradx != 0) {
        p1[i][j] = (grady/gradx) * grads[i][j-1] + (gradx - grady)*grads[i+1][j]/gradx;
        p2[i][j] = (grady/gradx) * grads[i][j+1] + (gradx - grady)*grads[i-1][j]/gradx;
      }
      else {
        p1[i][j] = 0;
        p2[i][j] = 0;
      }
    }
  }
}


double ** calculGradiant(double ** ims, int dimx, int dimy, double** masque, double** masque2, double** phis){
  double** grads = alloue_image_double(dimx, dimy);
  printf("%d %d\n", dimx, dimy);
  for(int i = 0; i < dimx; i++) {
    for(int j = 0; j < dimy; j++) {
      double gradx = 0;
      double grady = 0;
      //printf("%d %d\n", i , j);
      gradx = gradientPixel(ims, i, j, masque, dimx, dimy);
      grady = gradientPixel(ims, i, j, masque2, dimx, dimy);
      grads[i][j] = sqrt(pow(gradx, 2) + pow(grady, 2));
      phis[i][j] = atan(grady/gradx)*180/PI;
    }
  }
  printf("Pixel : %.0f\n", ims[1][1]);
  printf("\n\n\nGradientx = %.0f\n\n\n", gradientPixel(ims, 1, 1, masque, dimx, dimy));
  return grads;
}

double gradientPixel(double** ims, int i, int j, double** masque, int dimx, int dimy) {
  double grad = 0;
  grad = grad + masque[0][0] * ims[(i-1+dimx)%dimx][(j-1+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[0][0], ims[(i-1+dimx)%dimx][(j-1+dimy)%dimy], grad);
  grad = grad + masque[1][0] * ims[(i+dimx)%dimx][(j-1+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[1][0], ims[(i+dimx)%dimx][(j-1+dimy)%dimy], grad);
  grad = grad + masque[2][0] * ims[(i+1+dimx)%dimx][(j-1+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[2][0], ims[(i+1+dimx)%dimx][(j-1+dimy)%dimy], grad);
  grad = grad + masque[0][1] * ims[(i-1+dimx)%dimx][(j+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[0][1], ims[(i-1+dimx)%dimx][(j+dimy)%dimy], grad);
  grad = grad + masque[1][1] * ims[(i+dimx)%dimx][(j+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[1][1], ims[(i+dimx)%dimx][(j+dimy)%dimy], grad);
  grad = grad + masque[2][1] * ims[(i+1+dimx)%dimx][(j+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[2][1], ims[(i+1+dimx)%dimx][(j+dimy)%dimy], grad);
  grad = grad + masque[0][2] * ims[(i-1+dimx)%dimx][(j+1+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[0][2], ims[(i-1+dimx)%dimx][(j+1+dimy)%dimy], grad);
  grad = grad + masque[1][2] * ims[(i+dimx)%dimx][(j+1+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[1][2], ims[(i+dimx)%dimx][(j+1+dimy)%dimy], grad);
  grad = grad + masque[2][2] * ims[(i+1+dimx)%dimx][(j+1+dimy)%dimy];
  //printf("facteur : %f * %f = %.0f\n", masque[2][2], ims[(i+1+dimx)%dimx][(j+1+dimy)%dimy]//, grad);
  //printf("FINI");
  return grad;
}

void PrintMat(double** mat, int dimx, int dimy){

  for(unsigned int  i = 0; i < dimy; i++){
    printf("\n");
    //printf(" |");
    for(unsigned int  j = 0; j < dimx; j++){
      //if(mat[i][j] != 0)
        printf(" %.5f |", mat[i][j]);
    }
    printf("\n");

  }
  printf("\n");

}
