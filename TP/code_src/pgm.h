#ifndef _PGM
#define _PGM 1
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#define O_BINARY 0
#define IO_LEN  (1<<30)

#ifndef FALSE
#define FALSE (0)
#define TRUE (1)
#endif

#define MAGIC_PGM       "P5\n"
#define MAGIC_PPM       "P6\n"

unsigned char** lectureimagepgm(char* fic, int* plignes, int* pcols);
void ecritureimagepgm(char* nom, unsigned char** im, int lignes, int cols);
double** lectureimagedoubleraw(char* fic, int plignes, int pcol);
void ecritureimagedoubleraw(char*, double**,int lignes, int cols);

unsigned char ** alloue_image(int nl, int nc);
double ** alloue_image_double(int nl, int nc);
void libere_image(unsigned char** im) ;

 double** imuchar2double(unsigned char **im, int nl, int nc);
unsigned char**imdouble2uchar(double** im,int nl, int nc);
char**imdouble2char(double** im,int nl, int nc);
double ** crop_double(double **im,int oi, int oj, int fi, int fj);
unsigned char** crop(unsigned char **im,int oi, int oj, int fi, int fj);

int fft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int nl, int nc);
int ifft( double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , int nl,int nc);

int fftshift( double** imsr, double** imsi, double** imdr, double** imdi, int nl, int nc );
int nextpow2( double num );
int ispowerof2(double num);
double** padimdforfft(double** im, int* pnl, int* pnc);
double** padimd(double** im, int nl, int nc, int anl, int anc);
double** padimucforfft(unsigned char ** im, int* pnl, int* pnc);
double eqm(unsigned char **im1, unsigned char **im2,  int nl, int nc);
double psnr(unsigned char **im1, unsigned char **im2,  int nl, int nc) ;
double psnr_double(double** r, double** i, int nl, int nc);
double eqm_double(double** r, double** i, int nl, int nc);

//TP1
double** fftFiltreGaussien(double variance, int nl, int nc);
void filtrageGaussienImg(double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag , double variance, int dimx, int dimy);
double** cgims(double** ims, int dimMask, double resLissage, int dimx, int dimy);
double cgimsLocal(double** ims, int posx, int posy, double resLissage, int dimMask, int dimx, int dimy);
double *PsnrFiltreGaussien(double** ims_reel, double** ims_imag, double** imd_reel, double** imd_imag, int dimx, int dimy);
double* eqmConv(double variance, double** ims, int dimx, int dimy);
void eqmTotal(double** ims1, double** ims2, int dimx, int dimy, int dimx2, int dimy2);
void printEqm(double* eqmData);
double ** calculGradiant(double ** ims, int dimx, int dimy, double** masque, double** masque2, double** phis);
double gradientPixel(double** ims, int i, int j, double** masque, int dimx, int dimy);
void PrintMat(double** mat, int dimx, int dimy);
void calculGradiantP1P2(double ** ims, int dimx, int dimy, double** grads, double** p1, double** p2, double** masque, double** masque2);
void compareEtSuppr(double ** ims, int dimx, int dimy, double** grads, double** p1, double** p2);
double** lissageGaussien(double variance);
double ** laplacienFiltre(double **ims, int dimx, int dimy);
double ** logMasque(double variance);
void compareEtSupprLapl(double ** ims, int dimx, int dimy, double** lapls);

//TP2
double estimBruit(double** ims, int dimx, int dimy, int tailleBlocs, double pourcentile);
void PrintVect( double* vect, int taille);
double* histogrammeCreat(int* taille, double* min, double** ecartsLocaux, int dimx, int dimy, int cumul, double pas);
double **nlmeans(double** ims, int dimx, int dimy, int sizeRegion, int sizePatch, double sigma, double h);
double calculD(double **ims, int X, int Y, int x, int y, int dimx, int dimy, int sizePatch);
double** adaptRecursif(double** ims, int dimx, int dimy, double k);
double ** bilateral( double ** ims, int dimx, int dimy, double sigma1, double sigma2);
double ** iterationS( double ** st, int dimx, int dimy, double ** w);
double** median(double** ims, int dimx, int dimy, int n);
int verifStationnaire(double** s, double** scopy, int dimx, int dimy);


double** norme(double** real, double** imag, int nl, int nc);
double** phase(double** real, double** imag, int nl, int nc);
double** hamming_double(double** im, double**res, int  nl, int nc);
double** hamming_uc(unsigned char** im, int  nl, int nc);
