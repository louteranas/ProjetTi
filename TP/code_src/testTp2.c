#include "pgm.h"

	/*
		Exemple de code avec Entrees Sortie et transformations simples d'images
		S'utilise sous la forme  "exemple tangram.pgm res.pgm"
 	*/
int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL;
  double** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;
  double** m1, ** m2, ** grads, **lapls, ** phis, **p1, **p2;

  if (ac < 3) {printf("Usage : %s entree sortie\n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
	/* Calcul de son inverse video */
  double**im3=imuchar2double(im1,nl,nc);
  oldnl=nl; oldnc=nc;
	/*  la fft demande des puissances de 2. On padde avec des 0, mais les dimensions nl et nc changent */
  im7=padimdforfft(im3,&nl,&nc);
  //PrintMat(im7, nl, nc);
  printf("nl, nc = %d, %d\n", oldnl, oldnc);
  im5 = nlmeans(im7, nl, nc, 10, 1, 10, 4);
  double psn = psnr_double(im7, im5, nl, nc);
  printf("psnr : %f\n", psn);
  //PrintMat(im5, nl, nc);
  ecritureimagepgm(av[2],crop(imdouble2uchar(im5,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
  return 1;
}
