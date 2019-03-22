#include "pgm.h"
	/*
		Exemple de code avec Entrees Sortie et transformations simples d'images
		S'utilise sous la forme  "exemple tangram.pgm res.pgm"
 	*/
main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc, nl2, nc2, oldnl2,oldnc2;
  unsigned char **im2=NULL,** im1=NULL;
  double** im5, ** im6;

  //if (ac < 3) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm("../imagestp/Archive/formes2g.pgm",&nl,&nc);
  im2=lectureimagepgm("../imagestp/radio1.pgm",&nl2,&nc2);

  if (im1==NULL)  { puts("Lecture image 1 impossible"); exit(1); }
  if (im2==NULL)  { puts("Lecture image 2 impossible"); exit(1); }

	/* Calcul de son inverse video */
  double**im3=imuchar2double(im1,nl,nc);
  double**im4=imuchar2double(im2,nl2,nc2);

  oldnl=nl; oldnc=nc;
  oldnl2=nl2; oldnc2=nc2;

  im5=padimdforfft(im3,&nl,&nc);
  im6=padimdforfft(im4,&nl2,&nc2);

  // double** testConv = cgims(im7, 5, 10, nl, nc);
  // ecritureimagepgm(av[2],crop(imdouble2uchar(testConv,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);

  eqmTotal(im5, im6, nl, nc, nl2, nc2);
}
