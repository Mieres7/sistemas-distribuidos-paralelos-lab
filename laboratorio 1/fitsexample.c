#include "fitsio.h" /* required by every program that uses CFITSIO */
#include "stdio.h"

int main(int argc, char *argv[])
{
  fitsfile *fptr, *gptr; /* pointer to the FITS file; defined in fitsio.h */
  int status;
  long fpixel=1, naxis = 2, nelements, exposure;
  long naxes[2];
 
  status = 0; /* initialize status before calling fitsio routines */
 
  fits_open_file(&fptr, argv[1], READONLY, &status);
  
  fits_get_img_size(fptr, 2, naxes, &status);
  printf("naxes[0] = %ld\n", naxes[0]);
  double *myimage = (double *) malloc(naxes[0]*naxes[1]*sizeof(double));
 
  printf("size = (%ld, %ld)\n", naxes[1], naxes[0]);

  fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], NULL, myimage, NULL, &status);
  printf("%.10e \n", myimage[1024*512+128]);
  
  // int i, j;
  // double max = -1.0e-6;
  // for (i=0; i < naxes[0]*naxes[1]; i++)

  //   // if(myimage[i] == 255){
  //   //   printf("%d",i);
  //   //   break;
  //   // }

  //   if (myimage[i] > max)
  //     max = myimage[i];
      
  // printf("max = %f\n", myimage[12]);
  
  exit(0);
  
  return( status );
}

// gcc -o fitsexample fitsexample.c -lcfitsio
// ./fitsexample AS209_edge.fits
