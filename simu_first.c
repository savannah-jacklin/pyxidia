#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fitsio.h> 

#define LENMAX 600

int ximlist(char *image,double *uu);
void writeimage(int sizex,int sizey,double *pval,char *filename);
double *array(int n);

void errore(char *str);
void pausa(void);


double TEXP = 60;
double JD = 7000;

//example main to read a fits file line by line in a single zig-zag array
//this version crops to a specific size defined by x0 and y0
int main(int argc, char *argv[])
{
  int a,*pa,ix,jy,x0,y0;
  int size;
  double *uu,*pu;
  char *ima,str[LENMAX];


  if(argc!=2) goto Usage;
  ima = argv[1];

  printf(" bla\n");
  
  a = 2;
  pa = &a;
  printf(" a = %d\n",a);
  printf(" *pa = %d\n",*pa);
  size = 100;
  uu = array(size*size);

  x0 = 800;
  y0 = 3000;
  sprintf(str,"%s[%d:%d,%d:%d]",ima,x0,x0+size-1,y0,y0+size-1);
  ximlist(str,uu);
  
  pu = uu;
  for(jy=0; jy<size; jy++) {
    for(ix=0; ix<size; ix++,pu++) {
      if(isnan(*pu)) *pu = 0;
      //printf(" %4d %4d %g %d\n",ix,jy,*pu,isnan(*pu));      
      //getchar();
    }
  }
  writeimage(size,size,uu,"name.fits");  
  



  exit(0);
 Usage:
  printf(" Usage: ./simu {image}\n");
  exit(1);

}


double *array(int n)
{
  double *a;
  a=(double *) calloc(1,(size_t)(n*sizeof(double)));
  return a;
}

//write a fits file
void writeimage(int sizex,int sizey,double *pval,char *filename)
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status, ii, jj;
  int ival;
  long  fpixel, nelements, exposure;
  double *array[sizey];
  
  /* initialize FITS image parameters */
  //  int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
  int bitpix = DOUBLE_IMG; /* 16-bit unsigned short pixel values       */
  long naxis = 2;  /* 2-dimensional image                            */    
  long naxes[2];// = { 300, 200 };   /* image is 300 pixels wide by 200 rows */


  //if(b==0) sprintf(filename,"bla2_%d.fits",NB_GRID);
  //else sprintf(filename,"ela2_%d.fits",NB_GRID);
  naxes[0] = sizex;
  naxes[1] = sizey;
  
  /* allocate memory for the whole image */ 
  //  array[0] = (unsigned short *)malloc( naxes[0] * naxes[1] * sizeof( unsigned short ) );
  array[0] = (double *)malloc( naxes[0] * naxes[1] * sizeof( double ) );

  /* initialize pointers to the start of each row of the image */
  for( ii=1; ii<naxes[1]; ii++ ) array[ii] = array[ii-1] + naxes[0];
  
  remove(filename);               /* Delete old file if it already exists */
  
  status = 0;         /* initialize status before calling fitsio routines */
  
  if (fits_create_file(&fptr, filename, &status)) errore("wri_ima"); /* create new FITS file */
   
  /* write the required keywords for the primary array image.     */
  /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
  /* a FITS image with BITPIX = 16 (signed short integers) with   */
  /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
  /* FITS uses to store unsigned integers.  Note that the BSCALE  */
  /* and BZERO keywords will be automatically written by cfitsio  */
  /* in this case.                                                */
  
  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) ) errore("wri_ima");
   
  /* initialize the values in the image with a linear ramp function */
  for (jj = 0; jj < naxes[1]; jj++) {   
    for (ii = 0; ii < naxes[0]; ii++) {
      ival = ii + jj * sizex;
      array[jj][ii] = pval[ival];
    }
  }
  
  fpixel = 1;                               /* first pixel to write      */
  nelements = naxes[0] * naxes[1];          /* number of pixels to write */
  
  /* write the array of unsigned integers to the FITS file */
  //if ( fits_write_img(fptr, TUSHORT, fpixel, nelements, array[0], &status) )
  if ( fits_write_img(fptr, TDOUBLE, fpixel, nelements, array[0], &status) ) errore("wri_ima");
  
  free( array[0] );  /* free previously allocated memory */
  
  /* write another optional keyword to the header */
  /* Note that the ADDRESS of the value is passed in the routine */
  exposure = TEXP;
  if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, 
		       "Total Exposure Time", &status) ) errore("wri_ima");

  if ( fits_update_key(fptr, TDOUBLE, "JD", &JD, "JD", &status) ) errore("wri_ima");

  if ( fits_close_file(fptr, &status) ) errore("wri_ima"); /* close the file */
  
  return;
}

//list parameters of fits file? ask sebastiano
int ximlist(char *image,double *uu)
{
  fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int bitpix, naxis, ii, anynul;
  long naxes[2] = {1,1}, fpixel[2] = {1,1};
  double *pixels;
  char format[20], hdformat[20];
  int kk = 0;
  //char image[200];

  //  sprintf(image,"%s[1:125,1:125]",str);
  
  /*
    if (argc != 2) {
    printf("Usage:  imlist filename[ext][section filter] \n");
    printf("\n");
    printf("List the the pixel values in a FITS image \n");
    printf("\n");
    printf("Example: \n");
    printf("  imlist image.fits                    - list the whole image\n");
    printf("  imlist image.fits[100:110,400:410]   - list a section\n");
    printf("  imlist table.fits[2][bin (x,y) = 32] - list the pixels in\n");
    printf("         an image constructed from a 2D histogram of X and Y\n");
    printf("         columns in a table with a binning factor = 32\n");
    return(0);
    }
  */
  
  if (!fits_open_file(&fptr, image, READONLY, &status))
    {
      if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) )
        {
          if (naxis > 2 || naxis == 0)
	    printf("Error: only 1D or 2D images are supported\n");
          else
	    {
	      /* get memory for 1 row */
	      pixels = (double *) malloc(naxes[0] * sizeof(double));
	      
	      if (pixels == NULL) {
                printf("Memory allocation error\n");
                return(1);
	      }
	      
	      if (bitpix > 0) {  /* set the default output format string */
		strcpy(hdformat, " %7d");
		strcpy(format,   " %7.0f");
	      } else {
		strcpy(hdformat, " %15d");
		strcpy(format,   " %15.5f");
	      }
	      
	      //	      printf("\n      ");          /* print column header */
	      //for (ii = 1; ii <= naxes[0]; ii++) printf(hdformat, ii);
	      //	      printf("\n");                /* terminate header line */
	      
	      /* loop over all the rows in the image, top to bottom */
	      //for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
	      // from mlt02.car: reading bottom to top 
	      for (fpixel[1] = 1; fpixel[1] <= naxes[1]; fpixel[1]++) {
		if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL,
				  pixels, NULL, &status) )  /* read row of pixels */
		  break;  /* jump out of loop on error */
		
		//		  printf(" %4d ",fpixel[1]);  /* print row number */
		for (ii = 0; ii < naxes[0]; ii++) {
		  //printf(format, pixels[ii]);   /* print each value  */
		  uu[kk] = pixels[ii];
		  kk++;
		}
		//		  printf("\n");                    /* terminate line */
	      }
	      free(pixels);
	    }
	}
      fits_close_file(fptr, &status);
    } 
  
  if (status) fits_report_error(stderr, status); /* print any error message */
  return(status);
}

//error message
void errore(char *str)
{
  printf(" errore: %s\n\n",str);
  exit(1);
}

//pause
void pausa(void)
{
  printf(" pausa\n");
  getchar();

  return;
}

