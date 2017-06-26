#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fitsio.h> 

//The point of this code is to read in an image fits file and a psf model
//file and apply the psf model to the image using oversampled bins per pixel
//to fix the centroid on the appropriate point and then move forward from there
//see description in notebook for a more throrough explanation with pictures
//related to event injection
#define LENMAX 600


#define NB_PSF_MAX 5 // maximum oversampling
double *PPSF_0[NB_PSF_MAX*NB_PSF_MAX]; // array of psf realizations (pixel)

double PPSF_XC[NB_PSF_MAX*NB_PSF_MAX];
double PPSF_YC[NB_PSF_MAX*NB_PSF_MAX];
int PPSF_IA[NB_PSF_MAX*NB_PSF_MAX];

int NB_KERN,NB_PSF,NB_KSUB;
int NB_KERN2,NB_PSF2;

int ximlist(char *image,double *uu);
void writeimage(int sizex,int sizey,double *pval,char *filename);


int copia(int nt,double *va,double *vb);
int *iarray(int n);
double *array(int n);

void errore(char *str);
void pausa(void);





double *loadp(char *file);
int psfo(char *file);
int psfoff(int size_in,int size_out,double *pa,double *pb,int dx,int dy);
int rebin(double *pv_a,double *pv_b,int size,int ngrid);
int twod(int nt,int *ia,double *xc,double *yc);
void sort2i(int n, double *ra,int *rb);
int eval_centroid(int size,double *ima,double *xc,double *yc);
double pfraz(double a);

double TEXP = 60;
double JD = 7000;

int main(int argc, char *argv[])
{
  int a,*pa,ix,jy,x0,y0;
  int size;
  double *uu,*pu;
  char *ima,str[LENMAX];


  if(argc!=3) goto Usage;




  NB_KERN = 15;
  NB_PSF = 5;
  NB_KSUB = NB_KERN * NB_PSF;

  psfo(argv[2]);
  ima = argv[1];

  printf(" bla\n");
  
  //little sample code about pointers in C
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
  printf(" Usage: ./simu {image} {psf}\n");
  exit(1);

}


double *loadp(char *file)
{
  int size;
  double *psf;

  // TBD: read the header to get NB_KERN and NB_PSF
  size = NB_KERN * NB_PSF;
  printf(" size = %d\n",size);
  //pausa();
  psf = array(size*size);
  ximlist(file,psf);

  return psf;
}



int eval_centroid(int size,double *ima,double *xc,double *yc)
{
  int ii,jj,ix,jy;
  double *pa,val,pp[2],moy[2],num,den,dd;
  double ff[2][size];
  
  for(ii=0; ii<2; ii++) {
    moy[ii] = 0;
    for(jj=0; jj<size; jj++) ff[ii][jj] = 0;
  }

  pa = ima;
  for(jy=0; jy<size; jy++) {
    for(ix=0; ix<size; ix++) {
      val = *(pa + ix + jy*size);
      ff[0][ix] += val;
      ff[1][jy] += val;
    }
  }

  for(ii=0; ii<2; ii++) {
    for(jj=0; jj<size; jj++) {
      moy[ii] += ff[ii][jj];
    }
    moy[ii] /= size; // mumble, sono evidentemente uguali...
  }
  
  
  for(ii=0; ii<2; ii++) {
    for(jj=0,num=0,den=0; jj<size; jj++) {
      dd = ff[ii][jj]-moy[ii]; 
      if( dd > 0. ) {
	// .............. july 08, 2016
	// "half pixel": does it really matter?
	// likely not: values ppsf_xc,yc enter in interp2d(...)
	// so that eval_centroid and fit stay self-coherent
	// per come ho in mente la geometria, tuttavia, 
	// should actually use half pixel (-> cdr test with gauss2d()!)
	// ovvero: valore psf e' quello al centro del pixel
	// di coordinate, coerentemente, centro del pixel
	// it does matter, on the other hand, in self-built psf
	//num += dd * jj; // old
	num += dd * (jj+0.5); // new: july 08, 2016 - half pixel 
	den += dd;
      }
    }
    pp[ii] = num/den;
  }
  //  printf(" ====> %g, %g\n",pp[0],pp[1]);

  *xc = pp[0];
  *yc = pp[1];
  
  return 0;
}


int psfoff(int size_in,int size_out,double *pa,double *pb,int dx,int dy)
{
  int ix,jy,kx,ky,ia;
  double *qa,*qb;
  int s2,t2,shift,delta;

  s2 = (size_in-1)/2;
  t2 = (size_out-1)/2;
  delta = s2 - t2;
  shift = s2 + size_in*s2; // center
  shift -= t2*(1+size_in); // bottom-left 

  qa = pa; // input 
  qb = pb; // output
  // jy, ix: indici pa
  for(jy=0; jy<size_out; jy++) {
    for(ix=0; ix<size_out; ix++, qb++) {
      // ix, jy: indici pb 
      // kx, ky: indici pa 
      kx = ix + delta + dx;
      //if(kx==-1 || kx==size) kx = ix; // |offset| max = 1
      if(kx<0) kx = 0;
      if(kx>=size_in) kx = size_in-1;
      ky = jy + delta + dy;
      //if(ky==-1 || ky==size) ky = jy; // |offset| max = 1
      if(ky<0) ky = 0;
      if(ky>=size_in) ky = size_in-1;
      //ib = ix+jy*size; // may as well do this and *(qb+ib) = ....
      ia = kx+ky*size_in;
      *qb = *(qa+ia);
      //printf(" ix, jy = %2d, %2d => kx, ky = %2d, %2d => ia = %3d (ib=%3d)\n",ix,jy,kx,ky,ia,ix+jy*size);
      //pausa();
    } 
  }

  return 0;
}


int rebin(double *pv_a,double *pv_b,int size,int ngrid)
{
  int ix,jy,kx,ky,ival,kval,ngrid2,ksub,shift,pval,iw;
  double *pa,*pb,fsum,fsub,r,*pw;

  ksub = size * ngrid;
  ngrid2 = (ngrid-1)/2;

  // loop sui pixel - partendo dal centro
  // ... matricione input
  pa = pv_a; 
  // ... matrice output
  pb = pv_b;
  for(jy=0; jy<size; jy++) {
    for(ix=0; ix<size; ix++) {
      ival = ix + jy * size;
      // loop sulla griglia sub-pixel
      shift = (ngrid2+ix*ngrid) + (ngrid2+jy*ngrid)*ksub; // pixel center
      shift -= (ngrid2*(1+ksub)); // bottom-left corner 
      fsum = 0; // somma valori sub-pixel 
      for(ky=0; ky<ngrid; ky++) {
	for(kx=0; kx<ngrid; kx++) {
	  kval = shift + kx + ky * ksub; // moving around sub-pixel
	  pval = kx + ky * ngrid;
	  fsub = *(pa+kval);
	  fsum += fsub;
	}
      }
      *(pb+ival) = fsum;
    }
  }
 
  return 0;
}


int psfo(char *file)
{
  int i,j,ix,jy,dx,dy,index,size;
  double *psf,*pa,*pb,*pc,*subpsf;
  double xc,yc,xp,yp;
  int ksub2;

  printf("\n ....... psf %s\n",file);

  psf = loadp(file); 

  printf(" NB_PSF = %d\n",NB_PSF);
  printf(" NB_KERN = %d\n",NB_KERN);

  subpsf = array(NB_KSUB*NB_KSUB); // subpixel
  pb = subpsf;

  NB_PSF2 = (NB_PSF-1)/2;
  NB_KERN2 = (NB_KERN-1)/2; 
  for(jy=0; jy<NB_PSF; jy++) {
    for(ix=0; ix<NB_PSF; ix++) {
      index = ix + jy * NB_PSF;
      dx = -NB_PSF2 + ix;
      dy = -NB_PSF2 + jy;
      PPSF_0[index] = array(NB_KERN*NB_KERN); // pixel-based
      pa = PPSF_0[index];

      psfoff(NB_KSUB,NB_KSUB,psf,pb,dx,dy);
      rebin(pb,pa,NB_KERN,NB_PSF); 

      // psf devono essere centrate tutte nel pixel centrale
      // if not, move again to recenter
      // (at most, 1 pixel - although in principle I may do more....)
      eval_centroid(NB_KERN,pa,&xc,&yc);
      //      printf(" a: xc,yc = %g, %g\n",xc,yc);
      dx = 0;
      dy = 0;
      if(xc<NB_KERN2-2 || xc>NB_KERN2+2) errore(" load_psf: xc");
      if(yc<NB_KERN2-2 || yc>NB_KERN2+2) errore(" load_psf: yc");
      if(xc<NB_KERN2) dx = -1;
      else if(xc>NB_KERN2+1) dx = 1;
      if(yc<NB_KERN2) dy = -1;
      else if(yc>NB_KERN2+1) dy = 1;
      if(dx != 0 || dy != 0) {
	pc = array(NB_KERN*NB_KERN);
	psfoff(NB_KERN,NB_KERN,pa,pc,dx,dy); // pa -> pc
	copia(NB_KERN*NB_KERN,pc,pa); // pa = pc
	free(pc);
	eval_centroid(NB_KERN,pa,&xc,&yc);
	//printf(" =====> b: %2d %2d %7.4f %7.4f  %2d %2d\n",ix,jy,xc,yc,dx,dy);
	if(xc<NB_KERN2 || xc>NB_KERN2+1) errore(" load_psf: b/xc");
	if(yc<NB_KERN2 || yc>NB_KERN2+1) errore(" load_psf: b/yc");
      }
      //      pausa();
      PPSF_XC[index] = pfraz(xc);
      PPSF_YC[index] = pfraz(yc);
    }
  }
  free(psf);
  free(subpsf);

  twod(NB_PSF,PPSF_IA,PPSF_XC,PPSF_YC);   
  for(i=0; i<NB_PSF*NB_PSF; i++) {
    j = PPSF_IA[i];
    pa = PPSF_0[j];
    eval_centroid(NB_KERN,pa,&xc,&yc);
    printf(" %2d %2d %7.5f %7.5f  %7.5f %7.5f  %7.5f %7.5f\n",i,PPSF_IA[i],PPSF_XC[i],PPSF_YC[i],PPSF_XC[j],PPSF_YC[j],xc,yc);
  }

  return 0;
}

int twod(int nt,int *ia,double *xc,double *yc)
{
  int i,j,nt2;
  int *ib;
  double *va,*vb;

  nt2 = nt*nt;
  va = array(nt2);
  
  for(i=0; i<nt2; i++) {
    ia[i] = i;
    va[i] = yc[i];
  }
  sort2i(nt2,va-1,ia-1); 
  
  ib = iarray(nt);
  vb = array(nt);

  for(j=0; j<nt; j++) {
    for(i=0; i<nt; i++) {
      ib[i] = ia[i+j*nt];
      vb[i] = xc[ia[i+j*nt]];
    }
    sort2i(nt,vb-1,ib-1); 
    for(i=0; i<nt; i++) ia[i+j*nt] = ib[i];
  }
  
  free(va);
  free(vb);
  free(ib);

  return 0;
}


// ............................................... numerical recipes
void sort2i(int n, double *ra,int *rb)
{
  int l,j,ir,i;
  double rra;
  int rrb;
  int xl;
  
  l=(n >> 1)+1;
  ir=n;

  for(;;) {
    if (l>1) {
      //rra=ra[--l];
      //      rrb=rb[--l]; // that's would be wrong!
      xl = --l;
      rra=ra[xl];
      rrb=rb[xl];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      rrb=rb[ir];
      rb[ir]=rb[1];
      if (--ir == 1) {
	ra[1]=rra;
	rb[1]=rrb;
	return;
      }
    }
    i=l;
    j=l << 1;
    while (j<=ir) {
      if (j<ir && ra[j]<ra[j+1]) ++j;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	rb[i]=rb[j];
	j += (i=j);
      } else j = ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
  }
}



int *iarray(int n)
{
  int *a;
  a=(int *) calloc(1,(size_t)(n*sizeof(int)));
  return a;
}


double *array(int n)
{
  double *a;
  a=(double *) calloc(1,(size_t)(n*sizeof(double)));
  return a;
}


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


int copia(int nt,double *va,double *vb)
{
  int i;
  double *pa,*pb;

  pa = va;
  pb = vb;
  for(i=0; i<nt; i++,pa++,pb++) *pb = *pa;

  return 0;
}


double pfraz(double a)
{
  return a - (int)a;
}


void errore(char *str)
{
  printf(" errore: %s\n\n",str);
  exit(1);
}

void pausa(void)
{
  printf(" pausa\n");
  getchar();

  return;
}

