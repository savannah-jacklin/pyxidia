#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fitsio.h> 

#define LENMAX 600
#define PI 3.141592653589793238462643

//#define GAIN ??

#define NB_PSF_MAX 5 // maximum oversampling
double *PPSF_0[NB_PSF_MAX*NB_PSF_MAX]; // array of psf realizations (pixel)

// MUST be coherent with psf in input 
int  NB_KERN = 15;
int  NB_PSF = 5;


double PPSF_XC[NB_PSF_MAX*NB_PSF_MAX];
double PPSF_YC[NB_PSF_MAX*NB_PSF_MAX];
int PPSF_IA[NB_PSF_MAX*NB_PSF_MAX];

// pointers to PPSF-related quantities
double **PT_PPSF; 
double *PT_PPSF_XC;
double *PT_PPSF_YC;
int *PT_PPSF_IA;


int NB_KERN,NB_PSF,NB_KSUB;
int NB_KERN2,NB_PSF2;

//iniitializing functions
//open and writing fits files
int ximlist(char *image,double *uu);
void writeimage(int sizex,int sizey,double *pval,char *filename);

//basic c functions
int copia(int nt,double *va,double *vb);
char *charray(int n);
int *iarray(int n);
double *array(int n);
void errore(char *str);
void pausa(void);

//functions specific to this program and cfitsio and numerical recipes
int header(char *file);
int adds(int size,double *pval,double xc,double yc,double flux);
double *loadp(char *file);
int psfo(char *file);
int psfoff(int size_in,int size_out,double *pa,double *pb,int dx,int dy);
int rebin(double *pv_a,double *pv_b,int size,int ngrid);
int twod(int nt,int *ia,double *xc,double *yc);
void sort2i(int n, double *ra,int *rb); // nr
double poidev(double xm); // nr
double gammln(double xx); // nr
int eval_centroid(int size,double *ima,double *xc,double *yc);
double interp2d(int kern,double xc,double yc,double *psf);
double random1(void);
double pfraz(double a);

double TEXP = 60;
double JD;

//arguments list:
//[1] original telescope image
//[2] psf model to be used
//[3] injected star x-position
//[4] injected star y-position
//[5] injected star flux
//[6] size of each side of the trimmed image
//[7] output image name with injected star
int main(int argc, char *argv[])
{
  int ix,jy,x0,y0;
  int size;
  double *uu,*pu;
  double xs,ys,fs;
  char *ima,str[LENMAX], *xpos_s, *ypos_s, *insj_flux_s, *trim_size, *imcrop_x_s, *imcrop_y_s;

  //display usage info if simu is not called correctly
  if(argc!=10) goto Usage; //one greater than then number of argv's

  NB_KSUB = NB_KERN * NB_PSF;

  //psf model
  psfo(argv[2]);

  //original image
  ima = argv[1];
  //injected star x-position
  double xpos = strtod(argv[3], &xpos_s);
  //injected star y-position
  double ypos = strtod(argv[4], &ypos_s);
  //injected star flux
  double injs_flux = strtod(argv[5], &insj_flux_s);
  size = (int)strtod(argv[6], &trim_size);
  //printf("\n%d", size);
  //Print the file name
  //printf("Out File: %s", argv[7]);
  double imcrop_x = strtod(argv[8], &imcrop_x_s);
  double imcrop_y = strtod(argv[9], &imcrop_y_s);

  //size of the trimmed image (square with sides length size)
  //size = 100;
  uu = array(size*size);

  //read header of input image fits file
  header(ima);

  x0 = imcrop_x; //x position of trimming of whole image (bottom left corner)
  y0 = imcrop_y; ///y position of trimming of whole image (bottom left corner)
  //print image trimming parameters to screen
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
  writeimage(size,size,uu,"trimmed_image.fits"); //original trimmed image

  PT_PPSF = PPSF_0; 
  PT_PPSF_XC = PPSF_XC;
  PT_PPSF_YC = PPSF_YC;
  PT_PPSF_IA = PPSF_IA;

  
  xs = xpos; //x-position of input star 60
  ys = ypos; //y-position of input star 80 //total image bounds on the order of 90x90?
  fs = injs_flux; //changing this changes the flux of the input star
  adds(size,uu,xs,ys,fs);
  //writeimage(size,size,uu,"star_added_trimmed_image.fits"); //new image name containing the star
  writeimage(size,size,uu, argv[7]);

  exit(0);
 Usage:
  printf(" Usage: ./simu {image} {psf}\n");
  exit(1);
}


int adds(int size,double *pval,double xc,double yc,double flux)
{
  int ix,jy,xp,yp,shift,jval;
  double xs,ys;
  double *pa,*ps,xm,fl;
  int x2,y2,*pm;
  double *psf; 
  int *masque; 

  // full integer part
  xp = xc;
  yp = yc;

  // fractional part
  // ... here i use pixel coordinates: they are psf identifiers
  // psf that has "embedded" the correct further centroid shift
  xs = pfraz(xc);
  ys = pfraz(yc);

  // ............ shift
  shift = (xp + size*yp); // center (full integer part)
  shift -= NB_KERN2*(1+size); // bottom-left 

  masque = iarray(NB_KERN*NB_KERN);
  psf = array(NB_KERN*NB_KERN);

  //printf(" size = %d / NB_KSUB = %d, NB_QSUB = %d (flux=%g)\n",size,NB_KSUB,NB_QSUB,flux);
  //if(NB_QSUB>1000) pausa();
  //printf(" xc,yc = %g,%g => xp,yp = %d,%d, xs,ys = %g,%g => shift = %d / flux = %g\n",xc,yc,xp,yp,xs,ys,shift,flux);

  pm = masque;
  for(jy=0; jy<NB_KERN; jy++) {
    for(ix=0; ix<NB_KERN; ix++,pm++) {
      *pm = 0;
      x2 = ix + (xp - NB_KERN2); 
      y2 = jy + (yp - NB_KERN2); 
      if( x2<0 || x2>=size || y2<0 || y2>=size) *pm=1;
    }
  }

  interp2d(NB_KERN,xs,ys,psf);

  /*
    pq = array(NB_QSUB*NB_QSUB/NB_GRID/NB_GRID);
    rebin(psf,pq,NB_QSUB/NB_GRID,NB_GRID);
    eval_centroid(NB_QSUB/NB_GRID,pq,&xc,&yc);
    printf(" ===================> xc, yc = %g, %g\n",xc,yc);
    ps = pq;
    for(jy=0; jy<NB_QSUB/NB_GRID; jy++) {
    for(ix=0; ix<NB_QSUB/NB_GRID; ix++,ps++) {
    //    printf(" 1 %3d %3d %g\n",ix,jy,*ps);
    }
    }
    pausa();
  */


  ps = psf;
  pm = masque;
  for(jy=0; jy<NB_KERN; jy++) {
    for(ix=0; ix<NB_KERN; ix++,ps++,pm++) {
      if(*pm) continue;
      jval = ix + jy * size; // moving within image
      pa = pval + shift + jval;
      xm = flux*(*ps);
      //fl = xm; // star-noise free !
      fl = poidev(xm);
      //      if(fl>0) printf(" %2d %2d %g %g %g\n",ix,jy,*ps,xm,fl);
      *pa += fl;
      //*pa += xm; // TEST!!!!!!!!!
      //*pa = xm; // TEST!!!!!!!!!
    }
  }  

  free(masque);
  free(psf);

  return 0;
}



double *loadp(char *file)
{
  int size;
  double *psf;

  // TBD: read the header to get NB_KERN and NB_PSF
  size = NB_KERN * NB_PSF;
  //
  //printf("Size:\t  %d\n",size);
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
    moy[ii] /= size; // mumble, sono evidentemente uguali... ("mumble, I'm obviously the same")
  }
  
  
  for(ii=0; ii<2; ii++) {
    for(jj=0,num=0,den=0; jj<size; jj++) {
      dd = ff[ii][jj]-moy[ii]; 
      if( dd > 0. ) {
	num += dd * (jj+0.5);
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


double interp2d(int kern,double xc,double yc,double *psf)
{
  /*
    rectangular grid
    interpolation in between looks fine
    outside? moving psf around !
    coordinate 0,0: bottom left corner of central pixel
    pixels have linear size equal to 1 
    input: xc,yc in range (-1,2)
    (grid of 3x3 pixels around central one)
    initial xc, yc values set initial dx, dy
    (0,0) for central pixel
    xc < 0 => dx = +1 (offset centroid < 0)
    xc > 1 => dx = -1 (offset centroid > 0)
    yc < 0 => dy = +1 (offset centroid < 0)
    yc > 1 => dy = -1 (offset centroid > 0)
    for xc,yc falling at the pixel edges, the offset 
    for the 4 psf at the corners may change
    from within the fit procedure xc,yc should never
    excced the (-1,2)x(-1,2) patch - 
    in fact, not even getting to the edges of it
    (which would require an offset of 2 pixel for interpolation,
    which is not currenly managed)
  */
  int i,j,ix,jy,index[4];
  double a,b,num,den,t,u,yy[4],*pa,*p0,*p1,*p2,*p3;
  double x0,y0;
  int dx_i,dy_i,ia[4],dx[4],dy[4];
  double *ppsf[4]; // aug20, 2016 [undersampling ... (?)]
  int shift_x,shift_y;
  int i_psfoff;

  // look for initial position
  x0 = xc;
  y0 = yc;
  dx_i = 0;
  dy_i = 0;
  /*
    input (-1,2) is for fitting purpose only
    when simulating stars xc,yc are pfraz(..) => \in (0,1)
  */
  // from xc, yc (-1,2) to x0,y0 \in(0,1)
  // ............ set initial offset
  if(xc < 0) dx_i = 1;
  else if(xc > 1) dx_i = -1;
  x0 = xc + dx_i; // in(0,1)

  if(yc < 0) dy_i = 1;
  else if(yc > 1) dy_i = -1;
  y0 = yc + dy_i; // in(0,1) 


  //  printf(" ======= %g, %g => %d, %d => %g, %g\n",xc,yc,dx_i,dy_i,x0,y0);
  for(ix=0; ix<NB_PSF; ix++) {
    //printf(" ix: %d %g\n",ix,PT_PPSF_XC[PT_PPSF_IA[ix]]); 
    if(x0 < PT_PPSF_XC[PT_PPSF_IA[ix]]) break;
  }

  for(jy=0; jy<NB_PSF; jy++) {
    //printf(" jy: %d %g\n",jy,PT_PPSF_YC[PT_PPSF_IA[jy*NB_PSF]]);
    if(y0 < PT_PPSF_YC[PT_PPSF_IA[jy*NB_PSF]]) break;
  }
  //printf(" => ix, jy = %d, %d\n",ix,jy);
  
  /*
    index[0,1,2,3]: corner (ipsf index) 
    from bottom-left to top-left (counter-clock-wise)
  */
  // ...................................... bottom left
  index[0] = (ix-1) + NB_PSF * (jy-1);  

  // ...................................... bottom right
  index[1] = index[0]+1;

  // ...................................... top right
  index[2] = index[1] + NB_PSF;

  // ...................................... top left
  index[3] = index[0] + NB_PSF;
  
  //printf("\n .......................... interp2dpsf: xc,yc = %g,%g\n",xc,yc);
  //printf(" x0,y0 = %g,%g => ix,jy = %d,%d => index: %2d, %2d, %2d, %2d [offset=%d,%d]\n",x0,y0,ix,jy,index[0],index[1],index[2],index[3],dx_i,dy_i);
  //  pausa();
  /*
    bilinear interpolation for a psf to be centered in xc,yc
    starting from the 4 nearby psf's (at the corners of a
    rectangle around xc, yc) [cfr numerical recipes]
  */
  /*
    recall:
    PSF's are 1-d array running from bottom-left to top-right, row by row
  */

  // function evaluation
  // get the suitable psf around the given position
  pa = psf;
  // ..................... default, to be amended at the edges
  for(i=0; i<4; i++) {
    dx[i] = dx_i;
    dy[i] = dy_i;
    ia[i] = index[i];
  }
  // .................... edge shift
  shift_x = NB_PSF-1;
  shift_y = NB_PSF*(NB_PSF-1); 
  if(ix>0 && ix<NB_PSF && jy>0 && jy<NB_PSF) {
    // right in the middle .... no shift needed
  } else if(ix>0 && ix<NB_PSF && jy==0) { // bottom row
    for(i=0; i<4; i++) {
      if(i==0 || i==1) dy[i] += 1;
    }
    ia[0] = index[3]+shift_y;
    ia[1] = index[2]+shift_y;
  } else if(ix>0 && ix<NB_PSF && jy==NB_PSF) { // top row
    for(i=0; i<4; i++) {
      if(i==2 || i==3) dy[i] -= 1;
    }
    ia[2] = index[1]-shift_y;
    ia[3] = index[0]-shift_y;
  } else if(ix==0 && jy>0 && jy<NB_PSF) { // left column
    for(i=0; i<4; i++) {
      if(i==0 || i==3) dx[i] += 1;
    }
    ia[0] = index[1]+shift_x;
    ia[3] = index[2]+shift_x;
  } else if(ix==NB_PSF && jy>0 && jy<NB_PSF) { // right column
    for(i=0; i<4; i++) {
      if(i==1 || i==2) dx[i] -= 1;
    }
    ia[1] = index[0]-shift_x;
    ia[2] = index[3]-shift_x;
  } else { // the 4 corners
    // psf indexes at every corner
    ia[0] = NB_PSF*NB_PSF-1;
    ia[1] = shift_y;
    ia[2] = 0;
    ia[3] = shift_x;
    if(ix==0 && jy==0) { // bottom-left
      //.......... p0
      dx[0] += 1;
      dy[0] += 1;
      //.......... p1
      dx[1] += 0;
      dy[1] += 1;
      //.......... p2
      dx[2] += 0;
      dy[2] += 0;
      //.......... p3
      dx[3] += 1;
      dy[3] += 0;
    } else if(ix==NB_PSF && jy==0) { // bottom-right
      //.......... p0
      dx[0] += 0;
      dy[0] += 1;
      //.......... p1
      dx[1] += -1;
      dy[1] += 1;
      //.......... p2
      dx[2] += -1;
      dy[2] += 0;
      //.......... p3
      dx[3] += 0;
      dy[3] += 0;
    } else if(ix==NB_PSF && jy==NB_PSF) { // top-right
      //.......... p0
      dx[0] += 0;
      dy[0] += 0;
      //.......... p1
      dx[1] += -1;
      dy[1] += 0;
      //.......... p2
      dx[2] += -1;
      dy[2] += -1;
      //.......... p3
      dx[3] += 0;
      dy[3] += -1;
    } else { // top-left
      //.......... p0
      dx[0] += 1;
      dy[0] += 0;
      //.......... p1
      dx[1] += 0;
      dy[1] += 0;
      //.......... p2
      dx[2] += 0;
      dy[2] += -1;
      //.......... p3
      dx[3] += 1;
      dy[3] += -1;
    }
  }

  //  for(i=0; i<4; i++) printf(" %d: offset = %2d,%2d / ia = %d\n",i,dx[i],dy[i],ia[i]);

  for(i=0; i<4; i++) ppsf[i] = array(kern*kern);
  
  for(i=0,i_psfoff=0; i<4; i++) {
    i_psfoff += psfoff(kern,kern,PT_PPSF[PT_PPSF_IA[ia[i]]],ppsf[i],dx[i],dy[i]);
    //eval_centroid(kern,ppsf[i],&px,&py);
    //printf(" ====> xc, yc = %g, %g\n",px,py);
  }
  if(i_psfoff) return 0;

  p0 = ppsf[0];
  p1 = ppsf[1];
  p2 = ppsf[2];
  p3 = ppsf[3];

  // evaluation t,u 
  // along x
  a = xc;
  b = PT_PPSF_XC[PT_PPSF_IA[ia[0]]];
  //printf(" b[%d] = %g\n",ia[0],b);
  b -= dx_i;
  if(ix==0) b -= 1;
  num = a-b;  
  //printf(" t/num: a,b = %g, %g\n",a,b);

  a = PT_PPSF_XC[PT_PPSF_IA[ia[1]]];
  //printf(" a[%d] = %g\n",ia[1],a);
  a -= dx_i;
  if(ix==NB_PSF) a += 1;
  den = a-b;
  //printf(" t/den: a,b = %g, %g\n",a,b);
  
  t = num/den;
  //printf(" t = %g\n",t);
  
  // along y
  a = yc;
  b = PT_PPSF_YC[PT_PPSF_IA[ia[0]]];
  //printf(" b[%d] = %g\n",ia[0],b);
  b -= dy_i;
  if(jy==0)  b -= 1;
  num = a-b;
  //printf(" u/num: a,b = %g, %g\n",a,b);

  a = PT_PPSF_YC[PT_PPSF_IA[ia[3]]];
  //printf(" a[%d] = %g\n",ia[3],a);
  a -= dy_i;
  if(jy==NB_PSF) a += 1;
  den = a-b;
  //printf(" u/den: a,b = %g, %g\n",a,b);
  
  u = num/den;
  //printf(" u = %g\n",u);


  for(j=0; j<kern; j++) {
    for(i=0; i<kern; i++, pa++,p0++,p1++,p2++,p3++) {
      // values of the 4 psfs 
      yy[0] = *p0;
      yy[1] = *p1;
      yy[2] = *p2;
      yy[3] = *p3;
      // eval psf at xc, yc
      *pa = (1-t)*(1-u)*yy[0] + t*(1-u)*yy[1] + t*u*yy[2] + (1-t)*u*yy[3];
      //      printf(" %2d %2d %g\n",i,j,*pa);
    }
  }

  //  eval_centroid(kern,psf,&px,&py);
  //printf(" ==========> xc, yc = %g, %g\n",px,py);
  for(i=0; i<4; i++) free(ppsf[i]);
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
      // ix, jy: indices pb 
      // kx, ky: indices pa 
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
  int ix,jy,kx,ky,ival,kval,ngrid2,ksub,shift,pval;
  double *pa,*pb,fsum,fsub;

  ksub = size * ngrid;
  ngrid2 = (ngrid-1)/2;

  // loop on pixels starting from the center
  // ... input matrices
  pa = pv_a; 
  // ... matrix output
  pb = pv_b;
  for(jy=0; jy<size; jy++) {
    for(ix=0; ix<size; ix++) {
      ival = ix + jy * size;
      // loop on the sub-pixel grid
      shift = (ngrid2+ix*ngrid) + (ngrid2+jy*ngrid)*ksub; // pixel center
      shift -= (ngrid2*(1+ksub)); // bottom-left corner 
      fsum = 0; // sums sub-pixel values
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
  int i,j,ix,jy,dx,dy,index;
  double *psf,*pa,*pb,*pc,*subpsf;
  double xc,yc;

  //printf("\nPSF Model File: %s\n",file);

  psf = loadp(file); 

  //printf("NB_PSF:   %d\n",NB_PSF);
  //printf("NB_KERN:  %d\n",NB_KERN);

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

      // psf should be centered all in the central pixel
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
    //print fits file parameters
    //printf(" %2d %2d %7.5f %7.5f  %7.5f %7.5f  %7.5f %7.5f\n",i,PPSF_IA[i],PPSF_XC[i],PPSF_YC[i],PPSF_XC[j],PPSF_YC[j],xc,yc);
  }

  return 0;
}

//sort a two dimensional array
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

//sort a two dimensional array
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

//create a pointer to an array of characters
char *charray(int n)
{
  char *a;
  a=(char *) calloc(1,(size_t)(n*sizeof(char)));
  return a;
}

//create a pointer to an array of ints
int *iarray(int n)
{
  int *a;
  a=(int *) calloc(1,(size_t)(n*sizeof(int)));
  return a;
}

//create a pointer to an array of doubles
double *array(int n)
{
  double *a;
  a=(double *) calloc(1,(size_t)(n*sizeof(double)));
  return a;
}



/* 
   ----- numerical recipes
   return as a floating-point number an integer value that is
   a random deviate drawn from a Poisson distribution of mean xm
*/
double poidev(double xm)
{
  static double sq,alxm,g,oldm=(-1.0);
  double em,t,y;

  if( xm < 12.0) {
    if (xm != oldm) {
      oldm = xm;
      g = exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      ++em;
      t *= random1();
    } while(t>g);
  } else {
    if (xm != oldm) {
      oldm = xm;
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g = xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y = tan(PI*random1()); // original numeric; note because tangent, error of infinities will occur at 0.5 radians
	em = sq*y+xm;
      } while (em < 0.0);
      em = floor(em);
      t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while(random1() > t);
  }
  return em;
}

// return the value of ln(Gamma(xx)) for xx>0
double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6] = {76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for(j=0; j<=5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}



#define NT_KEY_MAX 10
char **PT_KEYWORD;
char *KEYWORD[] = {"NAXIS1","NAXIS2","MJD-OBS","FIN_"};
int nkey(void);
int listhead2(char *image, char *header[],int *flag);
int gethead(char *card, char *head);
int header(char *file)
{
  int nt;
  int i,*flag;
  char *hheader[NT_KEY_MAX];

  PT_KEYWORD = KEYWORD;

  nt = nkey();
  if(nt>NT_KEY_MAX) errore(" header: NT_KEY_MAX");

  flag = iarray(nt);
  for(i=0; i<nt; i++) hheader[i] = charray(80);

  listhead2(file,hheader,flag);

  for(i=0; i<nt; i++) {
    //printf(" %s %s\n",KEYWORD[i],hheader[i]);
    //printf("\n%s %s",PT_KEYWORD[i],hheader[i]);
    if(i==2) JD = atof(hheader[i]);
  }

  free(flag);
  for(i=0; i<nt; i++) free(hheader[i]);

  return 0;
}

int nkey(void)
{
  int i = 0;
  //while(strcmp(KEYWORD[i],"FIN_")) i++;
  while(strcmp(PT_KEYWORD[i],"FIN_")) i++;
  return i;
}

//grab fits file header
int gethead(char *card, char *head)
{
  int i,j,n,flag,istop;

  n = strlen(card);
  //printf("card : %s (n=%d) \n",card,n);
  i=0;
  j=0;
  flag=0;
  istop = 1;
  while(istop) {
    if(card[i]=='=') flag++;
    if(flag && card[i] != '=' && card[i] != ' ' && card[i] != '\0' && card[i] != '/') {
      //printf("card[%d] = %c\n",i,card[i]);
      head[j] = card[i];
      j++;
    }
    if(card[i]=='/' || i==n) istop = 0;
    i++;
  }
  //  printf("head %s [%d]\n",head,strlen(head));
  return 0;
}



//list headers contained in a fits file
// ........................................ cfitsio
int listhead2(char *image, char *header[],int *flag)
{
  int count;
  fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
  char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int single = 0, hdupos, nkeys, ii;

  if (!fits_open_file(&fptr, image, READONLY, &status))
    {
      fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
      
      /* List only a single header if a specific extension was given */ 
      if (hdupos != 1 || strchr(image, '[')) single = 1;
      
      for (; !status; hdupos++)  /* Main loop through each extension */
	{
	  fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of KEYWORD */
	  
	  //printf("Header listing for HDU #%d:\n", hdupos);
	  
	  for (ii = 1; ii <= nkeys; ii++) { /* Read and print each KEYWORD */
	    
	    if (fits_read_record(fptr, ii, card, &status))break;	    
	    count = 0;
	    //while(strcmp(KEYWORD[count],"FIN_")) {
	    //if(strncmp(card,KEYWORD[count],strlen(KEYWORD[count])) == 0) {
	    //printf("===== %s %lu\n",KEYWORD[count],strlen(KEYWORD[count]));
	    while(strcmp(PT_KEYWORD[count],"FIN_")) {
	      if(strncmp(card,PT_KEYWORD[count],strlen(PT_KEYWORD[count]))==0) {
		//		printf("===== %s %lu\n",PT_KEYWORD[count],strlen(PT_KEYWORD[count]));
		gethead(card,header[count]);
		flag[count]++;
		//printf("%2d : flag[%d] = %d\n",ii,count,flag[count]);
	      }
	      count++;
	      
	    }
	  }
	  //printf("END\n\n");  /* terminate listing with END */
	  
	  if (single) break;  /* quit if only listing a single header */
	  
	  fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
	}
      
      if (status == END_OF_FILE)  status = 0; /* Reset after normal error */
      
      fits_close_file(fptr, &status);
    }
  
  if (status) fits_report_error(stderr, status); /* print any error message */
  //  pausa();
  return status;
}

//write the new image
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




//gets fits file parameters and header information
int ximlist(char *image,double *uu)
{
  fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int bitpix, naxis, ii;
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

//random number generator from numerical perscriptions
double random1(void)
{
  return drand48(); // most likely NOT a good choice ... (eg 1005.4117) maybe not a good enough random number generator
}

//function to copy the variable
int copia(int nt,double *va,double *vb)
{
  int i;
  double *pa,*pb;

  pa = va;
  pb = vb;
  for(i=0; i<nt; i++,pa++,pb++) *pb = *pa;

  return 0;
}

//returns the fractional portion of a float
double pfraz(double a)
{
  return a - (int)a;
}

//error function
void errore(char *str)
{
  printf(" errore: %s\n\n",str);
  exit(1);
}

//pause function
void pausa(void)
{
  printf(" pausa\n");
  getchar();

  return;
}

