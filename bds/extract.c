#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <tcl.h>
#include <tk.h>

typedef
struct{double x,y,y2;}
spline_tab_entry;

typedef struct{int n,xscal,yscal;spline_tab_entry *tab;} SPLINE;

#define SPLINE_NMAX 1000

/* can be avoided for gcc */

void spline_init(double x[],int xscal,double y[],int yscal,int n,
		 SPLINE *spline)
{
    int i;
    double u[SPLINE_NMAX],sig,p;

    if (n>SPLINE_NMAX) {
	fprintf(stderr,"Error: to many values in spline_init\n");
	exit(1);
    }
    spline->n=n;
    spline->xscal=xscal;
    spline->yscal=yscal;
    spline->tab=(spline_tab_entry*)malloc(sizeof(spline_tab_entry)*n);
    spline->tab[0].y2=0.0;
    for (i=0;i<n;i++){
	switch (xscal){
	case 0:
	    spline->tab[i].x=x[i];
	    break;
	case 1:
	    spline->tab[i].x=log(x[i]);
	    break;
	}
	switch (yscal){
	case 0:
	    spline->tab[i].y=y[i];
	    break;
	case 1:
	    spline->tab[i].y=log(y[i]);
	    break;
	}	    
    }
    u[0]=0.0;
    for (i=1;i<n-1;i++){
	sig=(spline->tab[i].x-spline->tab[i-1].x)
	    /(spline->tab[i+1].x-spline->tab[i-1].x);
	p=1.0/(sig*spline->tab[i-1].y2+2.0);
	spline->tab[i].y2=(sig-1.0)*p;
	u[i]=(6.0*((spline->tab[i+1].y-spline->tab[i].y)
		   /(spline->tab[i+1].x-spline->tab[i].x)
		   -(spline->tab[i].y-spline->tab[i-1].y)
		   /(spline->tab[i].x-spline->tab[i-1].x))
	      /(spline->tab[i+1].x-spline->tab[i-1].x)-sig*u[i-1])*p;
    }
    spline->tab[n-1].y2=0.0;
    for (i=n-2;i>=0;i--){
	spline->tab[i].y2=spline->tab[i].y2*spline->tab[i+1].y2+u[i];
    }
}

double spline_int(SPLINE *spline,double x)
{
    int kmin,kmax,kpoint;
    double a,b,w;

    kmin=0;
    kmax=spline->n-1;
    switch(spline->xscal){
    case 0:	
	break;
    case 1:
	x=log(x);
	break;
    }
    if (x>spline->tab[kmax].x){
	if (spline->yscal){
	    return exp(spline->tab[kmax].y);
	}
	else{
	    return spline->tab[kmax].y;
	}
    }
    if (x<spline->tab[0].x){
      if (spline->yscal) {
	return exp(spline->tab[0].y);
      }
      else {
	return spline->tab[0].y;
      }
    }
    while (kmax-kmin>1){
	kpoint=(kmax+kmin)/2;
	if (spline->tab[kpoint].x>x){
	    kmax=kpoint;
	}
	else{
	    kmin=kpoint;
	}
    }
    w=spline->tab[kmax].x-spline->tab[kmin].x;
    a=(spline->tab[kmax].x-x)/w;
    b=(x-spline->tab[kmin].x)/w;
    x=a*spline->tab[kmin].y+b*spline->tab[kmax].y+
	(a*(a*a-1.0)*spline->tab[kmin].y2
	 +b*(b*b-1.0)*spline->tab[kmax].y2)*w*w/6.0;
    switch (spline->yscal){
    case 0:
	return x;
    case 1:
	return exp(x);
    }
    return x;
}


typedef struct{int n,xscal,yscal,nval;double *x,*y,*y2;} MSPLINE;

void mspline_init(double x[],int xscal,double y[],int yscal,int n,int nval,
	     MSPLINE *spline)
{
    int i,j;
    double u[SPLINE_NMAX],sig,p;

    if (n>SPLINE_NMAX) {
	fprintf(stderr,"Error: to many values in mspline_init\n");
	exit(1);
    }
    spline->n=n;
    spline->nval=nval;
    spline->xscal=xscal;
    spline->yscal=yscal;
    spline->x=(double*)malloc(sizeof(double)*n);
    spline->y=(double*)malloc(sizeof(double)*n*nval);
    spline->y2=(double*)malloc(sizeof(double)*n*nval);
    for (j=0;j<nval;j++){
	(spline->y2)[j]=0.0;
    }
    for (i=0;i<n;i++){
	switch (xscal){
	case 0:
	    (spline->x)[i]=x[i];
	    break;
	case 1:
	    (spline->x)[i]=log(x[i]);
	    break;
	}
	for (j=0;j<nval;j++){
	    switch (yscal){
	    case 0:
		(spline->y)[i*nval+j]=y[i*nval+j];
		break;
	    case 1:
		(spline->y)[i*nval+j]=log(y[i*nval+j]);
		break;
	    }
	}	    
    }
    for (j=0;j<nval;j++){
	u[0]=0.0;
	for (i=1;i<n-1;i++){
	    sig=((spline->x)[i]-(spline->x)[i-1])
		/((spline->x)[i+1]-(spline->x)[i-1]);
	    p=1.0/(sig*(spline->y2)[(i-1)*nval+j]+2.0);
	    (spline->y2)[i*nval+j]=(sig-1.0)*p;
	    u[i]=(6.0*(((spline->y)[(i+1)*nval+j]-(spline->y)[i*nval+j])
		       /((spline->x)[i+1]-(spline->x)[i])
		       -((spline->y)[i*nval+j]-(spline->y)[(i-1)*nval+j])
		       /((spline->x)[i]-(spline->x)[i-1]))
		  /((spline->x)[i+1]-(spline->x)[i-1])-sig*u[i-1])*p;
	}
	(spline->y2)[(n-1)*nval+j]=0.0;
	for (i=n-2;i>=0;i--){
	    (spline->y2)[i*nval+j]=
		(spline->y2)[i*nval+j]*(spline->y2)[(i+1)*nval+j]+u[i];
	}
    }
}

void mspline_int(MSPLINE *spline,double x,double y[])
{
    int kmin,kmax,kpoint,j,nval;
    double a,b,w,tmp;

    nval=spline->nval;
    kmin=0;
    kmax=spline->n-1;
    switch(spline->xscal){
    case 0:	
	break;
    case 1:
	x=log(x);
	break;
    }
    if (x>(spline->x)[kmax]){
	if (spline->yscal){
	    for (j=0;j<nval;j++) {y[j]=exp((spline->y)[kmax*nval+j]);}
	}
	else{
	    for (j=0;j<nval;j++) {y[j]=(spline->y)[kmax*nval+j];}
	}
	return;
    }
    if (x<(spline->x)[0]){
	if (spline->yscal){
	    for (j=0;j<nval;j++) {y[j]=exp((spline->y)[j]);}
	}
	else{
	    for (j=0;j<nval;j++) {y[j]=(spline->y)[j];}
	}
	return;
	for (j=0;j<nval;j++) {y[j]=0.0;}
	return;
    }
    while (kmax-kmin>1){
	kpoint=(kmax+kmin)/2;
	if ((spline->x)[kpoint]>x){
	    kmax=kpoint;
	}
	else{
	    kmin=kpoint;
	}
    }
    w=(spline->x)[kmax]-spline->x[kmin];
    a=((spline->x)[kmax]-x)/w;
    b=(x-(spline->x)[kmin])/w;
    for (j=0;j<nval;j++){
	tmp=a*(spline->y)[kmin*nval+j]+b*(spline->y)[kmax*nval+j]+
	    (a*(a*a-1.0)*(spline->y2)[kmin*nval+j]
	     +b*(b*b-1.0)*(spline->y2)[kmax*nval+j])*w*w/6.0;
	switch (spline->yscal){
	case 0:
	    y[j]=tmp;
	    break;
	case 1:
	    y[j]=exp(tmp);
	    break;
	}
    }
}

void spline_delete(SPLINE *spline)
{
    free(spline->tab);
    spline->tab=NULL;
}

void mspline_delete(MSPLINE *spline)
{
    free(spline->x);
    free(spline->y);
    free(spline->y2);
    spline->x=NULL;
    spline->y=NULL;
    spline->y2=NULL;
}

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef max
#define max(a,b) (((a)<(b))?(b):(a))
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif

#define RNDM_EPS 6e-8

static struct
{
  float u[97],c,cd,cm;
  int i,j;
} rndm5_store;

void rndmst5(int na1,int na2,int na3, int nb1)
{
  int i,j,nat;
  float s,t;
  rndm5_store.i=96;
  rndm5_store.j=32;
  for (i=0;i<97;i++)
    {
      s=0.0;
      t=0.5;
      for (j=0;j<24;j++)
	{
	  nat=(((na1*na2) % 179)*na3) % 179;
	  na1=na2;
	  na2=na3;
	  na3=nat;
	  nb1=(53*nb1+1) % 169;
	  if ((nb1*nat) % 64 >= 32)
	    {
	      s+=t;
	    }
	  t*=0.5;
	}
      rndm5_store.u[i]=s;
    }
  rndm5_store.c=    362436.0/16777216.0;
  rndm5_store.cd=  7654321.0/16777216.0;
  rndm5_store.cm= 16777213.0/16777216.0;
}

float rndm5()
{
  float temp;

  temp=rndm5_store.u[rndm5_store.i]-rndm5_store.u[rndm5_store.j];
  if (temp<0.0){
    temp+=1.0;
  }
  rndm5_store.u[rndm5_store.i]=temp;
  if (--rndm5_store.i<0) rndm5_store.i=96;
  if (--rndm5_store.j<0) rndm5_store.j=96;
  rndm5_store.c-=rndm5_store.cd;
  if (rndm5_store.c<0.0){
    rndm5_store.c+=rndm5_store.cm;
  }
  temp-=rndm5_store.c;
  if (temp<0.0){
    return temp+1.0;
  }
  else{
    return temp;
  }
}

static struct{
  int iset;
  float v1,v2;
} gasdev_data;

float gasdev()
{
  float r;
  if (gasdev_data.iset==0){
    for (;;){
      gasdev_data.v1=2.0*rndm5()-1.0;
      gasdev_data.v2=2.0*rndm5()-1.0;
      r=gasdev_data.v1*gasdev_data.v1+gasdev_data.v2*gasdev_data.v2;
      if ((r<=1.0) && (r!=0)){
	break;
      }
    }
    gasdev_data.iset=1;
    r=sqrt(-2.0*log((double)r)/r);
    gasdev_data.v1*=r;
    gasdev_data.v2*=r;
    return gasdev_data.v1;
  }
  else{
    gasdev_data.iset=0;
    return gasdev_data.v2;
  }
}

#undef RNDM_EPS

typedef struct{
  double energy,wgt,y,yp;
  double x,xp;
} PARTICLE;

typedef struct{
  double r11,r12,r21,r22;
} R_MATRIX;

class BEAM {
 public:
  int slices_per_bunch,macroparticles,bunches;
  int slices,n_field,which_field;
  double factor,transv_factor;
  PARTICLE *particle;
  R_MATRIX *sigma,*sigma_xx,*sigma_xy;
  double *z_position;
};

BEAM *bunch_make(int n_bunch,int n_slice,int n_macro,int n_field)
{
  BEAM *bunch;
  bunch=(BEAM*)malloc(sizeof(BEAM));
  bunch->bunches=n_bunch;
  bunch->slices_per_bunch=n_slice;
  bunch->macroparticles=n_macro;
  n_slice*=n_macro*n_bunch;
  bunch->slices=n_slice;
  bunch->particle=(PARTICLE*)malloc(sizeof(PARTICLE)*n_slice);
  bunch->sigma=(R_MATRIX*)malloc(sizeof(R_MATRIX)*n_slice);
  bunch->sigma_xx=(R_MATRIX*)malloc(sizeof(R_MATRIX)*n_slice);
  bunch->sigma_xy=(R_MATRIX*)malloc(sizeof(R_MATRIX)*n_slice);
  bunch->z_position=(double*)malloc(sizeof(double)*n_slice);
  return bunch;
}

BEAM *read_beam_file(char *name,int ibunch0)
{
  FILE *file;
  char buffer[1000],*point;
  int nbunch,nslice,nmacro,ibunch,islice,imacro,k=0;
  BEAM *beam;

  file=fopen(name,"r");
  point=fgets(buffer,1000,file);
  nbunch=strtol(point,&point,10);
  nslice=strtol(point,&point,10);
  nmacro=strtol(point,&point,10);
  if (nmacro==0) {
    nmacro=1;
  }
  for (ibunch=0;ibunch<ibunch0;ibunch++){
    for (islice=0;islice<nslice;islice++){
      for (imacro=0;imacro<nmacro;imacro++){
	point=fgets(buffer,1000,file);
      }
    }
    point=fgets(buffer,1000,file);
  }
  beam=bunch_make(1,nslice,nmacro,0);
  for (islice=0;islice<nslice;islice++){
    for (imacro=0;imacro<nmacro;imacro++){
      point=fgets(buffer,1000,file);
      beam->z_position[islice]=strtod(point,&point);
      beam->particle[k].wgt=strtod(point,&point);
      beam->particle[k].energy=strtod(point,&point);
      beam->particle[k].x=strtod(point,&point);
      beam->particle[k].xp=strtod(point,&point);
      beam->particle[k].y=strtod(point,&point);
      beam->particle[k].yp=strtod(point,&point);
      beam->sigma_xx[k].r11=strtod(point,&point);
      beam->sigma_xx[k].r12=strtod(point,&point);
      beam->sigma_xx[k].r22=strtod(point,&point);
      beam->sigma_xx[k].r21=beam->sigma_xx[k].r12;
      beam->sigma[k].r11=strtod(point,&point);
      beam->sigma[k].r12=strtod(point,&point);
      beam->sigma[k].r22=strtod(point,&point);
      beam->sigma[k].r21=beam->sigma[k].r12;
      beam->sigma_xy[k].r11=strtod(point,&point);
      beam->sigma_xy[k].r12=strtod(point,&point);
      beam->sigma_xy[k].r21=strtod(point,&point);
      beam->sigma_xy[k].r22=strtod(point,&point);
      k++;
    }
  }
  fclose(file);
  return beam;
}

void
data_MAD_track(BEAM *beam,char *name,int n_particles,int axis,
		    int bunch_number,double energy0,int which)
{
  FILE *file;
  int i,j,n,nb,np=0,nm;
  int slices;
  double z_rndm, x_rndm, xp_rndm, y_rndm, yp_rndm, e_rndm;
  double z_pos, energy_permille;
  double sigma_z;
  double z_cut=3.5,e_cut=3.5;
  int n_off,n_count;
  double *xvalue,*yvalue,*ytransv;
  SPLINE *s_e,*s_x,*s_xp,*s_y,*s_yp,*s_x_x,*s_x_xp,*s_xp_xp,*s_y_y,*s_y_yp,
      *s_yp_yp,*s_test;
  MSPLINE *s_transv;
  slices = beam->slices;

  nb = beam->bunches;
  n = beam->slices_per_bunch;
  nm=beam->macroparticles;
  if (bunch_number<0) {
    n_off=bunch_number*n*nm;
    n_count=n*nm;
  }
  else {
    n_off=0;
    n_count=nb*n*nm;
  }

  if (name){
    file=fopen(name,"w");
  }
  else{
    file=stdout;
  }

  s_transv=(MSPLINE*)alloca(sizeof(MSPLINE));
  s_test=(SPLINE*)alloca(sizeof(SPLINE));
  s_e=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_x=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_xp=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_y=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_yp=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_x_x=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_x_xp=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_xp_xp=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_y_y=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_y_yp=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  s_yp_yp=(SPLINE*)alloca(sizeof(SPLINE)*nm);
  xvalue=(double*)alloca(sizeof(double)*max(n,nm));
  yvalue=(double*)alloca(sizeof(double)*max(n,nm));
  ytransv=(double*)alloca(sizeof(double)*nm*11);

  for (i=0;i<n;i++) {
      xvalue[i]=beam->z_position[n_off/nm+i];
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->particle[n_off+i*nm+j].energy;
      }
      spline_init(xvalue,0,yvalue,0,n,s_e+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->particle[n_off+i*nm+j].x;
      }
      spline_init(xvalue,0,yvalue,0,n,s_x+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->particle[n_off+i*nm+j].xp;
      }
      spline_init(xvalue,0,yvalue,0,n,s_xp+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->particle[n_off+i*nm+j].y;
      }
      spline_init(xvalue,0,yvalue,0,n,s_y+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->particle[n_off+i*nm+j].yp;
      }
      spline_init(xvalue,0,yvalue,0,n,s_yp+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->sigma_xx[n_off+i*nm+j].r11;
      }
      spline_init(xvalue,0,yvalue,0,n,s_x_x+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->sigma_xx[n_off+i*nm+j].r12;
      }
      spline_init(xvalue,0,yvalue,0,n,s_x_xp+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->sigma_xx[n_off+i*nm+j].r22;
      }
      spline_init(xvalue,0,yvalue,0,n,s_xp_xp+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->sigma[n_off+i*nm+j].r11;
      }
      spline_init(xvalue,0,yvalue,0,n,s_y_y+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->sigma[n_off+i*nm+j].r12;
      }
      spline_init(xvalue,0,yvalue,0,n,s_y_yp+j);
  }
  for (j=0;j<nm;j++) {
      for (i=0;i<n;i++) {
	  yvalue[i]=beam->sigma[n_off+i*nm+j].r22;
      }
      spline_init(xvalue,0,yvalue,0,n,s_yp_yp+j);
  }
  for (i=0;i<nm;i++) {
      xvalue[i]=2.0*e_cut*(i+0.5)/(double)nm-e_cut;
  }
  sigma_z=(beam->z_position[n_off/nm+1]-beam->z_position[n_off/nm])*n
    /(2.0*z_cut);
  for (np=0;np<n_particles;np++){
    while (fabs(z_rndm=gasdev())>3.0) {}
    while (fabs(e_rndm=gasdev())>3.0) {}
    z_pos=sigma_z*z_rndm;
    for (j=0;j<nm;j++){
	ytransv[j*11]=spline_int(s_e+j,z_pos);
	ytransv[j*11+1]=spline_int(s_x+j,z_pos);
	ytransv[j*11+2]=spline_int(s_xp+j,z_pos);
	ytransv[j*11+3]=spline_int(s_y+j,z_pos);
	ytransv[j*11+4]=spline_int(s_yp+j,z_pos);
	ytransv[j*11+5]=spline_int(s_x_x+j,z_pos);
	ytransv[j*11+6]=spline_int(s_x_xp+j,z_pos);
	ytransv[j*11+7]=spline_int(s_xp_xp+j,z_pos);
	ytransv[j*11+8]=spline_int(s_y_y+j,z_pos);
	ytransv[j*11+9]=spline_int(s_y_yp+j,z_pos);
	ytransv[j*11+10]=spline_int(s_yp_yp+j,z_pos);
    }
    mspline_init(xvalue,0,ytransv,0,nm,11,s_transv);
    mspline_int(s_transv,e_rndm,yvalue);
    mspline_delete(s_transv);

    x_rndm=gasdev()*sqrt(yvalue[5]);
    xp_rndm=gasdev()*sqrt(yvalue[7]-yvalue[6]*yvalue[6]/yvalue[5])
      +x_rndm*yvalue[6]/yvalue[5];
    y_rndm=gasdev()*sqrt(yvalue[8]);
    yp_rndm=gasdev()*sqrt(yvalue[10]-yvalue[9]*yvalue[9]/yvalue[8])
      +y_rndm*yvalue[9]/yvalue[8];
    
    x_rndm+=yvalue[1];
    xp_rndm+=yvalue[2];
    y_rndm+=yvalue[3];
    yp_rndm+=yvalue[4];
    e_rndm=yvalue[0];

    switch (which) {
      case 1:
	fprintf(file,"%g %g %g %g %g %g \n",
		x_rndm*.001, xp_rndm*.001, y_rndm*.001, yp_rndm*.001, energy_permille, z_pos*.001);
	break;
    case 2:
      fprintf(file,"%g %g %g %g %g %g \n",
	      e_rndm,x_rndm,y_rndm,
	      z_pos,xp_rndm,yp_rndm);
      break;
    }
  }
  if (name){
    fclose(file);
  }  
}

main(int argc,char *argv[])
{
  rndmst5(12,34,56,78);
  data_MAD_track(read_beam_file(argv[1],strtol(argv[2],NULL,10)),
		 "particles.out",strtol(argv[3],NULL,10),1,0,-1.0,2);
  exit(0);
}
