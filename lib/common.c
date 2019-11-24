#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "define.h"
#include "common.h"

//  Timing variables
#ifdef _HAVE_OMP
#include <omp.h>
static double relbeg,relend,absbeg,absend;
#else //_HAVE_OMP
#include <time.h>
static time_t relbeg,relend,absbeg,absend;
#endif //_HAVE_OMP

///////////////////////////
//General purpose functions

static histo_t my_abs(histo_t x)
{
#ifdef _FLOAT32
	return fabsf(x);
#else
	return fabs(x);
#endif //_FLOAT32
}

void unravel_index(size_t ind,size_t* shape,size_t n_dim,size_t* index)
{
	index[0] = ind%shape[0];
	size_t ii;
	for (ii=1;ii<n_dim;ii++) {
		ind = (ind-index[ii-1])/shape[ii-1];
		index[ii] = ind%shape[ii];
	}
}

static size_t get_dichotomy_index(histo_t x,histo_t *tab,size_t min,size_t max)
{
	if (min == max-1) return min;
	else {
		size_t ind = (min+max)/2;
		if (x < tab[ind]) return get_dichotomy_index(x,tab,min,ind);
		else if (x > tab[ind]) return get_dichotomy_index(x,tab,ind,max);
		return ind;
	}
}

static histo_t interpol_lin(histo_t x,histo_t xmin,histo_t xmax,histo_t ymin,histo_t ymax) 
{
	return (x-xmin)/(xmax-xmin)*(ymax-ymin) + ymin;
}

histo_t interpol_bilin(histo_t x,histo_t y,histo_t *fx,histo_t *fy,histo_t* f) 
{
	histo_t dx = (x-fx[0])/(fx[1]-fx[0]);
	histo_t dy = (y-fy[0])/(fy[1]-fy[0]);
	histo_t dfx = f[1]-f[0];
	histo_t dfy = f[2]-f[0];
	histo_t dfxy = f[3]+f[0]-f[1]-f[2];
	return dfx*dx + dfy*dy + dfxy*dx*dy + f[0];
}

static histo_t interpol_poly(histo_t x,histo_t *tabx,histo_t *taby,size_t nx) 
{
	size_t ix,jx,ixclose = 0;
	histo_t* c = (histo_t*) malloc(nx*sizeof(histo_t));
	histo_t* d = (histo_t*) malloc(nx*sizeof(histo_t));
	histo_t diff = my_abs(x-tabx[ixclose]);
	for (ix=0;ix<nx;ix++) {
		histo_t tmpdiff = my_abs(x-tabx[ix]);
        if (tmpdiff<diff) {
        	ixclose = ix;
        	diff = tmpdiff;
        }
        c[ix] = taby[ix];
        d[ix] = taby[ix];
	}
	//printf("%d %.3f\n",(int) ixclose,taby[ixclose]);
	histo_t y = taby[ixclose];
	ixclose -= 1;
	for (jx=1;jx<nx;jx++) {
		for (ix=0;ix<nx-jx;ix++) {
			histo_t den = (c[ix+1]-d[ix])/(tabx[ix]-tabx[ix+jx]);
			c[ix] = (tabx[ix]-x)*den;
			d[ix] = (tabx[ix+jx]-x)*den;
			//printf("%.3f %.3f %.3f %.3f\n",c[ix+1]-d[ix],tabx[ix]-tabx[ix+jx],c[ix],d[ix]);
		}
		histo_t dy;
		if (2*(ixclose+1)<nx-jx) dy = c[ixclose+1];
		else {
			dy = d[ixclose];
			ixclose -= 1;			
		}
		y += dy;
		//printf("%d %.3f\n",(int) ixclose,y);
	}
	free(c);
	free(d);
	return y;
}

histo_t find_selection_1d(Selection selection,histo_t x)
{
	size_t end = selection.size;
	histo_t xstart = selection.x[0];
	histo_t xend = selection.x[end-1];
	
	if ((x<xstart)||(x>xend)) return 0.;
	size_t ind = get_dichotomy_index(x,selection.x,0,end);
	if (selection.interpol == POLY) {
		size_t ixmin = MAX(0,((long) ind)-2);
		size_t ixmax = MIN(end-1,ind+2);
		return interpol_poly(x,&(selection.x[ixmin]),&(selection.y[ixmin]),ixmax-ixmin+1);
	}
	return interpol_lin(x,selection.x[ind],selection.x[ind+1],selection.y[ind],selection.y[ind+1]);
}

histo_t find_selection_2d(Selection selection,histo_t x,histo_t y)
{
	size_t endx = selection.shape[0];
	size_t sizey = selection.shape[1];
	size_t endy = sizey+endx;
	histo_t xstart = selection.x[0];
	histo_t xend = selection.x[endx-1];
	histo_t ystart = selection.x[endx];
	histo_t yend = selection.x[endy-1];
	
	//printf("%.2f %.2f %.2f  %.2f %.2f %.2f\n",xstart,x,xend,ystart,y,yend);
	if ((x<xstart)||(x>xend)||(y<ystart)||(y>yend)) return 0.;
	size_t indx = get_dichotomy_index(x,selection.x,0,endx);
	size_t indy = get_dichotomy_index(y,selection.x,endx,endy);
	//printf("%zu %zu\n",indx,indy);
	histo_t fx[2] = {selection.x[indx],selection.x[indx+1]};
	histo_t fy[2] = {selection.x[indy],selection.x[indy+1]};
	indy -= endx;
	histo_t f[4] = {selection.y[indy+sizey*indx],selection.y[indy+sizey*(indx+1)],selection.y[indy+1+sizey*indx],selection.y[indy+1+sizey*(indx+1)]};
	return interpol_bilin(x,y,fx,fy,f);
}

void find_selection_range_1d(Selection selection,histo_t *min,histo_t *max,char *type)
{
	histo_t min_=1.,max_=-1.;
	size_t sizex = selection.shape[0];
	size_t ix;
	for (ix=0;ix<sizex;ix++) {
		if (selection.y[ix] != 0.) {
			histo_t x = selection.x[ix];
			if (x < min_) min_ = x;
			if (x > max_) max_ = x;
		}
	}
	if (min_ > *min) *min = min_;
	if (max_ < *max) *max = max_;
	if (verbose == INFO) printf(" - using %s bounds %.4f - %.4f\n",type,*min,*max);
}

void find_selection_range_2d(Selection selection,histo_t *min,histo_t *max,char *type)
{
	histo_t min_=1.,max_=-1.;
	size_t sizex = selection.shape[0];
	size_t sizey = selection.shape[1];
	size_t ix;
	for (ix=0;ix<sizex;ix++) {
		size_t iy;
		for (iy=0;iy<sizey;iy++) {
			if (selection.y[iy+sizey*ix] != 0.) {
				histo_t x = selection.x[ix];
				if (x < min_) min_ = x;
				if (x > max_) max_ = x;
				histo_t y = selection.x[sizex+iy];
				if (y < min_) min_ = y;
				if (y > max_) max_ = y;
			}
		}
	}
	if (min_ > *min) *min = min_;
	if (max_ < *max) *max = max_;
	if (verbose == INFO) printf(" - using %s bounds %.4f - %.4f\n",type,*min,*max);
}

static _Bool set_integration(INTEGRATION *current,char *new)
{
	if (!strcmp(new,"trapz")) {
		*current = TRAPZ;
		return 1;
	}
	else if (!strcmp(new,"gauleg")) {
		*current = GAULEG;
		return 1;
	}
	return 0;
}

_Bool set_precision_integration(Precision *precision,size_t n,histo_t min,histo_t max,char *integration,const Precision *precision_default)
{
	_Bool change = 0;
	if (n>0) {
		precision->n = n;
		change = 1;
	}
	if (max>=min) {
		precision->min = min;
		precision->max = max;
		change = 1;
	}
	if (set_integration(&(precision->integration),integration)) {
		change = 1;
	}
	if (!change) *precision = *precision_default;
	return change;
}

void nodes_weights_gauss_legendre(histo_t xmin,histo_t xmax,histo_t *x,histo_t *w,size_t n) {

	/*
	xmin: minimum x
	xmax: maximum x
	x: array of nodes
	w: array of weights
	n: size of x,w
	*/

	size_t ii,jj,mid = (n+1)/2;
	histo_t xmid = (xmax+xmin)/2.;
	histo_t xl = (xmax-xmin)/2.;
	histo_t p1=0.,p2=0.,p3=0.,pp=0.,z=0.,z1=0.;
	
	for (ii=0;ii<mid;ii++) {
		z = cos(M_PI*(ii+.75)/(n+.5));
		z1 = z + EPS + 1.;
		while (my_abs(z-z1)>EPS) {
			p1 = 1.;
			p2 = 0.;
			for (jj=0;jj<n;jj++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.*jj+1.)*z*p2-jj*p3)/(jj+1.);
			}
			pp=n*(z*p1-p2)/(z*z-1.);
			z1=z;
			z=z1-p1/pp;
		}
		x[ii] = xmid-xl*z;
		x[n-ii-1] = xmid+xl*z;
		w[ii] = 2.*xl/((1.-z*z)*pp*pp);
 		w[n-ii-1] = w[ii];
	}
}

void nodes_weights_trapz(histo_t xmin,histo_t xmax,histo_t *x,histo_t *w,size_t n) {

	/*
	xmin: minimum x
	xmax: maximum x
	x: array of nodes
	w: array of weights
	n: size of x,w
	*/

	histo_t step = (xmax-xmin)/(n-1);
	size_t ii;
	for (ii=0;ii<n;ii++) x[ii] = ii*step + xmin;
	for (ii=1;ii<n-1;ii++) w[ii] = (x[ii+1]-x[ii-1])/2.;
	w[0] = (x[1]-x[0])/2.;
	w[n-1] = (x[n-1]-x[n-2])/2.;
}

void init_integration(Integration *integration,Precision *precision)
{
	histo_t min=precision->min,max=precision->max;
	size_t size = sizeof(histo_t);
	size_t n = precision->n;
	integration->n = n;
	integration->x = (histo_t *) calloc(n,size);
	integration->w = (histo_t *) calloc(n,size);
	integration->type = precision->integration;
	
	if (integration->type == TRAPZ) {
		nodes_weights_trapz(min,max,integration->x,integration->w,integration->n);
	}
	if (integration->type == GAULEG) {
		nodes_weights_gauss_legendre(min,max,integration->x,integration->w,integration->n);
	}
}

void free_integration(Integration *integration)
{
	free(integration->x);
	free(integration->w);
}

static void legendre_all(histo_t mu,histo_t mu2,histo_t leg[])
{
	histo_t mu3 = mu2*mu;
	histo_t mu4 = mu2*mu2;
	histo_t mu5 = mu4*mu;
	histo_t mu6 = mu4*mu2;
	histo_t mu7 = mu6*mu;
	histo_t mu8 = mu6*mu2;
	//histo_t mu10 = mu8*mu2;
	//histo_t mu12 = mu10*mu2;
	leg[0] = 1.;
	leg[1] = mu;
	leg[2] = 0.5*(3.*mu2-1.);
	leg[3] = 0.5*(5.*mu3-3.*mu);
	leg[4] = 1./8.*(35.*mu4-30.*mu2+3.);
	leg[5] = 1./8.*(63.*mu5-70.*mu3+15.*mu);
	leg[6] = 1./16.*(231.*mu6-315.*mu4+105*mu2-5.);
	leg[7] = 1./16.*(429.*mu7-693.*mu5+315*mu3-35.*mu);
	leg[8] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
	//leg[5] = 1./256.*(46189.*mu10-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
	//leg[6] = 1./1024.*(676039.*mu12-1939938.*mu10+2078505.*mu8-1021020.*mu6+225225.*mu4-18018.*mu2+231.);
}

static void legendre_even(histo_t mu2,histo_t leg[])
{
	histo_t mu4 = mu2*mu2;
	histo_t mu6 = mu4*mu2;
	histo_t mu8 = mu6*mu2;
	histo_t mu10 = mu8*mu2;
	histo_t mu12 = mu10*mu2;
	leg[0] = 1.;
	leg[1] = 0.5*(3.*mu2-1.);
	leg[2] = 1./8.*(35.*mu4-30.*mu2+3.);
	leg[3] = 1./16.*(231.*mu6-315.*mu4+105*mu2-5.);
	leg[4] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
	leg[5] = 1./256.*(46189.*mu10-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
	leg[6] = 1./1024.*(676039.*mu12-1939938.*mu10+2078505.*mu8-1021020.*mu6+225225.*mu4-18018.*mu2+231.);
}

static void legendre_odd(histo_t mu,histo_t mu2,histo_t leg[])
{
	histo_t mu3 = mu2*mu;
	histo_t mu5 = mu3*mu2;
	histo_t mu7 = mu5*mu2;
	histo_t mu9 = mu7*mu2;
	histo_t mu11 = mu9*mu2;
	leg[0] = mu;
	leg[1] = 0.5*(5.*mu3-3.*mu);
	leg[2] = 1./8.*(63.*mu5-70.*mu3+15.*mu);
	leg[3] = 1./16.*(429.*mu7-693.*mu5+315*mu3-35.*mu);
	leg[4] = 1./128.*(12155.*mu9-25740.*mu7+18018*mu5-4620.*mu3+315.*mu);
	leg[5] = 1./256.*(88179.*mu11-230945.*mu9+218790*mu7-90090.*mu5+15015.*mu3-693.*mu);
}

void legendre(histo_t dist,histo_t leg[],MULTI_TYPE type) {

	if (type==MULTI_ALL) {
		legendre_all(dist,dist*dist,leg);
	}
	else if (type==MULTI_EVEN) {
		legendre_even(dist*dist,leg);
	}
	else if (type==MULTI_ODD) {
		legendre_odd(dist,dist*dist,leg);
	}
}

void timer(size_t i)
{
	/////
	// Timing routine
	// timer(0) -> initialize relative clock
	// timer(1) -> read relative clock
	// timer(2) -> read relative clock and initialize it afterwards
	// timer(4) -> initialize absolute clock
	// timer(5) -> read absolute clock
#ifdef _HAVE_OMP
	if(i==0)
		relbeg=omp_get_wtime();
	else if(i==1) {
		relend=omp_get_wtime();
		printf(" - relative time ellapsed %.1f ms\n",1000*(relend-relbeg));
	}    
	else if(i==2) {
		relend=omp_get_wtime();
		printf(" - relative time ellapsed %.1f ms\n",1000*(relend-relbeg));
		relbeg=omp_get_wtime();
	}
	else if(i==4)
		absbeg=omp_get_wtime();
	else if(i==5) {
		absend=omp_get_wtime();
		printf(" - total time ellapsed %.1f ms \n",1000*(absend-absbeg));
	}
#else //_HAVE_OMP
	int diff;
	
	if(i==0)
		relbeg=time(NULL);
	else if(i==1) {
		relend=time(NULL);
		diff=(int)(difftime(relend,relbeg));
		printf(" - relative time ellapsed %02d:%02d:%02d \n",
		 diff/3600,(diff/60)%60,diff%60);
	}    
	else if(i==2) {
		relend=time(NULL);
		diff=(size_t)(difftime(relend,relbeg));
		printf(" - relative time ellapsed %02d:%02d:%02d \n",
		 diff/3600,(diff/60)%60,diff%60);
		relbeg=time(NULL);
	}
	else if(i==4)
		absbeg=time(NULL);
	else if(i==5) {
		absend=time(NULL);
		diff=(size_t)(difftime(absend,absbeg));
		printf(" - total time ellapsed %02d:%02d:%02d \n",
		 diff/3600,(diff/60)%60,diff%60);
	}
#endif //_HAVE_OMP
}
