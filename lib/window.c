#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "define.h"
#include "common.h"

static Selection radials[MAX_RADIALS] = {{.size=0},{.size=0},{.size=0},{.size=0}};
static Selection angular;
static Selection window;
static Pole poles[MAX_POLES] = {{.n_ells=0},{.n_ells=0}};
static LOS los[MAX_LOS] = {{.type=LOS_MIDPOINT},{.type=LOS_MIDPOINT}};

void print_num_threads()
{
	//Calculate number of threads
	size_t num_threads=0;
#pragma omp parallel
	{
#pragma omp atomic
		num_threads++;
	}
	printf(" - using %zu threads\n",num_threads);
}

void set_num_threads(size_t num_threads)
{
	if (num_threads>0) omp_set_num_threads(num_threads);
	if (verbose == INFO) print_num_threads();
}

void set_verbosity(char* mode)
{
	if (!strcmp(mode,"quiet")) verbose = QUIET;
	if (!strcmp(mode,"info")) verbose = INFO;
	if (!strcmp(mode,"debug")) verbose = DEBUG;
}

void print_pole(Pole pole)
{
	printf("*** Multipoles\n");
	size_t ill;
	if (pole.type==MULTI_ALL) {
		printf(" - multi-type: all\n");
		printf(" - ells:");
		for (ill=0;ill<pole.n_ells;ill++) printf(" %zu",ill);
		printf("\n");
	}
	if (pole.type==MULTI_EVEN) {
		printf(" - multi-type: even\n");
		printf(" - ells:");
		for (ill=0;ill<pole.n_ells;ill++) printf(" %zu",2*ill);
		printf("\n");
	}
	if (pole.type==MULTI_ODD) {
		printf(" - multi-type: odd\n");
		printf(" - ells:");
		for (ill=0;ill<pole.n_ells;ill++) printf(" %zu",2*ill+1);
		printf("\n");
	}
}

void set_pole(size_t num,char *type,size_t n_ells)
{
	Pole pole;
	if (!strcmp(type,"all")) pole.type=MULTI_ALL;
	else if (!strcmp(type,"even")) pole.type=MULTI_EVEN;
	else if (!strcmp(type,"odd")) pole.type=MULTI_ODD;
	else {
		pole.type=MULTI_ALL;
		fprintf(stderr," - invalid multipole type. Choices: all, even or odd.\n");
		fprintf(stderr," - I choose all.\n");
	}
	pole.n_ells = n_ells;
	poles[num-1] = pole;
	if (verbose == INFO) print_pole(pole);
}

void clear_poles()
{
	size_t ipole;
	for (ipole=0;ipole<MAX_POLES;ipole++) poles[ipole].n_ells = 0;
}

void print_los(LOS l)
{
	printf("*** Line-of-sight\n");
	if (l.type==LOS_MIDPOINT) printf(" - los-type: midpoint\n");
	else if (l.type==LOS_ENDPOINT) printf(" - los-type: endpoint\n");
	else if (l.type==LOS_FIRSTPOINT) printf(" - los-type: firstpoint\n");
	printf(" - los-n: %zu\n",l.n);
}

void set_los(size_t num,char* type,size_t n)
{
	LOS l;
	if (!strcmp(type,"midpoint")) l.type = LOS_MIDPOINT;
	else if (!strcmp(type,"endpoint")) l.type = LOS_ENDPOINT;
	else if (!strcmp(type,"firstpoint")) l.type = LOS_FIRSTPOINT;
	else {
		l.type = LOS_MIDPOINT;
		fprintf(stderr," - invalid los type. Choices: midpoint, endpoint or firstpoint.\n");
		fprintf(stderr," - I choose midpoint.\n");
	}
	l.n = n;
	los[num-1] = l;
	if (verbose == INFO) print_los(l);
}

void print_window()
{
	size_t idim,start=0,end=0;
	for (idim=0;idim<window.n_dim;idim++) {
		end += window.shape[idim];
		printf(" - window %zu: %zu points from %.3f to %.3f\n",idim+1,window.shape[idim],window.x[start],window.x[end-1]);
		start = end;
	}
}

void set_window(histo_t* x,histo_t* y,size_t* s,size_t n)
{
	window.x = x;
	window.y = y;
	window.shape = s;
	window.n_dim = n;
	window.size = 1;
	size_t idim;
	for (idim=0;idim<window.n_dim;idim++) window.size *= window.shape[idim];
}

void print_angular_selection()
{
	size_t idim,start=0,end=0;
	for (idim=0;idim<angular.n_dim;idim++) {
		end += angular.shape[idim];
		printf(" - angular %zu: %zu points from %.3f to %.3f\n",idim+1,angular.shape[idim],angular.x[start],angular.x[end-1]);
		start = end;
	}
}

void set_angular_selection(histo_t* x,histo_t* y,size_t* s,size_t n)
{
	angular.x = x;
	angular.y = y;
	angular.shape = s;
	angular.n_dim = n;
	angular.size = 1;
	size_t idim;
	for (idim=0;idim<angular.n_dim;idim++) angular.size *= angular.shape[idim];
}

histo_t find_angular_selection_1d(histo_t x,INTERPOL interpol)
{
	return find_selection_1d(angular,x,interpol);
}

histo_t find_angular_selection_2d(histo_t x,histo_t y)
{
	return find_selection_2d(angular,x,y);
}

void clear_radial_selections()
{
	size_t irad;
	for (irad=0;irad<MAX_RADIALS;irad++) radials[irad].size = 0;
}

void print_radial_selections()
{
	size_t irad;
	for (irad=0;irad<MAX_RADIALS;irad++) {
		if (radials[irad].size>0) printf(" - radial %zu: %zu points from %.3f to %.3f\n",irad+1,radials[irad].size,radials[irad].x[0],radials[irad].x[radials[irad].size-1]);
	}
}

void set_radial_selection(size_t num,histo_t* x,histo_t* y,size_t n)
{
	Selection radial;
	radial.x = x;
	radial.y = y;
	radial.size = n;
	radials[num-1] = radial;
}

void run_2pcf_multi(char* type,size_t num_threads)
{
	timer(0);
	if (verbose == INFO) {
		printf("*** 2-point %s correlation function multipoles\n",type);
		print_angular_selection();
		print_radial_selections();
		print_window();
	}
	set_num_threads(num_threads);
	if (!strcmp(type,"global")) cross_2pcf_multi(angular,radials,window,poles[0],los[0]);
	if (!strcmp(type,"radial")) cross_2pcf_multi_radial(angular,radials,window,poles[0],los[0]);
	if (!strcmp(type,"angular")) cross_2pcf_multi_angular(radials,window,poles[0],los[0]);
	if (verbose == INFO) timer(1);
}

void run_3pcf_multi(char* type,size_t num_threads)
{
	timer(0);
	if (verbose == INFO) {
		printf("*** 3-point %s correlation function multipoles\n",type);
		print_angular_selection();
		print_radial_selections();
		print_window();
	}
	set_num_threads(num_threads);
	if (!strcmp(type,"global")) cross_3pcf_multi(angular,radials,window,poles,los);
	if (!strcmp(type,"radial")) cross_3pcf_multi_radial(angular,radials,window,poles,los);
	if (!strcmp(type,"angular")) cross_3pcf_multi_angular(angular,radials,window,poles,los);
	if (verbose == INFO) timer(1);
}

void run_4pcf_multi(char* type,size_t num_threads)
{
	timer(0);
	if (verbose == INFO) {
		printf("*** 4-point %s correlation function multipoles\n",type);
		print_angular_selection();
		print_radial_selections();
		print_window();
	}
	set_num_threads(num_threads);
	if (!strcmp(type,"global-global")) cross_4pcf_multi(angular,radials,window,poles,los);
	if (!strcmp(type,"radial-radial")) cross_4pcf_multi_radial_radial(angular,radials,window,poles,los);
	if (!strcmp(type,"radial-global")) cross_4pcf_multi_radial_global(angular,radials,window,poles,los);
	if (!strcmp(type,"angular-angular")) cross_4pcf_multi_angular_angular(angular,radials,window,poles,los);
	if (!strcmp(type,"angular-global")) cross_4pcf_multi_angular_global(angular,radials,window,poles,los);
	if (!strcmp(type,"angular-radial")) cross_4pcf_multi_angular_radial(angular,radials,window,poles,los);
	if (verbose == INFO) timer(1);
}
