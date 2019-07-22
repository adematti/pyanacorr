#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "define.h"
#include "common.h"

static histo_t my_sqrt(histo_t x)
{
#ifdef _FLOAT32
	return sqrtf(x);
#else
	return sqrt(x);
#endif	
}

static histo_t my_abs(histo_t x)
{
#ifdef _FLOAT32
	return fabsf(x);
#else
	return fabs(x);
#endif //_FLOAT32
}

static histo_t power(histo_t base,size_t exp)
{
    histo_t result = 1.;
    while (exp)
    {
        if (exp & 1) result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

static histo_t get_distance_losn(histo_t dist,size_t losn)
{
	return power(dist,losn);
}

static histo_t get_distance_mu(histo_t x1,histo_t x2,histo_t x12,LOS_TYPE type,histo_t* dist_los)
{
	if (type==LOS_ENDPOINT) {
		*dist_los = x1;
		return (x2*x2-x1*x1-x12*x12)/(2.*x1*x12);
	}
	else if (type==LOS_FIRSTPOINT) {
		*dist_los = x2;
		return (x2*x2-x1*x1+x12*x12)/(2.*x2*x12);
	}
	else { //MIDPOINT
		histo_t mus = (x2*x2-x1*x1-x12*x12)/(2.*x1*x12);
		histo_t x2 = my_sqrt(x1*x1 + mus*x1*x12 + x12*x12/4.);
		return get_distance_mu(x1,x2,x12/2.,LOS_FIRSTPOINT,dist_los);  
	}
}
/*
void cross_2pcf_multi(Selection angular,Selection* radials,Selection window,Pole pole,LOS los)
{
	histo_t *thready;
	size_t iy;
	size_t nw_tot = window.size*pole.n_ells;

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,pole,los) private(thready)
	{
		size_t ix,iy;
		MULTI_TYPE multi_type = pole.type;
		size_t n_ells = pole.n_ells;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		size_t n_s = window.size;

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			histo_t px = x*x*radials[0].y[ix];
			size_t is;
			for (is=0;is<n_s;is++) {
				histo_t s = window.x[is];
				size_t imu;
				const size_t nmu = 1000;
				for (imu=0;imu<nmu;imu++) {
					histo_t mu = -1. + imu*2./(nmu-1);
					histo_t xs2 = x*x + 2.*x*s*mu + s*s;
					if (xs2<0.) continue;
					histo_t xs = my_sqrt(xs2);
					histo_t costheta = (xs2+x*x-s*s)/(2.*x*xs);
					//printf("%.4f %.4f %.4f %.4f %.4f\n",costheta,s,x,xs,get_distance_mu(x,xs,s,LOS_MIDPOINT,&dist_los));
					histo_t weight = find_selection_1d(angular,costheta,LIN)*px*find_selection_1d(radials[1],xs,LIN);
					if (weight == 0.) continue;
					histo_t dist_los;
					histo_t leg[MAX_ELLS];
					legendre(get_distance_mu(x,xs,s,los.type,&dist_los),leg,multi_type);
					//histo_t mu = get_distance_mu(x,xs,s,LOS_MIDPOINT,&dist_los);
					//if ((mu<-1.) || (mu>1.)) printf("%.4e\n",mu);
					weight /= get_distance_losn(dist_los,los.n);
					size_t ill;
					for (ill=0;ill<n_ells;ill++) thready[ill+is*n_ells] += weight*leg[ill];
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			free(thready);
		}
	} //end omp parallel
}
*/
/*
void cross_2pcf_multi(Selection angular,Selection* radials,Selection window,Pole pole,LOS los)
{
	histo_t *thready;
	size_t iy;
	size_t nw_tot = window.size*pole.n_ells;

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,pole,los) private(thready)
	{
		size_t ix,iy;
		MULTI_TYPE multi_type = pole.type;
		size_t n_ells = pole.n_ells;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		size_t n_s = window.size;

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			histo_t px = x*x*radials[0].y[ix];
			size_t ixs;
			for (ixs=0;ixs<radials[1].size;ixs++) {
				histo_t xs = radials[1].x[ixs];
				size_t is;
				for (is=0;is<n_s;is++) {
					histo_t s = window.x[is];
					histo_t costheta = (xs*xs+x*x-s*s)/(2.*x*xs);
					if (costheta < angular.x[0]) break;
					//printf("%.4f %.4f %.4f %.4f %.4f\n",costheta,s,x,xs,get_distance_mu(x,xs,s,LOS_MIDPOINT,&dist_los));
					histo_t weight = xs*x/s*find_selection_1d(angular,costheta,LIN)*px*radials[1].y[ixs];
					if (weight == 0.) continue;
					histo_t dist_los;
					histo_t leg[MAX_ELLS];
					legendre(get_distance_mu(x,xs,s,los.type,&dist_los),leg,multi_type);
					//histo_t mu = get_distance_mu(x,xs,s,LOS_MIDPOINT,&dist_los);
					//if ((mu<-1.) || (mu>1.)) printf("%.4e\n",mu);
					weight /= get_distance_losn(dist_los,los.n);
					size_t ill;
					for (ill=0;ill<n_ells;ill++) thready[ill+is*n_ells] += weight*leg[ill];
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			free(thready);
		}
	} //end omp parallel
}
*/

void cross_2pcf_multi(Selection angular,Selection* radials,Selection window,Pole pole,LOS los)
{
	histo_t *thready;
	size_t iy;
	size_t nw_tot = window.size*pole.n_ells;

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,pole,los) private(thready)
	{
		size_t icos,iy;
		const histo_t sign[2] = {-1.,1.};
		MULTI_TYPE multi_type = pole.type;
		size_t n_ells = pole.n_ells;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		size_t n_s = window.size;

#pragma omp for nowait schedule(dynamic)
		for (icos=0;icos<angular.size;icos++) {
			histo_t costheta = angular.x[icos];
			histo_t pang = angular.y[icos];
			if (pang == 0.) continue;
			size_t ix;
			for (ix=0;ix<radials[0].size;ix++) {
				histo_t x = radials[0].x[ix];
				histo_t px = x*x*radials[0].y[ix];
				int is;
				for (is=n_s-1;is>=0;is--) {
					histo_t s = window.x[is];
					histo_t delta = x*x*(costheta*costheta-1.) + s*s;
					if (delta < 0.) break;
					delta = my_sqrt(delta);
					histo_t dist_los,leg[MAX_ELLS];
					size_t isign;
					for (isign=0;isign<2;isign++) {
						histo_t xs = x*costheta + sign[isign]*delta;
						histo_t weight = pang*px*find_selection_1d(radials[1],xs,LIN);
						weight *= my_abs(2.*xs*xs*xs/s/(xs*xs-x*x+s*s));
						if (weight == 0.) continue;
						legendre(get_distance_mu(x,xs,s,los.type,&dist_los),leg,multi_type);
						weight /= get_distance_losn(dist_los,los.n);
						size_t ill;
						for (ill=0;ill<n_ells;ill++) thready[ill+n_ells*is] += weight*leg[ill];
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			free(thready);
		}
	} //end omp parallel
}


void cross_2pcf_multi_radial(Selection angular,Selection* radials,Selection window,Pole pole,LOS los)
{
	histo_t *thready;
	size_t iy;
	size_t nw_tot = window.size*pole.n_ells;

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,pole,los) private(thready)
	{
		size_t ix,iy;
		MULTI_TYPE multi_type = pole.type;
		size_t n_ells = pole.n_ells;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		size_t n_s = window.size;

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			histo_t px = radials[0].y[ix]*radials[1].y[ix];
			size_t is;
			for (is=0;is<n_s;is++) {
				histo_t s = window.x[is];
				histo_t costheta = 1.-s*s/(2.*x*x);
				if (costheta < angular.x[0]) break;
				histo_t weight = find_selection_1d(angular,costheta,LIN)*px;
				if (weight == 0.) continue;
				histo_t dist_los;
				histo_t leg[MAX_ELLS];
				legendre(get_distance_mu(x,x,s,los.type,&dist_los),leg,multi_type);
				weight /= get_distance_losn(dist_los,los.n);
				size_t ill;
				for (ill=0;ill<n_ells;ill++) thready[ill+is*n_ells] += weight*leg[ill];
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			free(thready);
		}
	} //end omp parallel
}

void cross_2pcf_multi_angular(Selection* radials,Selection window,Pole pole,LOS los)
{
	histo_t *thready;
	size_t iy;
	size_t nw_tot = window.size*pole.n_ells;

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(radials,window,nw_tot,pole,los) private(thready)
	{
		size_t ix,iy;
		const histo_t sign[2] = {-1.,1.};
		MULTI_TYPE multi_type = pole.type;
		size_t n_ells = pole.n_ells;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		size_t n_s = window.size;

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			histo_t px = radials[0].y[ix];
			size_t is;
			for (is=0;is<n_s;is++) {
				histo_t s = window.x[is];
				size_t isign;
				for (isign=0;isign<2;isign++) {
					histo_t xs = x + sign[isign]*s;
					histo_t weight = px*find_selection_1d(radials[1],xs,LIN);
					if (weight == 0.) continue;
					histo_t dist_los;
					histo_t leg[MAX_ELLS];
					legendre(get_distance_mu(x,xs,s,los.type,&dist_los),leg,multi_type);
					weight /= get_distance_losn(dist_los,los.n);
					size_t ill;
					for (ill=0;ill<n_ells;ill++) thready[ill+is*n_ells] += weight*leg[ill];
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			free(thready);
		}
	} //end omp parallel
}

/*
void cross_3pcf_multi(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];
	histo_t* ysec[n_sec];
	histo_t* costhetasec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec,costhetasec)
	{
		size_t ix,iy,isec;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) {
			ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));
			costhetasec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));
		}

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			size_t ixs;
			for (isec=0;isec<n_sec;isec++) {
				histo_t *y2 = ysec[isec];
				for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
				histo_t *costheta2 = costhetasec[isec];
				MULTI_TYPE multi_type2 = poles[isec].type;
				size_t n_ells2 = poles[isec].n_ells;
				LOS los2 = los[isec];
				size_t n_s2 = window.shape[isec];
				histo_t dist_los,leg[MAX_ELLS];
				size_t start_win = (isec==1) ? window.shape[0] : 0;
				size_t start_ang = (isec==1) ? angular.shape[0] : 0;
				histo_t costheta_min = angular.x[start_ang];
				histo_t costheta_max = angular.x[start_ang+angular.shape[isec]-1];
				for (ixs=0;ixs<radials[isec].size;ixs++) {
					histo_t xs = radials[isec].x[ixs];
					histo_t px2 = radials[isec].y[ixs];
					size_t is;
					for (is=0;is<n_s2;is++) {
						histo_t s = window.x[is+start_win];
						histo_t costheta = (xs*xs+x*x-s*s)/(2.*x*xs);
						if (costheta < costheta_min) break;
						if (costheta > costheta_max) continue;
						legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
						histo_t weight = px2/get_distance_losn(dist_los,los2.n);
						size_t ill;
						for (ill=0;ill<n_ells2;ill++) y2[ill+is*n_ells2] += weight*leg[ill];
						costheta2[is] = costheta;
					}
				}
			}
			histo_t px = radials[0].y[ix];
			histo_t *y2 = ysec[0];
			histo_t *y3 = ysec[1];
			histo_t *costheta2 = costhetasec[0];
			histo_t *costheta3 = costhetasec[1];
			size_t n_ells2 = poles[0].n_ells;
			size_t n_ells3 = poles[1].n_ells;
			size_t n_s2 = window.shape[0];
			size_t n_s3 = window.shape[1];
			size_t is2;
			for (is2=0;is2<n_s2;is2++) {
				if (y2[is2*n_ells2] != 0.) {
					size_t is3;
					for (is3=0;is3<n_s3;is3++) {
						if (y3[is3*n_ells3] != 0.) {
							histo_t weight1 = px*find_selection_2d(angular,costheta2[is2],costheta3[is3]);
							size_t ill2,ill3;
							for (ill2=0;ill2<n_ells2;ill2++) {
								for (ill3=0;ill3<n_ells3;ill3++) {
									thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) {
				free(ysec[isec]);
				free(costhetasec[isec]);
			}
			free(thready);
		}
	} //end omp parallel
}
*/


void cross_3pcf_multi(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {angular.shape[0]*window.shape[0]*poles[0].n_ells,angular.shape[1]*window.shape[1]*poles[1].n_ells};
	size_t nw_tot = window.shape[0]*poles[0].n_ells*window.shape[1]*poles[1].n_ells;
	histo_t* ysec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec)
	{
		size_t ix,iy,isec;
		const histo_t sign[2] = {-1.,1.};
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			for (isec=0;isec<n_sec;isec++) {
				histo_t *y2 = ysec[isec];
				for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
				Selection radialsec = radials[isec+1];
				MULTI_TYPE multi_type2 = poles[isec].type;
				size_t n_ells2 = poles[isec].n_ells;
				LOS los2 = los[isec];
				size_t n_s2 = window.shape[isec];
				size_t start_win = (isec==1) ? window.shape[0] : 0;
				size_t start_ang = (isec==1) ? angular.shape[0] : 0;
				size_t icos;
				for (icos=0;icos<angular.shape[isec];icos++) {
					histo_t costheta = angular.x[icos+start_ang];
					int is;
					for (is=n_s2-1;is>=0;is--) {
						histo_t s = window.x[is+start_win];
						histo_t delta = x*x*(costheta*costheta-1.) + s*s;
						if (delta < 0.) break;
						delta = my_sqrt(delta);
						histo_t dist_los,leg[MAX_ELLS];
						size_t isign;
						for (isign=0;isign<2;isign++) {
							histo_t xs = x*costheta + sign[isign]*delta;
							histo_t weight = find_selection_1d(radialsec,xs,LIN);
							if (weight == 0.) continue;
							legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
							weight /= get_distance_losn(dist_los,los2.n);
							size_t ill;
							for (ill=0;ill<n_ells2;ill++) y2[ill+n_ells2*(is+n_s2*icos)] += weight*leg[ill];
						}
					}
				}
			}
			histo_t px = radials[0].y[ix];
			histo_t *y2 = ysec[0];
			histo_t *y3 = ysec[1];
			size_t n_ells2 = poles[0].n_ells;
			size_t n_ells3 = poles[1].n_ells;
			size_t n_s2 = window.shape[0];
			size_t n_s3 = window.shape[1];
			size_t n_cos2 = angular.shape[0];
			size_t n_cos3 = angular.shape[1];
			size_t icos2;
			for (icos2=0;icos2<n_cos2;icos2++) {
				size_t icos3;
				for (icos3=0;icos3<n_cos3;icos3++) {
					histo_t weight1 = px*angular.y[icos3+n_cos2*icos2];
					if (weight1 == 0.) continue;
					size_t is2;
					for (is2=0;is2<n_s2;is2++) {
						size_t ind2 = n_ells2*(is2+n_s2*icos2);
						if (y2[ind2] != 0.) {
							size_t is3;
							for (is3=0;is3<n_s3;is3++) {
								size_t ind3 = n_ells3*(is3+n_s3*icos3);
								if (y3[ind3] != 0.) {
									size_t ill2,ill3;
									for (ill2=0;ill2<n_ells2;ill2++) {
										for (ill3=0;ill3<n_ells3;ill3++) {
											thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+ind2]*y3[ill3+ind3];
										}
									}
								}
							}
						}	
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
			free(thready);
		}
	} //end omp parallel
}

void cross_3pcf_multi_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];
	histo_t* ysec[n_sec];
	histo_t* costhetasec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec,costhetasec)
	{
		size_t ix,iy,isec;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) {
			ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));
			costhetasec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));
		}

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			size_t ixs;
			for (ixs=0;ixs<radials[1].size;ixs++) {
				histo_t xs = radials[1].x[ixs];
				for (isec=0;isec<n_sec;isec++) {
					histo_t *y2 = ysec[isec];
					for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
					histo_t *costheta2 = costhetasec[isec];
					MULTI_TYPE multi_type2 = poles[isec].type;
					size_t n_ells2 = poles[isec].n_ells;
					LOS los2 = los[isec];
					size_t n_s2 = window.shape[isec];
					histo_t dist_los,leg[MAX_ELLS];
					size_t start_win = (isec==1) ? window.shape[0] : 0;
					size_t start_ang = (isec==1) ? angular.shape[0] : 0;
					histo_t costheta_min = angular.x[start_ang];
					histo_t costheta_max = angular.x[start_ang+angular.shape[isec]-1];
					size_t is;
					for (is=0;is<n_s2;is++) {
						histo_t s = window.x[is+start_win];
						histo_t costheta = (xs*xs+x*x-s*s)/(2.*x*xs);
						if (costheta < costheta_min) break;
						if (costheta > costheta_max) continue;
						legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
						histo_t weight = 1./get_distance_losn(dist_los,los2.n);
						size_t ill;
						for (ill=0;ill<n_ells2;ill++) y2[ill+is*n_ells2] = weight*leg[ill];
						costheta2[is] = costheta;
					}
				}
				histo_t pxxs = radials[0].y[ix]*radials[1].y[ixs]*radials[2].y[ixs];
				histo_t *y2 = ysec[0];
				histo_t *y3 = ysec[1];
				histo_t *costheta2 = costhetasec[0];
				histo_t *costheta3 = costhetasec[1];
				size_t n_ells2 = poles[0].n_ells;
				size_t n_ells3 = poles[1].n_ells;
				size_t n_s2 = window.shape[0];
				size_t n_s3 = window.shape[1];
				size_t is2;
				for (is2=0;is2<n_s2;is2++) {
					if (y2[is2*n_ells2] != 0.) {
						size_t is3;
						for (is3=0;is3<n_s3;is3++) {
							if (y3[is3*n_ells3] != 0.) {
								histo_t weight1 = pxxs*find_selection_2d(angular,costheta2[is2],costheta3[is3]);
								if (weight1 == 0.) continue;
								size_t ill2,ill3;
								for (ill2=0;ill2<n_ells2;ill2++) {
									for (ill3=0;ill3<n_ells3;ill3++) {
										thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
									}
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) {
				free(ysec[isec]);
				free(costhetasec[isec]);
			}
			free(thready);
		}
	} //end omp parallel
}


void cross_3pcf_multi_angular(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];
	histo_t* ysec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec)
	{
		size_t icos,iy,isec;
		const histo_t sign[2] = {-1.,1.};
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));

#pragma omp for nowait schedule(dynamic)
		for (icos=0;icos<angular.size;icos++) {
			histo_t costheta = angular.x[icos];
			histo_t pang = angular.y[icos];
			if (pang == 0.) continue;
			size_t ix;
			for (ix=0;ix<radials[0].size;ix++) {
				histo_t x = radials[0].x[ix];
				histo_t px = radials[0].y[ix];
				for (isec=0;isec<n_sec;isec++) {
					histo_t *y2 = ysec[isec];
					for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
					Selection radialsec = radials[isec+1];
					MULTI_TYPE multi_type2 = poles[isec].type;
					size_t n_ells2 = poles[isec].n_ells;
					LOS los2 = los[isec];
					size_t n_s2 = window.shape[isec];
					histo_t dist_los,leg[MAX_ELLS];
					size_t start_win = (isec==1) ? window.shape[0] : 0;
					int is;
					for (is=n_s2-1;is>=0;is--) {
						histo_t s = window.x[is+start_win];
						histo_t delta = x*x*(costheta*costheta-1.) + s*s;
						if (delta < 0.) break;
						delta = my_sqrt(delta);
						size_t isign;
						for (isign=0;isign<2;isign++) {
							histo_t xs = x*costheta + sign[isign]*delta;
							histo_t weight = find_selection_1d(radialsec,xs,LIN);
							if (weight == 0.) continue;
							legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
							//histo_t mu = get_distance_mu(x,xs,s,LOS_MIDPOINT,&dist_los);
							//printf("%.4f %.4f %.4f %.4f %.4f\n",costheta,s,x,xs,get_distance_mu(x,xs,s,los2.type,&dist_los));
							//if ((mu<-1.) || (mu>1.)) printf("%.9e\n",mu);
							weight /= get_distance_losn(dist_los,los2.n);
							size_t ill;
							for (ill=0;ill<n_ells2;ill++) y2[ill+is*n_ells2] += weight*leg[ill];
						}
					}
				}
				histo_t weight1 = pang*px;
				histo_t *y2 = ysec[0];
				histo_t *y3 = ysec[1];
				size_t n_ells2 = poles[0].n_ells;
				size_t n_ells3 = poles[1].n_ells;
				size_t n_s2 = window.shape[0];
				size_t n_s3 = window.shape[1];
				size_t is2;
				for (is2=0;is2<n_s2;is2++) {
					if (y2[is2*n_ells2] != 0.) {
						size_t is3;
						for (is3=0;is3<n_s3;is3++) {
							if (y3[is3*n_ells3] != 0.) {
								size_t ill2,ill3;
								for (ill2=0;ill2<n_ells2;ill2++) {
									for (ill3=0;ill3<n_ells3;ill3++) {
										thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
									}
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
			free(thready);
		}
	} //end omp parallel
}

void cross_4pcf_multi(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;
	
	histo_t *ysec[2];
	size_t isec;
	for (isec=0;isec<n_sec;isec++) {
		ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));
		Selection windowsec;
		windowsec.size = window.shape[isec];
		windowsec.x = &(window.x[(isec==1) ? window.shape[0] : 0]);
		windowsec.y = ysec[isec];
		cross_2pcf_multi(angular,&(radials[isec*2]),windowsec,poles[isec],los[isec]);
	}
	histo_t *y2 = ysec[0];
	histo_t *y3 = ysec[1];
	size_t n_ells2 = poles[0].n_ells;
	size_t n_ells3 = poles[1].n_ells;
	size_t n_s2 = window.shape[0];
	size_t n_s3 = window.shape[1];
	size_t is2;
	for (is2=0;is2<n_s2;is2++) {
		if (y2[is2*n_ells2] != 0.) {
			size_t is3;
			for (is3=0;is3<n_s3;is3++) {
				if (y3[is3*n_ells3] != 0.) {
					size_t ill2,ill3;
					for (ill2=0;ill2<n_ells2;ill2++) {
						for (ill3=0;ill3<n_ells3;ill3++) {
							window.y[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
						}
					}
				}
			}
		}
	}
	
	for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
}

void cross_4pcf_multi_radial_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];
	histo_t* ysec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec)
	{
		size_t ix,iy,isec;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			histo_t px = radials[0].y[ix]*radials[2].y[ix];
			size_t ixs;
			for (ixs=0;ixs<radials[1].size;ixs++) {
				histo_t xs = radials[1].x[ixs];
				for (isec=0;isec<n_sec;isec++) {
					histo_t *y2 = ysec[isec];
					for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
					MULTI_TYPE multi_type2 = poles[isec].type;
					size_t n_ells2 = poles[isec].n_ells;
					LOS los2 = los[isec];
					size_t n_s2 = window.shape[isec];
					histo_t dist_los,leg[MAX_ELLS];
					size_t start_win = (isec==1) ? window.shape[0] : 0;
					histo_t costheta_min = angular.x[0];
					histo_t costheta_max = angular.x[angular.size-1];
					size_t is;
					for (is=0;is<n_s2;is++) {
						histo_t s = window.x[is+start_win];
						histo_t costheta = (xs*xs+x*x-s*s)/(2.*x*xs);
						if (costheta < costheta_min) break;
						if (costheta > costheta_max) continue;
						legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
						histo_t weight = find_selection_1d(angular,costheta,LIN)/get_distance_losn(dist_los,los2.n);
						size_t ill;
						for (ill=0;ill<n_ells2;ill++) y2[ill+is*n_ells2] = weight*leg[ill];
					}
				}
				histo_t weight1 = px*radials[1].y[ixs]*radials[3].y[ixs];
				histo_t *y2 = ysec[0];
				histo_t *y3 = ysec[1];
				size_t n_ells2 = poles[0].n_ells;
				size_t n_ells3 = poles[1].n_ells;
				size_t n_s2 = window.shape[0];
				size_t n_s3 = window.shape[1];
				size_t is2;
				for (is2=0;is2<n_s2;is2++) {
					if (y2[is2*n_ells2] != 0.) {
						size_t is3;
						for (is3=0;is3<n_s3;is3++) {
							if (y3[is3*n_ells3] != 0.) {
								size_t ill2,ill3;
								for (ill2=0;ill2<n_ells2;ill2++) {
									for (ill3=0;ill3<n_ells3;ill3++) {
										thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
									}
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
			free(thready);
		}
	} //end omp parallel
}

void cross_4pcf_multi_radial_global(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];
	histo_t* ysec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec)
	{
		size_t ix,iy,isec;
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));

#pragma omp for nowait schedule(dynamic)
		for (ix=0;ix<radials[0].size;ix++) {
			histo_t x = radials[0].x[ix];
			for (isec=0;isec<n_sec;isec++) {
				histo_t *y2 = ysec[isec];
				for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
				Selection radialsec = radials[isec*2+1];
				MULTI_TYPE multi_type2 = poles[isec].type;
				size_t n_ells2 = poles[isec].n_ells;
				LOS los2 = los[isec];
				size_t n_s2 = window.shape[isec];
				histo_t dist_los,leg[MAX_ELLS];
				size_t start_win = (isec==1) ? window.shape[0] : 0;
				histo_t costheta_min = angular.x[0];
				histo_t costheta_max = angular.x[angular.size-1];
				size_t ixs;
				for (ixs=0;ixs<radialsec.size;ixs++) {
					histo_t xs = radialsec.x[ixs];
					histo_t pxs = radialsec.y[ixs];
					size_t is;
					for (is=0;is<n_s2;is++) {
						histo_t s = window.x[is+start_win];
						histo_t costheta = (xs*xs+x*x-s*s)/(2.*x*xs);
						if (costheta < costheta_min) break;
						if (costheta > costheta_max) continue;
						legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
						histo_t weight = find_selection_1d(angular,costheta,LIN)*pxs/get_distance_losn(dist_los,los2.n);
						size_t ill;
						for (ill=0;ill<n_ells2;ill++) y2[ill+is*n_ells2] += weight*leg[ill];
					}
				}
			}
			histo_t weight1 = radials[0].y[ix]*radials[2].y[ix];
			histo_t *y2 = ysec[0];
			histo_t *y3 = ysec[1];
			size_t n_ells2 = poles[0].n_ells;
			size_t n_ells3 = poles[1].n_ells;
			size_t n_s2 = window.shape[0];
			size_t n_s3 = window.shape[1];
			size_t is2;
			for (is2=0;is2<n_s2;is2++) {
				if (y2[is2*n_ells2] != 0.) {
					size_t is3;
					for (is3=0;is3<n_s3;is3++) {
						if (y3[is3*n_ells3] != 0.) {
							size_t ill2,ill3;
							for (ill2=0;ill2<n_ells2;ill2++) {
								for (ill3=0;ill3<n_ells3;ill3++) {
									thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
			free(thready);
		}
	} //end omp parallel
}

void cross_4pcf_multi_angular_angular(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {window.shape[0]*poles[0].n_ells,window.shape[1]*poles[1].n_ells};
	size_t nw_tot = nw_sec[0]*nw_sec[1];
	histo_t* ysec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec)
	{
		size_t icos,iy,isec;
		const histo_t sign[2] = {-1.,1.};
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));

#pragma omp for nowait schedule(dynamic)
		for (icos=0;icos<angular.size;icos++) {
			histo_t costheta = angular.x[icos];
			histo_t pang = angular.y[icos];
			if (pang == 0.) continue;
			size_t ix;
			for (isec=0;isec<n_sec;isec++) {
				histo_t *y2 = ysec[isec];
				for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
				Selection* radialsec = &(radials[isec*2]);
				MULTI_TYPE multi_type2 = poles[isec].type;
				size_t n_ells2 = poles[isec].n_ells;
				LOS los2 = los[isec];
				size_t n_s2 = window.shape[isec];
				histo_t dist_los,leg[MAX_ELLS];
				size_t start_win = (isec==1) ? window.shape[0] : 0;
				for (ix=0;ix<radialsec[0].size;ix++) {
					histo_t x = radialsec[0].x[ix];
					histo_t px = radialsec[0].y[ix];
					int is;
					for (is=n_s2-1;is>=0;is--) {
						histo_t s = window.x[is+start_win];
						histo_t delta = x*x*(costheta*costheta-1.) + s*s;
						if (delta < 0.) break;
						delta = my_sqrt(delta);
						size_t isign;
						for (isign=0;isign<2;isign++) {
							histo_t xs = x*costheta + sign[isign]*delta;
							histo_t weight = px*find_selection_1d(radialsec[1],xs,LIN);
							if (weight == 0.) continue;
							legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
							//histo_t mu = get_distance_mu(x,xs,s,LOS_MIDPOINT,&dist_los);
							//printf("%.4f %.4f %.4f %.4f %.4f\n",costheta,s,x,xs,get_distance_mu(x,xs,s,los2.type,&dist_los));
							//if ((mu<-1.) || (mu>1.)) printf("%.9e\n",mu);
							weight /= get_distance_losn(dist_los,los2.n);
							size_t ill;
							for (ill=0;ill<n_ells2;ill++) y2[ill+is*n_ells2] += weight*leg[ill];
						}
					}
				}
			}
			histo_t *y2 = ysec[0];
			histo_t *y3 = ysec[1];
			size_t n_ells2 = poles[0].n_ells;
			size_t n_ells3 = poles[1].n_ells;
			size_t n_s2 = window.shape[0];
			size_t n_s3 = window.shape[1];
			size_t is2;
			for (is2=0;is2<n_s2;is2++) {
				if (y2[is2*n_ells2] != 0.) {
					size_t is3;
					for (is3=0;is3<n_s3;is3++) {
						if (y3[is3*n_ells3] != 0.) {
							size_t ill2,ill3;
							for (ill2=0;ill2<n_ells2;ill2++) {
								for (ill3=0;ill3<n_ells3;ill3++) {
									thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += pang*y2[ill2+is2*n_ells2]*y3[ill3+is3*n_ells3];
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
			free(thready);
		}
	} //end omp parallel
}

void cross_4pcf_multi_angular_global(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	size_t iy,isec;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {angular.shape[0]*window.shape[0]*poles[0].n_ells,angular.shape[1]*window.shape[1]*poles[1].n_ells};
	size_t nw_tot = window.shape[0]*poles[0].n_ells*window.shape[1]*poles[1].n_ells;
	histo_t* ysec[n_sec];
	histo_t* thready;

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;
	
	for (isec=0;isec<n_sec;isec++) {
		ysec[isec] = (histo_t*) calloc(nw_sec[isec],sizeof(histo_t));
#pragma omp parallel default(none) shared(angular,radials,window,ysec,nw_tot,nw_sec,poles,los,isec) private(thready)
		{
			size_t icos,iy;
			const histo_t sign[2] = {-1.,1.};
			histo_t *y2 = ysec[isec];
			Selection* radialsec = &(radials[isec*2]);
			MULTI_TYPE multi_type2 = poles[isec].type;
			size_t n_ells2 = poles[isec].n_ells;
			LOS los2 = los[isec];
			size_t n_s2 = window.shape[isec];
			size_t start_win = (isec==1) ? window.shape[0] : 0;
			size_t start_ang = (isec==1) ? angular.shape[0] : 0;
			thready = (histo_t *) calloc(nw_sec[isec],sizeof(histo_t));
#pragma omp for nowait schedule(dynamic)
			for (icos=0;icos<angular.shape[isec];icos++) {
				histo_t costheta = angular.x[icos+start_ang];
				size_t ix;
				for (ix=0;ix<radialsec[0].size;ix++) {
					histo_t x = radialsec[0].x[ix];
					histo_t px = radialsec[0].y[ix];
					int is;
					for (is=n_s2-1;is>=0;is--) {
						histo_t s = window.x[is+start_win];
						histo_t delta = x*x*(costheta*costheta-1.) + s*s;
						if (delta < 0.) break;
						delta = my_sqrt(delta);
						size_t isign;
						histo_t dist_los,leg[MAX_ELLS];
						for (isign=0;isign<2;isign++) {
							histo_t xs = x*costheta + sign[isign]*delta;
							histo_t weight = px*find_selection_1d(radialsec[1],xs,LIN);
							if (weight == 0.) continue;
							legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
							weight /= get_distance_losn(dist_los,los2.n);
							size_t ill;
							for (ill=0;ill<n_ells2;ill++) thready[ill+n_ells2*(is+n_s2*icos)] += weight*leg[ill];
						}
					}
				}
			}
#pragma omp critical
			{
			for (iy=0;iy<nw_sec[isec];iy++) y2[iy] += thready[iy];
			free(thready);
			}
		}
	}
#pragma omp parallel default(none) shared(angular,radials,window,ysec,nw_tot,nw_sec,poles,los) private(thready)
	{
		size_t iy;
		histo_t *y2 = ysec[0];
		histo_t *y3 = ysec[1];
		size_t n_ells2 = poles[0].n_ells;
		size_t n_ells3 = poles[1].n_ells;
		size_t n_s2 = window.shape[0];
		size_t n_s3 = window.shape[1];
		size_t n_cos2 = angular.shape[0];
		size_t n_cos3 = angular.shape[1];
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		size_t icos2;
#pragma omp for nowait schedule(dynamic)
		for (icos2=0;icos2<n_cos2;icos2++) {
			size_t icos3;
			for (icos3=0;icos3<n_cos3;icos3++) {
				histo_t weight1 = angular.y[icos3+n_cos2*icos2];
				if (weight1 == 0.) continue;
				size_t is2;
				for (is2=0;is2<n_s2;is2++) {
					size_t ind2 = n_ells2*(is2+n_s2*icos2);
					if (y2[ind2] != 0.) {
						size_t is3;
						for (is3=0;is3<n_s3;is3++) {
							size_t ind3 = n_ells3*(is3+n_s3*icos3);
							if (y3[ind3] != 0.) {
								size_t ill2,ill3;
								for (ill2=0;ill2<n_ells2;ill2++) {
									for (ill3=0;ill3<n_ells3;ill3++) {
										thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+ind2]*y3[ill3+ind3];
									}
								}
							}
						}
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			free(thready);
		}
	} //end omp parallel
	for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
}


void cross_4pcf_multi_angular_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los)
{
	histo_t *thready;
	size_t iy;
	const size_t n_sec = 2;
	size_t nw_sec[2] = {angular.shape[0]*window.shape[0]*poles[0].n_ells,angular.shape[1]*window.shape[1]*poles[1].n_ells};
	size_t nw_tot = window.shape[0]*poles[0].n_ells*window.shape[1]*poles[1].n_ells;
	histo_t* ysec[n_sec];

	for (iy=0;iy<nw_tot;iy++) window.y[iy] = 0.;

#pragma omp parallel default(none) shared(angular,radials,window,nw_tot,nw_sec,poles,los) private(thready,ysec)
	{
		size_t ixs,iy,isec;
		const histo_t sign[2] = {-1.,1.};
		thready = (histo_t *) calloc(nw_tot,sizeof(histo_t));
		for (isec=0;isec<n_sec;isec++) ysec[isec] = (histo_t*) malloc(nw_sec[isec]*sizeof(histo_t));

#pragma omp for nowait schedule(dynamic)
		for (ixs=0;ixs<radials[1].size;ixs++) {
			histo_t xs = radials[1].x[ixs];
			for (isec=0;isec<n_sec;isec++) {
				histo_t *y2 = ysec[isec];
				for (iy=0;iy<nw_sec[isec];iy++) y2[iy] = 0.;
				Selection radialsec = radials[isec*2];
				MULTI_TYPE multi_type2 = poles[isec].type;
				size_t n_ells2 = poles[isec].n_ells;
				LOS los2 = los[isec];
				size_t n_s2 = window.shape[isec];
				size_t start_win = (isec==1) ? window.shape[0] : 0;
				size_t start_ang = (isec==1) ? angular.shape[0] : 0;
				size_t icos;
				for (icos=0;icos<angular.shape[isec];icos++) {
					histo_t costheta = angular.x[icos+start_ang];
					int is;
					for (is=n_s2-1;is>=0;is--) {
						histo_t s = window.x[is+start_win];
						histo_t delta = xs*xs*(costheta*costheta-1.) + s*s;
						if (delta < 0.) break;
						delta = my_sqrt(delta);
						histo_t dist_los,leg[MAX_ELLS];
						size_t isign;
						for (isign=0;isign<2;isign++) {
							histo_t x = xs*costheta + sign[isign]*delta;
							histo_t weight = find_selection_1d(radialsec,xs,LIN);
							if (weight == 0.) continue;
							legendre(get_distance_mu(x,xs,s,los2.type,&dist_los),leg,multi_type2);
							weight /= get_distance_losn(dist_los,los2.n);
							size_t ill;
							for (ill=0;ill<n_ells2;ill++) y2[ill+n_ells2*(is+n_s2*icos)] += weight*leg[ill];
						}
					}
				}
			}
			histo_t pxs = radials[1].y[ixs]*radials[3].y[ixs];
			histo_t *y2 = ysec[0];
			histo_t *y3 = ysec[1];
			size_t n_ells2 = poles[0].n_ells;
			size_t n_ells3 = poles[1].n_ells;
			size_t n_s2 = window.shape[0];
			size_t n_s3 = window.shape[1];
			size_t n_cos2 = angular.shape[0];
			size_t n_cos3 = angular.shape[1];
			size_t icos2;
			for (icos2=0;icos2<n_cos2;icos2++) {
				size_t icos3;
				for (icos3=0;icos3<n_cos3;icos3++) {
					histo_t weight1 = angular.y[icos3+n_cos2*icos2]*pxs;
					if (weight1 == 0.) continue;
					size_t is2;
					for (is2=0;is2<n_s2;is2++) {
						size_t ind2 = n_ells2*(is2+n_s2*icos2);
						if (y2[ind2] != 0.) {
							size_t is3;
							for (is3=0;is3<n_s3;is3++) {
								size_t ind3 = n_ells3*(is3+n_s3*icos3);
								if (y3[ind3] != 0.) {
									size_t ill2,ill3;
									for (ill2=0;ill2<n_ells2;ill2++) {
										for (ill3=0;ill3<n_ells3;ill3++) {
											thready[ill3+n_ells3*(ill2+n_ells2*(is3+n_s3*is2))] += weight1*y2[ill2+ind2]*y3[ill3+ind3];
										}
									}
								}
							}
						}	
					}
				}
			}
		} // end omp for
	
#pragma omp critical
		{
			for (iy=0;iy<nw_tot;iy++) window.y[iy] += thready[iy];
			for (isec=0;isec<n_sec;isec++) free(ysec[isec]);
			free(thready);
		}
	} //end omp parallel
}
