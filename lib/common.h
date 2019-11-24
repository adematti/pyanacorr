#ifndef _WINDOW_COMMON_
#define _WINDOW_COMMON_

#define EPS 3.e-14

//General-purpose functions

//void unravel_index(size_t ind,size_t* shape,size_t n_dim,size_t* index);

histo_t find_selection_1d(Selection selection,histo_t x);

histo_t find_selection_2d(Selection selection,histo_t x,histo_t y);

void find_selection_range_1d(Selection selection,histo_t *min,histo_t *max,char *type);

void find_selection_range_2d(Selection selection,histo_t *min,histo_t *max,char *type);

_Bool set_precision_integration(Precision *precision,size_t n,histo_t min,histo_t max,char *integration,const Precision *precision_default);

void init_integration(Integration *integration,Precision *precision);

void free_integration(Integration *integration);

void legendre(histo_t dist,histo_t leg[],MULTI_TYPE type);

void timer(size_t i);

//Correlators

void set_precision(char* type,size_t n_,histo_t min_,histo_t max_,char* integration_);

void cross_2pcf_multi(Selection angular,Selection* radials,Selection window,Pole pole,LOS los);

void cross_3pcf_multi_double_los(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_2pcf_multi_radial(Selection angular,Selection* radials,Selection window,Pole pole,LOS los);

void cross_3pcf_multi_radial_double_los(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_radial_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_2pcf_multi_angular(Selection* radials,Selection window,Pole pole,LOS los);

void cross_3pcf_multi_angular_double_los(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_angular_angular(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_radial_global(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_angular_global(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_angular_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

#endif //_WINDOW_COMMON_
