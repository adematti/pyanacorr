#ifndef _WINDOW_COMMON_
#define _WINDOW_COMMON_

//General-purpose functions

//void unravel_index(size_t ind,size_t* shape,size_t n_dim,size_t* index);

histo_t find_selection_1d(Selection selection,histo_t x,INTERPOL interpol);

histo_t find_selection_2d(Selection selection,histo_t x,histo_t y);

void legendre(histo_t dist,histo_t leg[],MULTI_TYPE type);

void timer(size_t i);

//Correlators

void cross_2pcf_multi(Selection angular,Selection* radials,Selection window,Pole pole,LOS los);

void cross_2pcf_multi_radial(Selection angular,Selection* radials,Selection window,Pole pole,LOS los);

void cross_2pcf_multi_angular(Selection* radials,Selection window,Pole pole,LOS los);

void cross_3pcf_multi(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_3pcf_multi_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_3pcf_multi_angular(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_radial_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_radial_global(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_angular_angular(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_angular_global(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

void cross_4pcf_multi_angular_radial(Selection angular,Selection* radials,Selection window,Pole* poles,LOS* los);

#endif //_WINDOW_COMMON_
