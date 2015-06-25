BEGIN {# read all of the cloud LWC file				   \
    getline; getline; nx=$1; ny=$2; nz=$3; getline; dx=$1; dy=$2;  \
    getline; for (k=1; k<=nz; k++) zc[k]=$k;			   \
    getline; for (k=1; k<=nz; k++) Tc[k]=$k;			   \
    while (getline > 0) {lwc[$1,$2,$3]=$4; re[$1,$2,$3]=$5;}	   \
    for (k=1; k<=nz; k++) {if (zc[k]<ztopcld) nzc=k};			\
    # read in aerosol profile						\
    while ((getline<"aero.t")>0) {na++; za[na]=$1; massa[na]=$2; reffa[na]=$3;}	\
    # add aerosol levels below cloud domain				\
    j=1;								\
    while (za[j] < zc[1]) {						\
	z[j]=za[j]; T[j]=Tc[1]+6.5*(zc[1]-z[j]);			\
	lwca[j]=massa[j]; rea[j]=reffa[j]; j++;}			\
    # interpolate aerosol properties to cloud grid			\
    i=j; kcb=j; kct=kcb+nzc-1;						\
    for (k=1; k<=nzc; k++) {						\
	if (zc[k]>=za[j+1]) j++; if (j>na-1) j=na-1;			\
	f=(zc[k]-za[j])/(za[j+1]-za[j]);				\
	z[i]=zc[k]; T[i]=Tc[k];						\
	lwca[i]=(1-f)*massa[j]+f*massa[j+1];				\
	rea[i]=(1-f)*reffa[j]+f*reffa[j+1];				\
	i++;}								\
    # add aerosol levels above cloud domain				\
    j1=j+1;								\
    for (j=j1; j<=na; j++) {						\
	z[i]=za[j]; T[i]=T[nzc]-6.5*(z[i]-z[nzc]);			\
	lwca[i]=massa[j]; rea[i]=reffa[j];				\
	i++;}  nz=i-1;							\
    # write header of particle file					\
    print "3   propgen particle file";					\
    print nx,ny,nz; print dx,dy;					\
    for (k=1; k<=nz; k++) {printf " %6.4f",z[k];} printf "\n";		\
    for (k=1; k<=nz; k++) {printf " %6.2f",T[k];} printf "\n";		\
    # write rest of particle file					\
    for (i=1; i<=nx; i++) for (j=1; j<=ny; j++) for (k=1; k<=nz; k++) {	\
		if ((k>=kcb)&&(k<=kct)) {kc=k-kcb+1;			\
		    if (lwc[i,j,kc]>0) printf "%3d %3d %3d %1d %1d %6.4f %5.2f %1d %7.5f %5.3f\n", \
					   i,j,k,2,1,lwc[i,j,kc],re[i,j,kc],2,lwca[k],rea[k]; }	\
		if (lwc[i,j,kc]==0) printf "%3d %3d %3d %1d %1d %7.5f %5.3f\n",	\
					i,j,k,1,2,lwca[k],rea[k];} }
