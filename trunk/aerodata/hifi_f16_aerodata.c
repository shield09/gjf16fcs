/*---------------------------------------------------------------------- */
/*                                                                       */
/* Hifi aerodata                                                         */
/* taken from Richard S. Russell's F-16 model                            */
/*                                                                       */
/*---------------------------------------------------------------------- */

double	*getALPHA1(){
FILE *fp = fopen("aerodata/ALPHA1.dat","r");
int i;
double *alpha1,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file ALPHA1.dat");

alpha1 = doubleVector(20);

for(i=0;i<20;i++){
	fscanf(fp,"%lf",&data);
	alpha1[i] = data;
	}
fclose(fp);
return(alpha1);
}

double	*getALPHA2(){
FILE *fp = fopen("aerodata/ALPHA2.dat","r");
int i;
double *alpha2,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file ALPHA2.dat");

alpha2 = doubleVector(14);

for(i=0;i<14;i++){
	fscanf(fp,"%lf",&data);
	alpha2[i] = data;
	}
fclose(fp);
return(alpha2);
}

double	*getBETA1(){
FILE *fp = fopen("aerodata/BETA1.dat","r");
int i;
double *beta1,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file BETA1.dat");

beta1 = doubleVector(19);

for(i=0;i<19;i++){
	fscanf(fp,"%lf",&data);
	beta1[i] = data;
	}
fclose(fp);
return(beta1);
}

double	*getDH1(){
FILE *fp = fopen("aerodata/DH1.dat","r");
int i;
double *dh1,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file DH1.dat");

dh1 = doubleVector(5);

for(i=0;i<5;i++){
	fscanf(fp,"%lf",&data);
	dh1[i] = data;
	}
fclose(fp);
return(dh1);
}

double	*getDH2(){
FILE *fp = fopen("aerodata/DH2.dat","r");
int i;
double *dh2,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file DH2.dat");

dh2 = doubleVector(3);

for(i=0;i<3;i++){
	fscanf(fp,"%lf",&data);
	dh2[i] = data;
	}
fclose(fp);
return(dh2);
}

double	*getDH3(){
FILE *fp = fopen("aerodata/DH3.dat","r");
int i;
double *dh3,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file DH3.dat");

dh3 = doubleVector(7);

for(i=0;i<7;i++){
	fscanf(fp,"%lf",&data);
	dh3[i] = data;
	}
fclose(fp);
return(dh3);
}

double _Cx(double alpha_in,double beta_in,double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1900;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		ndinfo.nPoints[2] = 5; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH1();
		fp = fopen("aerodata/f16CX.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CX.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
	x[2] = dele;

    return interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cz(double alpha_in,double beta_in, double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; /* alpha_in,beta_in,dele */
	double x[3];	/* Number of dimension */

	FILESIZE = 1900;	/* There are 1900 elements in the 20x19x5 3D array */

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); /* There are 1900 elements */
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	/* alpha_in npoints */
		ndinfo.nPoints[1] = 19; /* beta_in npoints  */
		ndinfo.nPoints[2] = 5;  /* dele npoints  */
		X = (double **) malloc(nDimension*sizeof(double*));
		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH1();

		fp = fopen("aerodata/f16CZ.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CZ.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}
	x[0] = alpha_in;
	x[1] = beta_in;
	x[2] = dele;
    return interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cm(double alpha_in,double beta_in,double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1900;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		ndinfo.nPoints[2] = 5; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH1();
		fp = fopen("aerodata/f16Cm.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cm.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
	x[2] = dele;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cy(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CY.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CY.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cn(double alpha_in, double beta_in, double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1140;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		ndinfo.nPoints[2] = 3;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH2();
		fp = fopen("aerodata/f16Cn.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cn.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
	x[2] = dele;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl(double alpha_in, double beta_in,double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1140;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		ndinfo.nPoints[2] = 3;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH2();
		fp = fopen("aerodata/f16Cl.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cl.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
	x[2] = dele;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cx_lef(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CX_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CX_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cz_lef(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CZ_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CZ_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cm_lef(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cm_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cm_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cy_lef(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CY_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CY_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cn_lef(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cn_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cn_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_lef(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; /* alpha_in,beta_in*/
	double x[2];	/* Number of dimension */
	FILESIZE = 266;	/* There are 266 elements in the 14x19 2D array */

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	/* alpha_in npoints */
		ndinfo.nPoints[1] = 19; /* beta_in npoints  */
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cl_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cl_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _CXq(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16CXq.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CXq.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CZq(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16CZq.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CZq.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CMq(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16Cmq.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cmq.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CYp(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16CYp.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CYp.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CYr(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16CYr.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CYr.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CNr(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16Cnr.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cnr.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CNp(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16Cnp.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cnp.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CLp(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16Clp.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Clp.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CLr(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16Clr.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Clr.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CXq_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCXq_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCXq_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CYr_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCYr_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCYr_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CYp_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCYp_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCYp_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CZq_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCZq_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCZq_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CLr_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dClr_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dClr_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CLp_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dClp_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dClp_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CMq_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCmq_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCmq_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CNr_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCnr_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCnr_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CNp_lef(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/f16dCnp_lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCnp_lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cy_r30(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CY_dr30.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CY_dr30.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cn_r30(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cn_dr30.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cn_dr30.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_r30(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cl_dr30.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cl_dr30.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cy_a20(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CY_da20.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CY_da20.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cy_a20_lef(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CY_da20lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CY_da20lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cn_a20(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cn_da20.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cn_da20.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cn_a20_lef(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cn_da20lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cn_da20lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_a20(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cl_da20.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cl_da20.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}
	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_a20_lef(double alpha_in, double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cl_da20lef.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cl_da20lef.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CNbeta(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16dCnbeta.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCnbeta.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CLbeta(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16dClbeta.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dClbeta.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_Cm(double alpha_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/f16dCm.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCm.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _eta_el(double el){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 5;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 5;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getDH1();
		fp = fopen("aerodata/ETA_DH1.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file ETA_DH1.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = el;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */

double _delta_Cm_ds(double alpha_in, double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 140;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 7;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getDH3();
		fp = fopen("aerodata/f16dCm_ds.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16dCm_ds.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}
	x[0] = alpha_in;
	x[1] = dele;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */

//EXTRA tables for control law

double _Cxde(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CXde.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CXde.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */

double _Czde(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16CZde.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16CZde.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */

double _Cmde(double alpha_in,double beta_in){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/f16Cmde.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file f16Cmde.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha_in;
	x[1] = beta_in;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


