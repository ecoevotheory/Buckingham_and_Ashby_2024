/***********************************************************************************************************
 * [t,SJ,SA,IJ,IA,EQFLAG] = LL_eco_dynamics_function(t_max,a0,g,q,c1a,c2a,rJ,rA,tolJ_current,alphaJ_current,f,bJ,bA,beta0,gamma,eqtol,init_pop,strain_totalH,strain_totalP)
 ***********************************************************************************************************/

/* Compile in Matlab using mex LL_eco_dynamics_function.c */

#include <mex.h>
#include <math.h>


/***********************************
 * Constant parameter values
 ***********************************/
#define MAXSTEPS 1e6 /* Maximum number of steps for ODE solver */
#define INTERVAL 1e2 /* Check if the system is close to equilibrium */
#define EPS 1e-6 /* ODE solver tolerance */
#define TINY 1e-30 /* Constant value for solver */
#define TINY2 1e-30 /* Constant value for solver to avoid tolerance issues */
/* RK solver parameters */
#define b21 0.2
#define b31 3.0/40.0
#define b32 9.0/40.0
#define b41 0.3
#define b42 -0.9
#define b43 1.2
#define b51 -11.0/54.0
#define b52 2.5
#define b53 -70.0/27.0
#define b54 35.0/27.0
#define b61 1631.0/55296
#define b62 175.0/512.0
#define b63 575.0/13824.0
#define b64 44275.0/110592
#define b65 253.0/4096.0
#define c1 37.0/378.0
#define c3 250.0/621.0
#define c4 125.0/594.0
#define c6 512.0/1771.0
#define dc5 -277.00/14336

/*************************************
 * Define structure for model parameters
 *************************************/
struct PARAM{
    double t_max;
    double a0;
    double g;
    double q;
    double c1a;
    double c2a; 
    double rJ;
    double rA;
    double f;
    double bJ;
    double bA;
    double beta0;
    double gamma;
    double eqtol;
    int strain_totalH;
    int strain_totalP;
};

/*************************************
 * Function prototypes
 *************************************/
int my_rungkut (double *T, double *SJ_out, double *SA_out, double *IJ_out, double *IA_out, double *EQFLAG, double *init_pop, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p);
void rkqs(double *SJ, double *SA, double *IJ, double *IA,  double *DSJDT, double *DSADT, double *DIJDT, double *DIADT,  double *h, double *hnext, double *SJ_SCALE, double *SA_SCALE, double *IJ_SCALE, double *IA_SCALE, double *tolJ,double *alphaJ, double *alphaJpower, struct PARAM *p);
void rkck(double *SJ, double *SA,  double *IJ, double *IA,  double *DSJDT, double *DSADT,  double *DIJDT, double *DIADT, double *SJout, double *SAout,  double *IJout, double *IAout, double *SJerr, double *SAerr,  double *IJerr, double *IAerr,  double h, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p);
void dynamic(double *SJ, double *SA, double *IJ, double *IA, double *DSJDT, double *DSADT, double *DIJDT, double *DIADT, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p);
double FMAX(double, double);
double FMIN(double, double);

/*************************************
 * Main function
 *************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *T, *SJ, *SA, *IJ, *IA, *EQFLAG, *init_pop, *tolJ, *alphaJ, *alphaJpower, *parameter;
    double *t_temp, *SJ_temp, *SA_temp, *IJ_temp, *IA_temp;
    int i, j, k, colLen, maxsteps;
    struct PARAM p;


   /* Allocate inputs */
    if(nrhs!=19){
       mexErrMsgTxt("Incorrect number of input arguments!\n");
    }

    else{
        parameter= mxGetPr(prhs[0]);
        p.t_max= *parameter;
        parameter= mxGetPr(prhs[1]);
        p.a0= *parameter;
        parameter= mxGetPr(prhs[2]);
        p.g= *parameter;
	parameter= mxGetPr(prhs[3]);
        p.q= *parameter;
        parameter= mxGetPr(prhs[4]);
        p.c1a= *parameter;
        parameter= mxGetPr(prhs[5]);        
        p.c2a= *parameter;
        parameter= mxGetPr(prhs[6]);
        p.rJ= *parameter;
	parameter= mxGetPr(prhs[7]);
        p.rA= *parameter;
        tolJ= mxGetPr(prhs[8]);
        alphaJ= mxGetPr(prhs[9]);
        parameter= mxGetPr(prhs[10]);
        p.f= *parameter;
        parameter= mxGetPr(prhs[11]);
        p.bJ= *parameter;
        parameter= mxGetPr(prhs[12]);
        p.bA= *parameter;
        parameter= mxGetPr(prhs[13]);
        p.beta0= *parameter;
	parameter= mxGetPr(prhs[14]);
        p.gamma= *parameter;
        parameter= mxGetPr(prhs[15]);
        p.eqtol= *parameter;
        init_pop= mxGetPr(prhs[16]);    
        parameter= mxGetPr(prhs[17]);
        p.strain_totalH= (int)*parameter;
        parameter= mxGetPr(prhs[18]);
        p.strain_totalP= (int)*parameter;
    }

    maxsteps = (int)MAXSTEPS;
    
    /* Allocate memory */
    t_temp = malloc(maxsteps*sizeof(double));
    SJ_temp = malloc(maxsteps*(p.strain_totalH)*(p.strain_totalP)*sizeof(double));
    IJ_temp = malloc(maxsteps*(p.strain_totalH)*(p.strain_totalP)*sizeof(double));
    SA_temp = malloc(maxsteps*(p.strain_totalH)*(p.strain_totalP)*sizeof(double));
    IA_temp = malloc(maxsteps*(p.strain_totalH)*(p.strain_totalP)*sizeof(double));
    
    /* Initialise this output */           
    plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    EQFLAG = mxGetPr(plhs[5]);
    EQFLAG[0] = 0;


    /* Call ODE solver */
    colLen = my_rungkut(t_temp, SJ_temp, SA_temp, IJ_temp, IA_temp,EQFLAG, init_pop, tolJ, alphaJ, alphaJpower, &p); 

    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, p.strain_totalH*p.strain_totalP, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(colLen, p.strain_totalH*p.strain_totalP, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(colLen, p.strain_totalH*p.strain_totalP, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(colLen, p.strain_totalH*p.strain_totalP, mxREAL);
    T = mxGetPr(plhs[0]);
    SJ = mxGetPr(plhs[1]);
    SA = mxGetPr(plhs[2]);
    IJ = mxGetPr(plhs[3]);
    IA = mxGetPr(plhs[4]);

    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        T[i] = t_temp[i];
        for (j=0;j<p.strain_totalH;j++) {
	    for (k=0;k<p.strain_totalP;k++) {
		SJ[i+j*colLen+k*colLen*p.strain_totalH]=SJ_temp[i+j*maxsteps+k*maxsteps*p.strain_totalH];
		SA[i+j*colLen+k*colLen*p.strain_totalH]=SA_temp[i+j*maxsteps+k*maxsteps*p.strain_totalH];
		IJ[i+j*colLen+k*colLen*p.strain_totalH]=IJ_temp[i+j*maxsteps+k*maxsteps*p.strain_totalH];
		IA[i+j*colLen+k*colLen*p.strain_totalH]=IA_temp[i+j*maxsteps+k*maxsteps*p.strain_totalH];
	    }
        }
    }
    
    /* Free memory */
    free(t_temp);
    free(SJ_temp);
    free(IJ_temp); 
    free(SA_temp);
    free(IA_temp); 
    
    return;

}

/*****************************************
 * ODE solver
 ****************************************/
int my_rungkut (double *T, double *SJ_out, double *SA_out, double *IJ_out, double *IA_out,  double *EQFLAG, double *init_pop, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p){
    
    double *SJ, *SA,  *IJ, *IA, *DSJDT, *DSADT, *DIJDT, *DIADT,  *SJ_SCALE, *SA_SCALE, *IJ_SCALE, *IA_SCALE;
    double *SJMIN, *SJMAX, *SAMIN, *SAMAX, *IJMIN, *IJMAX, *IAMIN, *IAMAX, hnext[1], h[1];
    double t, nextcheck;
    int i, j, exitflag, count, maxsteps;
    
    /* Allocate memory */
    SJ = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SA = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJ = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IA = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    DSJDT = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    DSADT = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    DIJDT = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    DIADT = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJ_SCALE = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SA_SCALE = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJ_SCALE = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IA_SCALE = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJMIN = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    SJMAX = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    SAMIN = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    SAMAX = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    IJMIN = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    IJMAX = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    IAMIN = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    IAMAX = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    
    /* Other parameters */
    exitflag = 1;
    count=0;
    h[0] = 1e-3;
    hnext[0] = 1e-3;
    t=0;
    nextcheck = INTERVAL;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations */
    for(i=0;i<p->strain_totalH;i++){
	for (j=0;j<p->strain_totalP;j++){
            SJ[i+j*p->strain_totalH] = init_pop[4*(i+j*p->strain_totalH)];
            SA[i+j*p->strain_totalH] = init_pop[4*(i+j*p->strain_totalH)+1];
            IJ[i+j*p->strain_totalH] = init_pop[4*(i+j*p->strain_totalH)+2];
            IA[i+j*p->strain_totalH] = init_pop[4*(i+j*p->strain_totalH)+3];
	}
    }
    
    /* Initialise equilibrium arrays */
    for(i=0;i<p->strain_totalH;i++){
	for (j=0;j<p->strain_totalP;j++){
            SJMIN[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH];
            SJMAX[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH];
            SAMIN[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH];
            SAMAX[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH];
            IJMIN[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH];
            IJMAX[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH];
            IAMIN[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH];
            IAMAX[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH];
	}
    }
    
    /* Update output */
    T[0]=t;
    for (i=0; i<p->strain_totalH; i++) {
	for (j=0; j<p->strain_totalP;j++){
            SJ_out[(i+j*p->strain_totalH)*maxsteps] = SJ[i+j*p->strain_totalH];
            SA_out[(i+j*p->strain_totalH)*maxsteps] = SA[i+j*p->strain_totalH];
            IJ_out[(i+j*p->strain_totalH)*maxsteps] = IJ[i+j*p->strain_totalH];
            IA_out[(i+j*p->strain_totalH)*maxsteps] = IA[i+j*p->strain_totalH];
	}
    }

    
    /* Main loop: */
    do{
        /* This ensures the final step lands us on the final time point */
        if(1.1*hnext[0]>(p->t_max-t)){
            hnext[0] = p->t_max-t;
            h[0] = p->t_max-t;
            t=p->t_max;
            exitflag=0;
        }
        else{
            h[0] = hnext[0];
            t+=h[0];
        }
        if(t>=p->t_max) {
            t=p->t_max;
            exitflag=0;
        }
        /* This is where the equations are first solved */

        dynamic(SJ, SA,  IJ, IA,  DSJDT, DSADT,  DIJDT, DIADT, tolJ, alphaJ, alphaJpower, p);

        
        /* Adjust the step size to maintain accuracy */
        for (i=0; i<p->strain_totalH; i++){
	    for (j=0; j<p->strain_totalP; j++){
                SJ[i+j*p->strain_totalH] = FMAX(SJ[i+j*p->strain_totalH],0);
                SJ_SCALE[i+j*p->strain_totalH]=fabs(SJ[i+j*p->strain_totalH])+fabs(DSJDT[i+j*p->strain_totalH]*(*h))+TINY;
                SA[i+j*p->strain_totalH] = FMAX(SA[i+j*p->strain_totalH],0);
                SA_SCALE[i+j*p->strain_totalH]=fabs(SA[i+j*p->strain_totalH])+fabs(DSADT[i+j*p->strain_totalH]*(*h))+TINY;
                IJ[i+j*p->strain_totalH] = FMAX(IJ[i+j*p->strain_totalH],0);
                IJ_SCALE[i+j*p->strain_totalH]=fabs(IJ[i+j*p->strain_totalH])+fabs(DIJDT[i+j*p->strain_totalH]*(*h))+TINY;
                IA[i+j*p->strain_totalH] = FMAX(IA[i+j*p->strain_totalH],0);
                IA_SCALE[i+j*p->strain_totalH]=fabs(IA[i+j*p->strain_totalH])+fabs(DIADT[i+j*p->strain_totalH]*(*h))+TINY;

	    }
        }
       

        /* RK solver & adaptive step-size */

        rkqs(SJ, SA, IJ, IA, DSJDT, DSADT, DIJDT, DIADT, h, hnext, SJ_SCALE, SA_SCALE, IJ_SCALE, IA_SCALE, tolJ, alphaJ, alphaJpower, p); 

        /* Make sure nothing has gone negative */
        for (i=0; i<p->strain_totalH; i++){
	    for (j=0; j<p->strain_totalP; j++){
                SJ[i+j*p->strain_totalH] = FMAX(SJ[i+j*p->strain_totalH],0);
                SA[i+j*p->strain_totalH] = FMAX(SA[i+j*p->strain_totalH],0);
                IJ[i+j*p->strain_totalH] = FMAX(IJ[i+j*p->strain_totalH],0);
                IA[i+j*p->strain_totalH] = FMAX(IA[i+j*p->strain_totalH],0);
	    }
        }
        
        /* Update output */
        count++;
        T[count] = t;
        for (i=0; i<p->strain_totalH; i++) {
	    for (j=0; j<p->strain_totalP; j++) {
                SJ_out[count + (i+j*p->strain_totalH)*maxsteps] = SJ[i+j*p->strain_totalH];
                SA_out[count + (i+j*p->strain_totalH)*maxsteps] = SA[i+j*p->strain_totalH];
                IJ_out[count + (i+j*p->strain_totalH)*maxsteps] = IJ[i+j*p->strain_totalH];
                IA_out[count + (i+j*p->strain_totalH)*maxsteps] = IA[i+j*p->strain_totalH];
	    }
        }


        /* For equilibrium check */
        for (i=0; i<p->strain_totalH; i++){
	    for (j=0; j<p->strain_totalP; j++) {
                SJMIN[i+j*p->strain_totalH] = FMIN(SJMIN[i+j*p->strain_totalH],SJ[i+j*p->strain_totalH]);
                SJMAX[i+j*p->strain_totalH] = FMAX(SJMAX[i+j*p->strain_totalH],SJ[i+j*p->strain_totalH]);
                SAMIN[i+j*p->strain_totalH] = FMIN(SAMIN[i+j*p->strain_totalH],SA[i+j*p->strain_totalH]);
                SAMAX[i+j*p->strain_totalH] = FMAX(SAMAX[i+j*p->strain_totalH],SA[i+j*p->strain_totalH]);
                IJMIN[i+j*p->strain_totalH] = FMIN(IJMIN[i+j*p->strain_totalH],IJ[i+j*p->strain_totalH]);
                IJMAX[i+j*p->strain_totalH] = FMAX(IJMAX[i+j*p->strain_totalH],IJ[i+j*p->strain_totalH]);
                IAMIN[i+j*p->strain_totalH] = FMIN(IAMIN[i+j*p->strain_totalH],IA[i+j*p->strain_totalH]);
                IAMAX[i+j*p->strain_totalH] = FMAX(IAMAX[i+j*p->strain_totalH],IA[i+j*p->strain_totalH]);
	    }

        }


        /* Check if we're close to equilibrium */
        if(t>nextcheck){
            exitflag = 0;
            for (i=0; i<p->strain_totalH; i++){
		for (j=0; j<p->strain_totalP; j++){
                    if(fabs(SJMAX[i+j*p->strain_totalH]-SJMIN[i+j*p->strain_totalH])>p->eqtol || fabs(SAMAX[i+j*p->strain_totalH]-SAMIN[i+j*p->strain_totalH])>p->eqtol || fabs(IJMAX[i+j*p->strain_totalH]-IJMIN[i+j*p->strain_totalH])>p->eqtol || fabs(IAMAX[i+j*p->strain_totalH]-IAMIN[i+j*p->strain_totalH])>p->eqtol ){
                        exitflag = 1;
                        break;
		    }
                }
            }
            /* If close to equilibrium, then break */
            if(exitflag==0){
                t=p->t_max; 
                T[count] = t; 
                EQFLAG[0] = 1; 
                break; 
            } 
            
            /* If not, then reset min/max values for each class */
            nextcheck+=INTERVAL;
            for (i=0; i<p->strain_totalH; i++){
		for (j=0; j<p->strain_totalP; j++){
                    SJMIN[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH];
                    SJMAX[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH];
                    SAMIN[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH];
                    SAMAX[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH];
                    IJMIN[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH];
                    IJMAX[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH];
                    IAMIN[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH];
                    IAMAX[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH];
		}
            }

        }
    }while(count<(maxsteps-1) && t<=p->t_max && exitflag);
    count++;

    
    
    /* Free memory */
    free(SJ);
    free(SA);
    free(IJ);
    free(IA);
    free(DSJDT);
    free(DSADT);
    free(DIJDT);
    free(DIADT);
    free(SJ_SCALE);
    free(SA_SCALE);
    free(IJ_SCALE);
    free(IA_SCALE);
    free(SJMIN);
    free(SAMIN);
    free(IJMIN);
    free(IAMIN);
    free(SJMAX);
    free(SAMAX);
    free(IJMAX);
    free(IAMAX);
    
    return count;
}

/***************************************
 * This generates the adaptive step-size
 **************************************/
void rkqs(double *SJ, double *SA, double *IJ, double *IA,   double *DSJDT, double *DSADT,  double *DIJDT, double *DIADT,  double *h, double *hnext, double *SJ_SCALE, double *SA_SCALE, double *IJ_SCALE, double *IA_SCALE, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p)
{
    double *SJ_temp, *SA_temp,  *IJ_temp, *IA_temp, *SJ_err, *SA_err, *IJ_err, *IA_err;
    double htemp, errmax;
    int i, j, count;
    
    /* Allocate memory */
    SJ_temp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SA_temp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJ_temp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IA_temp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJ_err = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SA_err = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJ_err = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IA_err = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    
    count = 0;
    for(;;)
    {
        rkck(SJ, SA, IJ, IA, DSJDT, DSADT, DIJDT, DIADT, SJ_temp, SA_temp,  IJ_temp, IA_temp, SJ_err, SA_err,  IJ_err, IA_err,*h, tolJ, alphaJ, alphaJpower, p); 
        
        errmax= 0.0;
        for(i=0;i<p->strain_totalH;i++){
	    for (j=0; j<p->strain_totalP; j++){

                errmax= FMAX(errmax, fabs(SJ_err[i+j*p->strain_totalH]/(SJ_SCALE[i+j*p->strain_totalH]))); 
                errmax= FMAX(errmax, fabs(SA_err[i+j*p->strain_totalH]/(SA_SCALE[i+j*p->strain_totalH])));   
                errmax= FMAX(errmax, fabs(IJ_err[i+j*p->strain_totalH]/(IJ_SCALE[i+j*p->strain_totalH])));   
                errmax= FMAX(errmax, fabs(IA_err[i+j*p->strain_totalH]/(IA_SCALE[i+j*p->strain_totalH])));     
	    }
        }

        errmax/= EPS;

        if(errmax<=1.0) break;
        htemp= 0.9*(*h)*pow(errmax, -0.25);
        *h= (*h>=0.0 ? FMAX(htemp, 0.1*(*h)) : FMIN(htemp, 0.1*(*h)));
        count++;
            

        if(count>1e4){
            printf("%f\n",errmax);
	    if(isnan(errmax)){
                mexErrMsgTxt("errmax is NAN\n");
            } 
            mexErrMsgTxt("stuck in loop!\n");
            break;
        }
    }    
    if(errmax > 1.89E-4) {
        *hnext= 0.9*(*h)*pow(errmax, -0.2);
    }
    else {
        *hnext= 5.0*(*h);
    }    
    *hnext = FMAX(*hnext, p->t_max/MAXSTEPS);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            SJ[i+j*p->strain_totalH] = SJ_temp[i+j*p->strain_totalH];
            SA[i+j*p->strain_totalH] = SA_temp[i+j*p->strain_totalH];
            IJ[i+j*p->strain_totalH] = IJ_temp[i+j*p->strain_totalH];
            IA[i+j*p->strain_totalH] = IA_temp[i+j*p->strain_totalH];
	}
    }
    
    /* Free memory */
    free(SJ_temp);
    free(IJ_temp);
    free(SJ_err);
    free(IJ_err);
    free(SA_temp);
    free(IA_temp);
    free(SA_err);
    free(IA_err);
}

/**************************************
 * Standard RK solver
 **************************************/
void rkck(double *SJ, double *SA, double *IJ, double *IA, double *DSJDT, double *DSADT, double *DIJDT, double *DIADT, double *SJout, double *SAout, double *IJout, double *IAout, double *SJerr, double *SAerr, double *IJerr, double *IAerr, double h, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p){
    int i, j;
    double *SJk1, *SJk2, *SJk3, *SJk4, *SJk5, *SJk6, *SJtemp;
    double *IJk1, *IJk2, *IJk3, *IJk4, *IJk5, *IJk6, *IJtemp;
    double *SAk1, *SAk2, *SAk3, *SAk4, *SAk5, *SAk6, *SAtemp;
    double *IAk1, *IAk2, *IAk3, *IAk4, *IAk5, *IAk6, *IAtemp;
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, 
            dc6=c6-0.25;
    
    /* Allocate memory */    
    SJk1 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJk2 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJk3 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJk4 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJk5 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJk6 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SJtemp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    
    SAk1 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SAk2 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SAk3 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SAk4 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SAk5 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SAk6 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SAtemp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));

    IJk1 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJk2 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJk3 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJk4 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJk5 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJk6 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IJtemp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));

    IAk1 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IAk2 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IAk3 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IAk4 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IAk5 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IAk6 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    IAtemp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            SJtemp[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH] + b21*h*DSJDT[i+j*p->strain_totalH];
            SAtemp[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH] + b21*h*DSADT[i+j*p->strain_totalH];
            IJtemp[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH] + b21*h*DIJDT[i+j*p->strain_totalH];
            IAtemp[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH] + b21*h*DIADT[i+j*p->strain_totalH];
	}
    }

    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk2, SAk2, IJk2, IAk2,  tolJ, alphaJ, alphaJpower, p);

    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            SJtemp[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH]+h*(b31*DSJDT[i+j*p->strain_totalH]+b32*SJk2[i+j*p->strain_totalH]);
	    SAtemp[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH]+h*(b31*DSADT[i+j*p->strain_totalH]+b32*SAk2[i+j*p->strain_totalH]);
            IJtemp[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH]+h*(b31*DIJDT[i+j*p->strain_totalH]+b32*IJk2[i+j*p->strain_totalH]);
            IAtemp[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH]+h*(b31*DIADT[i+j*p->strain_totalH]+b32*IAk2[i+j*p->strain_totalH]);
	}
    }    

    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk3, SAk3, IJk3, IAk3,  tolJ, alphaJ, alphaJpower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            SJtemp[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH]+h*(b41*DSJDT[i+j*p->strain_totalH]+b42*SJk2[i+j*p->strain_totalH]+b43*SJk3[i+j*p->strain_totalH]);
	    SAtemp[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH]+h*(b41*DSADT[i+j*p->strain_totalH]+b42*SAk2[i+j*p->strain_totalH]+b43*SAk3[i+j*p->strain_totalH]);
            IJtemp[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH]+h*(b41*DIJDT[i+j*p->strain_totalH]+b42*IJk2[i+j*p->strain_totalH]+b43*IJk3[i+j*p->strain_totalH]);
            IAtemp[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH]+h*(b41*DIADT[i+j*p->strain_totalH]+b42*IAk2[i+j*p->strain_totalH]+b43*IAk3[i+j*p->strain_totalH]);
	}
    }

    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk4, SAk4, IJk4, IAk4,  tolJ, alphaJ, alphaJpower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            SJtemp[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH]+h*(b51*DSJDT[i+j*p->strain_totalH]+b52*SJk2[i+j*p->strain_totalH]+b53*SJk3[i+j*p->strain_totalH]+b54*SJk4[i+j*p->strain_totalH]);
	    SAtemp[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH]+h*(b51*DSADT[i+j*p->strain_totalH]+b52*SAk2[i+j*p->strain_totalH]+b53*SAk3[i+j*p->strain_totalH]+b54*SAk4[i+j*p->strain_totalH]);
            IJtemp[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH]+h*(b51*DIJDT[i+j*p->strain_totalH]+b52*IJk2[i+j*p->strain_totalH]+b53*IJk3[i+j*p->strain_totalH]+b54*IJk4[i+j*p->strain_totalH]);
            IAtemp[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH]+h*(b51*DIADT[i+j*p->strain_totalH]+b52*IAk2[i+j*p->strain_totalH]+b53*IAk3[i+j*p->strain_totalH]+b54*IAk4[i+j*p->strain_totalH]);
	}
    }

    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk5, SAk5, IJk5, IAk5,  tolJ, alphaJ, alphaJpower, p); 
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){ 
            SJtemp[i+j*p->strain_totalH] = SJ[i+j*p->strain_totalH]+h*(b61*DSJDT[i+j*p->strain_totalH]+b62*SJk2[i+j*p->strain_totalH]+b63*SJk3[i+j*p->strain_totalH]+b64*SJk4[i+j*p->strain_totalH]+b65*SJk5[i+j*p->strain_totalH]);
	    SAtemp[i+j*p->strain_totalH] = SA[i+j*p->strain_totalH]+h*(b61*DSADT[i+j*p->strain_totalH]+b62*SAk2[i+j*p->strain_totalH]+b63*SAk3[i+j*p->strain_totalH]+b64*SAk4[i+j*p->strain_totalH]+b65*SAk5[i+j*p->strain_totalH]);
            IJtemp[i+j*p->strain_totalH] = IJ[i+j*p->strain_totalH]+h*(b61*DIJDT[i+j*p->strain_totalH]+b62*IJk2[i+j*p->strain_totalH]+b63*IJk3[i+j*p->strain_totalH]+b64*IJk4[i+j*p->strain_totalH]+b65*IJk5[i+j*p->strain_totalH]);
            IAtemp[i+j*p->strain_totalH] = IA[i+j*p->strain_totalH]+h*(b61*DIADT[i+j*p->strain_totalH]+b62*IAk2[i+j*p->strain_totalH]+b63*IAk3[i+j*p->strain_totalH]+b64*IAk4[i+j*p->strain_totalH]+b65*IAk5[i+j*p->strain_totalH]);
	}
    }

    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk6, SAk6, IJk6, IAk6,  tolJ, alphaJ, alphaJpower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            SJout[i+j*p->strain_totalH]= SJ[i+j*p->strain_totalH]+h*(c1*DSJDT[i+j*p->strain_totalH]+c3*SJk3[i+j*p->strain_totalH]+c4*SJk4[i+j*p->strain_totalH]+c6*SJk6[i+j*p->strain_totalH]);
            SJerr[i+j*p->strain_totalH]= h*(dc1*DSJDT[i+j*p->strain_totalH]+dc3*SJk3[i+j*p->strain_totalH]+dc4*SJk4[i+j*p->strain_totalH]+dc5*SJk5[i+j*p->strain_totalH]+dc6*SJk6[i+j*p->strain_totalH]);
            SAout[i+j*p->strain_totalH]= SA[i+j*p->strain_totalH]+h*(c1*DSADT[i+j*p->strain_totalH]+c3*SAk3[i+j*p->strain_totalH]+c4*SAk4[i+j*p->strain_totalH]+c6*SAk6[i+j*p->strain_totalH]);
            SAerr[i+j*p->strain_totalH]= h*(dc1*DSADT[i+j*p->strain_totalH]+dc3*SAk3[i+j*p->strain_totalH]+dc4*SAk4[i+j*p->strain_totalH]+dc5*SAk5[i+j*p->strain_totalH]+dc6*SAk6[i+j*p->strain_totalH]);
            IJout[i+j*p->strain_totalH]= IJ[i+j*p->strain_totalH]+h*(c1*DIJDT[i+j*p->strain_totalH]+c3*IJk3[i+j*p->strain_totalH]+c4*IJk4[i+j*p->strain_totalH]+c6*IJk6[i+j*p->strain_totalH]);
            IJerr[i+j*p->strain_totalH]= h*(dc1*DIJDT[i+j*p->strain_totalH]+dc3*IJk3[i+j*p->strain_totalH]+dc4*IJk4[i+j*p->strain_totalH]+dc5*IJk5[i+j*p->strain_totalH]+dc6*IJk6[i+j*p->strain_totalH]);
            IAout[i+j*p->strain_totalH]= IA[i+j*p->strain_totalH]+h*(c1*DIADT[i+j*p->strain_totalH]+c3*IAk3[i+j*p->strain_totalH]+c4*IAk4[i+j*p->strain_totalH]+c6*IAk6[i+j*p->strain_totalH]);
            IAerr[i+j*p->strain_totalH]= h*(dc1*DIADT[i+j*p->strain_totalH]+dc3*IAk3[i+j*p->strain_totalH]+dc4*IAk4[i+j*p->strain_totalH]+dc5*IAk5[i+j*p->strain_totalH]+dc6*IAk6[i+j*p->strain_totalH]);
	}
    }

    /* Free memory */
    free(SJk1);
    free(SJk2);
    free(SJk3);
    free(SJk4);
    free(SJk5);
    free(SJk6);
    free(SJtemp);
    free(IJk1);
    free(IJk2);
    free(IJk3);
    free(IJk4);
    free(IJk5);
    free(IJk6);
    free(IJtemp);
    free(SAk1);
    free(SAk2);
    free(SAk3);
    free(SAk4);
    free(SAk5);
    free(SAk6);
    free(SAtemp);
    free(IAk1);
    free(IAk2);
    free(IAk3);
    free(IAk4);
    free(IAk5);
    free(IAk6);
    free(IAtemp);
}

/**************************************
 * Population and evolutionary dynamics
 **************************************/
void dynamic(double *SJ, double *SA, double *IJ, double *IA,  double *DSJDT, double *DSADT, double *DIJDT, double *DIADT, double *tolJ, double *alphaJ, double *alphaJpower, struct PARAM *p){
    

    int i, j;

    double *IAallparasites, *IJallparasites, *totalparasites, N, betaallinfecteds;

    IAallparasites = malloc(p->strain_totalH*sizeof(double));
    IJallparasites = malloc(p->strain_totalH*sizeof(double));
    totalparasites = malloc(p->strain_totalP*sizeof(double));

    /* Population sums */
    N = 0;
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){
	    N = N + IJ[i+j*p->strain_totalH] + SJ[i+j*p->strain_totalH]+IA[i+j*p->strain_totalH] + SA[i+j*p->strain_totalH];
	}
    }

    betaallinfecteds = 0;
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){

	    betaallinfecteds = betaallinfecteds + sqrt(alphaJ[j])*IJ[i+j*p->strain_totalH] + sqrt(alphaJ[j])*IA[i+j*p->strain_totalH];

	}
    }

    for(i=0;i<p->strain_totalH;i++){ 
	    IAallparasites[i]=0;
    }
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){
	    IAallparasites[i]=IAallparasites[i]+IA[i+j*p->strain_totalH];
	}
    }

    for(i=0;i<p->strain_totalH;i++){ 
	    IJallparasites[i]=0;
    }
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){
	    IJallparasites[i]=IJallparasites[i]+IJ[i+j*p->strain_totalH];
	}
    }

    for(j=0;j<p->strain_totalP;j++){ 
	    totalparasites[j]=0;
    }
    for(j=0;j<p->strain_totalP;j++){ 
        for (i=0; i<p->strain_totalH; i++){
	    totalparasites[j]=totalparasites[j]+sqrt(alphaJ[j])*IJ[i+j*p->strain_totalH]+sqrt(alphaJ[j])*IA[i+j*p->strain_totalH];
	}
    }

    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
	    DSJDT[i+j*p->strain_totalH] = 0;
            DSADT[i+j*p->strain_totalH] = 0;
	    DIJDT[i+j*p->strain_totalH] = p->beta0*(1-p->rJ)*totalparasites[j]*SJ[i] -(p->g+p->gamma+p->bJ+alphaJ[j]*(1-tolJ[i]))*IJ[i+j*p->strain_totalH];
	    DIADT[i+j*p->strain_totalH] = p->g*IJ[i+j*p->strain_totalH]+p->beta0*(1-p->rA)*totalparasites[j]*SA[i]-(alphaJ[j]*(1-tolJ[i])+p->gamma+p->bA)*IA[i+j*p->strain_totalH];
	}

 	DSJDT[i] = p->a0*(1-p->c1a*(1-exp(p->c2a*tolJ[i]))/(1-exp(p->c2a)))*(1-p->q*N)*(SA[i]+p->f*IAallparasites[i])-(p->bJ+p->g+p->beta0*(1-p->rJ)*betaallinfecteds)*SJ[i]+p->gamma*IJallparasites[i];
        DSADT[i] = p->g*SJ[i]-(p->bA+p->beta0*(1-p->rA)*betaallinfecteds)*SA[i]+p->gamma*IAallparasites[i];

	  
    }
}


/***************************************
 * Return maximum of two inputs
 ***************************************/
double FMAX(double l, double r)
{
    if(l>r)return l;
    else   return r;
}

/***************************************
 * Return minimum of two inputs
 ***************************************/
double FMIN(double l, double r)
{
    if(l<r)return l;
    else   return r;
}
