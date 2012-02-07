/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*               6 DOF F-16 FIGHTER AIRCRAFT DYNAMICS                    */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*                                                                       */
/* Based on the F-16 model created by R. S. Russel in                    */
/* "Nonlinear F-16 simulation using Simulink and Matlab"                 */   
/* 2003 University of Minnesota and on the F-16 model                    */
/* created by Ying Huo in "Model of F-16 Fighter Aircraft"               */
/*                                                                       */
/* Aerodynamic data and the engine model have been obtained from         */
/* "NASA Technical Paper 1538" by Nguyen et al. 1979                     */   
/*                                                                       */
/* File "F16_dyn.c"                                                      */
/* Version 1.2 by E.R. van Oort & L. Sonneveldt                          */
/* Created with MATLAB 7.5 R2007B                                        */
/* Feb, 2009                                                             */
/*                                                                       */
/* Notes:                                                                */
/* -All units are SI.                                                    */
/* -Quaternion transformations are used.                                 */
/* -Flag is used to select between hifi and lofi aerodynamic model.      */
/* -Hifi aerodata is obtained from "aerodata/hifi_f16_aerodata.c"        */
/* -Lofi aerodata is obtained from "aerodata/lofi_f16_aerodata.c"        */
/* -"aerodata/mexndinterp.c" is used for interpolation of the data.      */                                                                    
/* -"aerodata/engine_model.c" contains the engine model.                 */  
/* -"aerodata/ISA_atmos.c" contains the ISA atmosphere model.            */
/*                                                                       */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*  Used variables:														 */
/*																		 */
/*  Input variables: (port 1)   										 */
/*      dth         throttle setting (between 0 and 1)          [-]      */
/*      de          elevator deflection                         [rad]    */
/*      da          aileron deflection                          [rad]	 */
/*      dr          rudder deflection                           [rad]	 */
/*																		 */
/*  Additonal input variables: (port 2)									 */
/*      dlef        leading edge flap deflection                [rad]    */
/*																		 */
/*  Fidelity flag: (port 3)											     */
/*		fi_flag     hifi/lofi aerodynamic model selection flag  [-]      */
/*																		 */
/*  State variables:													 */
/*      Vt          total airspeed                              [m/s]	 */
/*      beta        angle of sideslip                           [rad]	 */
/*      alpha       angle of attack                             [rad]	 */
/*      q0          quaternion component                        [-]      */
/*      q1          quaternion component                        [-]      */
/*      q2          quaternion component                        [-]      */
/*      q3          quaternion component                        [-]      */
/*      p_body      roll angular rate                           [rad/s]	 */
/*      q_body      pitch angular rate                          [rad/s]	 */
/*      r_body      yaw angular rate                            [rad/s]	 */
/*      x_earth     x position                                  [m]		 */
/*      y_earth     y position                                  [m]		 */
/*      z_earth     z position                                  [m]		 */
/*      power       power level (0-100%)                        [-]      */
/*																		 */
/*  State derivatives:													 */
/*      Vt_dot      rate of change in total airspeed            [m/s^2]	 */
/*      beta_dot    rate of change in angle of sideslip         [rad/s]	 */
/*      alpha_dot   rate of change in angle of attack           [rad/s]	 */
/*      q0_dot      rate of change in quaternion component      [1/s]    */
/*      q1_dot      rate of change in quaternion component      [1/s]    */
/*      q2_dot      rate of change in quaternion component      [1/s]    */
/*      q3_dot      rate of change in quaternion component      [1/s]    */
/*      p_body_dot  rate of change in roll angular rate         [rad/s^2]*/
/*      q_body_dot  rate of change in pitch angular rate        [rad/s^2]*/
/*      r_body_dot  rate of change in yaw angular rate          [rad/s^2]*/
/*      x_earth_dot rate of change in x position                [m/s]	 */
/*      y_earth_dot rate of change in y position                [m/s]	 */
/*      z_earth_dot rate of change in z position                [m/s]    */
/*      pow_dot     rate of change in power level               [1/s]    */
/*																		 */
/*  Output variables:												     */
/*      Vt          total airspeed                              [m/s]	 */
/*      beta        angle of sideslip                           [rad]	 */
/*      alpha       angle of attack                             [rad]	 */
/*      q0          quaternion component                        [-]      */
/*      q1          quaternion component                        [-]      */
/*      q2          quaternion component                        [-]      */
/*      q3          quaternion component                        [-]      */
/*      p_body      roll angular rate                           [rad/s]	 */
/*      q_body      pitch angular rate                          [rad/s]	 */
/*      r_body      yaw angular rate                            [rad/s]	 */
/*      x_earth     x position                                  [m]		 */
/*      y_earth     y position                                  [m]		 */
/*      z_earth     z position                                  [m]		 */
/*      power       power level (0-100%)                        [-]      */
/*																		 */
/*	Additional Output:  												 */
/*		ny          normal acceleration y-axis                  [-]		 */
/*		nz          normal acceleration z-axis                  [-]		 */
/*																		 */
/*  Real work vector: (Persistant memory)                    			 */
/*      C1          coefficients used in moment equations       [-]      */   
/*      C2                                                               */      
/*      C3                                                               */              
/*      C4                                                               */              
/*      C5                                                               */             
/*      C6                                                               */              
/*      C7                                                               */              
/*      C8                                                               */              
/*      C9                                                               */              
/*      Xbar         total force in body fixed x-axis           [N]      */   
/*      Ybar         total force in body fixed y-axis           [N]      */  
/*      Zbar         total force in body fixed z-axis           [N]      */  
/*      Lbar         total moment in body fixed x-axis          [Nm]     */  
/*      Mbar         total moment in body fixed y-axis          [Nm]     */
/*      Nbar         total moment in body fixed z-axis          [Nm]     */
/*      u_body       velocity in body fixed x-axis              [m/s]	 */
/*      v_body       velocity in body fixed y-axis              [m/s]	 */
/*      w_body       velocity in body fixed z-axis              [m/s]	 */
/*      u_body_dot   rate of change in velocity x-axis          [m/s^2]  */
/*      v_body_dot   rate of change in velocity y-axis          [m/s^2]  */
/*      w_body_dot   rate of change in velocity z-axis          [m/s^2]  */
/*      qbar         dynamic pressure                           [N/m]    */
/*      Mach         Mach number                                [-]      */
/*      Thrust       Total engine thrust                        [-]      */
/*      ny           normal acceleration y-axis                 [-]		 */
/*		nz           normal acceleration z-axis                 [-]		 */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define S_FUNCTION_NAME  F16_dyn
#define S_FUNCTION_LEVEL 2
 
/* include files */
#include <math.h>
#include "simstruc.h"

#include "aerodata/mexndinterp.c"
#include "aerodata/hifi_f16_aerodata.c"     /* hifi lookup tables */
#include "aerodata/lofi_f16_aerodata.c"     /* lofi lookup tables */

#include "aerodata/ISA_atmos.c"             /* ISA atmosphere model */
#include "aerodata/engine_model.c"          /* engine model */

/* input port 1: control inputs */
#define dth             (*u[0])
#define de              (*u[1])
#define da              (*u[2])
#define dr              (*u[3])
/* input port 2: leading edge flap deflection */
#define dlef            (*u2[0])
/* input port 3: fidelity flag, 0 = lofi model, 1 = hifi model */
#define fi_flag         (*u3[0])

/* 14 States */
#define Vt          x[0]
#define beta        x[1]
#define alpha       x[2]

#define q0          x[3]
#define q1          x[4]
#define q2          x[5]
#define q3          x[6]

#define p_body      x[7]
#define q_body      x[8]
#define r_body      x[9]

#define x_earth     x[10]
#define y_earth     x[11]
#define z_earth     x[12]

#define power       x[13]

/* State derivatives */
#define Vt_dot          dx[0]
#define beta_dot        dx[1]
#define alpha_dot       dx[2]

#define q0_dot          dx[3]
#define q1_dot          dx[4]
#define q2_dot          dx[5]
#define q3_dot          dx[6]

#define p_body_dot      dx[7]
#define q_body_dot      dx[8]
#define r_body_dot      dx[9]

#define x_earth_dot     dx[10]
#define y_earth_dot     dx[11]
#define z_earth_dot     dx[12]
#define pow_dot         dx[13]

/* Work Variables */
#define C1              ssGetRWork(S)[0]
#define C2              ssGetRWork(S)[1]
#define C3              ssGetRWork(S)[2]
#define C4              ssGetRWork(S)[3]
#define C5              ssGetRWork(S)[4]
#define C6              ssGetRWork(S)[5]
#define C7              ssGetRWork(S)[6]
#define C8              ssGetRWork(S)[7]
#define C9              ssGetRWork(S)[8]
#define Xbar            ssGetRWork(S)[9]
#define Ybar            ssGetRWork(S)[10]
#define Zbar            ssGetRWork(S)[11]
#define Lbar            ssGetRWork(S)[12]
#define Mbar            ssGetRWork(S)[13]
#define Nbar            ssGetRWork(S)[14]
#define u_body          ssGetRWork(S)[15]
#define v_body          ssGetRWork(S)[16]
#define w_body          ssGetRWork(S)[17]
#define u_body_dot      ssGetRWork(S)[18]
#define v_body_dot      ssGetRWork(S)[19]
#define w_body_dot      ssGetRWork(S)[20]
#define qbar            ssGetRWork(S)[21]
#define Mach            ssGetRWork(S)[22]
#define Thrust          ssGetRWork(S)[23]
#define ny              ssGetRWork(S)[24]
#define nz              ssGetRWork(S)[25]

/* Aircraft Parameters */
#define mass        9295.44 /*assumed fixed*/
#define Ixx         12874.8
#define Iyy         75673.6
#define Izz         85552.1
#define Ixz         1331.4
#define Sref        27.87
#define bref        9.144
#define cref        3.45
#define xcg         0.3
#define xcgr        0.35
#define heng        216.9 /* engine angular momentum, assumed fixed */

/* Additional parameters */
#define rtd           57.29577951
#define dtr           0.017453293
#define Pi			  3.141592654

/*=============================*/
/* Function: mdlInitalizeSizes */
/*=============================*/
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 1);  /* Number of expected parameters */

    ssSetNumContStates(S, 14);
    ssSetNumDiscStates(S, 0);

    ssSetNumInputPorts(S, 3);
    ssSetInputPortWidth(S, 0, 4);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);

    ssSetNumOutputPorts(S, 2);
    ssSetOutputPortWidth(S, 0, 14);
    ssSetOutputPortWidth(S, 1, 2);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 26);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}

/*===================================*/
/* Function: mdlInitalizeSampleTimes */
/*===================================*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);

}

/*==================================*/
/* Function: mdlInitalizeConditions */
/*==================================*/

#define MDL_INITIALIZE_CONDITIONS  
#if defined(MDL_INITIALIZE_CONDITIONS)

static void mdlInitializeConditions(SimStruct *S)
{
}
#endif 

/*====================*/
/* Function: mdlStart */
/*====================*/
#define MDL_START 
#if defined(MDL_START) 

static void mdlStart(SimStruct *S)
{
    real_T *x = ssGetContStates(S);
	int i;
        
	/* Initialize the state vector */
	for (i = 0; i < ssGetNumContStates(S); i++)
	{
        x[i] = mxGetPr(ssGetSFcnParam(S, 0))[i];    
	}
}
#endif 



/*======================*/
/* Function: mdlOutputs */
/*======================*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
	real_T *x = ssGetContStates(S);
    real_T *y = ssGetOutputPortRealSignal(S, 0);
    real_T *y1 = ssGetOutputPortRealSignal(S, 1);
    
	int i;
	 
	for (i = 0; i < ssGetNumContStates(S); i++)
	{
        y[i] = x[i]; /* outputs are the states */
	}
    
    y1[0] = ny;
    y1[1] = nz;
}


/*=====================*/
/* Function: mdlUpdate */
/*=====================*/
#undef MDL_UPDATE  
#if defined(MDL_UPDATE)
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif


/*==========================*/
/* Function: mdlDerivatives */
/*==========================*/
#define MDL_DERIVATIVES  
#if defined(MDL_DERIVATIVES)
  static void mdlDerivatives(SimStruct *S)
  {
    real_T *x  = ssGetContStates(S);
    real_T *dx = ssGetdX(S);
    
    InputRealPtrsType u = ssGetInputPortRealSignalPtrs(S, 0);
    InputRealPtrsType u2 = ssGetInputPortRealSignalPtrs(S, 1);
    InputRealPtrsType u3 = ssGetInputPortRealSignalPtrs(S, 2);
    
    /* Declare a lot of variables */
    
    double *temp;
    double Gamma;
    double CX_tot, CY_tot, CZ_tot, Cl_tot, Cm_tot, Cn_tot;
	double Cx, Cz, Cm, Cy, Cn, Cl;
    double Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp;
    double delta_Cx_lef, delta_Cz_lef,delta_Cm_lef;
    double delta_Cy_lef,delta_Cn_lef,delta_Cl_lef;
    double delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef; 
    double delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef;
    double delta_Cnp_lef;
    double delta_Cy_r30, delta_Cn_r30, delta_Cl_r30;
    double delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef;
    double delta_Cl_a20, delta_Cl_a20_lef;
    double delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el, delta_Cm_ds;
	double da_norm, dr_norm, dlef_norm;
	double ALPHA, BETA, DE, DA, DR, DLEF, alpha_lef;
    double cpow, g, dq;
      
    temp = (double *)malloc(9*sizeof(double)); 
    
    /* Usefull notation of moments of inertia,
     see e.g. Lewis and Stevens "Aircraft control and simulation" */
    Gamma = Ixx * Izz - (Ixz  * Ixz);
    C1 = ((Iyy - Izz)  * Izz  - (Ixz * Ixz))/ Gamma;
    C2 = ((Ixx - Iyy + Izz ) * Ixz ) / Gamma;
    C3 = Izz / Gamma;
    C4 = Ixz / Gamma;
    C5 = (Izz - Ixx) / Iyy;
    C6 = Ixz / Iyy ;
    C7 = 1 / Iyy;
    C8 = (Ixx * (Ixx - Iyy ) + Ixz * Ixz) / Gamma;
    C9 = Ixx / Gamma;  
    
    /* ISA atmosphere */
    atmos(-z_earth, Vt, temp);
        Mach = temp[0];
        qbar = temp[1];  
        g    = temp[2];
        
    /* Engine model, based on Ying Huo's m-files */
    cpow = tgear (dth);
    pow_dot = pdot ( power, cpow );
    Thrust = thrst ( power, -z_earth, Mach ); 
    
    /* Body velocity components */
    u_body = Vt * cos(alpha) * cos(beta);
    v_body = Vt * sin(beta);
    w_body = Vt * sin(alpha) * cos(beta);
	 
    /* Transformation rad to deg for the lookup tables */
    ALPHA = alpha * rtd;
    BETA = beta * rtd;
	DE = de * rtd;
	DA = da * rtd;
	DR = dr * rtd;	 
	DLEF = dlef * rtd;	 
    
	/* LEF tables are only valid up until alpha = 45 degrees*/
	if (ALPHA > 45)
	{
	alpha_lef = 45;
	}
	else
	{
	alpha_lef = ALPHA;
	}
    
    /* Limit of elevator deflection is 25 degrees*/
	if (DE > 25)
	{
	DE = 25;
	}
	if (DE < -25)
	{
	DE = -25;
	}	 
    
    /* Normalizing the control deflections */
	da_norm = DA/21.5;
	dr_norm = DR/30.0;
      
    if (fi_flag == 1)          /* hifi lookup tables */
    {   
        
        dlef_norm = (1 - DLEF/25.0);
        
        Cx = _Cx(ALPHA, BETA, DE);
        Cy = _Cy(ALPHA, BETA);
        Cz = _Cz(ALPHA, BETA, DE);
        Cl = _Cl(ALPHA, BETA, DE);
        Cm = _Cm(ALPHA, BETA, DE);
        Cn = _Cn(ALPHA, BETA, DE);
	 
        Cxq = _CXq(ALPHA);
        Cyp = _CYp(ALPHA);
        Cyr = _CYr(ALPHA);
        Czq = _CZq(ALPHA);
        Clp = _CLp(ALPHA);
        Clr = _CLr(ALPHA);
        Cmq = _CMq(ALPHA);
        Cnp = _CNp(ALPHA);
        Cnr = _CNr(ALPHA);
	 
        delta_Cx_lef = _Cx_lef(alpha_lef, BETA) - _Cx(ALPHA, BETA, 0);
        delta_Cy_lef = _Cy_lef(alpha_lef, BETA) - _Cy(ALPHA, BETA);
        delta_Cz_lef = _Cz_lef(alpha_lef, BETA) - _Cz(ALPHA, BETA, 0);
        delta_Cl_lef = _Cl_lef(alpha_lef, BETA) - _Cl(ALPHA, BETA, 0);
        delta_Cm_lef = _Cm_lef(alpha_lef, BETA) - _Cm(ALPHA, BETA, 0);
        delta_Cn_lef = _Cn_lef(alpha_lef, BETA) - _Cn(ALPHA, BETA, 0);
	 
        delta_Cxq_lef = _delta_CXq_lef(alpha_lef);
        delta_Cyp_lef = _delta_CYp_lef(alpha_lef);
        delta_Cyr_lef = _delta_CYr_lef(alpha_lef);
        delta_Czq_lef = _delta_CZq_lef(alpha_lef);
        delta_Clp_lef = _delta_CLp_lef(alpha_lef);
        delta_Clr_lef = _delta_CLr_lef(alpha_lef);
        delta_Cmq_lef = _delta_CMq_lef(alpha_lef);
        delta_Cnp_lef = _delta_CNp_lef(alpha_lef);
        delta_Cnr_lef = _delta_CNr_lef(alpha_lef);
	 
        delta_Cy_r30 = _Cy_r30(ALPHA, BETA) - _Cy(ALPHA, BETA);
        delta_Cl_r30 = _Cl_r30(ALPHA, BETA) - _Cl(ALPHA, BETA, 0);
        delta_Cn_r30 = _Cn_r30(ALPHA, BETA) - _Cn(ALPHA, BETA, 0);
	 
        delta_Cy_a20     = _Cy_a20(ALPHA, BETA) - _Cy(ALPHA, BETA);
        delta_Cy_a20_lef = _Cy_a20_lef(alpha_lef, BETA) - 
            _Cy_lef(alpha_lef, BETA) - delta_Cy_a20;
        delta_Cn_a20     = _Cn_a20(ALPHA, BETA) - _Cn(ALPHA, BETA, 0);
        delta_Cn_a20_lef = _Cn_a20_lef(alpha_lef, BETA) - 
            _Cn_lef(alpha_lef, BETA) - delta_Cn_a20;
        delta_Cl_a20     = _Cl_a20(ALPHA, BETA) - _Cl(ALPHA, BETA, 0);
        delta_Cl_a20_lef = _Cl_a20_lef(alpha_lef, BETA) - 
            _Cl_lef(alpha_lef, BETA) - delta_Cl_a20;

        delta_Cnbeta = _delta_CNbeta(ALPHA);
        delta_Clbeta = _delta_CLbeta(ALPHA);
        delta_Cm     = _delta_Cm(ALPHA);
        eta_el       = _eta_el(DE);
        delta_Cm_ds  = _delta_Cm_ds(ALPHA,DE);
    }

else if (fi_flag == 0)  /* lofi lookup tables (do not include dlef) */
    {  
        dlef_norm = 0.0;

        damping(ALPHA, temp);
            Cxq = temp[0];
            Cyr = temp[1];
            Cyp = temp[2];
            Czq = temp[3];
            Clr = temp[4];
            Clp = temp[5];
            Cmq = temp[6];
            Cnr = temp[7];
            Cnp = temp[8];
        
        dmomdcon(ALPHA, BETA, temp);
            delta_Cl_a20 = temp[0];    
            delta_Cl_r30 = temp[1];    
            delta_Cn_a20 = temp[2];    
            delta_Cn_r30 = temp[3];    

        clcn(ALPHA, BETA, temp);
            Cl = temp[0];
            Cn = temp[1];

        cxcm(ALPHA, DE, temp);
            Cx = temp[0];
            Cm = temp[1];

            Cy = -.02*BETA + .021*da_norm + .086*dr_norm;

        cz(ALPHA, BETA, DE, temp);
            Cz = temp[0];  
        
        /* All other coeffcients are zero for the lofi model */
        delta_Cx_lef    = 0.0;
        delta_Cz_lef    = 0.0;
        delta_Cm_lef    = 0.0;
        delta_Cy_lef    = 0.0;
        delta_Cn_lef    = 0.0;
        delta_Cl_lef    = 0.0;
        delta_Cxq_lef   = 0.0;
        delta_Cyr_lef   = 0.0;
        delta_Cyp_lef   = 0.0;
        delta_Czq_lef   = 0.0;
        delta_Clr_lef   = 0.0;
        delta_Clp_lef   = 0.0;
        delta_Cmq_lef   = 0.0;
        delta_Cnr_lef   = 0.0;
        delta_Cnp_lef   = 0.0;
        delta_Cy_r30    = 0.0;
        delta_Cy_a20    = 0.0;
        delta_Cy_a20_lef= 0.0;
        delta_Cn_a20_lef= 0.0;
        delta_Cl_a20_lef= 0.0;
        delta_Cnbeta    = 0.0;
        delta_Clbeta    = 0.0;
        delta_Cm        = 0.0;
        eta_el          = 1.0; 
        delta_Cm_ds     = 0.0;
    }

    /* Total force coefficients */
    
    /* Cx_tot */
	CX_tot = Cx + delta_Cx_lef * dlef_norm 
     + (cref/(2*Vt))*(Cxq + delta_Cxq_lef * dlef_norm) * q_body;
    /* Cy_tot */
    CY_tot = Cy + delta_Cy_lef * dlef_norm 
     + (delta_Cy_a20 + delta_Cy_a20_lef * dlef_norm) * da_norm
     + delta_Cy_r30 * dr_norm 
     + (bref / (2*Vt))*(Cyr + delta_Cyr_lef * dlef_norm) * r_body 
     + (bref/(2*Vt))*(Cyp + delta_Cyp_lef * dlef_norm) * p_body;
    /* Cz_tot */
    CZ_tot = Cz + delta_Cz_lef * dlef_norm 
     + (cref/(2*Vt))*(Czq + delta_Czq_lef * dlef_norm) * q_body;
    
    /* Total moment coefficients */
    
    /* Cl_tot */
    Cl_tot = Cl + delta_Cl_lef * dlef_norm 
     + (delta_Cl_a20 + delta_Cl_a20_lef * dlef_norm) * da_norm
     + delta_Cl_r30 * dr_norm 
     + (bref / ((2*Vt))*(Clr + delta_Clr_lef * dlef_norm)) * r_body 
     + ((bref / (2*Vt)) * (Clp + delta_Clp_lef * dlef_norm)) * p_body 
     + delta_Clbeta * beta * rtd;
    /* Cm_tot */
    Cm_tot = Cm * eta_el + CZ_tot * (xcgr - xcg) + delta_Cm_lef * dlef_norm 
     + (cref / (2*Vt))*(Cmq + delta_Cmq_lef * dlef_norm) * q_body 
     + delta_Cm + delta_Cm_ds;
    /* Cn_tot */
    Cn_tot = Cn + delta_Cn_lef * dlef_norm 
     - CY_tot * (xcgr - xcg)*(cref/bref) 
     + (delta_Cn_a20 + delta_Cn_a20_lef * dlef_norm) * da_norm 
     + ((bref / (2*Vt)) * (Cnr + delta_Cnr_lef * dlef_norm))* r_body
     + ((bref / (2*Vt)) * (Cnp + delta_Cnp_lef * dlef_norm)) * p_body 
     + delta_Cn_r30 * dr_norm + delta_Cnbeta * beta * rtd;
	
	/* Total forces */
    Xbar = qbar * Sref * CX_tot;
    Ybar = qbar * Sref * CY_tot;
    Zbar = qbar * Sref * CZ_tot;
    
    /* Total moments */
    Lbar = Cl_tot * qbar * Sref * bref;
    Mbar = Cm_tot * qbar * Sref * cref;
    Nbar = Cn_tot * qbar * Sref * bref;
    
    /* Derivatives */
    u_body_dot = r_body * v_body - q_body * w_body 
     + (Xbar + Thrust) / mass + 2*(q1*q3 - q0*q2)*g;
    v_body_dot = p_body * w_body - r_body * u_body 
     + Ybar / mass + 2*(q2*q3 + q0*q1)*g;
    w_body_dot = q_body * u_body - p_body * v_body 
     + Zbar / mass + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*g;
    
    Vt_dot      = (u_body * u_body_dot + v_body * v_body_dot 
     + w_body * w_body_dot) / Vt;
    beta_dot    = (v_body_dot * Vt - v_body * Vt_dot) 
     / (Vt * Vt * cos(beta));
    alpha_dot   = (u_body * w_body_dot - w_body * u_body_dot) 
     / (u_body * u_body + w_body * w_body);
     
    q0_dot      = 0.5 * (-p_body * q1 - q_body * q2 - r_body * q3);
    q1_dot      = 0.5 * ( p_body * q0 + r_body * q2 - q_body * q3);
    q2_dot      = 0.5 * ( q_body * q0 - r_body * q1 + p_body * q3);
    q3_dot      = 0.5 * ( r_body * q0 + q_body * q1 - p_body * q2);
     
    /* correction term from Moldy User’s Manual by K. Refson */
    dq = q0 * q0_dot + q1 * q1_dot + q2 * q2_dot + q3 * q3_dot;
    
    q0_dot -= dq * q0;
    q1_dot -= dq * q1;
    q2_dot -= dq * q2;
    q3_dot -= dq * q3;
      
    p_body_dot  = (C1 * r_body + C2 * p_body) * q_body 
     + C3 * Lbar + C4 * (Nbar + q_body * heng);
    q_body_dot  = C5 * p_body * r_body 
     - C6 * (p_body * p_body - r_body * r_body) 
     + C7 * (Mbar - heng * r_body);
    r_body_dot  = (C8 * p_body - C2 * r_body) * q_body 
     + C4 * Lbar + C9 * (Nbar + q_body * heng);
  
    x_earth_dot = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * u_body 
     + 2*(q1*q2 - q0*q3) * v_body 
     + 2*(q1*q3 + q0*q2) * w_body;
    y_earth_dot = 2*(q1*q2 + q0*q3) * u_body 
     + (q0*q0 - q1*q1 + q2*q2 - q3*q3) * v_body 
     + 2*(q2*q3 - q0*q1) * w_body;
    z_earth_dot = 2*(q1*q3 - q0*q2) * u_body 
     + 2*(q2*q3 + q0*q1) * v_body 
     + (q0*q0 - q1*q1 - q2*q2 + q3*q3) * w_body;
    
    /* normal accelerations */ 
    ny = Ybar/mass/g;
    nz = -Zbar/mass/g;    
    
    free(temp);
    
}   /* end mdlDerivatives */

#endif 


/*===================================*/
/* Function: mdlTerminate */
/*===================================*/
static void mdlTerminate(SimStruct *S)
{
}

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
