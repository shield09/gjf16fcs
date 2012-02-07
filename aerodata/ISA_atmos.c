/*---------------------------------------------------------------------- */
/*                                                                       */
/* ISA atmosphere model taken from "Elements of airplane performance"    */
/* by G .J. J. Ruijgrok, Delft university press/ 1996.                   */
/*                                                                       */
/* File "ISA_atmos.c"                                                    */
/* by L. Sonneveldt                                                      */
/* May, 2006                                                             */
/*                                                                       */
/*---------------------------------------------------------------------- */

void atmos(double alt, double vt, double *coeff )
{
    double rho0 = 1.22505;
    double Re = 6371000.0;
    double R = 287.04;
    double g0 = 9.80665;
    double pres0 = 101325.0;
    double temp0 = 288.16;
    double gamma = 1.403;
    double grad = -0.0065;
    double temp, rho, mach, qbar, grav, pres;

    temp = temp0 + grad * alt;   
    pres = pres0*(pow((temp/temp0),(-g0/(R*grad))));
    rho  = pres/(R*temp);
    
    if (alt >= 11000.0) {
       temp = temp0 + grad * 11000;
       pres  = pres0*(pow((temp/temp0),(-g0/(R*grad))))*pow(10.0,(-g0*(alt-11000)/(2.3026*R*temp)));
       rho  = pres0*(pow((temp/temp0),(-g0/(R*grad))))/(R*temp)*pow(10.0,(-g0*(alt-11000)/(2.3026*R*temp)));
    }

    mach = vt / sqrt(gamma * R * temp);
    qbar = .5 * rho * vt * vt;
    grav = g0*(Re*Re/((Re+alt)*(Re+alt)));
 
    coeff[0] = mach;
    coeff[1] = qbar;
    coeff[2] = grav;
}
