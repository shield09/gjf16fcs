/*---------------------------------------------------------------------- */
/*                                                                       */
/* tail and wing failure aerodata                                        */
/* based on Henri's AVL results                                          */
/* C_dmg = C_nom*C_fail                                                  */
/*                                                                       */
/* Note that for the wing dmg, some coefficients do not exist 
   for the nom. model, hence there C_dmg = C_fail;                       */
/*---------------------------------------------------------------------- */
void CYtail(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cyb00, cyb08, cybtail; 
    double cyp00, cyp08, cyptail; 
    double cyr00, cyr08, cyrtail; 
    double CYb_tail00[] = {-0.369080, -0.462490, -0.556678, -0.650650,
    -0.743378, -0.833921, -0.921499, -1.005694, -1.085914, -1.161698, -1.232879};    
    double CYb_tail08[] = {-0.386424, -0.488190, -0.592532, -0.698316,
    -0.804302, -0.909242, -1.012020, -1.111710, -1.207410, -1.298389, -1.384139};
    double CYp_tail00[] = {0.025770,  0.025561,  0.022990,  0.017692,
    0.009446, -0.001823, -0.016065, -0.033124, -0.052757, -0.074659, -0.098460}; 
    double CYp_tail08[] = {0.029440,  0.030270,  0.028635,  0.024045,
    0.016158,  0.004795, -0.010059, -0.028258, -0.049556, -0.073594, -0.099945};
    double CYr_tail00[] = {0.059898,  0.167651,  0.276877,  0.386836,
    0.496608,  0.605263,  0.711910,  0.816162,  0.917118,  1.014156,  1.106902};
    double CYr_tail08[] = {-0.043880,  0.076206,  0.199662,  0.325616,
    0.452951,  0.580444,  0.706903,  0.831300,  0.952499,  1.069533,  1.181622};    
    int k, L;
 
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from closest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cyb00 = CYb_tail00[k-1] + fabs(da)*(CYb_tail00[L-1]-CYb_tail00[k-1]);
    cyb08 = CYb_tail08[k-1] + fabs(da)*(CYb_tail08[L-1]-CYb_tail08[k-1]);
    cyp00 = CYp_tail00[k-1] + fabs(da)*(CYp_tail00[L-1]-CYp_tail00[k-1]);
    cyp08 = CYp_tail08[k-1] + fabs(da)*(CYp_tail08[L-1]-CYp_tail08[k-1]);    
    cyr00 = CYr_tail00[k-1] + fabs(da)*(CYr_tail00[L-1]-CYr_tail00[k-1]);
    cyr08 = CYr_tail08[k-1] + fabs(da)*(CYr_tail08[L-1]-CYr_tail08[k-1]);
    
    cybtail = ((mach/0.8)*cyb08+(0.8-mach)/0.8*cyb00)-((mach/0.8)*CYb_tail08[10]+(0.8-mach)/0.8*CYb_tail00[10]);
    cyptail = ((mach/0.8)*cyp08+(0.8-mach)/0.8*cyp00)/((mach/0.8)*CYp_tail08[10]+(0.8-mach)/0.8*CYp_tail00[10]);
    cyrtail = ((mach/0.8)*cyr08+(0.8-mach)/0.8*cyr00)/((mach/0.8)*CYr_tail08[10]+(0.8-mach)/0.8*CYr_tail00[10]);
    
    coeff[0] = cybtail;
    coeff[1] = cyptail;
    coeff[2] = cyrtail;
}

void Cltail(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double clb00, clb08, clbtail; 
    double clp00, clp08, clptail; 
    double clr00, clr08, clrtail; 
    double Clb_tail00[] = {0.033162,  0.031233,  0.026970,  0.020350,
    0.011435,  0.000352, -0.012716, -0.027579, -0.043969, -0.061610, -0.080208};
    double Clb_tail08[] = {0.036422,  0.034754,  0.030557,  0.023740,
    0.014294,  0.002281, -0.012167, -0.028868, -0.047565, -0.067948, -0.089677};
    double Clp_tail00[] = {-0.275899, -0.275873, -0.275965, -0.276315,
    -0.277084, -0.278436, -0.280528, -0.283501, -0.287465, -0.292501, -0.298645}; 
    double Clp_tail08[] = {-0.293351, -0.293301, -0.293335, -0.293602,
    -0.294276, -0.295544, -0.297595, -0.300600, -0.304708, -0.310027, -0.316623};
    double Clr_tail00[] = {-0.037148, -0.034531, -0.029254, -0.021220,
    -0.010418,  0.003093,  0.019178,  0.037715,  0.058440,  0.081090,  0.105354}; 
    double Clr_tail08[] = {-0.041606, -0.039197, -0.033840, -0.025365,
    -0.013690,  0.001197,  0.019227,  0.040276,  0.064118,  0.090455,  0.118931};
    int k, L;    
 
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from closest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    clb00 = Clb_tail00[k-1] + fabs(da)*(Clb_tail00[L-1]-Clb_tail00[k-1]);
    clb08 = Clb_tail08[k-1] + fabs(da)*(Clb_tail08[L-1]-Clb_tail08[k-1]);
    clp00 = Clp_tail00[k-1] + fabs(da)*(Clp_tail00[L-1]-Clp_tail00[k-1]);
    clp08 = Clp_tail08[k-1] + fabs(da)*(Clp_tail08[L-1]-Clp_tail08[k-1]);    
    clr00 = Clr_tail00[k-1] + fabs(da)*(Clr_tail00[L-1]-Clr_tail00[k-1]);
    clr08 = Clr_tail08[k-1] + fabs(da)*(Clr_tail08[L-1]-Clr_tail08[k-1]);
    
    clbtail = ((mach/0.8)*clb08+(0.8-mach)/0.8*clb00)-((mach/0.8)*Clb_tail08[10]+(0.8-mach)/0.8*Clb_tail00[10]);
    clptail = ((mach/0.8)*clp08+(0.8-mach)/0.8*clp00)/((mach/0.8)*Clp_tail08[10]+(0.8-mach)/0.8*Clp_tail00[10]);
    clrtail = ((mach/0.8)*clr08+(0.8-mach)/0.8*clr00)/((mach/0.8)*Clr_tail08[10]+(0.8-mach)/0.8*Clr_tail00[10]);
    
    coeff[0] = clbtail;
    coeff[1] = clptail;
    coeff[2] = clrtail;
}

void Cntail(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cnb00, cnb08, cnbtail; 
    double cnp00, cnp08, cnptail; 
    double cnr00, cnr08, cnrtail; 
    double Cnb_tail00[] = {-0.067245, -0.026510,  0.016379,  0.060873,
    0.106444,  0.152590,  0.198844,  0.244921,  0.290392,  0.334900,  0.378187};
    double Cnb_tail08[] = {-0.127900, -0.083384, -0.035790,  0.014357,
    0.066494,  0.120012,  0.174309,  0.228851,  0.283064,  0.336430,  0.388515};
    double Cnp_tail00[] = {-0.013665, -0.013744, -0.012700, -0.010292,
    -0.006322, -0.000650,  0.006807,  0.016081,  0.027149,  0.039945,  0.054348}; 
    double Cnp_tail08[] = {-0.015444, -0.016039, -0.015470, -0.013431,
    -0.009658, -0.003951,  0.003822,  0.013713,  0.025714,  0.039748,  0.055682};
    double Cnr_tail00[] = {-0.170364, -0.219297, -0.270493, -0.323679,
    -0.378499, -0.434554, -0.491403, -0.548861, -0.606387, -0.663591, -0.720136};
    double Cnr_tail08[] = {-0.182085, -0.237099, -0.295336, -0.356547,
    -0.420360, -0.486294, -0.553811, -0.622412, -0.691478, -0.760423, -0.828709};  
    int k, L;    
 
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from closest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cnb00 = Cnb_tail00[k-1] + fabs(da)*(Cnb_tail00[L-1]-Cnb_tail00[k-1]);
    cnb08 = Cnb_tail08[k-1] + fabs(da)*(Cnb_tail08[L-1]-Cnb_tail08[k-1]);
    cnp00 = Cnp_tail00[k-1] + fabs(da)*(Cnp_tail00[L-1]-Cnp_tail00[k-1]);
    cnp08 = Cnp_tail08[k-1] + fabs(da)*(Cnp_tail08[L-1]-Cnp_tail08[k-1]);    
    cnr00 = Cnr_tail00[k-1] + fabs(da)*(Cnr_tail00[L-1]-Cnr_tail00[k-1]);
    cnr08 = Cnr_tail08[k-1] + fabs(da)*(Cnr_tail08[L-1]-Cnr_tail08[k-1]);
    
    cnbtail = ((mach/0.8)*cnb08+(0.8-mach)/0.8*cnb00)-((mach/0.8)*Cnb_tail08[10]+(0.8-mach)/0.8*Cnb_tail00[10]);
    cnptail = ((mach/0.8)*cnp08+(0.8-mach)/0.8*cnp00)/((mach/0.8)*Cnp_tail08[10]+(0.8-mach)/0.8*Cnp_tail00[10]);
    cnrtail = ((mach/0.8)*cnr08+(0.8-mach)/0.8*cnr00)/((mach/0.8)*Cnr_tail08[10]+(0.8-mach)/0.8*Cnr_tail00[10]);
    
    coeff[0] = cnbtail;
    coeff[1] = cnptail;
    coeff[2] = cnrtail;
}

void CLwing(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cla00, cla08, clawing; 
    double clb00, clb08, clbwing; 
    double clp00, clp08, clpwing; 
    double clq00, clq08, clqwing; 
    double clr00, clr08, clrwing; 
    int k, L;     
    double Cla_wing00[] = {2.947413,  3.210292,  3.469855,  3.713924,  3.951623,
    4.194496,  4.455260,  4.743130,  5.046612,  5.342305,  5.612236};
    double Cla_wing08[] = {3.249808,  3.522976,  3.792655,  4.047637,  4.299032,
    4.560887,  4.849544,  5.177845,  5.533919,  5.890846,  6.227154};
    double Clb_wing00[] = {0.038887,  0.041362,  0.040033,  0.036817,  0.032629,
    0.027821,  0.022455,  0.016558,  0.010526,  0.004897,  0.000000};
    double Clb_wing08[] = {0.033207,  0.035356,  0.034414,  0.031930,  0.028646,
    0.024796,  0.020361,  0.015302,  0.009936,  0.004734,  0.000000};
    double Clp_wing00[] = {-0.471152, -0.490875, -0.500697, -0.499948, -0.487386,
    -0.460553, -0.415024, -0.344747, -0.248063, -0.130327,  0.000000};
    double Clp_wing08[] = {-0.510608, -0.532997, -0.545157, -0.546503, -0.535710,
    -0.509923, -0.463770, -0.389557, -0.283993, -0.151458,  0.000000};
    double Clq_wing00[] = {6.469810,  6.792623,  7.110568,  7.413656,  7.713890,
    8.026368,  8.368623,  8.753923,  9.167457,  9.578146,  9.961749};
    double Clq_wing08[] = {7.096896,  7.435208,  7.769571,  8.090858,  8.412883,
    8.753890, 9.136510,  9.579045, 10.065853, 10.561135, 11.036746};  
    double Clr_wing00[] = {-0.039063, -0.041330, -0.040132, -0.037184, -0.033238,
    -0.028579, -0.023250, -0.017271, -0.011055, -0.005174,  0.000000};
    double Clr_wing08[] = {-0.032099, -0.033947, -0.033137, -0.030992, -0.028074,
    -0.024537, -0.020336, -0.015420, -0.010097, -0.004847,  0.000000};      
    
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from closest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cla00 = Cla_wing00[k-1] + fabs(da)*(Cla_wing00[L-1]-Cla_wing00[k-1]);
    cla08 = Cla_wing08[k-1] + fabs(da)*(Cla_wing08[L-1]-Cla_wing08[k-1]);    
    clb00 = Clb_wing00[k-1] + fabs(da)*(Clb_wing00[L-1]-Clb_wing00[k-1]);
    clb08 = Clb_wing08[k-1] + fabs(da)*(Clb_wing08[L-1]-Clb_wing08[k-1]);
    clp00 = Clp_wing00[k-1] + fabs(da)*(Clp_wing00[L-1]-Clp_wing00[k-1]);
    clp08 = Clp_wing08[k-1] + fabs(da)*(Clp_wing08[L-1]-Clp_wing08[k-1]);
    clq00 = Clq_wing00[k-1] + fabs(da)*(Clq_wing00[L-1]-Clq_wing00[k-1]);
    clq08 = Clq_wing08[k-1] + fabs(da)*(Clq_wing08[L-1]-Clq_wing08[k-1]);    
    clr00 = Clr_wing00[k-1] + fabs(da)*(Clr_wing00[L-1]-Clr_wing00[k-1]);
    clr08 = Clr_wing08[k-1] + fabs(da)*(Clr_wing08[L-1]-Clr_wing08[k-1]);
    
    clawing = ((mach/0.8)*cla08+(0.8-mach)/0.8*cla00)-((mach/0.8)*Cla_wing08[10]+(0.8-mach)/0.8*Cla_wing00[10]);
    clbwing = ((mach/0.8)*clb08+(0.8-mach)/0.8*clb00);
    clpwing = ((mach/0.8)*clp08+(0.8-mach)/0.8*clp00);
    clqwing = ((mach/0.8)*clq08+(0.8-mach)/0.8*clq00)/((mach/0.8)*Clq_wing08[10]+(0.8-mach)/0.8*Clq_wing00[10]);
    clrwing = ((mach/0.8)*clr08+(0.8-mach)/0.8*clr00);
    
    coeff[0] = clawing;
    coeff[1] = clbwing;
    coeff[2] = clpwing;
    coeff[3] = clqwing;
    coeff[4] = clrwing;
}

void Cmwing(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cma00, cma08, cmawing; 
    double cmb00, cmb08, cmbwing; 
    double cmp00, cmp08, cmpwing; 
    double cmq00, cmq08, cmqwing; 
    double cmr00, cmr08, cmrwing; 
    int k, L;     
    double Cma_wing00[] = {-1.758322, -1.521810, -1.273889, -1.050288, -0.863123,
    -0.723782, -0.649314, -0.655475, -0.729471, -0.841405, -0.969287};
    double Cma_wing08[] = {-1.888201, -1.619678, -1.333673, -1.069286, -0.841054,
    -0.664221, -0.562210, -0.557914, -0.639430, -0.772952, -0.933073};
    double Cmb_wing00[] = {-0.007075, -0.005797, -0.005238, -0.005564, -0.006410,
    -0.007220, -0.007371, -0.006393, -0.004509, -0.002265,  0.000000};
    double Cmb_wing08[] = {-0.007805, -0.006812, -0.006283, -0.006413, -0.006953,
    -0.007469, -0.007421, -0.006372, -0.004506, -0.002291,  0.000000};
    double Cmp_wing00[] = {0.094016, 0.080642,  0.075086,  0.078988,  0.091651,
    0.109481,  0.124643,  0.125127,  0.103403,  0.060424,  0.000000};
    double Cmp_wing08[] = {0.103199,  0.087065,  0.079462,  0.082584,  0.096251,
    0.117000,  0.136212,  0.139807, 0.118029,  0.070470,  0.000000};
    double Cmq_wing00[] = {-8.279084, -8.187114, -8.063527, -7.924902, -7.792587,
    -7.690441, -7.649189, -7.697075, -7.822787, -7.991774, -8.179813};
    double Cmq_wing08[] = {-9.560942, -9.463394, -9.328935, -9.173249, -9.019029,
    -8.894692, -8.839551, -8.891730, -9.040734, -9.247063, -9.483506};  
    double Cmr_wing00[] = {0.009364,  0.008284,  0.007500,  0.007342,  0.007695,
    0.008134,  0.008040,  0.006888,  0.004848,  0.002440,  0.000000};
    double Cmr_wing08[] = {0.010041,  0.009322,  0.008644,  0.008308,  0.008301,
    0.008353,  0.007982,  0.006732,  0.004732,  0.002406,  0.000000};  
    
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from cmosest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cma00 = Cma_wing00[k-1] + fabs(da)*(Cma_wing00[L-1]-Cma_wing00[k-1]);
    cma08 = Cma_wing08[k-1] + fabs(da)*(Cma_wing08[L-1]-Cma_wing08[k-1]);    
    cmb00 = Cmb_wing00[k-1] + fabs(da)*(Cmb_wing00[L-1]-Cmb_wing00[k-1]);
    cmb08 = Cmb_wing08[k-1] + fabs(da)*(Cmb_wing08[L-1]-Cmb_wing08[k-1]);
    cmp00 = Cmp_wing00[k-1] + fabs(da)*(Cmp_wing00[L-1]-Cmp_wing00[k-1]);
    cmp08 = Cmp_wing08[k-1] + fabs(da)*(Cmp_wing08[L-1]-Cmp_wing08[k-1]);
    cmq00 = Cmq_wing00[k-1] + fabs(da)*(Cmq_wing00[L-1]-Cmq_wing00[k-1]);
    cmq08 = Cmq_wing08[k-1] + fabs(da)*(Cmq_wing08[L-1]-Cmq_wing08[k-1]);    
    cmr00 = Cmr_wing00[k-1] + fabs(da)*(Cmr_wing00[L-1]-Cmr_wing00[k-1]);
    cmr08 = Cmr_wing08[k-1] + fabs(da)*(Cmr_wing08[L-1]-Cmr_wing08[k-1]);
    
    cmawing = ((mach/0.8)*cma08+(0.8-mach)/0.8*cma00)-((mach/0.8)*Cma_wing08[10]+(0.8-mach)/0.8*Cma_wing00[10]);
    cmbwing = ((mach/0.8)*cmb08+(0.8-mach)/0.8*cmb00);
    cmpwing = ((mach/0.8)*cmp08+(0.8-mach)/0.8*cmp00);
    cmqwing = ((mach/0.8)*cmq08+(0.8-mach)/0.8*cmq00)/((mach/0.8)*Cmq_wing08[10]+(0.8-mach)/0.8*Cmq_wing00[10]);
    cmrwing = ((mach/0.8)*cmr08+(0.8-mach)/0.8*cmr00);
    
    coeff[0] = cmawing;
    coeff[1] = cmbwing;
    coeff[2] = cmpwing;
    coeff[3] = cmqwing;
    coeff[4] = cmrwing;
}

void CYwing(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cya00, cya08, cyawing; 
    double cyb00, cyb08, cybwing; 
    double cyp00, cyp08, cypwing; 
    double cyq00, cyq08, cyqwing; 
    double cyr00, cyr08, cyrwing; 
    int k, L;     
    double Cya_wing00[] = {-0.076776, -0.084241, -0.090536, -0.093901,
    -0.092935, -0.086748, -0.074611, -0.056493, -0.035568, -0.016007,  0.000000};
    double Cya_wing08[] = {-0.084542, -0.092588, -0.101216, -0.108137,
    -0.110842, -0.107200, -0.095291, -0.074381, -0.048276, -0.022429,  0.000000};
    double Cyb_wing00[] = {-1.231026, -1.231093, -1.231069, -1.231038,
    -1.231071, -1.231210, -1.231468, -1.231839, -1.232245, -1.232603, -1.232879};
    double Cyb_wing08[] = {-1.382503, -1.382560, -1.382534, -1.382475,
    -1.382458, -1.382535, -1.382740, -1.383076, -1.383469, -1.383837, -1.384139};
    double Cyp_wing00[] = {-0.128366, -0.127813, -0.127580, -0.127595,
    -0.127557, -0.126899, -0.124820, -0.120454, -0.113856, -0.106133, -0.098460};
    double Cyp_wing08[] = {-0.134896, -0.134242, -0.133852, -0.133814,
    -0.133940, -0.133610, -0.131760, -0.127114, -0.119473, -0.109969, -0.099945};
    double Cyq_wing00[] = {-0.111842, -0.120191, -0.127365, -0.130946,
    -0.128948, -0.120142, -0.103439, -0.078623, -0.049857, -0.022676,  0.000000};
    double Cyq_wing08[] = {-0.126008, -0.135334, -0.145845, -0.154191,
    -0.156788, -0.150773, -0.133579, -0.104189, -0.067796, -0.031701,  0.000000};  
    double Cyr_wing00[] = {1.104979,  1.105040,  1.105021,  1.104994,
    1.105030,  1.105169,  1.105430,  1.105810,  1.106231,  1.106608,  1.106902};
    double Cyr_wing08[] = {1.179962,  1.180010,  1.179989,  1.179941,
    1.179929,  1.180007,  1.180209,  1.180542,  1.180936,  1.181310,  1.181622};      
 
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from closest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cya00 = Cya_wing00[k-1] + fabs(da)*(Cya_wing00[L-1]-Cya_wing00[k-1]);
    cya08 = Cya_wing08[k-1] + fabs(da)*(Cya_wing08[L-1]-Cya_wing08[k-1]);    
    cyb00 = Cyb_wing00[k-1] + fabs(da)*(Cyb_wing00[L-1]-Cyb_wing00[k-1]);
    cyb08 = Cyb_wing08[k-1] + fabs(da)*(Cyb_wing08[L-1]-Cyb_wing08[k-1]);
    cyp00 = Cyp_wing00[k-1] + fabs(da)*(Cyp_wing00[L-1]-Cyp_wing00[k-1]);
    cyp08 = Cyp_wing08[k-1] + fabs(da)*(Cyp_wing08[L-1]-Cyp_wing08[k-1]);
    cyq00 = Cyq_wing00[k-1] + fabs(da)*(Cyq_wing00[L-1]-Cyq_wing00[k-1]);
    cyq08 = Cyq_wing08[k-1] + fabs(da)*(Cyq_wing08[L-1]-Cyq_wing08[k-1]);    
    cyr00 = Cyr_wing00[k-1] + fabs(da)*(Cyr_wing00[L-1]-Cyr_wing00[k-1]);
    cyr08 = Cyr_wing08[k-1] + fabs(da)*(Cyr_wing08[L-1]-Cyr_wing08[k-1]);
    
    cyawing = ((mach/0.8)*cya08+(0.8-mach)/0.8*cya00);
    cybwing = ((mach/0.8)*cyb08+(0.8-mach)/0.8*cyb00)-((mach/0.8)*Cyb_wing08[10]+(0.8-mach)/0.8*Cyb_wing00[10]);
    cypwing = ((mach/0.8)*cyp08+(0.8-mach)/0.8*cyp00)/((mach/0.8)*Cyp_wing08[10]+(0.8-mach)/0.8*Cyp_wing00[10]);
    cyqwing = ((mach/0.8)*cyq08+(0.8-mach)/0.8*cyq00);
    cyrwing = ((mach/0.8)*cyr08+(0.8-mach)/0.8*cyr00)/((mach/0.8)*Cyr_wing08[10]+(0.8-mach)/0.8*Cyr_wing00[10]);
    
    coeff[0] = cyawing;
    coeff[1] = cybwing;
    coeff[2] = cypwing;
    coeff[3] = cyqwing;
    coeff[4] = cyrwing;
}

void Clwing(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cla00, cla08, clawing; 
    double clb00, clb08, clbwing; 
    double clp00, clp08, clpwing; 
    double clq00, clq08, clqwing; 
    double clr00, clr08, clrwing; 
    int k, L;     
    double Cla_wing00[] = {0.265610,  0.290678,  0.306349,  0.311709,  0.306856,
    0.290862,  0.261301,  0.215190,  0.153274,  0.079863,  0.000000};
    double Cla_wing08[] = {0.281961,  0.310648,  0.329925,  0.338659,  0.336748,
    0.322855,  0.293755,  0.245311,  0.177428,  0.094039,  0.000000};
    double Clb_wing00[] = {-0.086593, -0.086350, -0.086368, -0.086352, -0.086180,
    -0.085794, -0.085142, -0.084182, -0.082964, -0.081602, -0.080208};
    double Clb_wing08[] = {-0.095109, -0.094880, -0.094891, -0.094893, -0.094778,
    -0.094493, -0.093984, -0.093200, -0.092167, -0.090966, -0.089677};
    double Clp_wing00[] = {-0.191398, -0.193281, -0.193838, -0.193756, -0.193919,
    -0.195569, -0.200591, -0.211690, -0.231244, -0.260290, -0.298645};
    double Clp_wing08[] = {-0.204242, -0.206629, -0.207473, -0.207442, -0.207407,
    -0.208624, -0.213100, -0.223853, -0.243767, -0.274508, -0.316623};
    double Clq_wing00[] = {0.395557,  0.423860,  0.439660,  0.442524,  0.432702,
    0.408820,  0.367235,  0.303307,  0.217316,  0.114227,  0.000000};
    double Clq_wing08[] = {0.425665,  0.459543,  0.480313,  0.487126,  0.480157,
    0.457720,  0.415283,  0.346794,  0.251569,  0.134130,  0.000000};  
    double Clr_wing00[] = {0.112016,  0.111795,  0.111799,  0.111765,  0.111582,
    0.111190,  0.110527,  0.109541,  0.108274,  0.106839,  0.105354};
    double Clr_wing08[] = {0.124432,  0.124237,  0.124236,  0.124218,  0.124093,
    0.123808,  0.123307,  0.122528,  0.121489,  0.120263,  0.118931};      
    
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from closest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cla00 = Cla_wing00[k-1] + fabs(da)*(Cla_wing00[L-1]-Cla_wing00[k-1]);
    cla08 = Cla_wing08[k-1] + fabs(da)*(Cla_wing08[L-1]-Cla_wing08[k-1]);    
    clb00 = Clb_wing00[k-1] + fabs(da)*(Clb_wing00[L-1]-Clb_wing00[k-1]);
    clb08 = Clb_wing08[k-1] + fabs(da)*(Clb_wing08[L-1]-Clb_wing08[k-1]);
    clp00 = Clp_wing00[k-1] + fabs(da)*(Clp_wing00[L-1]-Clp_wing00[k-1]);
    clp08 = Clp_wing08[k-1] + fabs(da)*(Clp_wing08[L-1]-Clp_wing08[k-1]);
    clq00 = Clq_wing00[k-1] + fabs(da)*(Clq_wing00[L-1]-Clq_wing00[k-1]);
    clq08 = Clq_wing08[k-1] + fabs(da)*(Clq_wing08[L-1]-Clq_wing08[k-1]);    
    clr00 = Clr_wing00[k-1] + fabs(da)*(Clr_wing00[L-1]-Clr_wing00[k-1]);
    clr08 = Clr_wing08[k-1] + fabs(da)*(Clr_wing08[L-1]-Clr_wing08[k-1]);
    
    clawing = ((mach/0.8)*cla08+(0.8-mach)/0.8*cla00);
    clbwing = ((mach/0.8)*clb08+(0.8-mach)/0.8*clb00)-((mach/0.8)*Clb_wing08[10]+(0.8-mach)/0.8*Clb_wing00[10]);
    clpwing = ((mach/0.8)*clp08+(0.8-mach)/0.8*clp00)/((mach/0.8)*Clp_wing08[10]+(0.8-mach)/0.8*Clp_wing00[10]);
    clqwing = ((mach/0.8)*clq08+(0.8-mach)/0.8*clq00);
    clrwing = ((mach/0.8)*clr08+(0.8-mach)/0.8*clr00)/((mach/0.8)*Clr_wing08[10]+(0.8-mach)/0.8*Clr_wing00[10]);
    
    coeff[0] = clawing;
    coeff[1] = clbwing;
    coeff[2] = clpwing;
    coeff[3] = clqwing;
    coeff[4] = clrwing;
}

void Cnwing(double mach, double dmg, double *coeff){

    int fix(double);
    double s, da;
    double cna00, cna08, cnawing; 
    double cnb00, cnb08, cnbwing; 
    double cnp00, cnp08, cnpwing; 
    double cnq00, cnq08, cnqwing; 
    double cnr00, cnr08, cnrwing; 
    int k, L;     
    double Cna_wing00[] = {0.033124,  0.035956,  0.039215,  0.041807,  0.042634,
    0.040944,  0.036090,  0.027812,  0.017700,  0.008021,  0.000000};
    double Cna_wing08[] = {0.036154,  0.039054,  0.043389,  0.047840,  0.050744,
    0.050660,  0.046239,  0.036760,  0.024123,  0.011287,  0.000000};
    double Cnb_wing00[] = {0.377385,  0.377409,  0.377393,  0.377362,  0.377355,
    0.377396,  0.377502,  0.377673,  0.377871,  0.378048,  0.378187};
    double Cnb_wing08[] = {0.387812,  0.387831,  0.387815,  0.387772,  0.387740,
    0.387752,  0.387833,  0.387987,  0.388179,  0.388362,  0.388515};
    double Cnp_wing00[] = {0.068487,  0.068281,  0.068157,  0.068163,  0.068207,
    0.068028,  0.067195,  0.065200,  0.062013,  0.058192,  0.054348};
    double Cnp_wing08[] = {0.072277,  0.072047,  0.071849,  0.071817,  0.071938,
    0.071935,  0.071248,  0.069141,  0.065444,  0.060726,  0.055682};
    double Cnq_wing00[] = {0.048559,  0.051713,  0.055659,  0.058782,  0.059558,
    0.057000,  0.050221,  0.038803,  0.024850,  0.011374,  0.000000};
    double Cnq_wing08[] = {0.054254,  0.057566,  0.063090,  0.068790,  0.072287,
    0.071651,  0.065092,  0.051648,  0.033946,  0.015975,  0.000000};  
    double Cnr_wing00[] = {-0.719297, -0.719318, -0.719304, -0.719277, -0.719271,
    -0.719314, -0.719421, -0.719597, -0.719802, -0.719989, -0.720136};
    double Cnr_wing08[] = {-0.827986, -0.828001, -0.827988, -0.827951, -0.827924,
    -0.827938, -0.828018, -0.828172, -0.828364, -0.828551, -0.828709};      
   
    s = 10*(dmg);
    k = fix(s);

    if (k <= 0) {               /*bounds of table for extrapolation*/
        k = 1;
    }
    else if (k >= 10){
        k = 9;
    }

    da = s - k;         /* amount from cnosest lower grid point*/

    /* L = k + fix(1.1*(da/fabs(da))); */ 
    L = k + fix(1.1*sign(da));

    k = k + 1;
    L = L + 1;

    cna00 = Cna_wing00[k-1] + fabs(da)*(Cna_wing00[L-1]-Cna_wing00[k-1]);
    cna08 = Cna_wing08[k-1] + fabs(da)*(Cna_wing08[L-1]-Cna_wing08[k-1]);    
    cnb00 = Cnb_wing00[k-1] + fabs(da)*(Cnb_wing00[L-1]-Cnb_wing00[k-1]);
    cnb08 = Cnb_wing08[k-1] + fabs(da)*(Cnb_wing08[L-1]-Cnb_wing08[k-1]);
    cnp00 = Cnp_wing00[k-1] + fabs(da)*(Cnp_wing00[L-1]-Cnp_wing00[k-1]);
    cnp08 = Cnp_wing08[k-1] + fabs(da)*(Cnp_wing08[L-1]-Cnp_wing08[k-1]);
    cnq00 = Cnq_wing00[k-1] + fabs(da)*(Cnq_wing00[L-1]-Cnq_wing00[k-1]);
    cnq08 = Cnq_wing08[k-1] + fabs(da)*(Cnq_wing08[L-1]-Cnq_wing08[k-1]);    
    cnr00 = Cnr_wing00[k-1] + fabs(da)*(Cnr_wing00[L-1]-Cnr_wing00[k-1]);
    cnr08 = Cnr_wing08[k-1] + fabs(da)*(Cnr_wing08[L-1]-Cnr_wing08[k-1]);
    
    cnawing = ((mach/0.8)*cna08+(0.8-mach)/0.8*cna00);
    cnbwing = ((mach/0.8)*cnb08+(0.8-mach)/0.8*cnb00)-((mach/0.8)*Cnb_wing08[10]+(0.8-mach)/0.8*Cnb_wing00[10]);
    cnpwing = ((mach/0.8)*cnp08+(0.8-mach)/0.8*cnp00)/((mach/0.8)*Cnp_wing08[10]+(0.8-mach)/0.8*Cnp_wing00[10]);
    cnqwing = ((mach/0.8)*cnq08+(0.8-mach)/0.8*cnq00);
    cnrwing = ((mach/0.8)*cnr08+(0.8-mach)/0.8*cnr00)/((mach/0.8)*Cnr_wing08[10]+(0.8-mach)/0.8*Cnr_wing00[10]);
    
    coeff[0] = cnawing;
    coeff[1] = cnbwing;
    coeff[2] = cnpwing;
    coeff[3] = cnqwing;
    coeff[4] = cnrwing;
}
