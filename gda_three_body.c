#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_deriv.h>
#include "gSL2cquad.c"
#include "gJmatrix.c"
#include "gHenon.c"
#include "gStack3.c"
#define PI 3.1415926535897932384626433832795028841971693993751
#define sign(x) x>0?1:-1
double A1,A2,B1,B2,C1,C2,D1,D2,E1,E2,F1,F2;
double J1p1,J1p2,J1p3,J1p4,J2p2,J2p3;
double e1fJ1, e2fJ2;
double Jz1p1,Jz1p2;
double coef1,coef2,coef3,coef4;
double a1,a2,J2,Jz2,omega2,Om2;
double sinI2,cosI2;
double cOm1,comega1;
double sOm1,somega1;
double dcoef1J,dcoef2J,dcoef3J,dcoef4J;
double dcoef1Jz,dcoef2Jz,dcoef3Jz,dcoef4Jz;
double dcoef1Om,dcoef2Om,dcoef3Om,dcoef4Om;

double omega1,Om1,J1,Jz1;
double abs_int_err, rel_int_err;

double G = 4*PI*PI; 
double m2 = 3e-5; 
double M = 1;

double omega2dot,Omega2dot;
double S,S2;

int dcur_sign_flag;

double updateCoefficients(const double *oarg)
{
	J1=oarg[0];
	Jz1=oarg[1];
	omega1=oarg[2];
	Om1=oarg[3];
	
	double cosI1=Jz1/J1;
	double sinI1=sqrt(1-(cosI1*cosI1));
	cOm1=cos(Om1);
	comega1=cos(omega1);
	somega1=sin(omega1);
	sOm1=sin(Om1);

	J1p1=J1;
	J1p2=J1*J1;
	J1p3=J1p2*J1;
	J1p4=J1p3*J1;	
	e1fJ1=sqrt(1-J1p2);
	Jz1p2=Jz1*Jz1;
	Jz1p1=Jz1;

	A1=cos(Om1)*cos(omega1)-(sin(omega1)*sin(Om1)*cosI1);
	B1=-(cos(Om1)*sin(omega1)+sin(Om1)*cos(omega1)*cosI1);
	C1=(sin(Om1)*cos(omega1)) + (cos(Om1)*sin(omega1)*cosI1);
	D1=-sin(omega1)*sin(Om1) + cos(omega1)*cos(Om1)*cosI1;
	E1=sin(omega1)*sinI1;
	F1=cos(omega1)*sinI1;
	
	coef1=A1*A2+ C1*C2+ E1*E2;
	coef2=A1*B2+ C1*D2+ E1*F2;
	coef3=A2*B1+ D1*C2+ F1*E2;
	coef4=B1*B2+ D1*D2+ F1*F2;

	//coefficients for Jz
	double dAJz1=-somega1*sOm1/J1p1;
	double dBJz1=-comega1*sOm1/J1p1;
	double dCJz1=somega1*cOm1/J1p1;
	double dDJz1=cOm1*comega1/J1p1;
	double dEJz1=-somega1*Jz1p1/J1p2/sqrt(1-Jz1p2/J1p2);
	double dFJz1=-comega1*Jz1p1/J1p2/sqrt(1-Jz1p2/J1p2);

	dcoef1Jz=dAJz1*A2+ dCJz1*C2+ dEJz1*E2;
	dcoef2Jz=dAJz1*B2+ dCJz1*D2+ dEJz1*F2;
	dcoef3Jz=A2*dBJz1+ dDJz1*C2+ dFJz1*E2;
	dcoef4Jz=dBJz1*B2+ dDJz1*D2+ dFJz1*F2;

	//coefficients for Om
	dcoef1Om=-C1*A2+ A1*C2;
	dcoef2Om=-C1*B2+ A1*D2;
	dcoef3Om=-A2*D1+ B1*C2;
	dcoef4Om=-D1*B2+ B1*D2;

	//coefficients for J
	double dA1J=somega1*sOm1*Jz1p1/J1p2;
	double dB1J=comega1*sOm1*Jz1p1/J1p2;
	double dC1J=-somega1*cOm1*Jz1p1/J1p2;
	double dD1J=-cOm1*comega1*Jz1p1/J1p2;
	double dE1J=somega1*Jz1p2/J1p3/sqrt(1-Jz1p2/J1p2);
	double dF1J=comega1*Jz1p2/J1p3/sqrt(1-Jz1p2/J1p2);

	dcoef1J=dA1J*A2+ dC1J*C2+ dE1J*E2;
	dcoef2J=dA1J*B2+ dC1J*D2+ dE1J*F2;
	dcoef3J=A2*dB1J+ dD1J*C2+ dF1J*E2;
	dcoef4J=dB1J*B2+ dD1J*D2+ dF1J*F2;
	
	
}

double hamiltorianInt(double f1,void *pf2)
{
	double *f2a=(double *)pf2;
	double f2=f2a[0];
	double cospsi=(coef1*cos(f1)*cos(f2)) + (coef2*cos(f1)*sin(f2)) +(coef3*cos(f2)*sin(f1)) + (coef4*sin(f1)*sin(f2));
	double r=a1*J1p2/(1+e1fJ1*cos(f1));
	double R=a2*J2p2/(1+e2fJ2*cos(f2));
	double dM1=J1p3/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1));
	double dM2=J2p3/(1+e2fJ2*cos(f2))/(1+e2fJ2*cos(f2));
	return dM1*dM2/sqrt(r*r + R*R - 2*r*R*cospsi);
}

double dHbdOmegaInt(double f1,void *pf2)
{
	double *f2a=(double *)pf2;
	double f2=f2a[0];
	double cospsi=(coef1*cos(f1)*cos(f2)) + (coef2*cos(f1)*sin(f2)) +(coef3*cos(f2)*sin(f1)) + (coef4*sin(f1)*sin(f2));
	double r=a1*J1p2/(1+e1fJ1*cos(f1));
	double R=a2*J2p2/(1+e2fJ2*cos(f2));
	double dM1=J1p3/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1));
	double dM2=J2p3/(1+e2fJ2*cos(f2))/(1+e2fJ2*cos(f2));
	double dcospsi=(dcoef1Om*cos(f1)*cos(f2)) + (dcoef2Om*cos(f1)*sin(f2)) +(dcoef3Om*cos(f2)*sin(f1)) + (dcoef4Om*sin(f1)*sin(f2));
	return dM1*dM2*r*R*dcospsi/pow(r*r + R*R - 2*r*R*cospsi,1.5);
}


double dHbdJInt(double f1,void *pf2)
{
	double *f2a=(double *)pf2;
	double f2=f2a[0];
	double cospsi=(coef1*cos(f1)*cos(f2)) + (coef2*cos(f1)*sin(f2)) +(coef3*cos(f2)*sin(f1)) + (coef4*sin(f1)*sin(f2));
	double r=a1*J1p2/(1+e1fJ1*cos(f1));
	double R=a2*J2p2/(1+e2fJ2*cos(f2));
	double dM1=J1p3/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1));
	double ddM11=3*J1p2/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1));
	double ddM12=2*J1p4*cos(f1)/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1))/e1fJ1;
	double dM2=J2p3/(1+e2fJ2*cos(f2))/(1+e2fJ2*cos(f2));
	double dr=((2*a1*J1p1)/(1+e1fJ1*cos(f1)))+(a1*J1p3*cos(f1)/e1fJ1/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1)));
	double dcospsi=(dcoef1J*cos(f1)*cos(f2)) + (dcoef2J*cos(f1)*sin(f2)) +(dcoef3J*cos(f2)*sin(f1)) + (dcoef4J*sin(f1)*sin(f2));
	double v1=(ddM11+ddM12)*dM2/sqrt(r*r + R*R - 2*r*R*cospsi);
	double v2=-dM1*dM2*(r*dr-r*R*dcospsi-dr*R*cospsi)/pow(r*r + R*R - 2*r*R*cospsi,1.5);
	return v1+v2;
}

double dHbdJzInt(double f1,void *pf2)
{
	double *f2a=(double *)pf2;
	double f2=f2a[0];
	double cospsi=(coef1*cos(f1)*cos(f2)) + (coef2*cos(f1)*sin(f2)) +(coef3*cos(f2)*sin(f1)) + (coef4*sin(f1)*sin(f2));
	double r=a1*J1p2/(1+e1fJ1*cos(f1));
	double R=a2*J2p2/(1+e2fJ2*cos(f2));
	double dM1=J1p3/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1));
	double dM2=J2p3/(1+e2fJ2*cos(f2))/(1+e2fJ2*cos(f2));
	double dcospsi=(dcoef1Jz*cos(f1)*cos(f2)) + (dcoef2Jz*cos(f1)*sin(f2)) +(dcoef3Jz*cos(f2)*sin(f1)) + (dcoef4Jz*sin(f1)*sin(f2));
	return dM1*dM2*r*R*dcospsi/pow(r*r + R*R - 2*r*R*cospsi,1.5);
}



double dHbdomegaInt(double f1,void *pf2)
{
	double *f2a=(double *)pf2;
	double f2=f2a[0];
	double cospsi=(coef1*cos(f1)*cos(f2)) + (coef2*cos(f1)*sin(f2)) +(coef3*cos(f2)*sin(f1)) + (coef4*sin(f1)*sin(f2));
	double r=a1*J1p2/(1+e1fJ1*cos(f1));
	double R=a2*J2p2/(1+e2fJ2*cos(f2));
	double dM1=J1p3/(1+e1fJ1*cos(f1))/(1+e1fJ1*cos(f1));
	double dM2=J2p3/(1+e2fJ2*cos(f2))/(1+e2fJ2*cos(f2));
	double dcospsi=(coef3*cos(f1)*cos(f2)) + (coef4*cos(f1)*sin(f2)) +(-coef1*cos(f2)*sin(f1)) + (-coef2*sin(f1)*sin(f2));
	return dM1*dM2*r*R*dcospsi/pow(r*r + R*R - 2*r*R*cospsi,1.5);
}


double hamiltorian()
{
	//double Havg=-(G*m2/4/PI/PI)*gauss_legendre_2D_cube(128,hamiltorianInt,NULL,0,2*PI,0,2*PI);
	double Havg=-(G*m2/4/PI/PI)*doubleIntegrateGSL(hamiltorianInt,0,2*PI,0,2*PI);
	double H9planets=-G*M*S*(3*Jz1p2 - J1p2)/a1/8/J1p4/J1p1;
	double Hautonomous=-omega2dot*sqrt(G*M*a1)*J1p1-Omega2dot*sqrt(G*M*a1)*Jz1p1;
	return Havg+H9planets+Hautonomous;
}

double dHbdomega()
{
	//return -(G*m2/4/PI/PI)*gauss_legendre_2D_cube(128,dHbdomegaInt,NULL,0,2*PI,0,2*PI);
	return -(G*m2/4/PI/PI)*doubleIntegrateGSL(dHbdomegaInt,0,2*PI,0,2*PI);
}

double dHbdOmega()
{
	//return -(G*m2/4/PI/PI)*gauss_legendre_2D_cube(128,dHbdOmegaInt,NULL,0,2*PI,0,2*PI);
	return -(G*m2/4/PI/PI)*doubleIntegrateGSL(dHbdOmegaInt,0,2*PI,0,2*PI);
}

double dHbdJ()
{
	//double HavgdJ=-(G*m2/4/PI/PI)*gauss_legendre_2D_cube(128,dHbdJInt,NULL,0,2*PI,0,2*PI);
	//printf("avgHdj:: %e\n",HavgdJ);
	double HavgdJ=-(G*m2/4/PI/PI)*doubleIntegrateGSL(dHbdJInt,0,2*PI,0,2*PI);
	double H9planetsdJ=5.0*G*M*S*(3*Jz1p2-J1p2)/8.0/a1/J1p4/J1p2 + G*M*S/a1/J1p4/4;
	double HautonomousdJ=-omega2dot*sqrt(G*M*a1);
	return HavgdJ+H9planetsdJ+HautonomousdJ;
}

double dHbdJz()
{
	//double HavgdJz=-(G*m2/4/PI/PI)*gauss_legendre_2D_cube(128,dHbdJzInt,NULL,0,2*PI,0,2*PI);
	double HavgdJz=-(G*m2/4/PI/PI)*doubleIntegrateGSL(dHbdJzInt,0,2*PI,0,2*PI);
	double H9planetsdJz=-0.75*G*M*S*(Jz1p1)/a1/J1p4/J1p1;
	double HautonomousdJz=-Omega2dot*sqrt(G*M*a1);
	return HavgdJz+H9planetsdJz+HautonomousdJz;
}

int initializePlanet9Terms()
{
	cosI2=Jz2/J2;
	sinI2=sqrt(1-(cosI2*cosI2));
	J2p2=J2*J2;
	J2p3=J2p2*J2;
	e2fJ2=sqrt(1-J2p2);	
	A2=cos(Om2)*cos(omega2)-(sin(omega2)*sin(Om2)*cosI2);
	B2=-(cos(Om2)*sin(omega2)+sin(Om2)*cos(omega2)*cosI2);	
	C2=(sin(Om2)*cos(omega2)) + (cos(Om2)*sin(omega2)*cosI2);	
	D2=-sin(omega2)*sin(Om2) + cos(omega2)*cos(Om2)*cosI2;	
	E2=sin(omega2)*sinI2;
	F2=cos(omega2)*sinI2;	
	//omega2, Omega2 dot
	//S= 0.0009543*pow(5.20,2)/pow(a1,2)+0.0002857*pow(9.54,2)/pow(a1,2)+0.00004365*pow(19.19,2)/pow(a1,2)+0.00005149*pow(30.110387,2)/pow(a1,2);
	//S2= 0.0009543*pow(5.20,2)/pow(a2,2) + 0.0002857*pow(9.54,2)/pow(a2,2) + 0.00004365*pow(19.19,2)/pow(a2,2) + 0.00005149*pow(30.110387,2)/pow(a2,2);
	//printf("S=%e; S2=%e\n",S,S2);
	S=0; S2=0;
	omega2dot=2*3*G*S2/4/a2/J2p3/J2/sqrt(G*M*a2);
	Omega2dot=-3*G*S2/4/a2/J2p3/J2/sqrt(G*M*a2);
	return 1;
}

double * dxbdt( double *org, double *rv)
{
	
	updateCoefficients(org);
	double dOmdt=dHbdJz();
	double domegadt=dHbdJ();
	double dJdt=-dHbdomega();
	double dJzdt=-dHbdOmega();
	rv[2]=domegadt/sqrt(G*M*a1);
	rv[3]=dOmdt/sqrt(G*M*a1);
	rv[0]=dJdt/sqrt(G*M*a1);
	rv[1]=dJzdt/sqrt(G*M*a1);
	//printf("continuous::	%e	%e	%e	%e\n",rv[0],rv[1],rv[2],rv[3]);	
	return rv;

}


int updateDSign(double *org)
{
	double rv[4];
	dxbdt(org,rv);
	
	dsign[0]=rv[0]>0?1:-1;
	dsign[1]=rv[1]>0?1:-1;
	dsign[2]=rv[2]>0?1:-1;
	dsign[3]=rv[3]>0?1:-1;
	return 1;
}

struct hParamStruct{
	int paramIndex;
	const double *org;
};

double dHdparam(double x, void * params)
{
	struct hParamStruct *dp=(struct hParamStruct *)(params); /* avoid unused parameter warning */
	const double *org=dp->org;
	double torg[4];
	torg[0]=org[0];torg[1]=org[1];torg[2]=org[2];torg[3]=org[3];
	torg[dp->paramIndex]=x;

	//updateCoefficients(torg);
	//double cH=hamiltorian();
	//double dcrossplus=dist_plus();
	//double dcrossminus=dist_minus();
	//printf("paramIndex:: %d\n",dp->paramIndex);				
	//printf ("dHdparam oe:%.15e	%.15e	%.15e	%.15e %.15e %.15e	%.15e\n", torg[0], torg[1],torg[2],torg[3],dcrossplus,dcrossminus,cH);
	updateCoefficients(torg);
	double H=hamiltorian();	
	return H;
}


double find_discrete_derivative(const double *org,int index)
{
	double h=1e-2;
	double prevResult;
	double prevError=1;
	double result,abserr;
	double error;
	
	gsl_function F;
	struct hParamStruct lparam;
	lparam.org=org;
	F.function = &dHdparam;
	F.params = &lparam;

	lparam.paramIndex=index;

	//printf("discrete derivative---calculation	index::%d\n",index);

	if(dsign[index]==-1)
		gsl_deriv_backward (&F, org[index], h, &result, &abserr);//dHbdJz();
	if(dsign[index]==1)
		gsl_deriv_forward (&F, org[index], h, &result, &abserr);//dHbdJz();
	error=fabs(result);
	//printf("#####  =++++++ error: %e prevError:%e\n",error,prevError);

	while(error<prevError)
	{
		prevResult=result;
		prevError=abserr;

		h/=2;
		if(dsign[index]==-1)
			gsl_deriv_backward (&F, org[index], h, &result, &abserr);//dHbdJz();
		if(dsign[index]==1)
			gsl_deriv_forward (&F, org[index], h, &result, &abserr);//dHbdJz();
		
		//printf("h: %e	result:%e	absrr:%e\n",h,result,abserr);

		error=abserr;//fabs(result-prevResult);
		
	}
	//printf("-> resultF/abserrF:: %e/%e (index: %d)\n",prevResult,prevError,index);
	return prevResult;

}

int
func_discrete (double t, const double org[], double f[],
      void *params)
{
  	(void)(t); /* avoid unused parameter warning */
  	//double mu = *(double *)params; /*unused parameter warning*/
	//printf ("--->>%.5e %.5e %.5e	%.5e	%.5e\n", t, org[0], org[1],org[2],org[3]);
	//printf("Func discrete called\n");
  	updateCoefficients(org);
	double dddJplus=dist_plus_dJ();
	double dddomegaplus=dist_plus_domega();
	double dddOmegaplus=dist_plus_dOmega();
	double dddJzplus=dist_plus_dJz();

	
	double dddJminus=dist_minus_dJ();
	double dddomegaminus=dist_minus_domega();
	double dddOmegaminus=dist_minus_dOmega();
	double dddJzminus=dist_minus_dJz();

	
	double dplus=dist_plus();
	double dminus=dist_minus();

	double dcur=fabs(dplus)<fabs(dminus)?dplus:dminus;
	int flag=fabs(dplus)<fabs(dminus)?1:2;
	dcur_sign_flag=sign(dcur);
	
	double dddJcur,dddomegacur,dddOmegacur,dddJzcur;
	
	if(flag==1)
	{
		dddJcur=dddJplus;
		dddomegacur=dddomegaplus;
		dddOmegacur=dddOmegaplus;
		dddJzcur=dddJzplus;
	}
	if(flag==2)
	{
		dddJcur=dddJminus;
		dddomegacur=dddomegaminus;
		dddOmegacur=dddOmegaminus;
		dddJzcur=dddJzminus;
	}

	if(dcur_sign_flag>0)
	{
		dsign[0]=sign(dddJcur);
		dsign[1]=sign(dddJzcur);
		dsign[2]=sign(dddomegacur);
		dsign[3]=sign(dddOmegacur);
	}
	else
	{
		dsign[0]=-sign(dddJcur);
		dsign[1]=-sign(dddJzcur);
		dsign[2]=-sign(dddomegacur);
		dsign[3]=-sign(dddOmegacur);
	}


	double dOmdt=find_discrete_derivative(org,1);
	//printf("dOmdt:: %e\n",dOmdt);	
	double domegadt=find_discrete_derivative(org,0);
	//printf("domegadt:: %e\n",domegadt);	
	double dJdt=-find_discrete_derivative(org,2);
	//printf("dJdt:: %e\n",dJdt);		
	double dJzdt=-find_discrete_derivative(org,3);
	//printf("dJzdt:: %e\n",dJzdt);

	f[2]=domegadt/sqrt(G*M*a1);
	f[3]=dOmdt/sqrt(G*M*a1);
	f[0]=dJdt/sqrt(G*M*a1);
	f[1]=dJzdt/sqrt(G*M*a1);
	return GSL_SUCCESS;
}

int
func (double t, const double org[], double f[],
      void *params)
{
  	(void)(t); /* avoid unused parameter warning */
  	//double mu = *(double *)params; /*unused parameter warning*/
	//printf ("--->>%.5e %.5e %.5e	%.5e	%.5e\n", t, org[0], org[1],org[2],org[3]);
  	updateCoefficients(org);
	
	double dplus=dist_plus();
	double dminus=dist_minus();

	double dcur=fabs(dplus)<fabs(dminus)?dplus:dminus;
	if(fabs(dcur)<0.01)
		return func_discrete(t,org,f,params);

	double dOmdt=dHbdJz();
	double domegadt=dHbdJ();
	double dJdt=-dHbdomega();
	double dJzdt=-dHbdOmega();

	f[2]=domegadt/sqrt(G*M*a1);
	f[3]=dOmdt/sqrt(G*M*a1);
	f[0]=dJdt/sqrt(G*M*a1);
	f[1]=dJzdt/sqrt(G*M*a1);
	return GSL_SUCCESS;
}



int jac (double t, const double y[], double *dfdy,double dfdt[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	//double mu = *(double *)params;
	gsl_matrix_view dfdy_mat
	= gsl_matrix_view_array (dfdy, 4, 4);
	
	double j11value=doubleIntegrateGSL(j11Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);;
	double j33value=-j11value;

	double j14value=doubleIntegrateGSL(j14Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);;
	double j23value=j14value;

	double j12value=doubleIntegrateGSL(j12Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j43value=-j12value;	

	double j22value=doubleIntegrateGSL(j22Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j44value=-j22value;

	double j21value=doubleIntegrateGSL(j21Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j34value=-j21value;

	double j32value=doubleIntegrateGSL(j32Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j41value=j32value;

	double j13value=doubleIntegrateGSL(j13Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j24value=doubleIntegrateGSL(j24Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j42value=doubleIntegrateGSL(j42Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);
	double j31value=doubleIntegrateGSL(j31Int,0,2*PI,0,2*PI)/sqrt(G*M*a1);

	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set (m, 0, 0, j11value);
	gsl_matrix_set (m, 0, 1, j12value);
	gsl_matrix_set (m, 0, 2, j13value);
	gsl_matrix_set (m, 0, 3, j14value);

	gsl_matrix_set (m, 1, 0, j21value);
	gsl_matrix_set (m, 1, 1, j22value);
	gsl_matrix_set (m, 1, 2, j23value);
	gsl_matrix_set (m, 1, 3, j24value);
	
	
	gsl_matrix_set (m, 2, 0, j31value);
	gsl_matrix_set (m, 2, 1, j32value);
	gsl_matrix_set (m, 2, 2, j33value);
	gsl_matrix_set (m, 2, 3, j34value);

	gsl_matrix_set (m, 3, 0, j41value);
	gsl_matrix_set (m, 3, 1, j42value);
	gsl_matrix_set (m, 3, 2, j43value);
	gsl_matrix_set (m, 3, 3, j44value);
	
	
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;

	return GSL_SUCCESS;
}



int dxbdt_henon_discrete( double t,const double org[], double rv[], void *params)
{
	updateCoefficients(org);
	double dddJplus=dist_plus_dJ();
	double dddomegaplus=dist_plus_domega();
	double dddOmegaplus=dist_plus_dOmega();
	double dddJzplus=dist_plus_dJz();

	
	double dddJminus=dist_minus_dJ();
	double dddomegaminus=dist_minus_domega();
	double dddOmegaminus=dist_minus_dOmega();
	double dddJzminus=dist_minus_dJz();

	
	double dplus=dist_plus();
	double dminus=dist_minus();

	double dcur=fabs(dplus)<fabs(dminus)?dplus:dminus;
	int flag=fabs(dplus)<fabs(dminus)?1:2;
	
	double dddJcur,dddomegacur,dddOmegacur,dddJzcur;
	
	if(flag==1)
	{
		dddJcur=dddJplus;
		dddomegacur=dddomegaplus;
		dddOmegacur=dddOmegaplus;
		dddJzcur=dddJzplus;
	}
	if(flag==2)
	{
		dddJcur=dddJminus;
		dddomegacur=dddomegaminus;
		dddOmegacur=dddOmegaminus;
		dddJzcur=dddJzminus;
	}

	if(dcur_sign_flag>0)
	{
		dsign[0]=sign(dddJcur);
		dsign[1]=sign(dddJzcur);
		dsign[2]=sign(dddomegacur);
		dsign[3]=sign(dddOmegacur);
	}
	else
	{
		dsign[0]=-sign(dddJcur);
		dsign[1]=-sign(dddJzcur);
		dsign[2]=-sign(dddomegacur);
		dsign[3]=-sign(dddOmegacur);
	}


	double dOmdt=find_discrete_derivative(org,1);
	//printf("dOmdt:: %e\n",dOmdt);	
	double domegadt=find_discrete_derivative(org,0);
	//printf("domegadt:: %e\n",domegadt);	
	double dJdt=-find_discrete_derivative(org,2);
	//printf("dJdt:: %e\n",dJdt);		
	double dJzdt=-find_discrete_derivative(org,3);
	//printf("dJzdt:: %e\n",dJzdt);

	int mode = *(int *)params;
	double K;

	rv[2]=domegadt/sqrt(G*M*a1);
	rv[3]=dOmdt/sqrt(G*M*a1);
	rv[0]=dJdt/sqrt(G*M*a1);
	rv[1]=dJzdt/sqrt(G*M*a1);

	
	//printf("org:: %e	%e	%e	%e	dplus:%e	dminus:%e\n",org[0],org[1],org[2],org[3],dplus,dminus);
	//printf("rv:: %e	%e	%e	%e\n",rv[0],rv[1],rv[2],rv[3]);
	rv[4]=((dddJplus*dJdt)+(dddJzplus*dJzdt)+(dddomegaplus*domegadt)+(dddOmegaplus*dOmdt))/sqrt(G*M*a1);
	rv[5]=((dddJminus*dJdt)+(dddJzminus*dJzdt)+(dddomegaminus*domegadt)+(dddOmegaminus*dOmdt))/sqrt(G*M*a1);

	rv[6]=1;
	K=1/rv[mode];
	rv[0]*=K;rv[1]*=K;rv[2]*=K;rv[3]*=K;rv[4]*=K;rv[5]*=K;rv[6]*=K;	

	return GSL_SUCCESS;

}

int dxbdt_henon( double t,const double org[], double rv[], void *params)
{
	(void)(t);
	updateCoefficients(org);
	
	double dplus=dist_plus();
	double dminus=dist_minus();
	//printf("dxdt_henon: org:: %e	%e	%e	%e	dplus:%e	dminus:%e\n",org[0],org[1],org[2],org[3],dplus,dminus);
	double dOmdt=dHbdJz();
	double domegadt=dHbdJ();
	double dJdt=-dHbdomega();
	double dJzdt=-dHbdOmega();
	int mode = *(int *)params;
	double K;

	rv[2]=domegadt/sqrt(G*M*a1);
	rv[3]=dOmdt/sqrt(G*M*a1);
	rv[0]=dJdt/sqrt(G*M*a1);
	rv[1]=dJzdt/sqrt(G*M*a1);
	
	double dddJplus=dist_plus_dJ();
	double dddomegaplus=dist_plus_domega();
	double dddOmegaplus=dist_plus_dOmega();
	double dddJzplus=dist_plus_dJz();

	rv[4]=((dddJplus*dJdt)+(dddJzplus*dJzdt)+(dddomegaplus*domegadt)+(dddOmegaplus*dOmdt))/sqrt(G*M*a1);
	
	
	double dddJminus=dist_minus_dJ();
	double dddomegaminus=dist_minus_domega();
	double dddOmegaminus=dist_minus_dOmega();
	double dddJzminus=dist_minus_dJz();

	rv[5]=((dddJminus*dJdt)+(dddJzminus*dJzdt)+(dddomegaminus*domegadt)+(dddOmegaminus*dOmdt))/sqrt(G*M*a1);

	rv[6]=1;
	K=1/rv[mode];
	//printf("a values:: %e	%e	%e	%e	%e	%e	%e\n",org[0],org[1],org[2],org[3],org[4],org[5],org[6]);
	//printf("rv values:: %e	%e	%e	%e	%e	%e	%e\n",rv[0],rv[1],rv[2],rv[3],rv[4],rv[5],rv[6]);
	
	rv[0]*=K;rv[1]*=K;rv[2]*=K;rv[3]*=K;rv[4]*=K;rv[5]*=K;rv[6]*=K;	
	//printf("rv values:: %e	%e	%e	%e	%e	%e	%e\n",rv[0],rv[1],rv[2],rv[3],rv[4],rv[5],rv[6]);	

	return GSL_SUCCESS;
}

int main (int argc, char* argv[])
{

	int i;
	
	if(argc!=4)
	{
		printf("number of arguments: %d/4\n",argc-1);
		printf("usage: ./a.out <IC file name> <outputfile name> <end time>\n");
		return 0;
	}

	char *inp_fname=argv[1];
	char *out_fname=argv[2];
	double end_time=atof(argv[3]);

	char line[512];
	abs_int_err=1e-10;rel_int_err=1e-5;

	double e2,i2,e1,i1;
	
	FILE* inp_fp = fopen(inp_fname, "r");
	FILE *out_fp=fopen(out_fname,"w");
	
	
	fgets(line,sizeof(line), inp_fp);
	fgets(line,sizeof(line), inp_fp);
	sscanf(line,"%lf	%lf	%lf	%lf	%lf	%lf\n",&a2,&e2,&i2,&omega2,&Om2,&m2);
	J2=sqrt(1-(e2*e2));
	Jz2=J2*cos(i2);

	fgets(line,sizeof(line),inp_fp);
	fgets(line,sizeof(line),inp_fp);
	double y[4],ty[4];
	double iH1;
	sscanf(line,"%lf	%lf	%lf	%lf	%lf",&a1,&e1,&i1,y+2,y+3);
	y[0]=sqrt(1-(e1*e1));	
	y[1]=sqrt(1-(e1*e1))*cos(i1);

	if(y[1]<0)
		y[1]*=-1;
	initializePlanet9Terms();
	initializeGSLstruct();
	gsl_odeiv2_system sys = {func, NULL, 4, NULL};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				1e6, 1e-6, 1e-10);
	

	//a1*=1000;
	//a2*=1000;
	double cH,dplus,dminus,dcur,dmin;
	double t = 0.0;
	double k=G*m2*a1*a1/8/a2/a2/a2/J2/J2/J2;
	double tk=sqrt(G*a1)/k;
	double tmax=end_time*tk;

	if(a1>a2)
	{
		double p2=pow(a1,1.5);
		tk=p2*(pow(M+m2,2)/M/m2)*pow(a1/a2,2);
		tmax=200.1*tk;
		printf("outer testpart\n");
	}
	double dt=tk/350;
	//printf("tK: %e tmax :%e dt: %e\n",tk,tmax,dt);
	double dttemp;
	double ti;
	int status;
	double h=1;
	double prevDistp,prevDistn,curDist,prevDist;
	

	updateCoefficients(y);
	cH=hamiltorian();
	double zerodot=dHbdJz();
	dplus=dist_plus();
	dminus=dist_minus();
	curDist=fabs(dplus)<fabs(dminus)?dplus:dminus;
	
	
	fprintf(out_fp,"a1	a2	e2	i2	omega2	Omega2	m2\n");
	fprintf(out_fp,"%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",a1,a2,sqrt(1-(J2*J2)),acos(Jz2/J2),omega2,Om2,m2);

	fprintf(out_fp,"time,e,i,omega,Omega,H,dmin\n");
	dmin=fabs(dplus)<fabs(dminus)?fabs(dplus):fabs(dminus);
	fprintf (out_fp,"%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n", t, sqrt(1-(y[0]*y[0])), acos(y[1]/y[0]),y[2],y[3],cH,dmin);
	i=0;
	double rv[4];

	while(t<end_time)
	{
		ti = t+ (h*dt);
		prevDistp=dplus;
		prevDistn=dminus;
		prevDist=curDist;
		ty[0]=y[0];ty[1]=y[1];ty[2]=y[2];ty[3]=y[3];
		status = gsl_odeiv2_driver_apply (d, &t, ti, ty);
		if (status != GSL_SUCCESS)
		  break;	
		updateCoefficients(ty);
		cH=hamiltorian();
		dplus=dist_plus();
		dminus=dist_minus();
		curDist=fabs(dplus)<fabs(dminus)?dplus:dminus;		

		if((prevDistp*dplus<0)||(prevDistn*dminus<0))
		{
			fflush(out_fp);
			t=(henon_step_stack(y,t));
			updateCoefficients(y);
			cH=hamiltorian();
			dplus=dist_plus();
			dminus=dist_minus();
			curDist=fabs(dplus)<fabs(dminus)?dplus:dminus;	
			gsl_odeiv2_driver_reset(d);
		}
		else
		{
			i++;
			y[0]=ty[0];y[1]=ty[1];y[2]=ty[2];y[3]=ty[3];
			
			dmin=fabs(dplus)<fabs(dminus)?fabs(dplus):fabs(dminus);	
			fprintf(out_fp,"%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n", t, sqrt(1-(y[0]*y[0])),acos(y[1]/y[0]),y[2],y[3],cH,dmin);
			h=1;
		}
		if(i%50==0)
			fflush(out_fp);

	}
	fclose(inp_fp);
	fclose(out_fp);
	freeGSLstruct();
  return 0;
}
