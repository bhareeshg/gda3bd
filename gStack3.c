#define STACK_SIZE 10000
double hamiltorian();
double updateCoefficients(const double *);
int dxbdt_henon( double ,const double *, double *, void *);
int dxbdt_henon_discrete( double ,const double *, double *, void *);
int updateDSign(double *);
double stack[STACK_SIZE];
int top=-1;
int dsign[4];
extern int dcur_sign_flag;
int push(double value)
{
	if(top<STACK_SIZE)
	{
		top++;
		stack[top]=value;
		return 1;
	}
	else
		return 0;
	
}

double pop()
{
	if(top >=0)
	{
		top--;
		return stack[top+1]; 
	}
}


int pullStackLimits(double * currentdistp,double *torg,gsl_odeiv2_driver * d)
{
	double ttorg[7];

	double to;
	double cH,prevH,prevdist,tol=1e-9;
	double dcrossplus,dcrossminus;
	int status;
	updateCoefficients(torg);
	cH=hamiltorian();	

	while(top>=0)
	{
		ttorg[0]=torg[0];
		ttorg[1]=torg[1];
		ttorg[2]=torg[2];
		ttorg[3]=torg[3];
		ttorg[4]=torg[4];
		ttorg[5]=torg[5];
		ttorg[6]=torg[6];

		to=pop();
		prevH=cH;
		prevdist=*currentdistp;

		//printf("from : %e to:%e\n",*currentdistp,to);
		fflush(stdout);
		status = gsl_odeiv2_driver_apply(d, currentdistp, to, ttorg);
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			  exit(0);
		}
	
		updateCoefficients(ttorg);
		cH=hamiltorian();
		if(fabs(prevH-cH)>tol)
		{
			//printf("try again adding:%e and: %e	error is:%e\n",to,(prevdist+*currentdistp)/2,fabs(prevH-cH));
			push(to);
			push((prevdist+*currentdistp)/2);
			*currentdistp=prevdist;
			continue;
		}
		else
		{
			torg[0]=ttorg[0];
			torg[1]=ttorg[1];
			torg[2]=ttorg[2];
			torg[3]=ttorg[3];
			torg[4]=ttorg[4];
			torg[5]=ttorg[5];
			torg[6]=ttorg[6];

			
			dcrossplus=dist_plus();
			dcrossminus=dist_minus();				
			//printf ("-->>>> %.5e	%.5e	%.5e	%.5e	%.5e %e %e %e	%e	%.15e\n", torg[6], torg[0], torg[1],torg[2],torg[3],torg[4],torg[5],dcrossplus,dcrossminus,cH);
		}
	}
	return 0;
}


double henon_step_stack(double *org,double t)
{
	double deltaxdt[7];
	int i;
	double cH;
	double torg[7];

	updateCoefficients(org);
	int mode;

	torg[0]=org[0];
	torg[1]=org[1];
	torg[2]=org[2];
	torg[3]=org[3];	
	torg[4]=dist_plus();
	torg[5]=dist_minus();
	torg[6]=t;

	if(fabs(torg[4])<fabs(torg[5]))
		mode=4;
	if(fabs(torg[4])>fabs(torg[5]))
		mode=5;

	double distFinal=0;
	double N1=10;
	double N2=5;
	double N3=3;
	double ddist=-torg[mode];
	double disti;
	if(torg[mode]<0)
		disti=-0;
	else
		disti=0;
	double ti;
	gsl_odeiv2_system sys = {dxbdt_henon_discrete, NULL, 7, &mode};

	gsl_odeiv2_driver * d =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
		                  ddist,1e-6, 1e-6);
	
	gsl_odeiv2_system sys1 = {dxbdt_henon, NULL, 7, &mode};
	gsl_odeiv2_driver * d1 =
	gsl_odeiv2_driver_alloc_y_new (&sys1, gsl_odeiv2_step_rkf45,
		                  ddist,1e-6, 1e-6);

	double rv[7];
	//printf("Henon stack step\n");
	//fflush(stdout);
	double curdist=torg[mode];
	int status;
	
	double prevdist=curdist;
	abs_int_err=1e-10;rel_int_err=1e-5;
	updateCoefficients(torg);
	cH=hamiltorian();
	double prevH=cH;
	double ttorg[7];
	int sign=curdist>0?1:-1;
	double tol=1e-12;
	double diststage=0.01;
	
	double dplus=dist_plus();
	double dminus=dist_minus();
	double dcur=fabs(dplus)<fabs(dminus)?dplus:dminus;
	dcur_sign_flag=dcur>0?1:-1;
	updateDSign(org);


	double r=pow(fabs(diststage/curdist),1/N1);
	//printf("r value is: %e\n",r);
	for(i=N1;i>0;i--)
	{
		//printf("pushed:: %e	top:%d\n",curdist-((curdist-sign*diststage)*i/N1),top);
		//push(curdist-((curdist-sign*diststage)*i/N1));
		//printf("pushed:: %e	top:%d\n",curdist*pow(r,i),top);
		push(curdist*pow(r,i));
	}
	pullStackLimits(&curdist,torg,d1);
	
	for(i=0;i<N2;i++)
	{
		//printf("pushed:: %e	top:%d\n",i*curdist/N2,top);
		push((curdist)*i/N2);
	}
	pullStackLimits(&curdist,torg,d);

	gsl_odeiv2_driver_reset(d);
	diststage=0.5;
	dcur_sign_flag*=-1;
	for(i=N3;i>=1;i--)
	{
		//printf("pushed:: %e	top:%d\n",curdist-((curdist+sign*diststage)*i/N1),top);
		push(curdist-((curdist+sign*diststage)*i/N1));
	}
	pullStackLimits(&curdist,torg,d);
	
	org[0]=torg[0];
	org[1]=torg[1];
	org[2]=torg[2];
	org[3]=torg[3];

	return torg[6];
}

