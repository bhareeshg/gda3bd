#include <gsl/gsl_integration.h>

struct params_struct_gsl{
	double (* integrand)(double,void *);
	double c;
	double d;
};
struct params_struct_gsl *inp_struct;

gsl_integration_cquad_workspace * w_outer;
gsl_integration_cquad_workspace * w_inner;
double *innerArgument;
size_t soCalls;
extern double abs_int_err, rel_int_err;
double f_inner(double , void *);

int initializeGSLstruct()
{
	inp_struct=(struct params_struct_gsl *)malloc(sizeof(struct params_struct_gsl));
	w_outer=gsl_integration_cquad_workspace_alloc (100);
	w_inner=gsl_integration_cquad_workspace_alloc (100);
	innerArgument=malloc(sizeof(double));
}

int freeGSLstruct()
{
	free(inp_struct);
	gsl_integration_cquad_workspace_free (w_outer);
	gsl_integration_cquad_workspace_free (w_inner);
	free(innerArgument);
}

double f_inner(double x, void * params) {
	gsl_integration_cquad_workspace * w= w_inner;
	innerArgument[0]=x;
	gsl_function F;
	F.function = inp_struct->integrand;
	F.params = innerArgument;
	double result, error;
	size_t nevals=100;
	gsl_integration_cquad (&F, inp_struct->c, inp_struct->d, abs_int_err, rel_int_err,w, &result, &error,&nevals);
	//soCalls+=nevals;
	//printf("nevals_inner::	%ld\n",nevals);
	return result;
}

double doubleIntegrateGSL(double (* integrand)(double,void *),double a,double b, double c,double d)
{
	soCalls=0;
	gsl_integration_cquad_workspace * w= w_outer;
	double result, error;
	inp_struct->integrand=integrand;
	inp_struct->c=c;
	inp_struct->d=d;
	gsl_function F;
	F.function = &f_inner;
	F.params = NULL;
	size_t nevals=100;
	gsl_integration_cquad (&F, a,b, abs_int_err, rel_int_err,w, &result, &error,&nevals);
	//printf("Result:	%e, Error: %e\n",result,error);
	//printf("nevals_outer::	%ld	soc::%ld\n",nevals,soCalls);
	return result;
}
