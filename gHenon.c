double dist_plus()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));
    return (a1*pow(J1,2))/(1 + (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) - 
   (a2*(1 - pow(e2,2)))/(1 + (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + 
           sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)));

}

double dist_minus()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));
	return (a1*pow(J1,2))/(1 - (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) - 
   (a2*(1 - pow(e2,2)))/(1 - (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + 
           sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)));
}


double dist_plus_dJ()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));
	return -((a1*pow(J1,2)*((sqrt(1 - pow(J1,2))*cos(omega1)*((pow(Jz1,2)*cos(I2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2))) + (Jz1*cos(Omega1 - Omega2)*sin(I2))/pow(J1,2)))/
           sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)) + 
          (sqrt(1 - pow(J1,2))*(-((Jz1*cos(I2))/pow(J1,2)) + (pow(Jz1,2)*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*
             ((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
             (cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5) - 
          (J1*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           (sqrt(1 - pow(J1,2))*sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)))))/
      pow(1 + (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)) + 
   (2*a1*J1)/(1 + (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) + 
   (a2*(1 - pow(e2,2))*((e2*(-(cos(omega2)*(-((pow(Jz1,2)*cos(I2)*cos(Omega1 - Omega2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))) - (Jz1*sin(I2))/pow(J1,2))) + 
             (pow(Jz1,2)*sin(omega2)*sin(Omega1 - Omega2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)) + 
        (e2*(-((Jz1*cos(I2))/pow(J1,2)) + (pow(Jz1,2)*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*
           ((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
           (-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
         pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
    pow(1 + (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2);
}

double dist_plus_dJz()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return -((a1*pow(J1,2)*((sqrt(1 - pow(J1,2))*cos(omega1)*(-((Jz1*cos(I2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))) - (cos(Omega1 - Omega2)*sin(I2))/J1))/
           sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)) + 
          (sqrt(1 - pow(J1,2))*(cos(I2)/J1 - (Jz1*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*
             ((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
             (cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
      pow(1 + (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)) + 
   (a2*(1 - pow(e2,2))*((e2*(-(cos(omega2)*((Jz1*cos(I2)*cos(Omega1 - Omega2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2))) + sin(I2)/J1)) - 
             (Jz1*sin(omega2)*sin(Omega1 - Omega2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))))/sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))\
         + (e2*(cos(I2)/J1 - (Jz1*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
           (-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
         pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
    pow(1 + (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2);
}


double dist_plus_domega()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return -((a1*pow(J1,2)*sqrt(1 - pow(J1,2))*(-((sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1)*sin(omega1)) + cos(omega1)*sin(I2)*sin(Omega1 - Omega2)))/
     (sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))*
       pow(1 + (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
          sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)));
}

double dist_plus_dOmega()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return -((a1*pow(J1,2)*((sqrt(1 - pow(J1,2))*(cos(Omega1 - Omega2)*sin(I2)*sin(omega1) + (Jz1*cos(omega1)*sin(I2)*sin(Omega1 - Omega2))/J1))/
           sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)) - 
          (sqrt(1 - pow(J1,2))*sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(I2)*((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*sin(Omega1 - Omega2)*
             (cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
      pow(1 + (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)) + 
   (a2*(1 - pow(e2,2))*((e2*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(omega2) - sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(omega2)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)) - 
        (e2*sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(I2)*((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*sin(Omega1 - Omega2)*
           (-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
         pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
    pow(1 + (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2);

}


double dist_minus_dJ()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return -((a1*pow(J1,2)*(-((sqrt(1 - pow(J1,2))*cos(omega1)*((pow(Jz1,2)*cos(I2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2))) + (Jz1*cos(Omega1 - Omega2)*sin(I2))/pow(J1,2)))/
             sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) - 
          (sqrt(1 - pow(J1,2))*(-((Jz1*cos(I2))/pow(J1,2)) + (pow(Jz1,2)*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*
             ((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
             (cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5) + 
          (J1*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           (sqrt(1 - pow(J1,2))*sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)))))/
      pow(1 - (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)) + 
   (2*a1*J1)/(1 - (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) + 
   (a2*(1 - pow(e2,2))*(-((e2*(-(cos(omega2)*(-((pow(Jz1,2)*cos(I2)*cos(Omega1 - Omega2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))) - (Jz1*sin(I2))/pow(J1,2))) + 
               (pow(Jz1,2)*sin(omega2)*sin(Omega1 - Omega2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))))/
           sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) - 
        (e2*(-((Jz1*cos(I2))/pow(J1,2)) + (pow(Jz1,2)*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,3)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*
           ((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
           (-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
         pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
    pow(1 - (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2);
}

double dist_minus_dJz()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return -((a1*pow(J1,2)*(-((sqrt(1 - pow(J1,2))*cos(omega1)*(-((Jz1*cos(I2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))) - (cos(Omega1 - Omega2)*sin(I2))/J1))/
             sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) - 
          (sqrt(1 - pow(J1,2))*(cos(I2)/J1 - (Jz1*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*
             ((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
             (cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
      pow(1 - (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)) + 
   (a2*(1 - pow(e2,2))*(-((e2*(-(cos(omega2)*((Jz1*cos(I2)*cos(Omega1 - Omega2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2))) + sin(I2)/J1)) - 
               (Jz1*sin(omega2)*sin(Omega1 - Omega2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2)))))/
           sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) - 
        (e2*(cos(I2)/J1 - (Jz1*cos(Omega1 - Omega2)*sin(I2))/(pow(J1,2)*sqrt(1 - pow(Jz1,2)/pow(J1,2))))*((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*
           (-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
         pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
    pow(1 - (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2);
}

double dist_minus_domega()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return 	(a1*pow(J1,2)*sqrt(1 - pow(J1,2))*(-((sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1)*sin(omega1)) + cos(omega1)*sin(I2)*sin(Omega1 - Omega2)))/
   (sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))*
     pow(1 - (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
        sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2));
}

double dist_minus_dOmega()
{
	double Omega1=Om1;
	double I2=acos(Jz2/J2),Omega2=Om2,e2=sqrt(1-(J2*J2));

	return -((a1*pow(J1,2)*(-((sqrt(1 - pow(J1,2))*(cos(Omega1 - Omega2)*sin(I2)*sin(omega1) + (Jz1*cos(omega1)*sin(I2)*sin(Omega1 - Omega2))/J1))/
             sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) + 
          (sqrt(1 - pow(J1,2))*sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(I2)*((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*sin(Omega1 - Omega2)*
             (cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
           pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
      pow(1 - (sqrt(1 - pow(J1,2))*(cos(omega1)*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2) - (Jz1*cos(Omega1 - Omega2)*sin(I2))/J1) + sin(I2)*sin(omega1)*sin(Omega1 - Omega2)))/
         sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2)) + 
   (a2*(1 - pow(e2,2))*(-((e2*(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(omega2) - sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(omega2)*sin(Omega1 - Omega2)))/
           sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2))) + 
        (e2*sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(I2)*((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2))*sin(Omega1 - Omega2)*
           (-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
         pow(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2),1.5)))/
    pow(1 - (e2*(-(cos(omega2)*(-(sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(I2)*cos(Omega1 - Omega2)) + (Jz1*sin(I2))/J1)) + sqrt(1 - pow(Jz1,2)/pow(J1,2))*sin(omega2)*sin(Omega1 - Omega2)))/
       sqrt(1 - pow((Jz1*cos(I2))/J1 + sqrt(1 - pow(Jz1,2)/pow(J1,2))*cos(Omega1 - Omega2)*sin(I2),2)),2);
}
