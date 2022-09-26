#include "DemirayUncoupled.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(DemirayUncoupled, FEUncoupledMaterial)
	ADD_PARAMETER(m_a, FE_RANGE_GREATER(0.0), "a");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER(0.0), "b");
END_FECORE_CLASS();

mat3ds DemirayUncoupled::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double a = m_a(mp);
	double b = m_b(mp);
	// deformation gradient
	double J = pt.m_J;
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();
	double I1 = B.tr();
	// calculate deviatoric stress
	return a*pow(J, -5.0/3.0)*(B - mat3dd(I1/3.0))*exp(0.5*b*(pow(J, -2.0/3.0)*I1-3.0));
}

tens4ds DemirayUncoupled::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double a = m_a(mp);
	double b = m_b(mp);
	// deformation gradient
	double J = pt.m_J;
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();
	// mat3ds C = pt.RightCauchyGreen();
	double I1 = B.tr();	
	// mat3ds Cm1 = C.inverse();


	mat3ds I(1,1,1,0,0,0);	// Identity
	tens4ds IxI = dyad1s(I);
	// tens4ds Cm1xCm1 = dyad1s(Cm1);
	// tens4ds Cm1xI = dyad1s(Cm1,I)/2.0; 
	// tens4ds IxCm1 = dyad1s(I,Cm1)/2.0;
	tens4ds I4  = dyad4s(I);
	// tens4ds I4th = pt.pull_back(I4)/J;

	tens4ds BxB = dyad1s(B);
	tens4ds IxB = dyad1s(I,B)/2.0; 
	tens4ds BxI = dyad1s(B,I)/2.0; 


	// double x = exp(0.5*b*(pow(J, -2.0/3.0)*I1-3.0));
	// mat3ds y = I-1.0/3.0*I1*Cm1;
	// double dxdI1 = b/2.0*pow(J, -2.0/3.0)*x;
	// double dxdI3 = -1.0/3.0*b/2.0*I1*pow(J, -8.0/3.0)*x;
	// mat3ds dydI1 = -1.0/3.0*Cm1;
	// tens4ds dydCm1 = -1.0/3.0*I1*I4;
	// mat3ds dSdI1 = a*pow(J, -2.0/3.0)*dxdI1*y+a*pow(J, -2.0/3.0)*x*dydI1;
	// mat3ds dSdI3 = -1.0/3.0*a*pow(J, -8.0/3.0)*x*y+a*pow(J, -2.0/3.0)*dxdI3*y;
	// tens4ds dsdCm1 = a*pow(J, -2.0/3.0)*x*dydCm1;

	// tens4ds tangent = 2.0*(0.5*dyad1s(dSdI1,I)+0.5*dyad1s(dSdI3,pow(J, 2.0)*Cm1)+0.5*ddots(dsdCm1,-I4th));
	//tens4ds tangent = 2.0*(dyad1(dSdI1,I)+dyad1(dSdI3,pow(J, 2.0)*Cm1)+ddot(tens4d(dsdCm1),tens4d(-I4th))).supersymm();

	// feLog("%f\t%f\t%f\n",y.xx(),y.xy(),y.xz());


	// feLog("y\n");
	// feLog("%f\t%f\t%f\n",y.xx(),y.xy(),y.xz());
	// feLog("%f\t%f\t%f\n",y.xy(),y.yy(),y.yz());
	// feLog("%f\t%f\t%f\n",y.xz(),y.yz(),y.zz());

	// feLog("dydI1\n");
	// feLog("%f\t%f\t%f\n",dydI1.xx(),dydI1.xy(),dydI1.xz());
	// feLog("%f\t%f\t%f\n",dydI1.xy(),dydI1.yy(),dydI1.yz());
	// feLog("%f\t%f\t%f\n",dydI1.xz(),dydI1.yz(),dydI1.zz());

	// feLog("dSdI1\n");
	// feLog("%f\t%f\t%f\n",dSdI1.xx(),dSdI1.xy(),dSdI1.xz());
	// feLog("%f\t%f\t%f\n",dSdI1.xy(),dSdI1.yy(),dSdI1.yz());
	// feLog("%f\t%f\t%f\n",dSdI1.xz(),dSdI1.yz(),dSdI1.zz());

	// feLog("dSdI3\n");
	// feLog("%f\t%f\t%f\n",dSdI3.xx(),dSdI3.xy(),dSdI3.xz());
	// feLog("%f\t%f\t%f\n",dSdI3.xy(),dSdI3.yy(),dSdI3.yz());
	// feLog("%f\t%f\t%f\n",dSdI3.xz(),dSdI3.yz(),dSdI3.zz());


	// feLog("dyad1s(dSdI1,I)\n");
	// feLog("%f\n",0.5*dyad1s(dSdI1,I)(0,0,0,0));


	// feLog("dyad1s(dSdI3,pow(J, 2.0)*Cm1)\n");
	// feLog("%f\n",0.5*dyad1s(dSdI3,pow(J, 2.0)*Cm1)(0,0,0,0));

	// feLog("dxdI1: %f\n",dxdI1);
	// feLog("dxdI3: %f\n",dxdI3);


	//return pt.push_forward(tangent);


	double cte = 2.0*a*pow(J, -2.0/3.0)*exp(0.5*b*(pow(J, -2.0/3.0)*I1-3.0));
	tens4ds tangent = cte*((pow(J, -2.0/3.0)*b/2.0*(BxB-1.0/3.0*I1*IxB)-1.0/3.0*IxB)-(1.0/3.0+0.5*b/3.0*I1*pow(J, -2.0/3.0))*(BxI-1.0/3.0*I1*IxI)+1.0/3.0*I1*I4);
	return tangent/J;

	// double cte = 2.0*a*pow(J, -2.0/3.0)*exp(0.5*b*(pow(J, -2.0/3.0)*I1-3.0));
	// tens4ds tangent = cte*((pow(J, -2.0/3.0)*b/2.0*(IxI-1.0/3.0*I1*Cm1xI)-1.0/3.0*Cm1xI)-(1.0/3.0+0.5*b/3.0*I1*pow(J, -2.0/3.0))*(IxCm1-1.0/3.0*I1*Cm1xCm1)+1.0/3.0*I1*I4th);
	// return pt.push_forward(tangent);
}
