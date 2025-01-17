/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* FLUID EVOLUTION :: BSSN evolution Class of COSMOS                                                     */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <algorithm>

//pressure with p=w\rho 
double Fmv0::pres(double rho)
{
	double w=fluidw*rho;
	
	return(w);
}

//dp/d\rho with p=w\rho
double Fmv0::dpres(double rho)
{
	double w=fluidw;
	
	return(w);
}

//for primitive variables from dynamical variables
//getting \rho and \Gamma(Lorentz factor) form E and p^2=p^\mu p_\mu with p=w\rho
void Fmv0::get_rhoGam(double Ene, double p2,double& rho,double& Gam)
{
	rho=0.5*((-1.+fluidw)*Ene+sqrt(pow(Ene*(1.-fluidw),2)+4.*fluidw*(Ene*Ene-p2)))/fluidw;
	Gam=sqrt((Ene+fluidw*rho)/((1.+fluidw)*rho));
	
	return;
}

//primitive variables from dynamical variables
void Fmv0::dyntoprim()
{
	#pragma omp parallel for 
	for(int l=ll[tab];l<=lu[tab];l++)
	{
		for(int k=kl[tab];k<=ku[tab];k++)
		{
			for(int j=jl[tab];j<=ju[tab];j++)
			{
				//proper: with respect to the fluid four-velocity u^\mu
				//Eulerian: with respect to the surface normal n^\mu
				////// fluid variables ///////////////////////////////////////////
				double Ene;							// Tnn
				double D;							// barion rest mass density by Eulerian
				double px;							// momentum density p_x by Eulerian
				double py;							// momentum density p_y by Eulerian
				double pz;							// momentum density p_z by Eulerian
				double pux;							// p^x
				double puy;							// p^y
				double puz;							// p^z
				double p2;							// p^\mu p_\mu
				double rho;							// proper energy density
				double Gam;							// Lorentz factor=-u^\mu n_\mu
				double EpP;							// Ene+Pressure
				////// fluid variables ///////////////////////////////////////////
				
				////// metric variables //////////////////////////////////////////
				double det;									// det \tilde \gamma
				double sqgam;								// sqrt(\gamma) !!not \tilde \gamma
				double gxx,gyy,gzz,gxy,gxz,gyz;				// metric components
				double gixx,giyy,gizz,gixy,gixz,giyz;		// inverse metric components
				double wa,ewa4i;							// conformal factor and exp(-4\psi)
				////// metric variables //////////////////////////////////////////

				////// substitution of metric variables //////////////////////////
				gxx=  get_bv(l,k,j, 7)+get_flat_df2x(j);						//
				gyy=  get_bv(l,k,j, 8)+get_flat_df2y(k);
				gzz=  get_bv(l,k,j, 9)+get_flat_df2z(l);
				gxy=  get_bv(l,k,j,10);
				gxz=  get_bv(l,k,j,11);
				gyz=  get_bv(l,k,j,12);
				wa=   get_bv(l,k,j,13);
				
				ewa4i=exp(-4.*wa);

				det=get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l);
				gixx= (gyy*gzz-gyz*gyz)/det;
				giyy= (gxx*gzz-gxz*gxz)/det;
				gizz= (gxx*gyy-gxy*gxy)/det;
				gixy=-(gxy*gzz-gyz*gxz)/det;
				gixz= (gxy*gyz-gyy*gxz)/det;
				giyz=-(gxx*gyz-gxy*gxz)/det;
				
				sqgam=exp(6.*wa)*sqrt(det);										//
				////// substitution of metric variables //////////////////////////
				
				
				////// substitution of dynamical fluid variables /////////////////
				Ene=get_bv(l,k,j,24)/sqgam;										//
				px=get_bv(l,k,j,25)/sqgam;
				py=get_bv(l,k,j,26)/sqgam;
				pz=get_bv(l,k,j,27)/sqgam;
				D=get_bv(l,k,j,28)/sqgam;
				
				pux=(gixx*px+gixy*py+gixz*pz)*ewa4i;
				puy=(gixy*px+giyy*py+giyz*pz)*ewa4i;
				puz=(gixz*px+giyz*py+gizz*pz)*ewa4i;
				
				p2=px*pux+py*puy+pz*puz;										//
				////// substitution of dynamical fluid variables /////////////////


				////// dynamical to primitive variables //////////////////////////
				get_rhoGam(Ene,p2,rho,Gam);										//

				// @@@@@@@@@@@@   primitive variables for fluid   @@@@@@@@@@@
				//  0:rho , 1:V^x  ,  2:V^y  ,  3:V^z   ,  4:epsilon
				
				set_primv(l,k,j,0)=rho;
				
				EpP=Ene+pres(rho);
				
				set_primv(l,k,j,1)=get_bv(l,k,j,0)*pux/EpP-get_bv(l,k,j,1);
				set_primv(l,k,j,2)=get_bv(l,k,j,0)*puy/EpP-get_bv(l,k,j,2);
				set_primv(l,k,j,3)=get_bv(l,k,j,0)*puz/EpP-get_bv(l,k,j,3);
				
				set_primv(l,k,j,4)=rho*Gam/D-1.;
				//set_primv(l,k,j,4)=rho/D-1.;									//
				////// dynamical to primitive variables //////////////////////////
				
			}
		}
	}
	//boundary_prim_quarter();
	
	return;
}

// for minmod function
double Fmv0::sign(double A){
    return (A>0)-(A<0);
}

// minmod function
double Fmv0::minmod(double a,double b)
{
	double w= sign(a)*max(0.,min(abs(a),sign(a)*b));
	
	return w;
}

