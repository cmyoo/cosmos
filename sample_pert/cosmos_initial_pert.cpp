/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* INITIAL DATA SETTING :: BSSN evolution Class of COSMOS                                                */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "../cosmos.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

using namespace std;
void Fmv::initial(double mu)
{
	//setting flatdet and Gamma
	set_flat();
	
	//double t0=2./(3.*Hb*(1.+fluidw));
	double kk=M_PI;
	//double w=fluidw;

	cout.precision(8);
	#pragma omp barrier

	for(int l=lli;l<=lui;l++)
	{
		double zz=get_z(l);
		double Z=funcf(zz);
		
		for(int k=kli;k<=kui;k++)
		{
			double yy=get_y(k);
			double Y=funcf(yy);
			
			for(int j=jli;j<=jui;j++)
			{
				double xx=get_x(j);
				double X=funcf(xx);
				
				// double dfx=df(xx);
				// double dfy=df(yy);
				// double dfz=df(zz);
				// double ddfx=ddf(xx);
				//double ddfy=ddf(yy);
				//double ddfz=ddf(zz);
				
				// double gxx=pow(dfx,2);
				// double gyy=pow(dfy,2);
				// double gzz=pow(dfz,2);

				double fxx=get_flat_df2x(j);
				double fyy=get_flat_df2y(k);
				double fzz=get_flat_df2z(l);

				double oafac=-mu;
				
				double g=(cos(X*kk)+cos(Y*kk)+cos(Z*kk));
				double Seed=oafac*g;
				
				double chi=0.;
   
				double xi=-0.181157064844222*Seed;
				
				double kap=0.178126451216426*Seed;

				double hpSk2=-0.0370205697380725*pow(kk,2);
				
				double hxx=hpSk2*(-oafac*cos(X*kk)+1./3.*Seed)*fxx;
				double hyy=hpSk2*(-oafac*cos(Y*kk)+1./3.*Seed)*fyy;
				double hzz=hpSk2*(-oafac*cos(Z*kk)+1./3.*Seed)*fzz;
				
				double ApSk2=0.07252724777633589*pow(kk,2);

				
				double Axx=ApSk2*(-oafac*cos(X*kk)+1./3.*Seed)*fxx;
				double Ayy=ApSk2*(-oafac*cos(Y*kk)+1./3.*Seed)*fyy;
				double Azz=ApSk2*(-oafac*cos(Z*kk)+1./3.*Seed)*fzz;
				
				double bx=0.;
				double by=0.;
				double bz=0.;
				
				//double det=gxx*gyy*gzz;
				
				double gxy=0.;
				double gxz=0.;
				double gyz=0.;
				
				double ek=-3.*Hb*(1.+kap);
				
				//double psii=grid_ipol(psi,jc,kc,lc,xx,yy,zz,2,deltax);
				double psii=1.+xi;
				//double psii=1.;

				double wa=log(psii);

				double alp=1.+chi;
				
				// W
				//set_bv(l,k,j,13)=wa;
				set_bv(l,k,j,13)=wa;
				// alpha
				set_bv(l,k,j,0)=alp;
				// beta^i
				set_bv(l,k,j,1)=bx;
				set_bv(l,k,j,2)=by;
				set_bv(l,k,j,3)=bz;
				// B^i
				set_bv(l,k,j,4)=0.;
				set_bv(l,k,j,5)=0.;
				set_bv(l,k,j,6)=0.;

				// g_{ij}
				set_bv(l,k,j, 7)=hxx;
				set_bv(l,k,j, 8)=hyy;
				set_bv(l,k,j, 9)=hzz;
				set_bv(l,k,j,10)=gxy;
				set_bv(l,k,j,11)=gxz;
				set_bv(l,k,j,12)=gyz;

				// A_{ij}
				set_bv(l,k,j,14)=Axx;
				set_bv(l,k,j,15)=Ayy;
				set_bv(l,k,j,16)=Azz;
				set_bv(l,k,j,17)=0.;
				set_bv(l,k,j,18)=0.;
				set_bv(l,k,j,19)=0.;

				// tr K
				set_bv(l,k,j,20)=ek;
				
				// Gamma^i
				//set_bv(l,k,j,21)=0.;
				//set_bv(l,k,j,22)=0.;
				//set_bv(l,k,j,23)=0.;

				set_bv(l,k,j,24)=0.;
			}
		}
	}
	
	#pragma omp barrier
	boundary_reflection();
	#pragma omp barrier
	enforce_const();							//det gamma=1 and tr A=0
	#pragma omp barrier
	boundary_reflection();
	#pragma omp barrier
	set_Gam();
	#pragma omp barrier
	boundary_reflection();
	#pragma omp barrier
	set_enemomini();
	#pragma omp barrier
	boundary_reflection_fluid();
	#pragma omp barrier
	dyntoprim();

	cout << "initial end" << endl;
		
	return;
}


