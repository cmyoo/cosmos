/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* INITIAL DATA SETTING :: BSSN evolution Class of COSMOS                                                */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

using namespace std;

//functions for superhorizon initial data with fluid
//void set_Psi_nonsph(double mu,double k,double xi2,double xi3);
//void initial_nonsph(double mu,double k,double xi2,double xi3);
	
//function to set fluid from geometrical variables by using constraint equations
//void set_enemomini()

//function to continue the evolution starting from the stored data
//void initial_continue(char *file_continue);

//functions for inhomogeneous grids
//double funcf(double X);
//double ifuncf(double X);
//double df(double X);
//double ddf(double X);
//double dddf(double X);
//void set_flat();

//function to set initial Gamma^i for evolution (Gamma^i will be treated as dynamical variables)
//void set_Gam();

//function to set initial scalar field profile
//void set_initial_scalar(double mu,double k,double xi2,double xi3);

//function to set initial fluid profile
//void set_initial_fluid(double mu,double k,double xi2,double xi3);

void Fmv::set_Psi_nonsph(double mu,double kk,double xi2,double xi3,double xit2,double xit3,double w)
{
	// parameters for the profile: kk, xi2, xi3, xit2, xit3, w(not used)
	
	double L=1.;						// box size
	double inr=0.6;						// Psi=1 for r>inr
	
	
	#pragma omp parallel for
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
				
				double r=sqrt(X*X+Y*Y+Z*Z);
				double a=0.25;

				if(r>L)
				set_psi(l,k,j)=1.;
				else if(r<inr)
				set_psi(l,k,j)=pow(M_E,(mu*pow(M_E,-(pow(a,3)*pow(kk,6)*pow(pow(X,2) + pow(Y,2) + pow(Z,2),3))/720.)*(-840*
         pow(840 + 10*pow(kk,2)*((3*xit2 + xit3)*pow(X,2) - 3*xit2*pow(Y,2) + xit3*(pow(Y,2) - 2*pow(Z,2)))*(pow(X,2) + pow(Y,2) + pow(Z,2)) + 7*pow(kk,4)*pow(pow(X,2) + pow(Y,2) + pow(Z,2),2),-1) + 
        32*pow(16 + (4*(6*X*xi2*Y + xi3*pow(X,2) + xi3*(pow(Y,2) - 2*pow(Z,2)) + pow(kk,2)*(pow(X,2) + pow(Y,2) + pow(Z,2))))/3. + 
           pow(6*X*xi2*Y + xi3*pow(X,2) + xi3*(pow(Y,2) - 2*pow(Z,2)) + pow(kk,2)*(pow(X,2) + pow(Y,2) + pow(Z,2)),2)/9.,-1)))/2.);
				else 
				set_psi(l,k,j)=pow(M_E,(mu*pow(M_E,-(pow(a,3)*pow(kk,6)*pow(pow(X,2) + pow(Y,2) + pow(Z,2),3))/720.)*(1 - pow(inr - L,-36)*pow(pow(inr - L,6) - pow(L - r,6),6))*
      (-840*pow(840 + 10*pow(kk,2)*((3*xit2 + xit3)*pow(X,2) - 3*xit2*pow(Y,2) + xit3*(pow(Y,2) - 2*pow(Z,2)))*(pow(X,2) + pow(Y,2) + pow(Z,2)) + 7*pow(kk,4)*pow(pow(X,2) + pow(Y,2) + pow(Z,2),2),-1) + 
        32*pow(16 + (4*(6*X*xi2*Y + xi3*pow(X,2) + xi3*(pow(Y,2) - 2*pow(Z,2)) + pow(kk,2)*(pow(X,2) + pow(Y,2) + pow(Z,2))))/3. + 
           pow(6*X*xi2*Y + xi3*pow(X,2) + xi3*(pow(Y,2) - 2*pow(Z,2)) + pow(kk,2)*(pow(X,2) + pow(Y,2) + pow(Z,2)),2)/9.,-1)))/2.);

			}
		}
	}
	#pragma omp barrier
	
	boundary_psi_initial();
	
	return;
}

void Fmv::set_Psi_nonsph(double mu,double kk,double xi2,double xi3)
{
	// parameters for the profile: kk, xi2, xi3

	double L=1.;						// box size
	double inr=0.8;						// Psi=1 for r>inr
	
	#pragma omp parallel for
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
				
				double r=sqrt(X*X+Y*Y+Z*Z);
				
				//see arXiv:2004.01042 for the functional form
				if(r>L)
				set_psi(l,k,j)=1.;
				else if(r<inr)
				set_psi(l,k,j)=pow(M_E,(pow(M_E,(-(pow(X,2)*(pow(kk,2) + 3*xi2 + xi3)) - 
				          (pow(kk,2) - 2*xi3)*pow(Y,2) - (pow(kk,2) - 3*xi2 + xi3)*pow(Z,2))/6.)*mu)/2.    );
				else 
				set_psi(l,k,j)=pow(M_E,(pow(M_E,(-(pow(X,2)*(pow(kk,2) + 3*xi2 + xi3)) - 
				          (pow(kk,2) - 2*xi3)*pow(Y,2) - (pow(kk,2) - 3*xi2 + xi3)*pow(Z,2))/6.)*mu*
					      (1 - pow(pow(inr - L,6) - pow(L - r,6),6)/pow(inr - L,36)))/2.);

			}
		}
	}
	#pragma omp barrier
	
	boundary_psi_initial();
	
	return;
}

void Fmv::initial_nonsph(double mu,double k,double xi2,double xi3,double xit2,double xit3,double w)
{
	//setting flatdet and Gamma
	set_flat();
	
	#pragma omp barrier
	set_Psi_nonsph(mu,k,xi2,xi3,xit2,xit3,w);

	#pragma omp barrier
	set_ini_from_Psi();
}

void Fmv::initial_nonsph(double mu,double k,double xi2,double xi3)
{
	//setting flatdet and Gamma
	set_flat();
	
	#pragma omp barrier
	set_Psi_nonsph(mu,k,xi2,xi3);

	#pragma omp barrier
	set_ini_from_Psi();
}

void Fmv::set_ini_from_Psi()
{
	double fGam=1.+fluidw;
	double Hbi2;
	double Hbi3;

	cout.precision(8);
	#pragma omp barrier
	
	double fac=2./(3.*fGam);

	double tt=fac/Hb;
	double HH=fac/tt;
	double aa=pow(tt/tini,fac);
	set_t(tt);

	// cout << "tt=" << tt 
	// << " dt0=" << dt0 
	// << " HH=" << HH 
	// << " aa=" << aa 
	// << " tini=" << tini 
	// << endl;

	// exit(0);

	Hbi2=pow(HH*aa,-2);
	Hbi3=pow(HH*aa,-3);
		
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		double zz=get_z(l);
		
		for(int k=kli;k<=kui;k++)
		{
			double yy=get_y(k);
			
			for(int j=jli;j<=jui;j++)
			{
				double xx=get_x(j);
				
				double fxx=get_flat_df2x(j);
				double fyy=get_flat_df2y(k);
				double fzz=get_flat_df2z(l);

				double dfx=df(xx);
				double dfy=df(yy);
				double dfz=df(zz);
				double ddfx=ddf(xx);
				double ddfy=ddf(yy);
				double ddfz=ddf(zz);
				
				double Psiv=get_psi(l,k,j);
				
				double Psi_x=(   -(get_psi(l,k,j+2)-get_psi(l,k,j-2))
						+8.*(get_psi(l,k,j+1)-get_psi(l,k,j-1))
						)*dxi12;
				double Psi_y=(   -(get_psi(l,k+2,j)-get_psi(l,k-2,j))
						+8.*(get_psi(l,k+1,j)-get_psi(l,k-1,j))
						)*dyi12;
				double Psi_z=(   -(get_psi(l+2,k,j)-get_psi(l-2,k,j))
						+8.*(get_psi(l+1,k,j)-get_psi(l-1,k,j))
						)*dzi12;
				
				double Psi_xx=(    -(     get_psi(l  ,k  ,j+2) + get_psi(l  ,k  ,j-2))
						+16.*( get_psi(l  ,k  ,j+1) + get_psi(l  ,k  ,j-1))
						-30.*  get_psi(l  ,k  ,j  )
						)*dxi*dxi12;
				double Psi_yy=(    -(     get_psi(l  ,k+2,j  ) + get_psi(l  ,k-2,j  ))
						+16.*( get_psi(l  ,k+1,j  ) + get_psi(l  ,k-1,j  ))
						-30.*  get_psi(l  ,k  ,j  )
						)*dyi*dyi12;
				double Psi_zz=(    -(     get_psi(l+2,k  ,j  ) + get_psi(l-2,k  ,j  ))
						+16.*( get_psi(l+1,k  ,j  ) + get_psi(l-1,k  ,j  ))
						-30.*  get_psi(l  ,k  ,j  )
						)*dzi*dzi12;
				
				double DDPsi_xy=( -( ( -(   get_psi(l  ,k+2,j+2) - get_psi(l  ,k+2,j-2))
						+8.*(get_psi(l  ,k+2,j+1) - get_psi(l  ,k+2,j-1)))
						-(   -(get_psi(l  ,k-2,j+2) - get_psi(l  ,k-2,j-2))
						+8.*(get_psi(l  ,k-2,j+1) - get_psi(l  ,k-2,j-1))))
						+8.*(  (-(get_psi(l  ,k+1,j+2) - get_psi(l  ,k+1,j-2))
						+8.*(get_psi(l  ,k+1,j+1) - get_psi(l  ,k+1,j-1)))
						-(-(get_psi(l  ,k-1,j+2) - get_psi(l  ,k-1,j-2))
						+8.*(get_psi(l  ,k-1,j+1) - get_psi(l  ,k-1,j-1))))
						)*dxi12*dyi12;
				double DDPsi_xz=( -( ( -(   get_psi(l+2,k  ,j+2) - get_psi(l+2,k  ,j-2))
						+8.*(get_psi(l+2,k  ,j+1) - get_psi(l+2,k  ,j-1)))
						-(   -(get_psi(l-2,k  ,j+2) - get_psi(l-2,k  ,j-2))
						+8.*(get_psi(l-2,k  ,j+1) - get_psi(l-2,k  ,j-1))))
						+8.*(  (-(get_psi(l+1,k  ,j+2) - get_psi(l+1,k  ,j-2))
						+8.*(get_psi(l+1,k  ,j+1) - get_psi(l+1,k  ,j-1)))
						-(-(get_psi(l-1,k  ,j+2) - get_psi(l-1,k  ,j-2))
						+8.*(get_psi(l-1,k  ,j+1) - get_psi(l-1,k  ,j-1))))
						)*dxi12*dzi12;
				double DDPsi_yz=( -( ( -(   get_psi(l+2,k+2,j  ) - get_psi(l+2,k-2,j  ))
						+8.*(get_psi(l+2,k+1,j  ) - get_psi(l+2,k-1,j  )))
						-(   -(get_psi(l-2,k+2,j  ) - get_psi(l-2,k-2,j  ))
						+8.*(get_psi(l-2,k+1,j  ) - get_psi(l-2,k-1,j  ))))
						+8.*(  (-(get_psi(l+1,k+2,j  ) - get_psi(l+1,k-2,j  ))
						+8.*(get_psi(l+1,k+1,j  ) - get_psi(l+1,k-1,j  )))
						-(-(get_psi(l-1,k+2,j  ) - get_psi(l-1,k-2,j  ))
						+8.*(get_psi(l-1,k+1,j  ) - get_psi(l-1,k-1,j  ))))
						)*dyi12*dzi12;
						
				double DDPsi_xx=Psi_xx-get_flat_Gamx(j)*Psi_x;
				double DDPsi_yy=Psi_yy-get_flat_Gamy(k)*Psi_y;
				double DDPsi_zz=Psi_zz-get_flat_Gamz(l)*Psi_z;

				double LapPsi=(Psi_xx-Psi_x*ddfx/dfx)*pow(dfx,-2)
								+(Psi_yy-Psi_y*ddfy/dfy)*pow(dfy,-2)
								+(Psi_zz-Psi_z*ddfz/dfz)*pow(dfz,-2);
				
				double DPDP=pow(Psi_x/dfx,2)+pow(Psi_y/dfy,2)+pow(Psi_z/dfz,2);
				
								
				double pxx=(-2.*Psiv*(DDPsi_xx-fxx*LapPsi/3.)+6.*(Psi_x*Psi_x-fxx*DPDP/3.))*pow(Psiv,-6);
				double pyy=(-2.*Psiv*(DDPsi_yy-fyy*LapPsi/3.)+6.*(Psi_y*Psi_y-fyy*DPDP/3.))*pow(Psiv,-6);
				double pzz=(-2.*Psiv*(DDPsi_zz-fzz*LapPsi/3.)+6.*(Psi_z*Psi_z-fzz*DPDP/3.))*pow(Psiv,-6);
				
				double pxy=(-2.*Psiv*DDPsi_xy+6.*Psi_x*Psi_y)*pow(Psiv,-6);
				double pxz=(-2.*Psiv*DDPsi_xz+6.*Psi_x*Psi_z)*pow(Psiv,-6);
				double pyz=(-2.*Psiv*DDPsi_yz+6.*Psi_y*Psi_z)*pow(Psiv,-6);
				
				double hxx=-4./(9.*pow(fGam,2)-4.)*pxx*Hbi2;
				double hyy=-4./(9.*pow(fGam,2)-4.)*pyy*Hbi2;
				double hzz=-4./(9.*pow(fGam,2)-4.)*pzz*Hbi2;
				double hxy=-4./(9.*pow(fGam,2)-4.)*pxy*Hbi2;
				double hxz=-4./(9.*pow(fGam,2)-4.)*pxz*Hbi2;
				double hyz=-4./(9.*pow(fGam,2)-4.)*pyz*Hbi2;
				
				
				double Axx=2./(3.*fGam+2)*pxx*Hbi2*HH;
				double Ayy=2./(3.*fGam+2)*pyy*Hbi2*HH;
				double Azz=2./(3.*fGam+2)*pzz*Hbi2*HH;
				double Axy=2./(3.*fGam+2)*pxy*Hbi2*HH;
				double Axz=2./(3.*fGam+2)*pxz*Hbi2*HH;
				double Ayz=2./(3.*fGam+2)*pyz*Hbi2*HH;
							
				double f=-4./3.*LapPsi*pow(Psiv,-5);
				
				double rhob=3./(8.*M_PI)*pow(HH,2);
				
				//CMC
				double dddfx=dddf(xx);
				double dddfy=dddf(yy);
				double dddfz=dddf(zz);
				double Psi_xxx=0.5*(-get_psi(l,k,j-2)+2.*get_psi(l,k,j-1)-2.*get_psi(l,k,j+1)+get_psi(l,k,j+2))*pow(dxi,3);
				double Psi_yyy=0.5*(-get_psi(l,k-2,j)+2.*get_psi(l,k-1,j)-2.*get_psi(l,k+1,j)+get_psi(l,k+2,j))*pow(dyi,3);
				double Psi_zzz=0.5*(-get_psi(l-2,k,j)+2.*get_psi(l-1,k,j)-2.*get_psi(l+1,k,j)+get_psi(l+2,k,j))*pow(dzi,3);
				
				double Psi_xxy=(    -(     get_psi(l  ,k-2  ,j+2) -8.*get_psi(l  ,k-1  ,j+2)+8.*get_psi(l  ,k+1  ,j+2)-get_psi(l  ,k+2  ,j+2)
						+ get_psi(l  ,k-2  ,j-2)-8.*get_psi(l  ,k-1  ,j-2)+8.*get_psi(l  ,k+1  ,j-2)-get_psi(l  ,k+2  ,j-2)		)
						+16.*( get_psi(l  ,k-2  ,j+1) -8.*get_psi(l  ,k-1  ,j+1)+8.*get_psi(l  ,k+1  ,j+1)-get_psi(l  ,k+2  ,j+1)
						+ get_psi(l  ,k-2  ,j-1)-8.*get_psi(l  ,k-1  ,j-1)+8.*get_psi(l  ,k+1  ,j-1)-get_psi(l  ,k+2  ,j-1)		)
						-30.*  (get_psi(l  ,k-2  ,j  )-8.*get_psi(l  ,k-1  ,j  )+8.*get_psi(l  ,k+1  ,j  )-get_psi(l  ,k+2  ,j  ))
						)*dxi*dxi12*dyi12;
				
				double Psi_xxz=(    -(     get_psi(l-2  ,k  ,j+2) -8.*get_psi(l-1  ,k  ,j+2)+8.*get_psi(l+1  ,k  ,j+2)-get_psi(l+2  ,k  ,j+2)
						+ get_psi(l-2  ,k  ,j-2)-8.*get_psi(l-1  ,k  ,j-2)+8.*get_psi(l+1  ,k  ,j-2)-get_psi(l+2  ,k  ,j-2)		)
						+16.*( get_psi(l-2  ,k  ,j+1) -8.*get_psi(l-1  ,k  ,j+1)+8.*get_psi(l+1  ,k  ,j+1)-get_psi(l+2  ,k  ,j+1)
						+ get_psi(l-2  ,k  ,j-1)-8.*get_psi(l-1  ,k  ,j-1)+8.*get_psi(l+1  ,k  ,j-1)-get_psi(l+2  ,k  ,j-1)		)
						-30.*  (get_psi(l-2  ,k  ,j  )-8.*get_psi(l-1  ,k  ,j  )+8.*get_psi(l+1  ,k  ,j  )-get_psi(l+2  ,k  ,j  ))
						)*dxi*dxi12*dzi12;
						
				double Psi_yyx=(    -(     get_psi(l  ,k+2  ,j-2) -8.*get_psi(l  ,k+2  ,j-1)+8.*get_psi(l  ,k+2  ,j+1)-get_psi(l  ,k+2  ,j+2)
						+ get_psi(l  ,k-2  ,j-2)-8.*get_psi(l  ,k-2  ,j-1)+8.*get_psi(l  ,k-2  ,j+1)-get_psi(l  ,k-2  ,j+2)		)
						+16.*( get_psi(l  ,k+1  ,j-2) -8.*get_psi(l  ,k+1  ,j-1)+8.*get_psi(l  ,k+1  ,j+1)-get_psi(l  ,k+1  ,j+2)
						+ get_psi(l  ,k-1  ,j-2)-8.*get_psi(l  ,k-1  ,j-1)+8.*get_psi(l  ,k-1  ,j+1)-get_psi(l  ,k-1  ,j+2)		)
						-30.*  (get_psi(l  ,k  ,j-2  )-8.*get_psi(l  ,k  ,j-1  )+8.*get_psi(l  ,k  ,j+1  )-get_psi(l  ,k  ,j+2  ))
						)*dyi*dyi12*dxi12;
						
				double Psi_yyz=(    -(     get_psi(l-2  ,k+2  ,j) -8.*get_psi(l-1  ,k+2  ,j)+8.*get_psi(l+1  ,k+2  ,j)-get_psi(l+2  ,k+2  ,j)
						+ get_psi(l-2  ,k-2  ,j)-8.*get_psi(l-1  ,k-2  ,j)+8.*get_psi(l+1  ,k-2  ,j)-get_psi(l+2  ,k-2  ,j)		)
						+16.*( get_psi(l-2  ,k+1  ,j) -8.*get_psi(l-1  ,k+1  ,j)+8.*get_psi(l+1  ,k+1  ,j)-get_psi(l+2  ,k+1  ,j)
						+ get_psi(l-2  ,k-1  ,j)-8.*get_psi(l-1  ,k-1  ,j)+8.*get_psi(l+1  ,k-1  ,j)-get_psi(l+2  ,k-1  ,j)		)
						-30.*  (get_psi(l-2  ,k  ,j  )-8.*get_psi(l-1  ,k  ,j  )+8.*get_psi(l+1  ,k  ,j  )-get_psi(l+2  ,k  ,j  ))
						)*dyi*dyi12*dzi12;
						
				double Psi_zzx=(    -(     get_psi(l+2  ,k  ,j-2) -8.*get_psi(l+2  ,k  ,j-1)+8.*get_psi(l+2  ,k  ,j+1)-get_psi(l+2  ,k  ,j+2)
						+ get_psi(l-2  ,k  ,j-2)-8.*get_psi(l-2  ,k  ,j-1)+8.*get_psi(l-2  ,k  ,j+1)-get_psi(l-2  ,k  ,j+2)		)
						+16.*( get_psi(l+1  ,k  ,j-2) -8.*get_psi(l+1  ,k  ,j-1)+8.*get_psi(l+1  ,k  ,j+1)-get_psi(l+1  ,k  ,j+2)
						+ get_psi(l-1  ,k  ,j-2)-8.*get_psi(l-1  ,k  ,j-1)+8.*get_psi(l-1  ,k  ,j+1)-get_psi(l-1  ,k  ,j+2)		)
						-30.*  (get_psi(l  ,k  ,j-2  )-8.*get_psi(l  ,k  ,j-1  )+8.*get_psi(l  ,k  ,j+1  )-get_psi(l  ,k  ,j+2  ))
						)*dzi*dzi12*dxi12;
						
				double Psi_zzy=(    -(     get_psi(l+2  ,k-2  ,j) -8.*get_psi(l+2  ,k -1 ,j)+8.*get_psi(l+2  ,k+1  ,j)-get_psi(l+2  ,k+2  ,j)
						+ get_psi(l-2  ,k-2  ,j)-8.*get_psi(l-2  ,k-1  ,j)+8.*get_psi(l-2  ,k+1  ,j)-get_psi(l-2  ,k+2  ,j)		)
						+16.*( get_psi(l+1  ,k-2  ,j) -8.*get_psi(l+1  ,k-1  ,j)+8.*get_psi(l+1  ,k+1  ,j)-get_psi(l+1  ,k+2  ,j)
						+ get_psi(l-1  ,k-2  ,j)-8.*get_psi(l-1  ,k-1  ,j)+8.*get_psi(l-1  ,k+1  ,j)-get_psi(l-1  ,k+2  ,j)		)
						-30.*  (get_psi(l  ,k-2  ,j  )-8.*get_psi(l  ,k-1  ,j  )+8.*get_psi(l  ,k+1  ,j  )-get_psi(l  ,k+2  ,j  ))
						)*dzi*dzi12*dyi12;
				double DLapPsi_x=-3.*ddfx*pow(dfx,-3)*Psi_xx+pow(dfx,-2)*Psi_xxx+(3.*pow(ddfx,2)*pow(dfx,-4)-dddfx*pow(dfx,-3))*Psi_x
								+pow(dfy,-2)*Psi_yyx-ddfy*pow(dfy,-3)*DDPsi_xy
								+pow(dfz,-2)*Psi_zzx-ddfz*pow(dfz,-3)*DDPsi_xz;
								
				double DLapPsi_y=-3.*ddfy*pow(dfy,-3)*Psi_yy+pow(dfy,-2)*Psi_yyy+(3.*pow(ddfy,2)*pow(dfy,-4)-dddfy*pow(dfy,-3))*Psi_y
								+pow(dfx,-2)*Psi_xxy-ddfx*pow(dfx,-3)*DDPsi_xy
								+pow(dfz,-2)*Psi_zzy-ddfz*pow(dfz,-3)*DDPsi_yz;
				
				double DLapPsi_z=-3.*ddfz*pow(dfz,-3)*Psi_zz+pow(dfz,-2)*Psi_zzz+(3.*pow(ddfz,2)*pow(dfz,-4)-dddfz*pow(dfz,-3))*Psi_z
								+pow(dfy,-2)*Psi_yyz-ddfy*pow(dfy,-3)*DDPsi_yz
								+pow(dfx,-2)*Psi_xxz-ddfx*pow(dfx,-3)*DDPsi_xz;
				double f_x=-4./3.*(DLapPsi_x*pow(Psiv,-5)-5.*LapPsi*Psi_x*pow(Psiv,-6));
				double f_y=-4./3.*(DLapPsi_y*pow(Psiv,-5)-5.*LapPsi*Psi_y*pow(Psiv,-6));
				double f_z=-4./3.*(DLapPsi_z*pow(Psiv,-5)-5.*LapPsi*Psi_z*pow(Psiv,-6));
				double kap=0.;
				double chi=-(3.*fGam-2.)/(3.*fGam)*f*Hbi2;
				double xi=-1./(6.*fGam)*f*Hbi2;
				double delta=f*Hbi2;
				double u_x=2./(3.*fGam*(3.*fGam+2.))*f_x*Hbi3*aa;
				double u_y=2./(3.*fGam*(3.*fGam+2.))*f_y*Hbi3*aa;
				double u_z=2./(3.*fGam*(3.*fGam+2.))*f_z*Hbi3*aa;
				
				//comoving
				//double kap=-1./(3.*fGam+2)*f*Hbi2;
				//double chi=-3.*(fGam-1.)/(3.*fGam+2.)*f*Hbi2;
				//double xi=-1./(2.*(3.*fGam+2))*f*Hbi2;
				//double delta=3.*fGam/(3.*fGam+2)*f*Hbi2;
				//double u_x=0.;
				//double u_y=0.;
				//double u_z=0.;
				
				double psii=Psiv*(1.+xi);

				double ek=-3.*HH*(1.+kap);
				double wa=log(psii*sqrt(aa));
				double alp=1.+chi;
				double rho=rhob*(1.+delta);
				
				// w
				set_bv(l,k,j,13)=wa;
				// alpha
				set_bv(l,k,j,0)=alp;
				// beta^i
				set_bv(l,k,j,1)=0.;
				set_bv(l,k,j,2)=0.;
				set_bv(l,k,j,3)=0.;
				// B^i
				set_bv(l,k,j,4)=0.;
				set_bv(l,k,j,5)=0.;
				set_bv(l,k,j,6)=0.;

				// g_{ij}
				set_bv(l,k,j, 7)=hxx;
				set_bv(l,k,j, 8)=hyy;
				set_bv(l,k,j, 9)=hzz;
				set_bv(l,k,j,10)=hxy;
				set_bv(l,k,j,11)=hxz;
				set_bv(l,k,j,12)=hyz;

				// A_{ij}
				set_bv(l,k,j,14)=Axx;
				set_bv(l,k,j,15)=Ayy;
				set_bv(l,k,j,16)=Azz;
				set_bv(l,k,j,17)=Axy;
				set_bv(l,k,j,18)=Axz;
				set_bv(l,k,j,19)=Ayz;

				// tr K
				set_bv(l,k,j,20)=ek;
									
				double sqgam=dfx*dfy*dfz*exp(6.*wa);
				
				set_bv(l,k,j,24)=sqgam*rho;
				
				double rhop=rho+pres(rho);
				
				set_bv(l,k,j,25)=sqgam*rhop*u_x;
				set_bv(l,k,j,26)=sqgam*rhop*u_y;
				set_bv(l,k,j,27)=sqgam*rhop*u_z;
				
				double eps=pow(rho/rhob,fluidw/(1.+fluidw))-1.;
				
				set_bv(l,k,j,28) = rho*sqgam/(1.+eps);
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


//setting initial Gamma^i for evolution (Gamma^i will be treated as dynamical variables)
void Fmv0::set_Gam()
{
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				double 
				gxx_p,gyy_p,gzz_p,
				gxy_p,gxz_p,gyz_p,
				gyx_p,gzx_p,gzy_p;
				
				
				double dxx_x_p,dyy_x_p,dzz_x_p,
				dxy_x_p,dxz_x_p,dyz_x_p,
				dyx_x_p,dzx_x_p,dzy_x_p;
				double dxx_y_p,dyy_y_p,dzz_y_p,
				dxy_y_p,dxz_y_p,dyz_y_p,
				dyx_y_p,dzx_y_p,dzy_y_p;
				double dxx_z_p,dyy_z_p,dzz_z_p,
				dxy_z_p,dxz_z_p,dyz_z_p,
				dyx_z_p,dzx_z_p,dzy_z_p;

				double det_p,det_pi;
				double gixx_p,giyy_p,gizz_p,
				gixy_p,gixz_p,giyz_p,
				giyx_p,gizx_p,gizy_p;

				double crdx_xx,crdx_yy,crdx_zz,
				crdx_xy,crdx_xz,crdx_yz;
				double crdy_xx,crdy_yy,crdy_zz,
				crdy_xy,crdy_xz,crdy_yz;
				double crdz_xx,crdz_yy,crdz_zz,
				crdz_xy,crdz_xz,crdz_yz;

				double crx_xx,crx_yy,crx_zz,
				crx_xy,crx_xz,crx_yz;
				double cry_xx,cry_yy,cry_zz,
				cry_xy,cry_xz,cry_yz;
				double crz_xx,crz_yy,crz_zz,
				crz_xy,crz_xz,crz_yz;

				
				double gamma0_x,gamma0_y,gamma0_z;
				
				double fxx,fyy,fzz,dfxx,dfyy,dfzz;
				double Gam_ux_xx,Gam_uy_yy,Gam_uz_zz;
				double Del_ux_xx,Del_ux_yy,Del_ux_zz,Del_ux_xy,Del_ux_xz,Del_ux_yz,
						Del_uy_xx,Del_uy_yy,Del_uy_zz,Del_uy_xy,Del_uy_xz,Del_uy_yz,
						Del_uz_xx,Del_uz_yy,Del_uz_zz,Del_uz_xy,Del_uz_xz,Del_uz_yz;
				
				
				//fij for inhomogeneous grid 
				fxx=get_flat_df2x(j);
				fyy=get_flat_df2y(k);
				fzz=get_flat_df2z(l);

				gxx_p=  get_bv(l,k,j, 7)+fxx;
				gyy_p=  get_bv(l,k,j, 8)+fyy;
				gzz_p=  get_bv(l,k,j, 9)+fzz;
				gxy_p=  get_bv(l,k,j,10);
				gxz_p=  get_bv(l,k,j,11);
				gyz_p=  get_bv(l,k,j,12);
				gyx_p=	gxy_p;
				gzx_p=	gxz_p;
				gzy_p=	gyz_p;
				
				
				//bar Gamma_ui_jk for inhomogeneous grid 
				Gam_ux_xx=get_flat_Gamx(j);
				Gam_uy_yy=get_flat_Gamy(k);
				Gam_uz_zz=get_flat_Gamz(l);
				
				dfxx=2.*fxx*Gam_ux_xx;
				dfyy=2.*fyy*Gam_uy_yy;
				dfzz=2.*fzz*Gam_uz_zz;

				// \del_x \tilde{gamma}_{ij} *0.5 
				dxx_x_p=(get_f_x(l,k,j,7)+dfxx)*0.5;
				dyy_x_p=get_f_x(l,k,j,8)*0.5;
				dzz_x_p=get_f_x(l,k,j,9)*0.5;
				dxy_x_p=get_f_x(l,k,j,10)*0.5;
				dxz_x_p=get_f_x(l,k,j,11)*0.5;
				dyz_x_p=get_f_x(l,k,j,12)*0.5;
				dyx_x_p=dxy_x_p;
				dzx_x_p=dxz_x_p;
				dzy_x_p=dyz_x_p;

				// \del_y \tilde \gamma_{ij} * 0.5 
				dxx_y_p=get_f_y(l,k,j,7)*0.5;
				dyy_y_p=(get_f_y(l,k,j,8)+dfyy)*0.5;
				dzz_y_p=get_f_y(l,k,j,9)*0.5;
				dxy_y_p=get_f_y(l,k,j,10)*0.5;
				dxz_y_p=get_f_y(l,k,j,11)*0.5;
				dyz_y_p=get_f_y(l,k,j,12)*0.5;
				dyx_y_p=dxy_y_p;
				dzx_y_p=dxz_y_p;
				dzy_y_p=dyz_y_p;

				// \del_z \tilde \gamma_{ij} * 0.5
				dxx_z_p=get_f_z(l,k,j,7)*0.5;
				dyy_z_p=get_f_z(l,k,j,8)*0.5;
				dzz_z_p=(get_f_z(l,k,j,9)+dfzz)*0.5;
				dxy_z_p=get_f_z(l,k,j,10)*0.5;
				dxz_z_p=get_f_z(l,k,j,11)*0.5;
				dyz_z_p=get_f_z(l,k,j,12)*0.5;
				dyx_z_p=dxy_z_p;
				dzx_z_p=dxz_z_p;
				dzy_z_p=dyz_z_p;
				
				// tilde{gamma}^{ij} -> gi
			//	det_p=gxx_p*gyy_p*gzz_p+gxy_p*gyz_p*gzx_p+gxz_p*gyx_p*gzy_p
			//			-gxz_p*gyy_p*gzx_p-gxy_p*gyx_p*gzz_p-gxx_p*gyz_p*gzy_p;
			//	if(det_p<1.e-16) 
			//	 det_p=1.e-16;
				det_p=get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l);
				det_pi=1./det_p;
				gixx_p= (gyy_p*gzz_p-gyz_p*gzy_p)*det_pi;
				giyy_p= (gxx_p*gzz_p-gxz_p*gzx_p)*det_pi;
				gizz_p= (gxx_p*gyy_p-gxy_p*gyx_p)*det_pi;
				gixy_p=-(gyx_p*gzz_p-gyz_p*gzx_p)*det_pi;
				gixz_p= (gyx_p*gzy_p-gyy_p*gzx_p)*det_pi;
				giyz_p=-(gxx_p*gzy_p-gxy_p*gzx_p)*det_pi;
				
				giyx_p=gixy_p;
				gizx_p=gixz_p;
				gizy_p=giyz_p;
				


				//  \tilde{\Gamma}_{ijk}=\tilde{\gamma}_{im}\tilde{\Gamma}^m_{jk} -> crd{i}_{jk}
				// NOTE:d.._._p is 0.5*derivatives
				crdx_xx=dxx_x_p+dxx_x_p -dxx_x_p;
				crdx_yy=dxy_y_p+dxy_y_p -dyy_x_p;
				crdx_zz=dxz_z_p+dxz_z_p -dzz_x_p;
				crdx_xy=dxx_y_p+dxy_x_p -dxy_x_p;
				crdx_xz=dxx_z_p+dxz_x_p -dxz_x_p;
				crdx_yz=dxy_z_p+dxz_y_p -dyz_x_p;
				
				crdy_xx=dyx_x_p+dyx_x_p -dxx_y_p;
				crdy_yy=dyy_y_p+dyy_y_p -dyy_y_p;
				crdy_zz=dyz_z_p+dyz_z_p -dzz_y_p;
				crdy_xy=dyx_y_p+dyy_x_p -dxy_y_p;
				crdy_xz=dyx_z_p+dyz_x_p -dxz_y_p;
				crdy_yz=dyy_z_p+dyz_y_p -dyz_y_p;
				
				crdz_xx=dzx_x_p+dzx_x_p -dxx_z_p;
				crdz_yy=dzy_y_p+dzy_y_p -dyy_z_p;
				crdz_zz=dzz_z_p+dzz_z_p -dzz_z_p;
				crdz_xy=dzx_y_p+dzy_x_p -dxy_z_p;
				crdz_xz=dzx_z_p+dzz_x_p -dxz_z_p;
				crdz_yz=dzy_z_p+dzz_y_p -dyz_z_p;
				
				//  \tilde{\Gamma}^i_{jk} -> cr{i}_{jk}
				crx_xx=gixx_p*crdx_xx +gixy_p*crdy_xx +gixz_p*crdz_xx;
				crx_yy=gixx_p*crdx_yy +gixy_p*crdy_yy +gixz_p*crdz_yy;
				crx_zz=gixx_p*crdx_zz +gixy_p*crdy_zz +gixz_p*crdz_zz;
				crx_xy=gixx_p*crdx_xy +gixy_p*crdy_xy +gixz_p*crdz_xy;
				crx_xz=gixx_p*crdx_xz +gixy_p*crdy_xz +gixz_p*crdz_xz;
				crx_yz=gixx_p*crdx_yz +gixy_p*crdy_yz +gixz_p*crdz_yz;
				
				cry_xx=giyx_p*crdx_xx +giyy_p*crdy_xx +giyz_p*crdz_xx;
				cry_yy=giyx_p*crdx_yy +giyy_p*crdy_yy +giyz_p*crdz_yy;
				cry_zz=giyx_p*crdx_zz +giyy_p*crdy_zz +giyz_p*crdz_zz;
				cry_xy=giyx_p*crdx_xy +giyy_p*crdy_xy +giyz_p*crdz_xy;
				cry_xz=giyx_p*crdx_xz +giyy_p*crdy_xz +giyz_p*crdz_xz;
				cry_yz=giyx_p*crdx_yz +giyy_p*crdy_yz +giyz_p*crdz_yz;
				
				crz_xx=gizx_p*crdx_xx +gizy_p*crdy_xx +gizz_p*crdz_xx;
				crz_yy=gizx_p*crdx_yy +gizy_p*crdy_yy +gizz_p*crdz_yy;
				crz_zz=gizx_p*crdx_zz +gizy_p*crdy_zz +gizz_p*crdz_zz;
				crz_xy=gizx_p*crdx_xy +gizy_p*crdy_xy +gizz_p*crdz_xy;
				crz_xz=gizx_p*crdx_xz +gizy_p*crdy_xz +gizz_p*crdz_xz;
				crz_yz=gizx_p*crdx_yz +gizy_p*crdy_yz +gizz_p*crdz_yz;
				
				//Delta_ui_jk=tilde Gamma - bar Gamma
				
				//Del_ux_xx=0.5*(gixx_p*Dgam_x_xx+gixy_p*(2.*Dgam_x_yx-Dgam_y_xx)+gixz_p*(2.*Dgam_x_zx-Dgam_z_xx));
				Del_ux_xx=crx_xx-Gam_ux_xx;
				Del_ux_yy=crx_yy;
				Del_ux_zz=crx_zz;
				Del_ux_xy=crx_xy;
				Del_ux_xz=crx_xz;
				Del_ux_yz=crx_yz;
				
				Del_uy_xx=cry_xx;
				//Del_uy_yy=0.5*(giyy_p*Dgam_y_yy+giyz_p*(2.*Dgam_y_zy-Dgam_z_yy)+giyx_p*(2.*Dgam_y_xy-Dgam_x_yy));
				Del_uy_yy=cry_yy-Gam_uy_yy;
				Del_uy_zz=cry_zz;
				Del_uy_xy=cry_xy;
				Del_uy_xz=cry_xz;
				Del_uy_yz=cry_yz;
				
				Del_uz_xx=crz_xx;
				Del_uz_yy=crz_yy;
				//Del_uz_zz=0.5*(gizz_p*Dgam_z_zz+gizx_p*(2.*Dgam_z_xz-Dgam_x_zz)+gizy_p*(2.*Dgam_z_yz-Dgam_y_zz));
				Del_uz_zz=crz_zz-Gam_uz_zz;
				Del_uz_xy=crz_xy;
				Del_uz_xz=crz_xz;
				Del_uz_yz=crz_yz;
				
				//tilde Gamma0^i =- cal D_j tilde gamma ^ji=tilde gamma^jk Delta^i_jk
				gamma0_x=gixx_p*Del_ux_xx+giyy_p*Del_ux_yy+gizz_p*Del_ux_zz
						+2.*(gixy_p*Del_ux_xy+gixz_p*Del_ux_xz+giyz_p*Del_ux_yz);
				gamma0_y=gixx_p*Del_uy_xx+giyy_p*Del_uy_yy+gizz_p*Del_uy_zz
						+2.*(gixy_p*Del_uy_xy+gixz_p*Del_uy_xz+giyz_p*Del_uy_yz);
				gamma0_z=gixx_p*Del_uz_xx+giyy_p*Del_uz_yy+gizz_p*Del_uz_zz
						+2.*(gixy_p*Del_uz_xy+gixz_p*Del_uz_xz+giyz_p*Del_uz_yz);
				
				// Gamma^i
				set_bv(l,k,j,21)=gamma0_x;
				set_bv(l,k,j,22)=gamma0_y;
				set_bv(l,k,j,23)=gamma0_z;
				
			}
		}
	}
	
	return;
}

//function to continue the evolution starting from the stored data
void Fmv0::initial_continue(ifstream& fcontinue)
{
	set_zero_all();
	#pragma omp barrier
	
	//setting flatdet and Gamma
	set_flat();

	string buf;
	
	// ifstream fcontinue(file_continue);
	// if(fcontinue.fail()){
	// 	cout << "Initial data file to continue is not found." << endl;
	// 	abort();
	// }

	cout.precision(8);
	
	double t_continue;
	
	getline(fcontinue, buf);
	cout << buf << endl;
	
	sscanf(buf.data(),"##time=%lf",
	       &t_continue);

	sscanf(buf.data(),"##time=%lf dtp=%lf dtpp=%lf",
	       &t_continue,&dtp,&dtpp);
	
	set_t(t_continue);
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				double gxx,gyy,gzz,gxy,gxz,gyz,wa;
				double akxx,akyy,akzz,akxy,akxz,akyz,ek;
				double zgx,zgy,zgz;
				double bx,by,bz,bbx,bby,bbz;
				double alp;
				double Ene,px,py,pz,D;
				double phii,Pi;

				getline(fcontinue, buf);
				
				if(fluidevo && scalarevo)
				{
					sscanf(buf.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
					&alp,&bx,&by,&bz,&bbx,&bby,&bbz,&gxx,&gyy,&gzz,
					&gxy,&gxz,&gyz,&wa,&akxx,&akyy,&akzz,&akxy,&akxz,&akyz,
					&ek,&zgx,&zgy,&zgz,&Ene,&px,&py,&pz,&D,&phii,&Pi);
				}
				else if(fluidevo)
				{
					sscanf(buf.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
					&alp,&bx,&by,&bz,&bbx,&bby,&bbz,&gxx,&gyy,&gzz,
					&gxy,&gxz,&gyz,&wa,&akxx,&akyy,&akzz,&akxy,&akxz,&akyz,
					&ek,&zgx,&zgy,&zgz,&Ene,&px,&py,&pz,&D);
				}
				else
				{
					sscanf(buf.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
						&alp,&bx,&by,&bz,&bbx,&bby,&bbz,&gxx,&gyy,&gzz,
						&gxy,&gxz,&gyz,&wa,&akxx,&akyy,&akzz,&akxy,&akxz,&akyz,
						&ek,&zgx,&zgy,&zgz);
				}	

				// w
				set_bv(l,k,j,13)=wa;
				// alpha
				set_bv(l,k,j,0)=alp;
				// beta^i
				set_bv(l,k,j,1)=bx;
				set_bv(l,k,j,2)=by;
				set_bv(l,k,j,3)=bz;
				// B^i
				set_bv(l,k,j,4)=bbx;
				set_bv(l,k,j,5)=bby;
				set_bv(l,k,j,6)=bbz;

				// g_{ij}
				set_bv(l,k,j, 7)=gxx;
				set_bv(l,k,j, 8)=gyy;
				set_bv(l,k,j, 9)=gzz;
				set_bv(l,k,j,10)=gxy;
				set_bv(l,k,j,11)=gxz;
				set_bv(l,k,j,12)=gyz;

				// A_{ij}
				set_bv(l,k,j,14)=akxx;
				set_bv(l,k,j,15)=akyy;
				set_bv(l,k,j,16)=akzz;
				set_bv(l,k,j,17)=akxy;
				set_bv(l,k,j,18)=akxz;
				set_bv(l,k,j,19)=akyz;

				// tr K
				set_bv(l,k,j,20)=ek;
				
				// Gamma^i
				set_bv(l,k,j,21)=zgx;
				set_bv(l,k,j,22)=zgy;
				set_bv(l,k,j,23)=zgz;
				
				if(fluidevo)
				{
					set_bv(l,k,j,24)=Ene;
					set_bv(l,k,j,25)=px;
					set_bv(l,k,j,26)=py;
					set_bv(l,k,j,27)=pz;
					set_bv(l,k,j,28)=D;
				}
				
				if(scalarevo)
				{
					set_bv(l,k,j,nsc)=phii;
					set_bv(l,k,j,nscp)=Pi;
				}
			}
			getline(fcontinue, buf);
		}
		getline(fcontinue, buf);
	}
	
	return;
}

double Fmv0::funcf(double X)
{
	double w=X-amp/(1.+amp)*boxL/M_PI*sin(M_PI*X/boxL);
	
	return(w);
}
double Fmv0::df(double X)
{
	double w=1.-amp/(1.+amp)*cos(M_PI*X/boxL);
	
	return(w);
}
double Fmv0::ddf(double X)
{
	
	double w=amp/(1.+amp)*M_PI/boxL*sin(M_PI*X/boxL);
	
	return(w);
}
double Fmv0::dddf(double X)
{
	double w=amp/(1.+amp)*pow(M_PI/boxL,2)*cos(M_PI*X/boxL);
	
	return(w);
}

void Fmv0::set_flat()
{
	double xx,yy,zz;
	double Gamu,delGam,dff,ddff,dddff;
	
	#pragma omp parallel for private(zz,Gamu,delGam,dff,ddff,dddff)
	for(int l=ll[tab];l<=lu[tab];l++)
	{
		zz=get_z(l);
		dff=df(zz);
		
		ddff=ddf(zz);
		dddff=dddf(zz);

		Gamu=ddff/dff;
	
		delGam=dddff/dff-pow(ddff/dff,2);

		set_flat_df2z(l)=pow(dff,2);
		set_flat_Gamz(l)=Gamu;
		set_flat_dGamz(l)=delGam;
	}
	
	#pragma omp parallel for private(yy,Gamu,delGam,dff,ddff,dddff)
	for(int k=kl[tab];k<=ku[tab];k++)
	{
		yy=get_y(k);
		dff=df(yy);
		
		ddff=ddf(yy);
		dddff=dddf(yy);

		Gamu=ddff/dff;
	
		delGam=dddff/dff-pow(ddff/dff,2);

		set_flat_df2y(k)=pow(dff,2);
		set_flat_Gamy(k)=Gamu;
		set_flat_dGamy(k)=delGam;

	}
	
	#pragma omp parallel for private(xx,Gamu,delGam,dff,ddff,dddff)
	for(int j=jl[tab];j<=ju[tab];j++)
	{
		xx=get_x(j);
		dff=df(xx);
		
		ddff=ddf(xx);
		dddff=dddf(xx);

		Gamu=ddff/dff;
	
		delGam=dddff/dff-pow(ddff/dff,2);

		set_flat_df2x(j)=pow(dff,2);
		set_flat_Gamx(j)=Gamu;
		set_flat_dGamx(j)=delGam;

	}
	
	return;
}

//X-amp/(1.+amp)*xu/M_PI*sin(M_PI*X/xu);
double Fmv0::ifuncf(double x)
{
	double w=x;
	int n;
	
	for(n=0;n<1000;n++)
	{
		//double w1=x+amp/(1.+amp)*xu/M_PI*sin(M_PI*w/xu);
		
		double w1=x-funcf(w)+w;
		if(abs(w-w1)<1.0e-14)
		{
			w=w1;
			break;
		}
		w=w1;
	}
	
	if(n==10000)
	 cout << "ifuncf not converge" << endl;
	
	return(w);
}

//function to set fluid from geometrical variables by using constraint equations
void Fmv0::set_enemomini()
{
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				double 
				gxx_p,gyy_p,gzz_p,
				gxy_p,gxz_p,gyz_p,
				gyx_p,gzx_p,gzy_p,wa_p,
				akxx_p,akyy_p,akzz_p,
				akxy_p,akxz_p,akyz_p,
				akyx_p,akzx_p,akzy_p,ek_p;

				double zgx_p,zgy_p,zgz_p;
				
				double dxx_x_p,dyy_x_p,dzz_x_p,
				dxy_x_p,dxz_x_p,dyz_x_p,
				dyx_x_p,dzx_x_p,dzy_x_p,wa_x_p,
				ek_x_p,
				zgx_x_p,zgy_x_p,zgz_x_p;
				double dxx_y_p,dyy_y_p,dzz_y_p,
				dxy_y_p,dxz_y_p,dyz_y_p,
				dyx_y_p,dzx_y_p,dzy_y_p,wa_y_p,
				ek_y_p,
				zgx_y_p,zgy_y_p,zgz_y_p;
				double dxx_z_p,dyy_z_p,dzz_z_p,
				dxy_z_p,dxz_z_p,dyz_z_p,
				dyx_z_p,dzx_z_p,dzy_z_p,wa_z_p,
				ek_z_p,
				zgx_z_p,zgy_z_p,zgz_z_p;

				double wa_xx,wa_yy,wa_zz,
				wa_xy,wa_xz,wa_yz;

				double gxx_xx,gxx_yy,gxx_zz,
				gxx_xy,gxx_xz,gxx_yz,
				gxx_yx,gxx_zx,gxx_zy;
				double gyy_xx,gyy_yy,gyy_zz,
				gyy_xy,gyy_xz,gyy_yz,
				gyy_yx,gyy_zx,gyy_zy;
				double gzz_xx,gzz_yy,gzz_zz,
				gzz_xy,gzz_xz,gzz_yz,
				gzz_yx,gzz_zx,gzz_zy;
				double gxy_xx,gxy_yy,gxy_zz,
				gxy_xy,gxy_xz,gxy_yz,
				gxy_yx,gxy_zx,gxy_zy;
				double gxz_xx,gxz_yy,gxz_zz,
				gxz_xy,gxz_xz,gxz_yz,
				gxz_yx,gxz_zx,gxz_zy;
				double gyz_xx,gyz_yy,gyz_zz,
				gyz_xy,gyz_xz,gyz_yz,
				gyz_yx,gyz_zx,gyz_zy;

				double det_p,det_pi;
				double gixx_p,giyy_p,gizz_p,
				gixy_p,gixz_p,giyz_p,
				giyx_p,gizx_p,gizy_p;

				double akx_ux_p,aky_uy_p,akz_uz_p,
				akx_uy_p,akx_uz_p,aky_uz_p,
				aky_ux_p,akz_ux_p,akz_uy_p;
				
				double crdx_xx,crdx_yy,crdx_zz,
				crdx_xy,crdx_xz,crdx_yz;
				double crdy_xx,crdy_yy,crdy_zz,
				crdy_xy,crdy_xz,crdy_yz;
				double crdz_xx,crdz_yy,crdz_zz,
				crdz_xy,crdz_xz,crdz_yz;

				double crx_xx,crx_yy,crx_zz,
				crx_xy,crx_xz,crx_yz;
				double cry_xx,cry_yy,cry_zz,
				cry_xy,cry_xz,cry_yz;
				double crz_xx,crz_yy,crz_zz,
				crz_xy,crz_xz,crz_yz;
				
				double gamma0_x,gamma0_y,gamma0_z;

				double wa_cdxx,wa_cdyy,wa_cdzz,
				wa_cdxy,wa_cdxz,wa_cdyz,
				wa_cdyx,wa_cdzx,wa_cdzy,
				wa_lap,wawa;

				double rc_xx_p,rc_yy_p,rc_zz_p,
				rc_xy_p,rc_xz_p,rc_yz_p,
				rc_yx_p,rc_zx_p,rc_zy_p;
				double aaaa,ricci;
				
				double fxx,fyy,fzz,dfxx,dfyy,dfzz;
				double Gam_ux_xx,Gam_uy_yy,Gam_uz_zz;
				double delG_x_ux_xx,delG_y_uy_yy,delG_z_uz_zz;
				double Dgam_x_xx,Dgam_x_yy,Dgam_x_zz,Dgam_x_xy,Dgam_x_xz,Dgam_x_yz,Dgam_x_yx,Dgam_x_zx,Dgam_x_zy,
						Dgam_y_xx,Dgam_y_yy,Dgam_y_zz,Dgam_y_xy,Dgam_y_xz,Dgam_y_yz,Dgam_y_yx,Dgam_y_zx,Dgam_y_zy,
						Dgam_z_xx,Dgam_z_yy,Dgam_z_zz,Dgam_z_xy,Dgam_z_xz,Dgam_z_yz,Dgam_z_yx,Dgam_z_zx,Dgam_z_zy;
				double Del_ux_xx,Del_ux_yy,Del_ux_zz,Del_ux_xy,Del_ux_xz,Del_ux_yz,Del_ux_yx,Del_ux_zx,Del_ux_zy,
						Del_uy_xx,Del_uy_yy,Del_uy_zz,Del_uy_xy,Del_uy_xz,Del_uy_yz,Del_uy_yx,Del_uy_zx,Del_uy_zy,
						Del_uz_xx,Del_uz_yy,Del_uz_zz,Del_uz_xy,Del_uz_xz,Del_uz_yz,Del_uz_yx,Del_uz_zx,Del_uz_zy;
				double Dgam_x_uxx,Dgam_x_uyy,Dgam_x_uzz,Dgam_x_uxy,Dgam_x_uxz,Dgam_x_uyz,
						Dgam_y_uxx,Dgam_y_uyy,Dgam_y_uzz,Dgam_y_uxy,Dgam_y_uxz,Dgam_y_uyz,
						Dgam_z_uxx,Dgam_z_uyy,Dgam_z_uzz,Dgam_z_uxy,Dgam_z_uxz,Dgam_z_uyz;
				double lapgam_xx,lapgam_yy,lapgam_zz,
						lapgam_xy,lapgam_xz,lapgam_yz;
				double gDDg_xx,gDDg_yy,gDDg_zz,
						gDDg_xy,gDDg_xz,gDDg_yz;
				double DGam_x_ux,DGam_y_uy,DGam_z_uz,
						DGam_x_uy,DGam_x_uz,DGam_y_uz,
						DGam_y_ux,DGam_z_ux,DGam_z_uy;
				
				//\tilde A_ij,k=\del_k \tilde A_ij
				double daxx_x,dayy_x,dazz_x,daxy_x,daxz_x,dayz_x,
					daxx_y,dayy_y,dazz_y,daxy_y,daxz_y,dayz_y,
					daxx_z,dayy_z,dazz_z,daxy_z,daxz_z,dayz_z;
				
				//D^i wa
				double Daxx_x,Dayy_x,Dazz_x,Daxy_x,Daxz_x,Dayz_x;
				double Daxx_y,Dayy_y,Dazz_y,Daxy_y,Daxz_y,Dayz_y;
				double Daxx_z,Dayy_z,Dazz_z,Daxy_z,Daxz_z,Dayz_z;
				
				double Da_x,Da_y,Da_z,M_x,M_y,M_z;
				
				double Ene=0.;
				double p_x=0.;
				double p_y=0.;
				double p_z=0.;
				double p2,EpP,M_ux,M_uy,M_uz;
				
				double rho,Gam,sqgam;
				double ewa4i;

				gxx_p=  get_bv(l,k,j, 7)+get_flat_df2x(j);
				gyy_p=  get_bv(l,k,j, 8)+get_flat_df2y(k);
				gzz_p=  get_bv(l,k,j, 9)+get_flat_df2z(l);
				gxy_p=  get_bv(l,k,j,10);
				gxz_p=  get_bv(l,k,j,11);
				gyz_p=  get_bv(l,k,j,12);
				gyx_p=	gxy_p;
				gzx_p=	gxz_p;
				gzy_p=	gyz_p;
				wa_p=   get_bv(l,k,j,13);

				ewa4i=exp(-4.*wa_p);

				akxx_p= get_bv(l,k,j,14);
				akyy_p= get_bv(l,k,j,15);
				akzz_p= get_bv(l,k,j,16);
				akxy_p= get_bv(l,k,j,17);
				akxz_p= get_bv(l,k,j,18);
				akyz_p= get_bv(l,k,j,19);
				akyx_p=	akxy_p;
				akzx_p=	akxz_p;
				akzy_p=	akyz_p;
				ek_p=   get_bv(l,k,j,20);
				
				zgx_p=  get_bv(l,k,j,21);
				zgy_p=  get_bv(l,k,j,22);
				zgz_p=  get_bv(l,k,j,23);

				//\del_x \Gamma^i
				zgx_x_p=get_f_x(l,k,j,21);
				zgy_x_p=get_f_x(l,k,j,22);
				zgz_x_p=get_f_x(l,k,j,23);

				// \del_y \Gamma^i
				zgx_y_p=get_f_y(l,k,j,21);
				zgy_y_p=get_f_y(l,k,j,22);
				zgz_y_p=get_f_y(l,k,j,23);

				// \del_z \Gamma^i
				zgx_z_p=get_f_z(l,k,j,21);
				zgy_z_p=get_f_z(l,k,j,22);
				zgz_z_p=get_f_z(l,k,j,23);

				//fij for inhomogeneous grid 
				fxx=get_flat_df2x(j);
				fyy=get_flat_df2y(k);
				fzz=get_flat_df2z(l);

				//bar Gamma_ui_jk for inhomogeneous grid 
				Gam_ux_xx=get_flat_Gamx(j);
				Gam_uy_yy=get_flat_Gamy(k);
				Gam_uz_zz=get_flat_Gamz(l);
				
				//calculation of derivatives of flat Christoffel
				//double delG_1_u2_34
				delG_x_ux_xx=get_flat_dGamx(j);
				delG_y_uy_yy=get_flat_dGamy(k);
				delG_z_uz_zz=get_flat_dGamz(l);

				dfxx=2.*fxx*Gam_ux_xx;
				dfyy=2.*fyy*Gam_uy_yy;
				dfzz=2.*fzz*Gam_uz_zz;
				
				// \del_x \tilde{gamma}_{ij} *0.5 
				dxx_x_p=(get_f_x(l,k,j,7)+dfxx)*0.5;
				dyy_x_p=get_f_x(l,k,j,8)*0.5;
				dzz_x_p=get_f_x(l,k,j,9)*0.5;
				dxy_x_p=get_f_x(l,k,j,10)*0.5;
				dxz_x_p=get_f_x(l,k,j,11)*0.5;
				dyz_x_p=get_f_x(l,k,j,12)*0.5;
				dyx_x_p=dxy_x_p;
				dzx_x_p=dxz_x_p;
				dzy_x_p=dyz_x_p;

				// \del_x \psi 
				wa_x_p=get_f_x(l,k,j,13);

				// \del_x K 
				ek_x_p=get_f_x(l,k,j,20);
				// \del_y \tilde \gamma_{ij} * 0.5 
				dxx_y_p=get_f_y(l,k,j,7)*0.5;
				dyy_y_p=(get_f_y(l,k,j,8)+dfyy)*0.5;
				dzz_y_p=get_f_y(l,k,j,9)*0.5;
				dxy_y_p=get_f_y(l,k,j,10)*0.5;
				dxz_y_p=get_f_y(l,k,j,11)*0.5;
				dyz_y_p=get_f_y(l,k,j,12)*0.5;
				dyx_y_p=dxy_y_p;
				dzx_y_p=dxz_y_p;
				dzy_y_p=dyz_y_p;

				// \del_y \psi
				wa_y_p=get_f_y(l,k,j,13);

				// \del_y trK
				ek_y_p=get_f_y(l,k,j,20);
				
				// \del_z \tilde \gamma_{ij} * 0.5
				dxx_z_p=get_f_z(l,k,j,7)*0.5;
				dyy_z_p=get_f_z(l,k,j,8)*0.5;
				dzz_z_p=(get_f_z(l,k,j,9)+dfzz)*0.5;
				dxy_z_p=get_f_z(l,k,j,10)*0.5;
				dxz_z_p=get_f_z(l,k,j,11)*0.5;
				dyz_z_p=get_f_z(l,k,j,12)*0.5;
				dyx_z_p=dxy_z_p;
				dzx_z_p=dxz_z_p;
				dzy_z_p=dyz_z_p;

				// \del_z \psi
				wa_z_p=get_f_z(l,k,j,13);
				
				// \del_z trK
				ek_z_p=get_f_z(l,k,j,20);

				//\del_x \del_x g_ij
				gxx_xx=get_f_xx(l,k,j,7)+2.*fxx*(delG_x_ux_xx+2.*pow(Gam_ux_xx,2));
				gyy_xx=get_f_xx(l,k,j,8);
				gzz_xx=get_f_xx(l,k,j,9);
				gxy_xx=get_f_xx(l,k,j,10);
				gxz_xx=get_f_xx(l,k,j,11);
				gyz_xx=get_f_xx(l,k,j,12);
				
				//\del_y \del_y g_ij
				gxx_yy=get_f_yy(l,k,j,7);
				gyy_yy=get_f_yy(l,k,j,8)+2.*fyy*(delG_y_uy_yy+2.*pow(Gam_uy_yy,2));
				gzz_yy=get_f_yy(l,k,j,9);
				gxy_yy=get_f_yy(l,k,j,10);
				gxz_yy=get_f_yy(l,k,j,11);
				gyz_yy=get_f_yy(l,k,j,12);
				
				//\del_z \del_z g_ij
				gxx_zz=get_f_zz(l,k,j,7);
				gyy_zz=get_f_zz(l,k,j,8);
				gzz_zz=get_f_zz(l,k,j,9)+2.*fzz*(delG_z_uz_zz+2.*pow(Gam_uz_zz,2));
				gxy_zz=get_f_zz(l,k,j,10);
				gxz_zz=get_f_zz(l,k,j,11);
				gyz_zz=get_f_zz(l,k,j,12);

				//\del_x\del_y g_ij
				gxx_xy=get_f_xy(l,k,j,7);
				gyy_xy=get_f_xy(l,k,j,8);
				gzz_xy=get_f_xy(l,k,j,9);
				gxy_xy=get_f_xy(l,k,j,10);
				gxz_xy=get_f_xy(l,k,j,11);
				gyz_xy=get_f_xy(l,k,j,12);
				
				gxx_yx=gxx_xy;
				gxy_yx=gxy_xy;
				gxz_yx=gxz_xy;
				gyy_yx=gyy_xy;
				gyz_yx=gyz_xy;
				gzz_yx=gzz_xy;

				//\del_x \del_z g_ij
				gxx_xz=get_f_xz(l,k,j,7);
				gyy_xz=get_f_xz(l,k,j,8);
				gzz_xz=get_f_xz(l,k,j,9);
				gxy_xz=get_f_xz(l,k,j,10);
				gxz_xz=get_f_xz(l,k,j,11);
				gyz_xz=get_f_xz(l,k,j,12);
				
				gxx_zx=gxx_xz;
				gxy_zx=gxy_xz;
				gxz_zx=gxz_xz;
				gyy_zx=gyy_xz;
				gyz_zx=gyz_xz;
				gzz_zx=gzz_xz;

				//\del_y\del_z g_ij
				gxx_yz=get_f_yz(l,k,j,7);
				gyy_yz=get_f_yz(l,k,j,8);
				gzz_yz=get_f_yz(l,k,j,9);
				gxy_yz=get_f_yz(l,k,j,10);
				gxz_yz=get_f_yz(l,k,j,11);
				gyz_yz=get_f_yz(l,k,j,12);
				
				gxx_zy=gxx_yz;
				gxy_zy=gxy_yz;
				gxz_zy=gxz_yz;
				gyy_zy=gyy_yz;
				gyz_zy=gyz_yz;
				gzz_zy=gzz_yz;

				//\del_i\del_j \psi
				wa_xx=get_f_xx(l,k,j,13);
				wa_yy=get_f_yy(l,k,j,13);
				wa_zz=get_f_zz(l,k,j,13);
				wa_xy=get_f_xy(l,k,j,13);
				wa_xz=get_f_xz(l,k,j,13);
				wa_yz=get_f_yz(l,k,j,13);

				// tilde{gamma}^{ij} -> gi
			//	det_p=gxx_p*gyy_p*gzz_p+gxy_p*gyz_p*gzx_p+gxz_p*gyx_p*gzy_p
			//			-gxz_p*gyy_p*gzx_p-gxy_p*gyx_p*gzz_p-gxx_p*gyz_p*gzy_p;
			//	if(det_p<1.e-16) 
			//	 det_p=1.e-16;
				det_p=get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l);
				det_pi=1./det_p;
				gixx_p= (gyy_p*gzz_p-gyz_p*gzy_p)*det_pi;
				giyy_p= (gxx_p*gzz_p-gxz_p*gzx_p)*det_pi;
				gizz_p= (gxx_p*gyy_p-gxy_p*gyx_p)*det_pi;
				gixy_p=-(gyx_p*gzz_p-gyz_p*gzx_p)*det_pi;
				gixz_p= (gyx_p*gzy_p-gyy_p*gzx_p)*det_pi;
				giyz_p=-(gxx_p*gzy_p-gxy_p*gzx_p)*det_pi;
				
				giyx_p=gixy_p;
				gizx_p=gixz_p;
				gizy_p=giyz_p;

				//  \tilde{A}_i^k=\tilde{A}_{ij}\tilde{\gamma}^{jk} -> ak{i}_u{k}
				akx_ux_p=akxx_p*gixx_p +akxy_p*giyx_p +akxz_p*gizx_p;
				aky_ux_p=akyx_p*gixx_p +akyy_p*giyx_p +akyz_p*gizx_p;
				akz_ux_p=akzx_p*gixx_p +akzy_p*giyx_p +akzz_p*gizx_p;
				akx_uy_p=akxx_p*gixy_p +akxy_p*giyy_p +akxz_p*gizy_p;
				aky_uy_p=akyx_p*gixy_p +akyy_p*giyy_p +akyz_p*gizy_p;
				akz_uy_p=akzx_p*gixy_p +akzy_p*giyy_p +akzz_p*gizy_p;
				akx_uz_p=akxx_p*gixz_p +akxy_p*giyz_p +akxz_p*gizz_p;
				aky_uz_p=akyx_p*gixz_p +akyy_p*giyz_p +akyz_p*gizz_p;
				akz_uz_p=akzx_p*gixz_p +akzy_p*giyz_p +akzz_p*gizz_p;

				//  \tilde{\Gamma}_{ijk}=\tilde{\gamma}_{im}\tilde{\Gamma}^m_{jk} -> crd{i}_{jk}
				// NOTE:d.._._p is 0.5*derivatives
				crdx_xx=dxx_x_p+dxx_x_p -dxx_x_p;
				crdx_yy=dxy_y_p+dxy_y_p -dyy_x_p;
				crdx_zz=dxz_z_p+dxz_z_p -dzz_x_p;
				crdx_xy=dxx_y_p+dxy_x_p -dxy_x_p;
				crdx_xz=dxx_z_p+dxz_x_p -dxz_x_p;
				crdx_yz=dxy_z_p+dxz_y_p -dyz_x_p;
				
				crdy_xx=dyx_x_p+dyx_x_p -dxx_y_p;
				crdy_yy=dyy_y_p+dyy_y_p -dyy_y_p;
				crdy_zz=dyz_z_p+dyz_z_p -dzz_y_p;
				crdy_xy=dyx_y_p+dyy_x_p -dxy_y_p;
				crdy_xz=dyx_z_p+dyz_x_p -dxz_y_p;
				crdy_yz=dyy_z_p+dyz_y_p -dyz_y_p;
				
				crdz_xx=dzx_x_p+dzx_x_p -dxx_z_p;
				crdz_yy=dzy_y_p+dzy_y_p -dyy_z_p;
				crdz_zz=dzz_z_p+dzz_z_p -dzz_z_p;
				crdz_xy=dzx_y_p+dzy_x_p -dxy_z_p;
				crdz_xz=dzx_z_p+dzz_x_p -dxz_z_p;
				crdz_yz=dzy_z_p+dzz_y_p -dyz_z_p;
				
				//  \tilde{\Gamma}^i_{jk} -> cr{i}_{jk}
				crx_xx=gixx_p*crdx_xx +gixy_p*crdy_xx +gixz_p*crdz_xx;
				crx_yy=gixx_p*crdx_yy +gixy_p*crdy_yy +gixz_p*crdz_yy;
				crx_zz=gixx_p*crdx_zz +gixy_p*crdy_zz +gixz_p*crdz_zz;
				crx_xy=gixx_p*crdx_xy +gixy_p*crdy_xy +gixz_p*crdz_xy;
				crx_xz=gixx_p*crdx_xz +gixy_p*crdy_xz +gixz_p*crdz_xz;
				crx_yz=gixx_p*crdx_yz +gixy_p*crdy_yz +gixz_p*crdz_yz;
				
				cry_xx=giyx_p*crdx_xx +giyy_p*crdy_xx +giyz_p*crdz_xx;
				cry_yy=giyx_p*crdx_yy +giyy_p*crdy_yy +giyz_p*crdz_yy;
				cry_zz=giyx_p*crdx_zz +giyy_p*crdy_zz +giyz_p*crdz_zz;
				cry_xy=giyx_p*crdx_xy +giyy_p*crdy_xy +giyz_p*crdz_xy;
				cry_xz=giyx_p*crdx_xz +giyy_p*crdy_xz +giyz_p*crdz_xz;
				cry_yz=giyx_p*crdx_yz +giyy_p*crdy_yz +giyz_p*crdz_yz;
				
				crz_xx=gizx_p*crdx_xx +gizy_p*crdy_xx +gizz_p*crdz_xx;
				crz_yy=gizx_p*crdx_yy +gizy_p*crdy_yy +gizz_p*crdz_yy;
				crz_zz=gizx_p*crdx_zz +gizy_p*crdy_zz +gizz_p*crdz_zz;
				crz_xy=gizx_p*crdx_xy +gizy_p*crdy_xy +gizz_p*crdz_xy;
				crz_xz=gizx_p*crdx_xz +gizy_p*crdy_xz +gizz_p*crdz_xz;
				crz_yz=gizx_p*crdx_yz +gizy_p*crdy_yz +gizz_p*crdz_yz;
				
				//cal D_i tilde gamma_jk
				
				//Dgam_1_11=2.*(d11_1_p-Gam_u1_11*g11_p);
				Dgam_x_xx=2.*(dxx_x_p-Gam_ux_xx*gxx_p);
				
				Dgam_x_yy=2.*dyy_x_p;
				Dgam_x_zz=2.*dzz_x_p;

				//Dgam_1_12=2.*d12_1_p-Gam_u1_11*g12_p;
				Dgam_x_xy=2.*dxy_x_p-Gam_ux_xx*gxy_p;
				Dgam_x_xz=2.*dxz_x_p-Gam_ux_xx*gxz_p;
				
				Dgam_x_yz=2.*dyz_x_p;
				
				Dgam_x_yx=Dgam_x_xy;
				Dgam_x_zx=Dgam_x_xz;
				Dgam_x_zy=Dgam_x_yz;
				
				Dgam_y_xx=2.*dxx_y_p;
				Dgam_y_yy=2.*(dyy_y_p-Gam_uy_yy*gyy_p);
				Dgam_y_zz=2.*dzz_y_p;
				Dgam_y_xy=2.*dyx_y_p-Gam_uy_yy*gyx_p;
				Dgam_y_xz=2.*dxz_y_p;
				Dgam_y_yz=2.*dyz_y_p-Gam_uy_yy*gyz_p;
				Dgam_y_yx=Dgam_y_xy;
				Dgam_y_zx=Dgam_y_xz;
				Dgam_y_zy=Dgam_y_yz;
				
				Dgam_z_xx=2.*dxx_z_p;
				Dgam_z_yy=2.*dyy_z_p;
				Dgam_z_zz=2.*(dzz_z_p-Gam_uz_zz*gzz_p);
				Dgam_z_xy=2.*dxy_z_p;
				Dgam_z_xz=2.*dzx_z_p-Gam_uz_zz*gzx_p;
				Dgam_z_yz=2.*dzy_z_p-Gam_uz_zz*gzy_p;
				Dgam_z_yx=Dgam_z_xy;
				Dgam_z_zx=Dgam_z_xz;
				Dgam_z_zy=Dgam_z_yz;
				
				//Delta_ui_jk=tilde Gamma - bar Gamma
				
				//Del_ux_xx=0.5*(gixx_p*Dgam_x_xx+gixy_p*(2.*Dgam_x_yx-Dgam_y_xx)+gixz_p*(2.*Dgam_x_zx-Dgam_z_xx));
				Del_ux_xx=crx_xx-Gam_ux_xx;
				Del_ux_yy=crx_yy;
				Del_ux_zz=crx_zz;
				Del_ux_xy=crx_xy;
				Del_ux_xz=crx_xz;
				Del_ux_yz=crx_yz;
				Del_ux_yx=crx_xy;
				Del_ux_zx=crx_xz;
				Del_ux_zy=crx_yz;
				
				Del_uy_xx=cry_xx;
				//Del_uy_yy=0.5*(giyy_p*Dgam_y_yy+giyz_p*(2.*Dgam_y_zy-Dgam_z_yy)+giyx_p*(2.*Dgam_y_xy-Dgam_x_yy));
				Del_uy_yy=cry_yy-Gam_uy_yy;
				Del_uy_zz=cry_zz;
				Del_uy_xy=cry_xy;
				Del_uy_xz=cry_xz;
				Del_uy_yz=cry_yz;
				Del_uy_yx=cry_xy;
				Del_uy_zx=cry_xz;
				Del_uy_zy=cry_yz;
				
				Del_uz_xx=crz_xx;
				Del_uz_yy=crz_yy;
				//Del_uz_zz=0.5*(gizz_p*Dgam_z_zz+gizx_p*(2.*Dgam_z_xz-Dgam_x_zz)+gizy_p*(2.*Dgam_z_yz-Dgam_y_zz));
				Del_uz_zz=crz_zz-Gam_uz_zz;
				Del_uz_xy=crz_xy;
				Del_uz_xz=crz_xz;
				Del_uz_yz=crz_yz;
				Del_uz_yx=crz_xy;
				Del_uz_zx=crz_xz;
				Del_uz_zy=crz_yz;
				
				//tilde Gamma0^i =- cal D_j tilde gamma ^ji=tilde gamma^jk Delta^i_jk
				gamma0_x=gixx_p*Del_ux_xx+giyy_p*Del_ux_yy+gizz_p*Del_ux_zz
						+2.*(gixy_p*Del_ux_xy+gixz_p*Del_ux_xz+giyz_p*Del_ux_yz);
				gamma0_y=gixx_p*Del_uy_xx+giyy_p*Del_uy_yy+gizz_p*Del_uy_zz
						+2.*(gixy_p*Del_uy_xy+gixz_p*Del_uy_xz+giyz_p*Del_uy_yz);
				gamma0_z=gixx_p*Del_uz_xx+giyy_p*Del_uz_yy+gizz_p*Del_uz_zz
						+2.*(gixy_p*Del_uz_xy+gixz_p*Del_uz_xz+giyz_p*Del_uz_yz);
				
				// Gamma^i
				set_bv(l,k,j,21)=gamma0_x;
				set_bv(l,k,j,22)=gamma0_y;
				set_bv(l,k,j,23)=gamma0_z;

				// \tilde D_i D_j \psi = \psi,ij -\tilde{\Gamma}^k_ij \psi_k
				// ok even for non-Cartesian
				wa_cdxx=wa_xx-(crx_xx*wa_x_p +cry_xx*wa_y_p +crz_xx*wa_z_p);
				wa_cdyy=wa_yy-(crx_yy*wa_x_p +cry_yy*wa_y_p +crz_yy*wa_z_p);
				wa_cdzz=wa_zz-(crx_zz*wa_x_p +cry_zz*wa_y_p +crz_zz*wa_z_p);
				wa_cdxy=wa_xy-(crx_xy*wa_x_p +cry_xy*wa_y_p +crz_xy*wa_z_p);
				wa_cdxz=wa_xz-(crx_xz*wa_x_p +cry_xz*wa_y_p +crz_xz*wa_z_p);
				wa_cdyz=wa_yz-(crx_yz*wa_x_p +cry_yz*wa_y_p +crz_yz*wa_z_p);
				wa_cdyx=wa_cdxy;
				wa_cdzx=wa_cdxz;
				wa_cdzy=wa_cdyz;
				
				// \tilde Laplacian \psi
				wa_lap=gixx_p*wa_cdxx +gixy_p*wa_cdxy +gixz_p*wa_cdxz
						+giyx_p*wa_cdyx +giyy_p*wa_cdyy +giyz_p*wa_cdyz
						+gizx_p*wa_cdzx +gizy_p*wa_cdzy +gizz_p*wa_cdzz;
				
				// (\tilde D_i \psi)(\tilde D^i \psi)
				wawa=gixx_p*wa_x_p*wa_x_p +gixy_p*wa_x_p*wa_y_p +gixz_p*wa_x_p*wa_z_p
					+giyx_p*wa_y_p*wa_x_p +giyy_p*wa_y_p*wa_y_p +giyz_p*wa_y_p*wa_z_p
					+gizx_p*wa_z_p*wa_x_p +gizy_p*wa_z_p*wa_y_p +gizz_p*wa_z_p*wa_z_p;

				//cal D_i tilde gamma^jk
				
				//Dgam_1_u23=-Del_u2_1x*gix3_p-Del_u2_1y*giy3_p-Del_u2_1z*giz3_p-Del_u3_1x*gix2_p-Del_u3_1y*giy2_p-Del_u3_1z*giz2_p;
				Dgam_x_uxx=-Del_ux_xx*gixx_p-Del_ux_xy*giyx_p-Del_ux_xz*gizx_p-Del_ux_xx*gixx_p-Del_ux_xy*giyx_p-Del_ux_xz*gizx_p;
				Dgam_x_uyy=-Del_uy_xx*gixy_p-Del_uy_xy*giyy_p-Del_uy_xz*gizy_p-Del_uy_xx*gixy_p-Del_uy_xy*giyy_p-Del_uy_xz*gizy_p;
				Dgam_x_uzz=-Del_uz_xx*gixz_p-Del_uz_xy*giyz_p-Del_uz_xz*gizz_p-Del_uz_xx*gixz_p-Del_uz_xy*giyz_p-Del_uz_xz*gizz_p;
				Dgam_x_uxy=-Del_ux_xx*gixy_p-Del_ux_xy*giyy_p-Del_ux_xz*gizy_p-Del_uy_xx*gixx_p-Del_uy_xy*giyx_p-Del_uy_xz*gizx_p;
				Dgam_x_uxz=-Del_ux_xx*gixz_p-Del_ux_xy*giyz_p-Del_ux_xz*gizz_p-Del_uz_xx*gixx_p-Del_uz_xy*giyx_p-Del_uz_xz*gizx_p;
				Dgam_x_uyz=-Del_uy_xx*gixz_p-Del_uy_xy*giyz_p-Del_uy_xz*gizz_p-Del_uz_xx*gixy_p-Del_uz_xy*giyy_p-Del_uz_xz*gizy_p;
				//Dgam_x_uyx=Dgam_x_uxy;
				//Dgam_x_uzx=Dgam_x_uxz;
				//Dgam_x_uzy=Dgam_x_uyz;
				
				Dgam_y_uxx=-Del_ux_yx*gixx_p-Del_ux_yy*giyx_p-Del_ux_yz*gizx_p-Del_ux_yx*gixx_p-Del_ux_yy*giyx_p-Del_ux_yz*gizx_p;
				Dgam_y_uyy=-Del_uy_yx*gixy_p-Del_uy_yy*giyy_p-Del_uy_yz*gizy_p-Del_uy_yx*gixy_p-Del_uy_yy*giyy_p-Del_uy_yz*gizy_p;
				Dgam_y_uzz=-Del_uz_yx*gixz_p-Del_uz_yy*giyz_p-Del_uz_yz*gizz_p-Del_uz_yx*gixz_p-Del_uz_yy*giyz_p-Del_uz_yz*gizz_p;
				Dgam_y_uxy=-Del_ux_yx*gixy_p-Del_ux_yy*giyy_p-Del_ux_yz*gizy_p-Del_uy_yx*gixx_p-Del_uy_yy*giyx_p-Del_uy_yz*gizx_p;
				Dgam_y_uxz=-Del_ux_yx*gixz_p-Del_ux_yy*giyz_p-Del_ux_yz*gizz_p-Del_uz_yx*gixx_p-Del_uz_yy*giyx_p-Del_uz_yz*gizx_p;
				Dgam_y_uyz=-Del_uy_yx*gixz_p-Del_uy_yy*giyz_p-Del_uy_yz*gizz_p-Del_uz_yx*gixy_p-Del_uz_yy*giyy_p-Del_uz_yz*gizy_p;
				//Dgam_y_uyx=Dgam_y_uxy;
				//Dgam_y_uzx=Dgam_y_uxz;
				//Dgam_y_uzy=Dgam_y_uyz;
				
				Dgam_z_uxx=-Del_ux_zx*gixx_p-Del_ux_zy*giyx_p-Del_ux_zz*gizx_p-Del_ux_zx*gixx_p-Del_ux_zy*giyx_p-Del_ux_zz*gizx_p;
				Dgam_z_uyy=-Del_uy_zx*gixy_p-Del_uy_zy*giyy_p-Del_uy_zz*gizy_p-Del_uy_zx*gixy_p-Del_uy_zy*giyy_p-Del_uy_zz*gizy_p;
				Dgam_z_uzz=-Del_uz_zx*gixz_p-Del_uz_zy*giyz_p-Del_uz_zz*gizz_p-Del_uz_zx*gixz_p-Del_uz_zy*giyz_p-Del_uz_zz*gizz_p;
				Dgam_z_uxy=-Del_ux_zx*gixy_p-Del_ux_zy*giyy_p-Del_ux_zz*gizy_p-Del_uy_zx*gixx_p-Del_uy_zy*giyx_p-Del_uy_zz*gizx_p;
				Dgam_z_uxz=-Del_ux_zx*gixz_p-Del_ux_zy*giyz_p-Del_ux_zz*gizz_p-Del_uz_zx*gixx_p-Del_uz_zy*giyx_p-Del_uz_zz*gizx_p;
				Dgam_z_uyz=-Del_uy_zx*gixz_p-Del_uy_zy*giyz_p-Del_uy_zz*gizz_p-Del_uz_zx*gixy_p-Del_uz_zy*giyy_p-Del_uz_zz*gizy_p;
				//Dgam_z_uyx=Dgam_z_uxy;
				//Dgam_z_uzx=Dgam_z_uxz;
				//Dgam_z_uzy=Dgam_z_uyz;
				

				//\tilde gamma^kl \del_k \del_l \tilde \gamma_ij
				
				lapgam_xx=gixx_p*gxx_xx +gixy_p*gxx_xy +gixz_p*gxx_xz
						+giyx_p*gxx_yx +giyy_p*gxx_yy +giyz_p*gxx_yz
						+gizx_p*gxx_zx +gizy_p*gxx_zy +gizz_p*gxx_zz;
				
				lapgam_yy=gixx_p*gyy_xx +gixy_p*gyy_xy +gixz_p*gyy_xz 
						+giyx_p*gyy_yx +giyy_p*gyy_yy +giyz_p*gyy_yz 
						+gizx_p*gyy_zx +gizy_p*gyy_zy +gizz_p*gyy_zz;
				
				lapgam_zz=gixx_p*gzz_xx +gixy_p*gzz_xy +gixz_p*gzz_xz 
						+giyx_p*gzz_yx +giyy_p*gzz_yy +giyz_p*gzz_yz 
						+gizx_p*gzz_zx +gizy_p*gzz_zy +gizz_p*gzz_zz;
				
				lapgam_xy=gixx_p*gxy_xx +gixy_p*gxy_xy +gixz_p*gxy_xz 
						+giyx_p*gxy_yx +giyy_p*gxy_yy +giyz_p*gxy_yz 
						+gizx_p*gxy_zx +gizy_p*gxy_zy +gizz_p*gxy_zz;
				
				lapgam_xz=gixx_p*gxz_xx +gixy_p*gxz_xy +gixz_p*gxz_xz 
						+giyx_p*gxz_yx +giyy_p*gxz_yy +giyz_p*gxz_yz 
						+gizx_p*gxz_zx +gizy_p*gxz_zy +gizz_p*gxz_zz;
				
				lapgam_yz=gixx_p*gyz_xx +gixy_p*gyz_xy +gixz_p*gyz_xz 
						+giyx_p*gyz_yx +giyy_p*gyz_yy +giyz_p*gyz_yz 
						+gizx_p*gyz_zx +gizy_p*gyz_zy +gizz_p*gyz_zz;
				
				
				//\tilde gamma^kl \cal D_k \cal D_l \tilde \gamma_ij
						
				//gDDg_11=lapgamma_11-2.*gi11_p*delG_1_u1_11*g11_p
				//			-8.*Gam_u1_11*(gix1_p*d11_x_p+giy1_p*d11_y_p+giz1_p*d11_z_p)+4.*gi11_p*pow(Gam_u1_11,2)*g11_p
				//			-(gixx_p*Gam_ux_xx*Dgam_x_11+giyy_p*Gam_uy_yy*Dgam_y_11+gizz_p*Gam_uz_zz*Dgam_z_11);
				gDDg_xx=lapgam_xx-2.*gixx_p*delG_x_ux_xx*gxx_p
							-8.*Gam_ux_xx*(gixx_p*dxx_x_p+giyx_p*dxx_y_p+gizx_p*dxx_z_p)+4.*gixx_p*pow(Gam_ux_xx,2)*gxx_p
							-(gixx_p*Gam_ux_xx*Dgam_x_xx+giyy_p*Gam_uy_yy*Dgam_y_xx+gizz_p*Gam_uz_zz*Dgam_z_xx);
				gDDg_yy=lapgam_yy-2.*giyy_p*delG_y_uy_yy*gyy_p
							-8.*Gam_uy_yy*(gixy_p*dyy_x_p+giyy_p*dyy_y_p+gizy_p*dyy_z_p)+4.*giyy_p*pow(Gam_uy_yy,2)*gyy_p
							-(gixx_p*Gam_ux_xx*Dgam_x_yy+giyy_p*Gam_uy_yy*Dgam_y_yy+gizz_p*Gam_uz_zz*Dgam_z_yy);
				gDDg_zz=lapgam_zz-2.*gizz_p*delG_z_uz_zz*gzz_p
							-8.*Gam_uz_zz*(gixz_p*dzz_x_p+giyz_p*dzz_y_p+gizz_p*dzz_z_p)+4.*gizz_p*pow(Gam_uz_zz,2)*gzz_p
							-(gixx_p*Gam_ux_xx*Dgam_x_zz+giyy_p*Gam_uy_yy*Dgam_y_zz+gizz_p*Gam_uz_zz*Dgam_z_zz);
				
				//gDDg_12=lapgam_12-gi11_p*delG_1_u1_11*g12_p-gi22_p*delG_2_u2_22*g12_p
				//			-4.*Gam_u1_11*(gix1_p*d12_x_p+giy1_p*d12_y_p+giz1_p*d12_z_p)-4.*Gam_u2_22*(gix2_p*d12_x_p+giy2_p*d12_y_p+giz2_p*d12_z_p)
				//			+(gi11_p*pow(Gam_u1_11,2)+2.*gi12_p*Gam_u1_11*Gam_u2_22+gi22_p*pow(Gam_u2_22,2))*g12_p
				//			-(gixx_p*Gam_ux_xx*Dgam_x_12+giyy_p*Gam_uy_yy*Dgam_y_12+gizz_p*Gam_uz_zz*Dgam_z_12);
				gDDg_xy=lapgam_xy-gixx_p*delG_x_ux_xx*gxy_p-giyy_p*delG_y_uy_yy*gxy_p
							-4.*Gam_ux_xx*(gixx_p*dxy_x_p+giyx_p*dxy_y_p+gizx_p*dxy_z_p)-4.*Gam_uy_yy*(gixy_p*dxy_x_p+giyy_p*dxy_y_p+gizy_p*dxy_z_p)
							+(gixx_p*pow(Gam_ux_xx,2)+2.*gixy_p*Gam_ux_xx*Gam_uy_yy+giyy_p*pow(Gam_uy_yy,2))*gxy_p
							-(gixx_p*Gam_ux_xx*Dgam_x_xy+giyy_p*Gam_uy_yy*Dgam_y_xy+gizz_p*Gam_uz_zz*Dgam_z_xy);
				gDDg_xz=lapgam_xz-gixx_p*delG_x_ux_xx*gxz_p-gizz_p*delG_z_uz_zz*gxz_p
							-4.*Gam_ux_xx*(gixx_p*dxz_x_p+giyx_p*dxz_y_p+gizx_p*dxz_z_p)-4.*Gam_uz_zz*(gixz_p*dxz_x_p+giyz_p*dxz_y_p+gizz_p*dxz_z_p)
							+(gixx_p*pow(Gam_ux_xx,2)+2.*gixz_p*Gam_ux_xx*Gam_uz_zz+gizz_p*pow(Gam_uz_zz,2))*gxz_p
							-(gixx_p*Gam_ux_xx*Dgam_x_xz+giyy_p*Gam_uy_yy*Dgam_y_xz+gizz_p*Gam_uz_zz*Dgam_z_xz);
				gDDg_yz=lapgam_yz-giyy_p*delG_y_uy_yy*gyz_p-gizz_p*delG_z_uz_zz*gyz_p
							-4.*Gam_uy_yy*(gixy_p*dyz_x_p+giyy_p*dyz_y_p+gizy_p*dyz_z_p)-4.*Gam_uz_zz*(gixz_p*dyz_x_p+giyz_p*dyz_y_p+gizz_p*dyz_z_p)
							+(giyy_p*pow(Gam_uy_yy,2)+2.*giyz_p*Gam_uy_yy*Gam_uz_zz+gizz_p*pow(Gam_uz_zz,2))*gyz_p
							-(gixx_p*Gam_ux_xx*Dgam_x_yz+giyy_p*Gam_uy_yy*Dgam_y_yz+gizz_p*Gam_uz_zz*Dgam_z_yz);
				
				//\cal D_i \tilde Gamma^j
				
				//DGam_1_u1=zg1_1_p+Gam_u1_11*zg1_p;
				DGam_x_ux=zgx_x_p+Gam_ux_xx*zgx_p;
				DGam_y_uy=zgy_y_p+Gam_uy_yy*zgy_p;
				DGam_z_uz=zgz_z_p+Gam_uz_zz*zgz_p;
				//DGam_1_u2=zg2_1_p;
				DGam_x_uy=zgy_x_p;
				DGam_x_uz=zgz_x_p;
				DGam_y_uz=zgz_y_p;
				DGam_y_ux=zgx_y_p;
				DGam_z_ux=zgx_z_p;
				DGam_z_uy=zgy_z_p;
				
				
				// exp(-4\psi)* R_{jk} -> rc_{ij};  !!! NOT R_{jk}
				//rc_12_p=0.5*(-gDDg_12
				//			+g1x_p*DGam_2_ux+g1y_p*DGam_2_uy+g1z_p*DGam_2_uz
				//			+g2x_p*DGam_1_ux+g2y_p*DGam_1_uy+g2z_p*DGam_1_uz)
				//		-0.5*(Dgam_x_x1*Dgam_2_uxx+Dgam_y_y1*Dgam_2_uyy+Dgam_z_z1*Dgam_2_uzz
				//			+Dgam_x_y1*Dgam_2_uxy+Dgam_x_z1*Dgam_2_uxz+Dgam_y_z1*Dgam_2_uyz
				//			+Dgam_y_x1*Dgam_2_uxy+Dgam_z_x1*Dgam_2_uxz+Dgam_z_y1*Dgam_2_uyz
				//			
				//			+Dgam_x_x2*Dgam_1_uxx+Dgam_y_y2*Dgam_1_uyy+Dgam_z_z2*Dgam_1_uzz
				//			+Dgam_x_y2*Dgam_1_uxy+Dgam_x_z2*Dgam_1_uxz+Dgam_y_z2*Dgam_1_uyz
				//			+Dgam_y_x2*Dgam_1_uxy+Dgam_z_x2*Dgam_1_uxz+Dgam_z_y2*Dgam_1_uyz
				//			
				//			-zgx_p*Dgam_x_12-zgy_p*Dgam_y_12-zgz_p*Dgam_z_12)
				//		
				//		-Del_ux_1x*Del_ux_x2-Del_uy_1y*Del_uy_y2-Del_uz_1z*Del_uz_z2
				//		-Del_ux_1y*Del_uy_x2-Del_ux_1z*Del_uz_x2-Del_uy_1z*Del_uz_y2
				//		-Del_uy_1x*Del_ux_y2-Del_uz_1x*Del_ux_z2-Del_uz_1y*Del_uy_z2;
				
				rc_xx_p=(0.5*(-gDDg_xx
							+gxx_p*DGam_x_ux+gxy_p*DGam_x_uy+gxz_p*DGam_x_uz
							+gxx_p*DGam_x_ux+gxy_p*DGam_x_uy+gxz_p*DGam_x_uz)
						-0.5*(Dgam_x_xx*Dgam_x_uxx+Dgam_y_yx*Dgam_x_uyy+Dgam_z_zx*Dgam_x_uzz
							+Dgam_x_yx*Dgam_x_uxy+Dgam_x_zx*Dgam_x_uxz+Dgam_y_zx*Dgam_x_uyz
							+Dgam_y_xx*Dgam_x_uxy+Dgam_z_xx*Dgam_x_uxz+Dgam_z_yx*Dgam_x_uyz
							
							+Dgam_x_xx*Dgam_x_uxx+Dgam_y_yx*Dgam_x_uyy+Dgam_z_zx*Dgam_x_uzz
							+Dgam_x_yx*Dgam_x_uxy+Dgam_x_zx*Dgam_x_uxz+Dgam_y_zx*Dgam_x_uyz
							+Dgam_y_xx*Dgam_x_uxy+Dgam_z_xx*Dgam_x_uxz+Dgam_z_yx*Dgam_x_uyz
							
							//-zgx_p*Dgam_x_xx-zgy_p*Dgam_y_xx-zgz_p*Dgam_z_xx)
							-gamma0_x*Dgam_x_xx-gamma0_y*Dgam_y_xx-gamma0_z*Dgam_z_xx)
						
						-Del_ux_xx*Del_ux_xx-Del_uy_xy*Del_uy_yx-Del_uz_xz*Del_uz_zx
						-Del_ux_xy*Del_uy_xx-Del_ux_xz*Del_uz_xx-Del_uy_xz*Del_uz_yx
						-Del_uy_xx*Del_ux_yx-Del_uz_xx*Del_ux_zx-Del_uz_xy*Del_uy_zx
						
						-2.*wa_cdxx -2.*wa_lap*gxx_p -4.*wawa*gxx_p +4.*wa_x_p*wa_x_p)*ewa4i;
				
				rc_yy_p=(0.5*(-gDDg_yy
							+gyx_p*DGam_y_ux+gyy_p*DGam_y_uy+gyz_p*DGam_y_uz
							+gyx_p*DGam_y_ux+gyy_p*DGam_y_uy+gyz_p*DGam_y_uz)
						-0.5*(Dgam_x_xy*Dgam_y_uxx+Dgam_y_yy*Dgam_y_uyy+Dgam_z_zy*Dgam_y_uzz
							+Dgam_x_yy*Dgam_y_uxy+Dgam_x_zy*Dgam_y_uxz+Dgam_y_zy*Dgam_y_uyz
							+Dgam_y_xy*Dgam_y_uxy+Dgam_z_xy*Dgam_y_uxz+Dgam_z_yy*Dgam_y_uyz
							
							+Dgam_x_xy*Dgam_y_uxx+Dgam_y_yy*Dgam_y_uyy+Dgam_z_zy*Dgam_y_uzz
							+Dgam_x_yy*Dgam_y_uxy+Dgam_x_zy*Dgam_y_uxz+Dgam_y_zy*Dgam_y_uyz
							+Dgam_y_xy*Dgam_y_uxy+Dgam_z_xy*Dgam_y_uxz+Dgam_z_yy*Dgam_y_uyz
							
							//-zgx_p*Dgam_x_yy-zgy_p*Dgam_y_yy-zgz_p*Dgam_z_yy)
							-gamma0_x*Dgam_x_yy-gamma0_y*Dgam_y_yy-gamma0_z*Dgam_z_yy)
						
						-Del_ux_yx*Del_ux_xy-Del_uy_yy*Del_uy_yy-Del_uz_yz*Del_uz_zy
						-Del_ux_yy*Del_uy_xy-Del_ux_yz*Del_uz_xy-Del_uy_yz*Del_uz_yy
						-Del_uy_yx*Del_ux_yy-Del_uz_yx*Del_ux_zy-Del_uz_yy*Del_uy_zy
						
						-2.*wa_cdyy -2.*wa_lap*gyy_p -4.*wawa*gyy_p +4.*wa_y_p*wa_y_p)*ewa4i;

				rc_zz_p=(0.5*(-gDDg_zz
							+gzx_p*DGam_z_ux+gzy_p*DGam_z_uy+gzz_p*DGam_z_uz
							+gzx_p*DGam_z_ux+gzy_p*DGam_z_uy+gzz_p*DGam_z_uz)
						-0.5*(Dgam_x_xz*Dgam_z_uxx+Dgam_y_yz*Dgam_z_uyy+Dgam_z_zz*Dgam_z_uzz
							+Dgam_x_yz*Dgam_z_uxy+Dgam_x_zz*Dgam_z_uxz+Dgam_y_zz*Dgam_z_uyz
							+Dgam_y_xz*Dgam_z_uxy+Dgam_z_xz*Dgam_z_uxz+Dgam_z_yz*Dgam_z_uyz
							
							+Dgam_x_xz*Dgam_z_uxx+Dgam_y_yz*Dgam_z_uyy+Dgam_z_zz*Dgam_z_uzz
							+Dgam_x_yz*Dgam_z_uxy+Dgam_x_zz*Dgam_z_uxz+Dgam_y_zz*Dgam_z_uyz
							+Dgam_y_xz*Dgam_z_uxy+Dgam_z_xz*Dgam_z_uxz+Dgam_z_yz*Dgam_z_uyz
							
							//-zgx_p*Dgam_x_zz-zgy_p*Dgam_y_zz-zgz_p*Dgam_z_zz)
							-gamma0_x*Dgam_x_zz-gamma0_y*Dgam_y_zz-gamma0_z*Dgam_z_zz)
						
						-Del_ux_zx*Del_ux_xz-Del_uy_zy*Del_uy_yz-Del_uz_zz*Del_uz_zz
						-Del_ux_zy*Del_uy_xz-Del_ux_zz*Del_uz_xz-Del_uy_zz*Del_uz_yz
						-Del_uy_zx*Del_ux_yz-Del_uz_zx*Del_ux_zz-Del_uz_zy*Del_uy_zz
						
						-2.*wa_cdzz -2.*wa_lap*gzz_p -4.*wawa*gzz_p +4.*wa_z_p*wa_z_p)*ewa4i;

				rc_xy_p=(0.5*(-gDDg_xy
							+gxx_p*DGam_y_ux+gxy_p*DGam_y_uy+gxz_p*DGam_y_uz
							+gyx_p*DGam_x_ux+gyy_p*DGam_x_uy+gyz_p*DGam_x_uz)
						-0.5*(Dgam_x_xx*Dgam_y_uxx+Dgam_y_yx*Dgam_y_uyy+Dgam_z_zx*Dgam_y_uzz
							+Dgam_x_yx*Dgam_y_uxy+Dgam_x_zx*Dgam_y_uxz+Dgam_y_zx*Dgam_y_uyz
							+Dgam_y_xx*Dgam_y_uxy+Dgam_z_xx*Dgam_y_uxz+Dgam_z_yx*Dgam_y_uyz
							
							+Dgam_x_xy*Dgam_x_uxx+Dgam_y_yy*Dgam_x_uyy+Dgam_z_zy*Dgam_x_uzz
							+Dgam_x_yy*Dgam_x_uxy+Dgam_x_zy*Dgam_x_uxz+Dgam_y_zy*Dgam_x_uyz
							+Dgam_y_xy*Dgam_x_uxy+Dgam_z_xy*Dgam_x_uxz+Dgam_z_yy*Dgam_x_uyz
							
							//-zgx_p*Dgam_x_xy-zgy_p*Dgam_y_xy-zgz_p*Dgam_z_xy)
							-gamma0_x*Dgam_x_xy-gamma0_y*Dgam_y_xy-gamma0_z*Dgam_z_xy)
						
						-Del_ux_xx*Del_ux_xy-Del_uy_xy*Del_uy_yy-Del_uz_xz*Del_uz_zy
						-Del_ux_xy*Del_uy_xy-Del_ux_xz*Del_uz_xy-Del_uy_xz*Del_uz_yy
						-Del_uy_xx*Del_ux_yy-Del_uz_xx*Del_ux_zy-Del_uz_xy*Del_uy_zy
						
						-2.*wa_cdxy -2.*wa_lap*gxy_p -4.*wawa*gxy_p +4.*wa_x_p*wa_y_p)*ewa4i;

				rc_xz_p=(0.5*(-gDDg_xz
							+gxx_p*DGam_z_ux+gxy_p*DGam_z_uy+gxz_p*DGam_z_uz
							+gzx_p*DGam_x_ux+gzy_p*DGam_x_uy+gzz_p*DGam_x_uz)
						-0.5*(Dgam_x_xx*Dgam_z_uxx+Dgam_y_yx*Dgam_z_uyy+Dgam_z_zx*Dgam_z_uzz
							+Dgam_x_yx*Dgam_z_uxy+Dgam_x_zx*Dgam_z_uxz+Dgam_y_zx*Dgam_z_uyz
							+Dgam_y_xx*Dgam_z_uxy+Dgam_z_xx*Dgam_z_uxz+Dgam_z_yx*Dgam_z_uyz
							
							+Dgam_x_xz*Dgam_x_uxx+Dgam_y_yz*Dgam_x_uyy+Dgam_z_zz*Dgam_x_uzz
							+Dgam_x_yz*Dgam_x_uxy+Dgam_x_zz*Dgam_x_uxz+Dgam_y_zz*Dgam_x_uyz
							+Dgam_y_xz*Dgam_x_uxy+Dgam_z_xz*Dgam_x_uxz+Dgam_z_yz*Dgam_x_uyz
							
							//-zgx_p*Dgam_x_xz-zgy_p*Dgam_y_xz-zgz_p*Dgam_z_xz)
							-gamma0_x*Dgam_x_xz-gamma0_y*Dgam_y_xz-gamma0_z*Dgam_z_xz)
						
						-Del_ux_xx*Del_ux_xz-Del_uy_xy*Del_uy_yz-Del_uz_xz*Del_uz_zz
						-Del_ux_xy*Del_uy_xz-Del_ux_xz*Del_uz_xz-Del_uy_xz*Del_uz_yz
						-Del_uy_xx*Del_ux_yz-Del_uz_xx*Del_ux_zz-Del_uz_xy*Del_uy_zz
						
						-2.*wa_cdxz -2.*wa_lap*gxz_p -4.*wawa*gxz_p +4.*wa_x_p*wa_z_p)*ewa4i;

				rc_yz_p=(0.5*(-gDDg_yz
							+gyx_p*DGam_z_ux+gyy_p*DGam_z_uy+gyz_p*DGam_z_uz
							+gzx_p*DGam_y_ux+gzy_p*DGam_y_uy+gzz_p*DGam_y_uz)
						-0.5*(Dgam_x_xy*Dgam_z_uxx+Dgam_y_yy*Dgam_z_uyy+Dgam_z_zy*Dgam_z_uzz
							+Dgam_x_yy*Dgam_z_uxy+Dgam_x_zy*Dgam_z_uxz+Dgam_y_zy*Dgam_z_uyz
							+Dgam_y_xy*Dgam_z_uxy+Dgam_z_xy*Dgam_z_uxz+Dgam_z_yy*Dgam_z_uyz
							
							+Dgam_x_xz*Dgam_y_uxx+Dgam_y_yz*Dgam_y_uyy+Dgam_z_zz*Dgam_y_uzz
							+Dgam_x_yz*Dgam_y_uxy+Dgam_x_zz*Dgam_y_uxz+Dgam_y_zz*Dgam_y_uyz
							+Dgam_y_xz*Dgam_y_uxy+Dgam_z_xz*Dgam_y_uxz+Dgam_z_yz*Dgam_y_uyz
							
							//-zgx_p*Dgam_x_yz-zgy_p*Dgam_y_yz-zgz_p*Dgam_z_yz)
							-gamma0_x*Dgam_x_yz-gamma0_y*Dgam_y_yz-gamma0_z*Dgam_z_yz)
						
						-Del_ux_yx*Del_ux_xz-Del_uy_yy*Del_uy_yz-Del_uz_yz*Del_uz_zz
						-Del_ux_yy*Del_uy_xz-Del_ux_yz*Del_uz_xz-Del_uy_yz*Del_uz_yz
						-Del_uy_yx*Del_ux_yz-Del_uz_yx*Del_ux_zz-Del_uz_yy*Del_uy_zz
						
						-2.*wa_cdyz -2.*wa_lap*gyz_p -4.*wawa*gyz_p +4.*wa_y_p*wa_z_p)*ewa4i;
				
				rc_yx_p=rc_xy_p;
				rc_zx_p=rc_xz_p;
				rc_zy_p=rc_yz_p;
				
				
				// R 
				ricci=gixx_p*rc_xx_p +gixy_p*rc_xy_p +gixz_p*rc_xz_p 
						+giyx_p*rc_yx_p +giyy_p*rc_yy_p +giyz_p*rc_yz_p 
						+gizx_p*rc_zx_p +gizy_p*rc_zy_p +gizz_p*rc_zz_p;
				

				//\tilde A_ij,k
				daxx_x=get_f_x(l,k,j,14);
				dayy_x=get_f_x(l,k,j,15);
				dazz_x=get_f_x(l,k,j,16);
				daxy_x=get_f_x(l,k,j,17);
				daxz_x=get_f_x(l,k,j,18);
				dayz_x=get_f_x(l,k,j,19);
				
				daxx_y=get_f_y(l,k,j,14);
				dayy_y=get_f_y(l,k,j,15);
				dazz_y=get_f_y(l,k,j,16);
				daxy_y=get_f_y(l,k,j,17);
				daxz_y=get_f_y(l,k,j,18);
				dayz_y=get_f_y(l,k,j,19);
				
				daxx_z=get_f_z(l,k,j,14);
				dayy_z=get_f_z(l,k,j,15);
				dazz_z=get_f_z(l,k,j,16);
				daxy_z=get_f_z(l,k,j,17);
				daxz_z=get_f_z(l,k,j,18);
				dayz_z=get_f_z(l,k,j,19);
				
				//D_3 \tilde A_12=Da12_3
				Daxx_x=daxx_x-2.*Gam_ux_xx*akxx_p;
				Dayy_x=dayy_x;
				Dazz_x=dazz_x;
				Daxy_x=daxy_x-Gam_ux_xx*akxy_p;
				Daxz_x=daxz_x-Gam_ux_xx*akxz_p;
				Dayz_x=dayz_x;
				
				Daxx_y=daxx_y;
				Dayy_y=dayy_y-2.*Gam_uy_yy*akyy_p;
				Dazz_y=dazz_y;
				Daxy_y=daxy_y-Gam_uy_yy*akxy_p;
				Daxz_y=daxz_y;
				Dayz_y=dayz_y-Gam_uy_yy*akyz_p;
				
				Daxx_z=daxx_z;
				Dayy_z=dayy_z;
				Dazz_z=dazz_z-2.*Gam_uz_zz*akzz_p;
				Daxy_z=daxy_z;
				Daxz_z=daxz_z-Gam_uz_zz*akxz_p;
				Dayz_z=dayz_z-Gam_uz_zz*akyz_p;
				
				//\tilde D^1 \tilde A_x1=Da_x
				Da_x=gixx_p*Daxx_x+gixy_p*Daxx_y+gixz_p*Daxx_z
					+giyx_p*Daxy_x+giyy_p*Daxy_y+giyz_p*Daxy_z
					+gizx_p*Daxz_x+gizy_p*Daxz_y+gizz_p*Daxz_z
					
					-Del_ux_xx*akx_ux_p-Del_ux_yx*akx_uy_p-Del_ux_zx*akx_uz_p
					-Del_uy_xx*aky_ux_p-Del_uy_yx*aky_uy_p-Del_uy_zx*aky_uz_p
					-Del_uz_xx*akz_ux_p-Del_uz_yx*akz_uy_p-Del_uz_zx*akz_uz_p
					
					-gamma0_x*akxx_p-gamma0_y*akxy_p-gamma0_z*akxz_p;
					
				
				Da_y=gixx_p*Daxy_x+gixy_p*Daxy_y+gixz_p*Daxy_z
					+giyx_p*Dayy_x+giyy_p*Dayy_y+giyz_p*Dayy_z
					+gizx_p*Dayz_x+gizy_p*Dayz_y+gizz_p*Dayz_z
					
					-Del_ux_xy*akx_ux_p-Del_ux_yy*akx_uy_p-Del_ux_zy*akx_uz_p
					-Del_uy_xy*aky_ux_p-Del_uy_yy*aky_uy_p-Del_uy_zy*aky_uz_p
					-Del_uz_xy*akz_ux_p-Del_uz_yy*akz_uy_p-Del_uz_zy*akz_uz_p
					
					-gamma0_x*akyx_p-gamma0_y*akyy_p-gamma0_z*akyz_p;
					
				Da_z=gixx_p*Daxz_x+gixy_p*Daxz_y+gixz_p*Daxz_z
					+giyx_p*Dayz_x+giyy_p*Dayz_y+giyz_p*Dayz_z
					+gizx_p*Dazz_x+gizy_p*Dazz_y+gizz_p*Dazz_z
				
					-Del_ux_xz*akx_ux_p-Del_ux_yz*akx_uy_p-Del_ux_zz*akx_uz_p
					-Del_uy_xz*aky_ux_p-Del_uy_yz*aky_uy_p-Del_uy_zz*aky_uz_p
					-Del_uz_xz*akz_ux_p-Del_uz_yz*akz_uy_p-Del_uz_zz*akz_uz_p
					
					-gamma0_x*akzx_p-gamma0_y*akzy_p-gamma0_z*akzz_p;
				
				if(scalarevo)
				{
					//scalar field and conjugate momentum, those derivatives
					double phi,Pi,phi_x,phi_y,phi_z,dphidphi;
					
					//scalar field potential and derivative
					double Vpot;
					
					//substitution start
					phi=get_bv(l,k,j,nsc);
					Pi=get_bv(l,k,j,nscp);
					
					phi_x=get_f_x(l,k,j,nsc);
					phi_y=get_f_y(l,k,j,nsc);
					phi_z=get_f_z(l,k,j,nsc);

					// (\tilde D_i \phi)(\tilde D^i \phi)
					dphidphi=gixx_p*phi_x*phi_x +gixy_p*phi_x*phi_y +gixz_p*phi_x*phi_z
						+giyx_p*phi_y*phi_x +giyy_p*phi_y*phi_y +giyz_p*phi_y*phi_z
						+gizx_p*phi_z*phi_x +gizy_p*phi_z*phi_y +gizz_p*phi_z*phi_z;

					Vpot=funcV(phi);
					
					
					//////////////////////////////////////////////////////
					// scalar field contribution to stress energy tensor
					//////////////////////////////////////////////////////
					
					Ene=0.5*(pow(Pi,2)+ewa4i*dphidphi)+Vpot;
					
					p_x=Pi*phi_x;
					p_y=Pi*phi_y;
					p_z=Pi*phi_z;
					
				}
				
				//M_x=Da_x+6.*(akx_ux_p*wa_x_p+akx_uy_p*wa_y_p+akx_uz_p*wa_z_p)-2./3.*ek_x_p-pi8*j_x;
				//M_y=Da_y+6.*(aky_ux_p*wa_x_p+aky_uy_p*wa_y_p+aky_uz_p*wa_z_p)-2./3.*ek_y_p-pi8*j_y;
				//M_z=Da_z+6.*(akz_ux_p*wa_x_p+akz_uy_p*wa_y_p+akz_uz_p*wa_z_p)-2./3.*ek_z_p-pi8*j_z;

				M_x=(Da_x+6.*(akx_ux_p*wa_x_p+akx_uy_p*wa_y_p+akx_uz_p*wa_z_p)-2./3.*ek_x_p-pi8*p_x)/pi8;
				M_y=(Da_y+6.*(aky_ux_p*wa_x_p+aky_uy_p*wa_y_p+aky_uz_p*wa_z_p)-2./3.*ek_y_p-pi8*p_y)/pi8;
				M_z=(Da_z+6.*(akz_ux_p*wa_x_p+akz_uy_p*wa_y_p+akz_uz_p*wa_z_p)-2./3.*ek_z_p-pi8*p_z)/pi8;
				
				sqgam=sqrt(det_p)*exp(6.*wa_p);
				
				set_bv(l,k,j,25)=M_x*sqgam;
				set_bv(l,k,j,26)=M_y*sqgam;
				set_bv(l,k,j,27)=M_z*sqgam;

				
				//  \tilde{A}_i^k=\tilde{A}_{ij}\tilde{\gamma}^{jk} -> ak{i}_u{k}
				akx_ux_p=akxx_p*gixx_p +akxy_p*giyx_p +akxz_p*gizx_p;
				aky_ux_p=akyx_p*gixx_p +akyy_p*giyx_p +akyz_p*gizx_p;
				akz_ux_p=akzx_p*gixx_p +akzy_p*giyx_p +akzz_p*gizx_p;
				akx_uy_p=akxx_p*gixy_p +akxy_p*giyy_p +akxz_p*gizy_p;
				aky_uy_p=akyx_p*gixy_p +akyy_p*giyy_p +akyz_p*gizy_p;
				akz_uy_p=akzx_p*gixy_p +akzy_p*giyy_p +akzz_p*gizy_p;
				akx_uz_p=akxx_p*gixz_p +akxy_p*giyz_p +akxz_p*gizz_p;
				aky_uz_p=akyx_p*gixz_p +akyy_p*giyz_p +akyz_p*gizz_p;
				akz_uz_p=akzx_p*gixz_p +akzy_p*giyz_p +akzz_p*gizz_p;

				//\tilde A_i^j \tilde A_j^i
				aaaa=akx_ux_p*akx_ux_p 
					+akx_uy_p*aky_ux_p 
					+akx_uz_p*akz_ux_p 
					+aky_ux_p*akx_uy_p  
					+aky_uy_p*aky_uy_p 
					+aky_uz_p*akz_uy_p 
					+akz_ux_p*akx_uz_p  
					+akz_uy_p*aky_uz_p 
					+akz_uz_p*akz_uz_p;
				
				//hamc = ricci + 2.*pow(ek_p,2)/3. - aaaa - 2.*lambda - pi16*get_enemom(l,k,j,0);
				//hamc = ricci - aaaa - 2.*lambda - pi16*get_enemom(l,k,j,0);
				Ene=(ricci +2.*pow(ek_p,2)/3.- aaaa - 2.*lambda- pi16*Ene)/pi16;
				set_bv(l,k,j,24) = Ene*sqgam;
				
				M_ux=(gixx_p*M_x+gixy_p*M_y+gixz_p*M_z)*ewa4i;
				M_uy=(giyx_p*M_x+giyy_p*M_y+giyz_p*M_z)*ewa4i;
				M_uz=(gizx_p*M_x+gizy_p*M_y+gizz_p*M_z)*ewa4i;
				
				p2=M_x*M_ux+M_y*M_uy+M_z*M_uz;
				
				get_rhoGam(Ene,p2,rho,Gam);

				set_primv(l,k,j,0)=rho;
				
				EpP=Ene+pres(rho);
				
				set_primv(l,k,j,1)=get_bv(l,k,j,0)*M_ux/EpP-get_bv(l,k,j,1);
				set_primv(l,k,j,2)=get_bv(l,k,j,0)*M_uy/EpP-get_bv(l,k,j,2);
				set_primv(l,k,j,3)=get_bv(l,k,j,0)*M_uz/EpP-get_bv(l,k,j,3);
				
				double rhob=pow(Hb,2)*3./pi8;
				double eps=pow(rho/rhob,fluidw/(1.+fluidw))-1.;
				
				set_primv(l,k,j,4)=eps;
				
				set_bv(l,k,j,28) = rho*Gam*sqgam/(1.+eps);
				
			}
		}
	}
	
	return;
}

void Fmv::set_initial_scalar(double mu,double kk,double xi2,double xi3)
{
	double L=1.;								// box size
	double inr=0.8;								// =0 for r>inr
	
	cout << "kks=" << kk << endl;				// wave number kk=\sqrt(\sum (k^i)^2) 
	cout << "mus=" << mu << endl;				// amplitude mus 
												// see arXiv:2004.01042 for the profile
	#pragma omp parallel for
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
				
				double r=sqrt(X*X+Y*Y+Z*Z);
				
				//see arXiv:2004.01042 for the functional form
				if(r>L)
				 set_bv(l,k,j,nsc)=0.;
				else if(r<inr)
				 set_bv(l,k,j,nsc)=pow(M_E,(pow(M_E,(-(pow(X,2)*(pow(kk,2) + 3*xi2 + xi3)) - 
				          (pow(kk,2) - 2*xi3)*pow(Y,2) - (pow(kk,2) - 3*xi2 + xi3)*pow(Z,2))/6.)*mu)/2.    )-1.;
				else 
				 set_bv(l,k,j,nsc)=pow(M_E,(pow(M_E,(-(pow(X,2)*(pow(kk,2) + 3*xi2 + xi3)) - 
				          (pow(kk,2) - 2*xi3)*pow(Y,2) - (pow(kk,2) - 3*xi2 + xi3)*pow(Z,2))/6.)*mu*
					      (1 - pow(pow(inr - L,6) - pow(L - r,6),6)/pow(inr - L,36)))/2.)-1.;
				set_bv(l,k,j,nscp)=0.;
			}
		}
	}
	#pragma omp barrier
	
	boundary_reflection_scalar();
	
	return;
}

void Fmv::set_initial_fluid(double mu,double kk,double xi2,double xi3)
{
	double L=1.;								// box size
	double inr=0.8;								// =0 for r>inr
	
	cout << "kk=" << kk << endl;				// wave number kk=\sqrt(\sum (k^i)^2) 
												// see arXiv:2004.01042 for the profile
	#pragma omp parallel for
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
				
				double r=sqrt(X*X+Y*Y+Z*Z);
				
				//see arXiv:2004.01042 for the functional form
				if(r>L)
				set_bv(l,k,j,24)=1.;
				else if(r<inr)
				set_bv(l,k,j,24)=pow(M_E,(pow(M_E,(-(pow(X,2)*(pow(kk,2) + 3*xi2 + xi3)) - 
				          (pow(kk,2) - 2*xi3)*pow(Y,2) - (pow(kk,2) - 3*xi2 + xi3)*pow(Z,2))/6.)*mu)/2.    );
				else 
				set_bv(l,k,j,24)=pow(M_E,(pow(M_E,(-(pow(X,2)*(pow(kk,2) + 3*xi2 + xi3)) - 
				          (pow(kk,2) - 2*xi3)*pow(Y,2) - (pow(kk,2) - 3*xi2 + xi3)*pow(Z,2))/6.)*mu*
					      (1 - pow(pow(inr - L,6) - pow(L - r,6),6)/pow(inr - L,36)))/2.);
				
				set_bv(l,k,j,25)=0.;
				set_bv(l,k,j,26)=0.;
				set_bv(l,k,j,27)=0.;
				set_bv(l,k,j,28)=get_bv(l,k,j,24);
			}
		}
	}
	#pragma omp barrier
	
	boundary_reflection_fluid();
	
	return;
}

