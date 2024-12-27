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
	double L=1.;						// box size
	double inr=0.8;						// Psi=1 for r>inr
	double kk=10.;
	
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
				
				//see arXiv:2112.12335 for the functional form
				//we use psi for temporary storage for phi
				if(r>L)
				set_psi(l,k,j)=0. ;
				else if(r<inr)
				set_psi(l,k,j)=mu*pow(M_E,-(pow(kk,2)*pow(r,2))/6.);
				else 
				set_psi(l,k,j)=mu*pow(M_E,-(pow(kk,2)*pow(r,2))/6.)*(1 - pow(inr - L,-36)*pow(pow(inr - L,6) - pow(L - r,6),6));
			}
		}
	}
	#pragma omp barrier
	
	boundary_psi_initial();

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
				
				double Phiv=get_psi(l,k,j);
				
				double Phi_x=(   -(get_psi(l,k,j+2)-get_psi(l,k,j-2))
						+8.*(get_psi(l,k,j+1)-get_psi(l,k,j-1))
						)*dxi12;
				double Phi_y=(   -(get_psi(l,k+2,j)-get_psi(l,k-2,j))
						+8.*(get_psi(l,k+1,j)-get_psi(l,k-1,j))
						)*dyi12;
				double Phi_z=(   -(get_psi(l+2,k,j)-get_psi(l-2,k,j))
						+8.*(get_psi(l+1,k,j)-get_psi(l-1,k,j))
						)*dzi12;
				
				double Phi_xx=(    -(     get_psi(l  ,k  ,j+2) + get_psi(l  ,k  ,j-2))
						+16.*( get_psi(l  ,k  ,j+1) + get_psi(l  ,k  ,j-1))
						-30.*  get_psi(l  ,k  ,j  )
						)*dxi*dxi12;
				double Phi_yy=(    -(     get_psi(l  ,k+2,j  ) + get_psi(l  ,k-2,j  ))
						+16.*( get_psi(l  ,k+1,j  ) + get_psi(l  ,k-1,j  ))
						-30.*  get_psi(l  ,k  ,j  )
						)*dyi*dyi12;
				double Phi_zz=(    -(     get_psi(l+2,k  ,j  ) + get_psi(l-2,k  ,j  ))
						+16.*( get_psi(l+1,k  ,j  ) + get_psi(l-1,k  ,j  ))
						-30.*  get_psi(l  ,k  ,j  )
						)*dzi*dzi12;
				
				double DDPhi_xy=( -( ( -(   get_psi(l  ,k+2,j+2) - get_psi(l  ,k+2,j-2))
						+8.*(get_psi(l  ,k+2,j+1) - get_psi(l  ,k+2,j-1)))
						-(   -(get_psi(l  ,k-2,j+2) - get_psi(l  ,k-2,j-2))
						+8.*(get_psi(l  ,k-2,j+1) - get_psi(l  ,k-2,j-1))))
						+8.*(  (-(get_psi(l  ,k+1,j+2) - get_psi(l  ,k+1,j-2))
						+8.*(get_psi(l  ,k+1,j+1) - get_psi(l  ,k+1,j-1)))
						-(-(get_psi(l  ,k-1,j+2) - get_psi(l  ,k-1,j-2))
						+8.*(get_psi(l  ,k-1,j+1) - get_psi(l  ,k-1,j-1))))
						)*dxi12*dyi12;
				double DDPhi_xz=( -( ( -(   get_psi(l+2,k  ,j+2) - get_psi(l+2,k  ,j-2))
						+8.*(get_psi(l+2,k  ,j+1) - get_psi(l+2,k  ,j-1)))
						-(   -(get_psi(l-2,k  ,j+2) - get_psi(l-2,k  ,j-2))
						+8.*(get_psi(l-2,k  ,j+1) - get_psi(l-2,k  ,j-1))))
						+8.*(  (-(get_psi(l+1,k  ,j+2) - get_psi(l+1,k  ,j-2))
						+8.*(get_psi(l+1,k  ,j+1) - get_psi(l+1,k  ,j-1)))
						-(-(get_psi(l-1,k  ,j+2) - get_psi(l-1,k  ,j-2))
						+8.*(get_psi(l-1,k  ,j+1) - get_psi(l-1,k  ,j-1))))
						)*dxi12*dzi12;
				double DDPhi_yz=( -( ( -(   get_psi(l+2,k+2,j  ) - get_psi(l+2,k-2,j  ))
						+8.*(get_psi(l+2,k+1,j  ) - get_psi(l+2,k-1,j  )))
						-(   -(get_psi(l-2,k+2,j  ) - get_psi(l-2,k-2,j  ))
						+8.*(get_psi(l-2,k+1,j  ) - get_psi(l-2,k-1,j  ))))
						+8.*(  (-(get_psi(l+1,k+2,j  ) - get_psi(l+1,k-2,j  ))
						+8.*(get_psi(l+1,k+1,j  ) - get_psi(l+1,k-1,j  )))
						-(-(get_psi(l-1,k+2,j  ) - get_psi(l-1,k-2,j  ))
						+8.*(get_psi(l-1,k+1,j  ) - get_psi(l-1,k-1,j  ))))
						)*dyi12*dzi12;
						
				double DDPhi_xx=Phi_xx-get_flat_Gamx(j)*Phi_x;
				double DDPhi_yy=Phi_yy-get_flat_Gamy(k)*Phi_y;
				double DDPhi_zz=Phi_zz-get_flat_Gamz(l)*Phi_z;

				double LapPhi=(Phi_xx-Phi_x*ddfx/dfx)*pow(dfx,-2)
								+(Phi_yy-Phi_y*ddfy/dfy)*pow(dfy,-2)
								+(Phi_zz-Phi_z*ddfz/dfz)*pow(dfz,-2);
				
				double DPDP=pow(Phi_x/dfx,2)+pow(Phi_y/dfy,2)+pow(Phi_z/dfz,2);
				
								
				double pxx=Phi_x*Phi_x-fxx*DPDP/3.;
				double pyy=Phi_y*Phi_y-fyy*DPDP/3.;
				double pzz=Phi_z*Phi_z-fzz*DPDP/3.;
				
				double pxy=Phi_x*Phi_y;
				double pxz=Phi_x*Phi_z;
				double pyz=Phi_y*Phi_z;

				fac=32.*M_PI/((5.+3.*fluidw)*(1.*3.*fluidw));
				
				double hxx=fac*pxx*Hbi2;
				double hyy=fac*pyy*Hbi2;
				double hzz=fac*pzz*Hbi2;
				double hxy=fac*pxy*Hbi2;
				double hxz=fac*pxz*Hbi2;
				double hyz=fac*pyz*Hbi2;
				
				fac=16.*M_PI/(5.+3.*fluidw);
				
				double Axx=fac*pxx*Hbi2*HH;
				double Ayy=fac*pyy*Hbi2*HH;
				double Azz=fac*pzz*Hbi2*HH;
				double Axy=fac*pxy*Hbi2*HH;
				double Axz=fac*pxz*Hbi2*HH;
				double Ayz=fac*pyz*Hbi2*HH;
							
				double rhob=3./(8.*M_PI)*pow(HH,2);
				
				//CMC
				double kap=0.;
				double delta=-4.*M_PI/3.*DPDP*Hbi2;
				double chi=-(1.+3.*fluidw)/(3.*(1.+fluidw))*delta;
				double xi=2.*M_PI/(9.*(1.+fluidw))*DPDP*Hbi2;

				double u_x=-16.*M_PI/(9.*fGam*(3.*fGam+2.))*Hbi3*aa
				*(pow(dfx,-2)*Phi_x*DDPhi_xx+pow(dfy,-2)*Phi_y*DDPhi_xy+pow(dfz,-2)*Phi_z*DDPhi_xz);
				
				double u_y=-16.*M_PI/(9.*fGam*(3.*fGam+2.))*Hbi3*aa
				*(pow(dfx,-2)*Phi_x*DDPhi_xy+pow(dfy,-2)*Phi_y*DDPhi_yy+pow(dfz,-2)*Phi_z*DDPhi_yz);
				
				double u_z=-16.*M_PI/(9.*fGam*(3.*fGam+2.))*Hbi3*aa
				*(pow(dfx,-2)*Phi_x*DDPhi_xz+pow(dfy,-2)*Phi_y*DDPhi_yz+pow(dfz,-2)*Phi_z*DDPhi_zz);
				
				double vpi=-2./(3.*fGam+2.)*Hbi2*HH*LapPhi;
				double lam=2./((3.*fGam+2.)*(3.*fGam-2.))*Hbi2*LapPhi;
				
				//comoving
				//double kap=-1./(3.*fGam+2)*f*Hbi2;
				//double chi=-3.*(fGam-1.)/(3.*fGam+2.)*f*Hbi2;
				//double xi=-1./(2.*(3.*fGam+2))*f*Hbi2;
				//double delta=3.*fGam/(3.*fGam+2)*f*Hbi2;
				//double u_x=0.;
				//double u_y=0.;
				//double u_z=0.;
				
				double psii=(1.+xi);

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

				set_bv(l,k,j,nsc)=Phiv+lam;
				set_bv(l,k,j,nscp)=vpi;
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

