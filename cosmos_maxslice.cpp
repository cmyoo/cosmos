/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* MAXIMAL SLICE EQUATION SOLVER :: BSSN evolution Class of COSMOS                                       */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*                             NOTE boundary conditions are restricted to the reflection boundary cond.  */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

//setting coefficients for maximal slice elliptic equation
void Fmv::coefficients()
{
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				////////// variables /////////////////////////////////////////////////////
				double gxx_p,gyy_p,gzz_p,												//
				gxy_p,gxz_p,gyz_p,
				gyx_p,gzx_p,gzy_p,wa_p,
				akxx_p,akyy_p,akzz_p,
				akxy_p,akxz_p,akyz_p,
				akyx_p,akzx_p,akzy_p,ek_p;
				
				double fxx,fyy,fzz,dfxx,dfyy,dfzz;
				
				double Gam_ux_xx,Gam_uy_yy,Gam_uz_zz;
				
				double dxx_x_p,dyy_x_p,dzz_x_p,
				dxy_x_p,dxz_x_p,dyz_x_p,
				dyx_x_p,dzx_x_p,dzy_x_p,wa_x_p;
				
				double dxx_y_p,dyy_y_p,dzz_y_p,
				dxy_y_p,dxz_y_p,dyz_y_p,
				dyx_y_p,dzx_y_p,dzy_y_p,wa_y_p;
				
				double dxx_z_p,dyy_z_p,dzz_z_p,
				dxy_z_p,dxz_z_p,dyz_z_p,
				dyx_z_p,dzx_z_p,dzy_z_p,wa_z_p;


				double det_p,det_pi,ewa4i;
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

				double aaaa;

				double trs=0.;
				
				double Ene=0.;															//
				////////// variables /////////////////////////////////////////////////////

				//fij for inhomogeneous grid 
				fxx=get_flat_df2x(j);
				fyy=get_flat_df2y(k);
				fzz=get_flat_df2z(l);

				////////// metric components /////////////////////////////////////////////////////
				gxx_p=  get_bv(l,k,j, 7)+fxx;													//
				gyy_p=  get_bv(l,k,j, 8)+fyy;
				gzz_p=  get_bv(l,k,j, 9)+fzz;
				gxy_p=  get_bv(l,k,j,10);
				gxz_p=  get_bv(l,k,j,11);
				gyz_p=  get_bv(l,k,j,12);
				gyx_p=	gxy_p;
				gzx_p=	gxz_p;
				gzy_p=	gyz_p;
				wa_p=   get_bv(l,k,j,13);
				ewa4i=exp(-4.*wa_p);															//
				////////// metric components /////////////////////////////////////////////////////

				////////// extrinsic curvature ///////////////////////////////////////////////////
				akxx_p= get_bv(l,k,j,14);														//
				akyy_p= get_bv(l,k,j,15);
				akzz_p= get_bv(l,k,j,16);
				akxy_p= get_bv(l,k,j,17);
				akxz_p= get_bv(l,k,j,18);
				akyz_p= get_bv(l,k,j,19);
				akyx_p=	akxy_p;
				akzx_p=	akxz_p;
				akzy_p=	akyz_p;
				ek_p=   get_bv(l,k,j,20);														//
				////////// extrinsic curvature ///////////////////////////////////////////////////
				
				////////// metric derivatives  ///////////////////////////////////////////////////
				//bar Gamma_ui_jk for inhomogeneous grid 										//
				Gam_ux_xx=get_flat_Gamx(j);
				Gam_uy_yy=get_flat_Gamy(k);
				Gam_uz_zz=get_flat_Gamz(l);
				
				// reference metric derivatives
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
				wa_z_p=get_f_z(l,k,j,13);														//
				////////// metric derivatives  ///////////////////////////////////////////////////

				////////// metric determinant and inverse metric /////////////////////////////////
				// tilde{gamma}^{ij} -> gi														//
				det_p=fxx*fyy*fzz;
				det_pi=1./det_p;
				gixx_p= (gyy_p*gzz_p-gyz_p*gzy_p)*det_pi;
				giyy_p= (gxx_p*gzz_p-gxz_p*gzx_p)*det_pi;
				gizz_p= (gxx_p*gyy_p-gxy_p*gyx_p)*det_pi;
				gixy_p=-(gyx_p*gzz_p-gyz_p*gzx_p)*det_pi;
				gixz_p= (gyx_p*gzy_p-gyy_p*gzx_p)*det_pi;
				giyz_p=-(gxx_p*gzy_p-gxy_p*gzx_p)*det_pi;										//
				////////// metric determinant and inverse metric /////////////////////////////////
				
				////////// coefficient of derivative terms ///////////////////////////////////////
				set_coef(l,k,j,0)=gixx_p*ewa4i;													// 
				set_coef(l,k,j,1)=giyy_p*ewa4i;
				set_coef(l,k,j,2)=gizz_p*ewa4i;
				set_coef(l,k,j,3)=gixy_p*ewa4i;
				set_coef(l,k,j,4)=gixz_p*ewa4i;
				set_coef(l,k,j,5)=giyz_p*ewa4i;													//
				////////// coefficient of derivative terms ///////////////////////////////////////
				
				////////// cals of coefficient of first derivative terms /////////////////////////
				giyx_p=gixy_p;																	//
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
				
				//maxslice mod
				//for maxslice 2
				set_coef(l,k,j,6)=(gixx_p*crx_xx+giyy_p*crx_yy+gizz_p*crx_zz
						+2.*(gixy_p*crx_xy+gixz_p*crx_xz+giyz_p*crx_yz)
						-2.*(gixx_p*wa_x_p+gixy_p*wa_y_p+gixz_p*wa_z_p))*ewa4i;
				set_coef(l,k,j,7)=(gixx_p*cry_xx+giyy_p*cry_yy+gizz_p*cry_zz
						+2.*(gixy_p*cry_xy+gixz_p*cry_xz+giyz_p*cry_yz)
						-2.*(giyx_p*wa_x_p+giyy_p*wa_y_p+giyz_p*wa_z_p))*ewa4i;
				set_coef(l,k,j,8)=(gixx_p*crz_xx+giyy_p*crz_yy+gizz_p*crz_zz
						+2.*(gixy_p*crz_xy+gixz_p*crz_xz+giyz_p*crz_yz)
						-2.*(gizx_p*wa_x_p+gizy_p*wa_y_p+gizz_p*wa_z_p))*ewa4i;					//
				////////// cals of coefficient of first derivative terms /////////////////////////
				
				////////// cals of source term ///////////////////////////////////////////////////
				////////////// for scalar field ///////////////////								//
				if(scalarevo)
				{
					//substitution start
					double phi=get_bv(l,k,j,nsc);
					double Pi=get_bv(l,k,j,nscp);
					
					double phi_x=get_f_x(l,k,j,nsc);
					double phi_y=get_f_y(l,k,j,nsc);
					double phi_z=get_f_z(l,k,j,nsc);

					// (\tilde D_i \phi)(\tilde D^i \phi)
					double dphidphi=gixx_p*	phi_x*phi_x +gixy_p*phi_x*phi_y +gixz_p*phi_x*phi_z
						+giyx_p*phi_y*phi_x +giyy_p*phi_y*phi_y +giyz_p*phi_y*phi_z
						+gizx_p*phi_z*phi_x +gizy_p*phi_z*phi_y +gizz_p*phi_z*phi_z;
					
					double Vpot=funcV(phi);
				
					Ene=0.5*(pow(Pi,2)+ewa4i*dphidphi)+Vpot;
					
					double sxx=phi_x*phi_x+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxx_p;
					double syy=phi_y*phi_y+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gyy_p;
					double szz=phi_z*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gzz_p;
					double sxy=phi_x*phi_y+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxy_p;
					double sxz=phi_x*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxz_p;
					double syz=phi_y*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gyz_p;
					
					trs=(gixx_p*sxx +gixy_p*sxy +gixz_p*sxz
						+giyx_p*sxy +giyy_p*syy +giyz_p*syz
						+gizx_p*sxz +gizy_p*syz +gizz_p*szz)*ewa4i;
				}
				
				//for fluid
				if(fluidevo)
				{
					double sqgam=sqrt(det_p)*exp(6.*wa_p);
					
					double fEne=get_bv(l,k,j,24)/sqgam;
					double fp_x=get_bv(l,k,j,25)/sqgam;
					double fp_y=get_bv(l,k,j,26)/sqgam;
					double fp_z=get_bv(l,k,j,27)/sqgam;
					
					double press=pres(get_primv(l,k,j,0));
					
					double EpP=fEne+press;
					
					press/=ewa4i;
					double fsxx=fp_x*fp_x/EpP+press*gxx_p;
					double fsyy=fp_y*fp_y/EpP+press*gyy_p;
					double fszz=fp_z*fp_z/EpP+press*gzz_p;
					double fsxy=fp_x*fp_y/EpP+press*gxy_p;
					double fsxz=fp_x*fp_z/EpP+press*gxz_p;
					double fsyz=fp_y*fp_z/EpP+press*gyz_p;

					double ftrs=(gixx_p*fsxx +gixy_p*fsxy +gixz_p*fsxz
						+giyx_p*fsxy +giyy_p*fsyy +giyz_p*fsyz
						+gizx_p*fsxz +gizy_p*fsyz +gizz_p*fszz)*ewa4i;
					
					Ene+=fEne;
					trs+=ftrs;
				}
				
				
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
				
				//source term 
				set_coef(l,k,j,9)=(aaaa + pow(ek_p,2)/3.+pi4*(Ene+trs));						//
				////////// cals of source term ///////////////////////////////////////////////////

				//initial data for the iteration 
				set_bv(l,k,j,0)=get_bv0(l,k,j,0);
				
			}
		}
	}
	
	return;
}

//iteration with Jacobi method for maximal slice cond.
//boundary condition is restricted to the reflection boundary cond.
void Fmv::maxslice(int msstep,double convcond,double acc)
{
	//for convergence cond.
	double errm=0.;						// max error
	double norm=1.;						// norm
	double normt=1.0e-6;				// temporary norm for parallelization
	bool conv=false;					// convergence bool
	
	dtk=0.;								// value of dtrK/dt
	
	//setting coefficients of elliptic eq.
	coefficients();

	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);
	
	//iteration steps
	for(int n=0;n<msstep;n++)
	{
		errm=0.;
		
		double maxalpha=1.0e-3;

		//reflection boundary for lapse(bv(i=0))
		boundary_reflection_even(0);
		
		//alpha=bv 
		set_alpha();
		
		double sum=0.;
		double vol=0.;
		
		#pragma omp parallel for default(shared) reduction(max:errm,normt)
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int j=jli;j<=jui;j++)
				{
					double dev=0.4*(
								get_coef(l,k,j,0)*
								(    -(    get_alpha(l  ,k  ,j+2) + get_alpha(l  ,k  ,j-2))
								+16.*( get_alpha(l  ,k  ,j+1) + get_alpha(l  ,k  ,j-1))
								//-30.*  get_alpha(l  ,k  ,j  )
								)/12.
						
								+get_coef(l,k,j,1)*
								(    -(    get_alpha(l  ,k+2,j  ) + get_alpha(l  ,k-2,j  ))
								+16.*( get_alpha(l  ,k+1,j  ) + get_alpha(l  ,k-1,j  ))
								//-30.*  get_alpha(l  ,k  ,j  )
								)/12.
								
								+get_coef(l,k,j,2)*
								(    -(    get_alpha(l+2,k  ,j  ) + get_alpha(l-2,k  ,j  ))
								+16.*( get_alpha(l+1,k  ,j  ) + get_alpha(l-1,k  ,j  ))
								//-30.*  get_alpha(l  ,k  ,j  )
								)/12.
								+2.0*(
								
								get_coef(l,k,j,3)*
								( -( ( -(      get_alpha(l  ,k+2,j+2) - get_alpha(l  ,k+2,j-2))
								+8.*(get_alpha(l  ,k+2,j+1) - get_alpha(l  ,k+2,j-1)))
								-(   -(get_alpha(l  ,k-2,j+2) - get_alpha(l  ,k-2,j-2))
								+8.*(get_alpha(l  ,k-2,j+1) - get_alpha(l  ,k-2,j-1))))
								+8.*(  (-(get_alpha(l  ,k+1,j+2) - get_alpha(l  ,k+1,j-2))
								+8.*(get_alpha(l  ,k+1,j+1) - get_alpha(l  ,k+1,j-1)))
								-(-(get_alpha(l  ,k-1,j+2) - get_alpha(l  ,k-1,j-2))
								+8.*(get_alpha(l  ,k-1,j+1) - get_alpha(l  ,k-1,j-1))))
								)/144.
								
								+get_coef(l,k,j,4)*
								( -( ( -(      get_alpha(l+2,k  ,j+2) - get_alpha(l+2,k  ,j-2))
								+8.*(get_alpha(l+2,k  ,j+1) - get_alpha(l+2,k  ,j-1)))
								-(   -(get_alpha(l-2,k  ,j+2) - get_alpha(l-2,k  ,j-2))
								+8.*(get_alpha(l-2,k  ,j+1) - get_alpha(l-2,k  ,j-1))))
								+8.*(  (-(get_alpha(l+1,k  ,j+2) - get_alpha(l+1,k  ,j-2))
								+8.*(get_alpha(l+1,k  ,j+1) - get_alpha(l+1,k  ,j-1)))
								-(-(get_alpha(l-1,k  ,j+2) - get_alpha(l-1,k  ,j-2))
								+8.*(get_alpha(l-1,k  ,j+1) - get_alpha(l-1,k  ,j-1))))
								)/144.

								+get_coef(l,k,j,5)*
								( -( ( -(      get_alpha(l+2,k+2,j  ) - get_alpha(l+2,k-2,j  ))
								+8.*(get_alpha(l+2,k+1,j  ) - get_alpha(l+2,k-1,j  )))
								-(   -(get_alpha(l-2,k+2,j  ) - get_alpha(l-2,k-2,j  ))
								+8.*(get_alpha(l-2,k+1,j  ) - get_alpha(l-2,k-1,j  ))))
								+8.*(  (-(get_alpha(l+1,k+2,j  ) - get_alpha(l+1,k-2,j  ))
								+8.*(get_alpha(l+1,k+1,j  ) - get_alpha(l+1,k-1,j  )))
								-(-(get_alpha(l-1,k+2,j  ) - get_alpha(l-1,k-2,j  ))
								+8.*(get_alpha(l-1,k+1,j  ) - get_alpha(l-1,k-1,j  ))))
								)/144.
								
								)
								-dx*(
								get_coef(l,k,j,6)*
								(    -(get_alpha(l,k,j+2)-get_alpha(l,k,j-2))
								+8.*(get_alpha(l,k,j+1)-get_alpha(l,k,j-1))
								)/12.
								
								+get_coef(l,k,j,7)*
								(   -(get_alpha(l,k+2,j)-get_alpha(l,k-2,j))
								+8.*(get_alpha(l,k+1,j)-get_alpha(l,k-1,j))
								)/12.

								+get_coef(l,k,j,8)*
								(   -(get_alpha(l+2,k,j)-get_alpha(l-2,k,j))
								+8.*(get_alpha(l+1,k,j)-get_alpha(l-1,k,j))
								)/12.
								
								)
								
								-dx*dx*(get_alpha(l,k,j)*get_coef(l,k,j,9)-dtk)
									)/(get_coef(l,k,j,0)+get_coef(l,k,j,1)+get_coef(l,k,j,2))
									-get_alpha(l,k,j);
								
					double absv=abs(get_alpha(l,k,j))+abs(get_coef(l,k,j,9)*dx*dx);
					
					double err=abs(dev)/norm;
					
					if(err>1.0e6)
					{
						cout << "large error in maxslice" << endl;
						cout << "j=" << j << " "
						<< "k=" << k << " " << "l=" << l << " "
						<< "dev=" << dev << " " << "norm=" << norm << endl;
						
						cout << " get_alpha(l,k,j)=" << get_alpha(l,k,j) 
							<< " get_alpha(l,k,j+1)=" << get_alpha(l,k,j+1) 
							<< " get_alpha(l,k,j-1)=" << get_alpha(l,k,j-1) 
							<< " get_alpha(l,k+1,j)=" << get_alpha(l,k+1,j) 
							<< " get_alpha(l,k-1,j)=" << get_alpha(l,k-1,j) 
							<< " get_alpha(l+1,k,j)=" << get_alpha(l+1,k,j) 
							<< " get_alpha(l-1,k,j)=" << get_alpha(l-1,k,j) 
							<< " get_alpha(l,k-1,j-1)=" << get_alpha(l,k-1,j-1) 
							<< " get_alpha(l,k+1,j+1)=" << get_alpha(l,k+1,j+1) 
							<< " get_coef(l,k,j,0)=" << get_coef(l,k,j,0) 
							<< " get_coef(l,k,j,1)=" << get_coef(l,k,j,1) 
							<< " get_coef(l,k,j,2)=" << get_coef(l,k,j,2) 
							<< " get_coef(l,k,j,3)=" << get_coef(l,k,j,3) 
							<< " get_coef(l,k,j,4)=" << get_coef(l,k,j,4) 
							<< " get_coef(l,k,j,5)=" << get_coef(l,k,j,5) 
							<< " get_coef(l,k,j,6)=" << get_coef(l,k,j,6) 
							<< " get_coef(l,k,j,7)=" << get_coef(l,k,j,7) 
							<< " get_coef(l,k,j,8)=" << get_coef(l,k,j,8) 
							<< " get_coef(l,k,j,9)=" << get_coef(l,k,j,9) 
							<< " dtk=" << dtk 
							<< endl;
					
						exit(1);
					}
					
					set_bv(l,k,j,0)=get_bv(l,k,j,0)+acc*dev;
					
					if (errm<err)
					errm=err;
					
					if (absv > normt)
					normt=absv;
				}
			}
		}

		#pragma omp barrier
		
		norm=normt;
		normt=1.0e-10;
		
		//cal of dtrk/dt as integrability cond.
		#pragma omp parallel for default(shared) reduction(+:sum,vol) reduction(max:maxalpha)
		for(int l=lli;l<=lui;l++)
		{
			double facl=1.;
			if(l==lui || l==lli)
			facl=0.5;
			
			for(int k=kli;k<=kui;k++)
			{
				double fack=1.;
				if(k==kui || k==kli)
				fack=0.5;

				for(int j=jli;j<=jui;j++)
				{
					double facj=1.;
					if(j==jui || j==jli)
					facj=0.5;
				
					if(maxalpha< get_bv(l,k,j,0))
					maxalpha=get_bv(l,k,j,0);
					
					double fac=facl*fack*facj*sqrt(get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l));
					sum+=(get_bv(l,k,j,0)*get_coef(l,k,j,9))*fac*exp(get_bv(l,k,j,13)*6.);
					vol+=fac*exp(get_bv(l,k,j,13)*6.);
				}
			}
		}
		
		dtk=sum/(vol*maxalpha);
		
		#pragma omp parallel for default(shared)
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int j=jli;j<=jui;j++)
				{
					set_bv(l,k,j,0)=get_bv(l,k,j,0)/maxalpha;
				}
			}
		}
		
		cout << "maxslice step=" << n
			<< " errm =" << errm 
			<< endl;
		
		if(errm<convcond)
		{
			cout << "maxslice step=" << n 
			<< "errm maxslice=" << errm << endl;

			conv=true;
			boundary_reflection_even(0);
			break;
			
		}
	}
	
	if(!conv)
	{
		cout << "not converge in maxslice step in multi" 
			<< " errm=" << errm << endl;

		boundary_reflection_even(0);
	}
	
}

//iteration with SOR for maximal slice cond. no parallelization
void Fmv::maxslice_sor(int msstep,double convcond,double acc)
{
	//for convergence cond.
	double errm=0.;						// max error
	double norm=1.;						// norm
	double normt=1.0e-6;				// temporary norm for parallelization
	bool conv=false;					// convergence bool
	
	dtk=0.;								// value of dtrK/dt
	
	//setting coefficients of elliptic eq.
	coefficients();
	
	//iteration steps
	for(int n=0;n<msstep;n++)
	{
		errm=0.;
		double maxalpha=1.0e-3;
		double sum=0.;
		double vol=0.;

		boundary_reflection_even(0);
		
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int j=jli;j<=jui;j++)
				{
					double dev=0.4*(
								get_coef(l,k,j,0)*
								(    -(    get_bv(l  ,k  ,j+2,0) + get_bv(l  ,k  ,j-2,0))
								+16.*( get_bv(l  ,k  ,j+1,0) + get_bv(l  ,k  ,j-1,0))
								//-30.*  get_bv(l  ,k  ,j  ,0)
								)/12.
						
								+get_coef(l,k,j,1)*
								(    -(    get_bv(l  ,k+2,j  ,0) + get_bv(l  ,k-2,j  ,0))
								+16.*( get_bv(l  ,k+1,j  ,0) + get_bv(l  ,k-1,j  ,0))
								//-30.*  get_bv(l  ,k  ,j  ,0)
								)/12.
								
								+get_coef(l,k,j,2)*
								(    -(    get_bv(l+2,k  ,j  ,0) + get_bv(l-2,k  ,j  ,0))
								+16.*( get_bv(l+1,k  ,j  ,0) + get_bv(l-1,k  ,j  ,0))
								//-30.*  get_bv(l  ,k  ,j  ,0)
								)/12.
								+2.0*(
								
								get_coef(l,k,j,3)*
								( -( ( -(      get_bv(l  ,k+2,j+2,0) - get_bv(l  ,k+2,j-2,0))
								+8.*(get_bv(l  ,k+2,j+1,0) - get_bv(l  ,k+2,j-1,0)))
								-(   -(get_bv(l  ,k-2,j+2,0) - get_bv(l  ,k-2,j-2,0))
								+8.*(get_bv(l  ,k-2,j+1,0) - get_bv(l  ,k-2,j-1,0))))
								+8.*(  (-(get_bv(l  ,k+1,j+2,0) - get_bv(l  ,k+1,j-2,0))
								+8.*(get_bv(l  ,k+1,j+1,0) - get_bv(l  ,k+1,j-1,0)))
								-(-(get_bv(l  ,k-1,j+2,0) - get_bv(l  ,k-1,j-2,0))
								+8.*(get_bv(l  ,k-1,j+1,0) - get_bv(l  ,k-1,j-1,0))))
								)/144.
								
								+get_coef(l,k,j,4)*
								( -( ( -(      get_bv(l+2,k  ,j+2,0) - get_bv(l+2,k  ,j-2,0))
								+8.*(get_bv(l+2,k  ,j+1,0) - get_bv(l+2,k  ,j-1,0)))
								-(   -(get_bv(l-2,k  ,j+2,0) - get_bv(l-2,k  ,j-2,0))
								+8.*(get_bv(l-2,k  ,j+1,0) - get_bv(l-2,k  ,j-1,0))))
								+8.*(  (-(get_bv(l+1,k  ,j+2,0) - get_bv(l+1,k  ,j-2,0))
								+8.*(get_bv(l+1,k  ,j+1,0) - get_bv(l+1,k  ,j-1,0)))
								-(-(get_bv(l-1,k  ,j+2,0) - get_bv(l-1,k  ,j-2,0))
								+8.*(get_bv(l-1,k  ,j+1,0) - get_bv(l-1,k  ,j-1,0))))
								)/144.

								+get_coef(l,k,j,5)*
								( -( ( -(      get_bv(l+2,k+2,j  ,0) - get_bv(l+2,k-2,j  ,0))
								+8.*(get_bv(l+2,k+1,j  ,0) - get_bv(l+2,k-1,j  ,0)))
								-(   -(get_bv(l-2,k+2,j  ,0) - get_bv(l-2,k-2,j  ,0))
								+8.*(get_bv(l-2,k+1,j  ,0) - get_bv(l-2,k-1,j  ,0))))
								+8.*(  (-(get_bv(l+1,k+2,j  ,0) - get_bv(l+1,k-2,j  ,0))
								+8.*(get_bv(l+1,k+1,j  ,0) - get_bv(l+1,k-1,j  ,0)))
								-(-(get_bv(l-1,k+2,j  ,0) - get_bv(l-1,k-2,j  ,0))
								+8.*(get_bv(l-1,k+1,j  ,0) - get_bv(l-1,k-1,j  ,0))))
								)/144.
								
								)
								-dx*(
								get_coef(l,k,j,6)*
								(    -(get_bv(l,k,j+2,0)-get_bv(l,k,j-2,0))
								+8.*(get_bv(l,k,j+1,0)-get_bv(l,k,j-1,0))
								)/12.
								
								+get_coef(l,k,j,7)*
								(   -(get_bv(l,k+2,j,0)-get_bv(l,k-2,j,0))
								+8.*(get_bv(l,k+1,j,0)-get_bv(l,k-1,j,0))
								)/12.

								+get_coef(l,k,j,8)*
								(   -(get_bv(l+2,k,j,0)-get_bv(l-2,k,j,0))
								+8.*(get_bv(l+1,k,j,0)-get_bv(l-1,k,j,0))
								)/12.
								
								)
								
								-dx*dx*(get_bv(l,k,j,0)*get_coef(l,k,j,9)-dtk)
									)/(get_coef(l,k,j,0)+get_coef(l,k,j,1)+get_coef(l,k,j,2))
									-get_bv(l,k,j,0);

								
					double absv=abs(get_bv(l,k,j,0))+abs(get_coef(l,k,j,9)*dx*dx);
					
					double err=abs(dev)/norm;
					
					if(err>1.0e6)
					{
						cout << "large error in maxslice" << endl;
						cout << "j=" << j << " "
						<< "k=" << k << " " << "l=" << l << " "
						<< "dev=" << dev << " " << "norm=" << norm << endl;
						
						cout << " get_bv(l,k,j,0)=" << get_bv(l,k,j,0) 
							<< " get_bv(l,k,j+1,0)=" << get_bv(l,k,j+1,0) 
							<< " get_bv(l,k,j-1,0)=" << get_bv(l,k,j-1,0) 
							<< " get_bv(l,k+1,j,0)=" << get_bv(l,k+1,j,0) 
							<< " get_bv(l,k-1,j,0)=" << get_bv(l,k-1,j,0) 
							<< " get_bv(l+1,k,j,0)=" << get_bv(l+1,k,j,0) 
							<< " get_bv(l-1,k,j,0)=" << get_bv(l-1,k,j,0) 
							<< " get_bv(l,k-1,j-1,0)=" << get_bv(l,k-1,j-1,0) 
							<< " get_bv(l,k+1,j+1,0)=" << get_bv(l,k+1,j+1,0) 
							<< " get_coef(l,k,j,0)=" << get_coef(l,k,j,0) 
							<< " get_coef(l,k,j,1)=" << get_coef(l,k,j,1) 
							<< " get_coef(l,k,j,2)=" << get_coef(l,k,j,2) 
							<< " get_coef(l,k,j,3)=" << get_coef(l,k,j,3) 
							<< " get_coef(l,k,j,4)=" << get_coef(l,k,j,4) 
							<< " get_coef(l,k,j,5)=" << get_coef(l,k,j,5) 
							<< " get_coef(l,k,j,6)=" << get_coef(l,k,j,6) 
							<< " get_coef(l,k,j,7)=" << get_coef(l,k,j,7) 
							<< " get_coef(l,k,j,8)=" << get_coef(l,k,j,8) 
							<< " get_coef(l,k,j,9)=" << get_coef(l,k,j,9) 
							<< " dtk=" << dtk 
							<< endl;
					
						exit(1);
					}
					
					set_bv(l,k,j,0)=get_bv(l,k,j,0)+acc*dev;
					
					if (errm<err)
					errm=err;
					
					if (absv > normt)
					normt=absv;
				}
			}
		}

		norm=normt;
		normt=1.0e-10;
		
		//cal of dtrk/dt as integrability cond.
		#pragma omp parallel for default(shared) reduction(+:sum,vol) reduction(max:maxalpha)
		for(int l=lli;l<=lui;l++)
		{
			double facl=1.;
			if(l==lui || l==lli)
			facl=0.5;
			
			for(int k=kli;k<=kui;k++)
			{
				double fack=1.;
				if(k==kui || k==kli)
				fack=0.5;

				for(int j=jli;j<=jui;j++)
				{
					double facj=1.;
					if(j==jui || j==jli)
					facj=0.5;
				
					if(maxalpha< get_bv(l,k,j,0))
					maxalpha=get_bv(l,k,j,0);
					
					double fac=facl*fack*facj*sqrt(get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l));
					sum+=(get_bv(l,k,j,0)*get_coef(l,k,j,9))*fac*exp(get_bv(l,k,j,13)*6.);
					vol+=fac*exp(get_bv(l,k,j,13)*6.);
				}
			}
		}
		
		dtk=sum/(vol*maxalpha);
				
		cout << "maxslice step=" << n
			<< " errm =" << errm 
			<< endl;

		#pragma omp parallel for default(shared) 
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int j=jli;j<=jui;j++)
				{
					set_bv(l,k,j,0)=get_bv(l,k,j,0)/maxalpha;
				}
			}
		}
		
		if(errm<convcond)
		{
			cout << "maxslice step=" << n 
			<< "errm maxslice=" << errm << endl;

			conv=true;
			boundary_reflection_even(0);
			break;
			
		}
		
	
	}
	
	if(!conv)
	{
		cout << "not converge in maxslice step sor" 
			 << " errm=" << errm << endl;

		boundary_reflection_even(0);
	}
	
}



