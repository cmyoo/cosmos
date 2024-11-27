/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* CONSTRAINTS SOLVER :: BSSN evolution Class of COSMOS                                                  */
/*-------------------------------------------------------------------------------------------------------*/
/* for Hamiltonian constraint, the conformal factor is the variable to be solved                         */
/*                                         see (7.94) in "3+1 Formalism in General Relativity"           */
/*                                                                by Eric Gourgoulhon                    */
/*                                       Eq.(6.111) in the arXiv version gr-qc/0703035                   */
/*-------------------------------------------------------------------------------------------------------*/
/* for Momentum constraints, vector w defined by \tilde A^ij =(\tilde L w)^ij + \tilde A^ij_TT           */
/*     the venctor w^i are the variables to be solved                                                    */
/*     for the definition of \tilde L                                                                    */
/*                see (9.11)[(8.11) in arXiv ver.] in "3+1 Formalism in General Relativity"              */
/*     a similar decomposition can be seen in (9.9)[(8.9)] but for \hat A^ij                             */
/*     \tilde A_ij will be updated in every iteration step => w^i is expected to be zero finally         */
/*-------------------------------------------------------------------------------------------------------*/
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/* NOTE: boundary conditions are restricted to the reflection boundary cond.                             */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

void Fmv::coefficients_solveham()
{
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				////////// variables /////////////////////////////////////////////////////////////////////
				double gxx_p,gyy_p,gzz_p,																//
				gxy_p,gxz_p,gyz_p,
				gyx_p,gzx_p,gzy_p,
				akxx_p,akyy_p,akzz_p,
				akxy_p,akxz_p,akyz_p,
				akyx_p,akzx_p,akzy_p;

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

				double aaaa;																			//
				////////// variables /////////////////////////////////////////////////////////////////////

				gxx_p=  get_bv(l,k,j, 7)+get_flat_df2x(j);
				gyy_p=  get_bv(l,k,j, 8)+get_flat_df2y(k);
				gzz_p=  get_bv(l,k,j, 9)+get_flat_df2z(l);
				gxy_p=  get_bv(l,k,j,10);
				gxz_p=  get_bv(l,k,j,11);
				gyz_p=  get_bv(l,k,j,12);
				gyx_p=	gxy_p;
				gzx_p=	gxz_p;
				gzy_p=	gyz_p;
				
				akxx_p= get_bv(l,k,j,14);
				akyy_p= get_bv(l,k,j,15);
				akzz_p= get_bv(l,k,j,16);
				akxy_p= get_bv(l,k,j,17);
				akxz_p= get_bv(l,k,j,18);
				akyz_p= get_bv(l,k,j,19);
				akyx_p=	akxy_p;
				akzx_p=	akxz_p;
				akzy_p=	akyz_p;
				
				double dfxxh=get_flat_df2x(j)*get_flat_Gamx(j);
				double dfyyh=get_flat_df2y(k)*get_flat_Gamy(k);
				double dfzzh=get_flat_df2z(l)*get_flat_Gamz(l);

				// \del_x \tilde{gamma}_{ij} *0.5 
				dxx_x_p=(   -(get_bv(l,k,j+2, 7)-get_bv(l,k,j-2, 7))
						+8.*(get_bv(l,k,j+1, 7)-get_bv(l,k,j-1, 7))
						)*dxi24+dfxxh;
				dyy_x_p=(   -(get_bv(l,k,j+2, 8)-get_bv(l,k,j-2, 8))
						+8.*(get_bv(l,k,j+1, 8)-get_bv(l,k,j-1, 8))
						)*dxi24;
				dzz_x_p=(   -(get_bv(l,k,j+2, 9)-get_bv(l,k,j-2, 9))
						+8.*(get_bv(l,k,j+1, 9)-get_bv(l,k,j-1, 9))
						)*dxi24;
				dxy_x_p=(   -(get_bv(l,k,j+2,10)-get_bv(l,k,j-2,10))
						+8.*(get_bv(l,k,j+1,10)-get_bv(l,k,j-1,10))
						)*dxi24;
				dxz_x_p=(   -(get_bv(l,k,j+2,11)-get_bv(l,k,j-2,11))
						+8.*(get_bv(l,k,j+1,11)-get_bv(l,k,j-1,11))
						)*dxi24;
				dyz_x_p=(   -(get_bv(l,k,j+2,12)-get_bv(l,k,j-2,12))
						+8.*(get_bv(l,k,j+1,12)-get_bv(l,k,j-1,12))
						)*dxi24;
				dyx_x_p=dxy_x_p;
				dzx_x_p=dxz_x_p;
				dzy_x_p=dyz_x_p;
				
				// \del_y \tilde \gamma_{ij} * 0.5 
				dxx_y_p=(   -(get_bv(l,k+2,j, 7)-get_bv(l,k-2,j, 7))
						+8.*(get_bv(l,k+1,j, 7)-get_bv(l,k-1,j, 7))
						)*dyi24;
				dyy_y_p=(   -(get_bv(l,k+2,j, 8)-get_bv(l,k-2,j, 8))
						+8.*(get_bv(l,k+1,j, 8)-get_bv(l,k-1,j, 8))
						)*dyi24+dfyyh;
				dzz_y_p=(   -(get_bv(l,k+2,j, 9)-get_bv(l,k-2,j, 9))
						+8.*(get_bv(l,k+1,j, 9)-get_bv(l,k-1,j, 9))
						)*dyi24;
				dxy_y_p=(   -(get_bv(l,k+2,j,10)-get_bv(l,k-2,j,10))
						+8.*(get_bv(l,k+1,j,10)-get_bv(l,k-1,j,10))
						)*dyi24;
				dxz_y_p=(   -(get_bv(l,k+2,j,11)-get_bv(l,k-2,j,11))
						+8.*(get_bv(l,k+1,j,11)-get_bv(l,k-1,j,11))
						)*dyi24;
				dyz_y_p=(   -(get_bv(l,k+2,j,12)-get_bv(l,k-2,j,12))
						+8.*(get_bv(l,k+1,j,12)-get_bv(l,k-1,j,12))
						)*dyi24;
				dyx_y_p=dxy_y_p;
				dzx_y_p=dxz_y_p;
				dzy_y_p=dyz_y_p;

				// \del_z \tilde \gamma_{ij} * 0.5
				dxx_z_p=(   -(get_bv(l+2,k,j, 7)-get_bv(l-2,k,j, 7))
						+8.*(get_bv(l+1,k,j, 7)-get_bv(l-1,k,j, 7))
						)*dzi24;
				dyy_z_p=(-   (get_bv(l+2,k,j, 8)-get_bv(l-2,k,j, 8))
						+8.*(get_bv(l+1,k,j, 8)-get_bv(l-1,k,j, 8))
						)*dzi24;
				dzz_z_p=(   -(get_bv(l+2,k,j, 9)-get_bv(l-2,k,j, 9))
						+8.*(get_bv(l+1,k,j, 9)-get_bv(l-1,k,j, 9))
						)*dzi24+dfzzh;
				dxy_z_p=(   -(get_bv(l+2,k,j,10)-get_bv(l-2,k,j,10))
						+8.*(get_bv(l+1,k,j,10)-get_bv(l-1,k,j,10))
						)*dzi24;
				dxz_z_p=(   -(get_bv(l+2,k,j,11)-get_bv(l-2,k,j,11))
						+8.*(get_bv(l+1,k,j,11)-get_bv(l-1,k,j,11))
						)*dzi24;
				dyz_z_p=(   -(get_bv(l+2,k,j,12)-get_bv(l-2,k,j,12))
						+8.*(get_bv(l+1,k,j,12)-get_bv(l-1,k,j,12))
						)*dzi24;
				dyx_z_p=dxy_z_p;
				dzx_z_p=dxz_z_p;
				dzy_z_p=dyz_z_p;
				
				// tilde{gamma}^{ij} -> gi
				det_p=gxx_p*gyy_p*gzz_p+gxy_p*gyz_p*gzx_p+gxz_p*gyx_p*gzy_p
						-gxz_p*gyy_p*gzx_p-gxy_p*gyx_p*gzz_p-gxx_p*gyz_p*gzy_p;
				if(det_p<1.e-16) 
				 det_p=1.e-16;
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
				
				//for psi
				set_coef(l,k,j,0)=gixx_p*crx_xx+giyy_p*crx_yy+gizz_p*crx_zz
									+2.*(gixy_p*crx_xy+gixz_p*crx_xz+giyz_p*crx_yz);
				set_coef(l,k,j,1)=gixx_p*cry_xx+giyy_p*cry_yy+gizz_p*cry_zz
									+2.*(gixy_p*cry_xy+gixz_p*cry_xz+giyz_p*cry_yz);
				set_coef(l,k,j,2)=gixx_p*crz_xx+giyy_p*crz_yy+gizz_p*crz_zz
									+2.*(gixy_p*crz_xy+gixz_p*crz_xz+giyz_p*crz_yz);
				
				//cal D_i tilde gamma_jk
				double Dgam_x_xx,Dgam_x_yy,Dgam_x_zz,Dgam_x_xy,Dgam_x_xz,Dgam_x_yz,Dgam_x_yx,Dgam_x_zx,Dgam_x_zy,
						Dgam_y_xx,Dgam_y_yy,Dgam_y_zz,Dgam_y_xy,Dgam_y_xz,Dgam_y_yz,Dgam_y_yx,Dgam_y_zx,Dgam_y_zy,
						Dgam_z_xx,Dgam_z_yy,Dgam_z_zz,Dgam_z_xy,Dgam_z_xz,Dgam_z_yz,Dgam_z_yx,Dgam_z_zx,Dgam_z_zy;

				//bar Gamma_ui_jk for inhomogeneous grid 
				double Gam_ux_xx,Gam_uy_yy,Gam_uz_zz;
				
				Gam_ux_xx=get_flat_Gamx(j);
				Gam_uy_yy=get_flat_Gamy(k);
				Gam_uz_zz=get_flat_Gamz(l);

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
				double Del_ux_xx,Del_ux_yy,Del_ux_zz,Del_ux_xy,Del_ux_xz,Del_ux_yz,Del_ux_yx,Del_ux_zx,Del_ux_zy,
						Del_uy_xx,Del_uy_yy,Del_uy_zz,Del_uy_xy,Del_uy_xz,Del_uy_yz,Del_uy_yx,Del_uy_zx,Del_uy_zy,
						Del_uz_xx,Del_uz_yy,Del_uz_zz,Del_uz_xy,Del_uz_xz,Del_uz_yz,Del_uz_yx,Del_uz_zx,Del_uz_zy;
				
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
				
				//cal D_i tilde gamma^jk
				double Dgam_x_uxx,Dgam_x_uyy,Dgam_x_uzz,Dgam_x_uxy,Dgam_x_uxz,Dgam_x_uyz,Dgam_x_uyx,Dgam_x_uzx,Dgam_x_uzy,
						Dgam_y_uxx,Dgam_y_uyy,Dgam_y_uzz,Dgam_y_uxy,Dgam_y_uxz,Dgam_y_uyz,Dgam_y_uyx,Dgam_y_uzx,Dgam_y_uzy,
						Dgam_z_uxx,Dgam_z_uyy,Dgam_z_uzz,Dgam_z_uxy,Dgam_z_uxz,Dgam_z_uyz,Dgam_z_uyx,Dgam_z_uzx,Dgam_z_uzy;
				
				//Dgam_1_u23=-Del_u2_1x*gix3_p-Del_u2_1y*giy3_p-Del_u2_1z*giz3_p-Del_u3_1x*gix2_p-Del_u3_1y*giy2_p-Del_u3_1z*giz2_p;
				Dgam_x_uxx=-Del_ux_xx*gixx_p-Del_ux_xy*giyx_p-Del_ux_xz*gizx_p-Del_ux_xx*gixx_p-Del_ux_xy*giyx_p-Del_ux_xz*gizx_p;
				Dgam_x_uyy=-Del_uy_xx*gixy_p-Del_uy_xy*giyy_p-Del_uy_xz*gizy_p-Del_uy_xx*gixy_p-Del_uy_xy*giyy_p-Del_uy_xz*gizy_p;
				Dgam_x_uzz=-Del_uz_xx*gixz_p-Del_uz_xy*giyz_p-Del_uz_xz*gizz_p-Del_uz_xx*gixz_p-Del_uz_xy*giyz_p-Del_uz_xz*gizz_p;
				Dgam_x_uxy=-Del_ux_xx*gixy_p-Del_ux_xy*giyy_p-Del_ux_xz*gizy_p-Del_uy_xx*gixx_p-Del_uy_xy*giyx_p-Del_uy_xz*gizx_p;
				Dgam_x_uxz=-Del_ux_xx*gixz_p-Del_ux_xy*giyz_p-Del_ux_xz*gizz_p-Del_uz_xx*gixx_p-Del_uz_xy*giyx_p-Del_uz_xz*gizx_p;
				Dgam_x_uyz=-Del_uy_xx*gixz_p-Del_uy_xy*giyz_p-Del_uy_xz*gizz_p-Del_uz_xx*gixy_p-Del_uz_xy*giyy_p-Del_uz_xz*gizy_p;
				Dgam_x_uyx=Dgam_x_uxy;
				Dgam_x_uzx=Dgam_x_uxz;
				Dgam_x_uzy=Dgam_x_uyz;
				
				Dgam_y_uxx=-Del_ux_yx*gixx_p-Del_ux_yy*giyx_p-Del_ux_yz*gizx_p-Del_ux_yx*gixx_p-Del_ux_yy*giyx_p-Del_ux_yz*gizx_p;
				Dgam_y_uyy=-Del_uy_yx*gixy_p-Del_uy_yy*giyy_p-Del_uy_yz*gizy_p-Del_uy_yx*gixy_p-Del_uy_yy*giyy_p-Del_uy_yz*gizy_p;
				Dgam_y_uzz=-Del_uz_yx*gixz_p-Del_uz_yy*giyz_p-Del_uz_yz*gizz_p-Del_uz_yx*gixz_p-Del_uz_yy*giyz_p-Del_uz_yz*gizz_p;
				Dgam_y_uxy=-Del_ux_yx*gixy_p-Del_ux_yy*giyy_p-Del_ux_yz*gizy_p-Del_uy_yx*gixx_p-Del_uy_yy*giyx_p-Del_uy_yz*gizx_p;
				Dgam_y_uxz=-Del_ux_yx*gixz_p-Del_ux_yy*giyz_p-Del_ux_yz*gizz_p-Del_uz_yx*gixx_p-Del_uz_yy*giyx_p-Del_uz_yz*gizx_p;
				Dgam_y_uyz=-Del_uy_yx*gixz_p-Del_uy_yy*giyz_p-Del_uy_yz*gizz_p-Del_uz_yx*gixy_p-Del_uz_yy*giyy_p-Del_uz_yz*gizy_p;
				Dgam_y_uyx=Dgam_y_uxy;
				Dgam_y_uzx=Dgam_y_uxz;
				Dgam_y_uzy=Dgam_y_uyz;
				
				Dgam_z_uxx=-Del_ux_zx*gixx_p-Del_ux_zy*giyx_p-Del_ux_zz*gizx_p-Del_ux_zx*gixx_p-Del_ux_zy*giyx_p-Del_ux_zz*gizx_p;
				Dgam_z_uyy=-Del_uy_zx*gixy_p-Del_uy_zy*giyy_p-Del_uy_zz*gizy_p-Del_uy_zx*gixy_p-Del_uy_zy*giyy_p-Del_uy_zz*gizy_p;
				Dgam_z_uzz=-Del_uz_zx*gixz_p-Del_uz_zy*giyz_p-Del_uz_zz*gizz_p-Del_uz_zx*gixz_p-Del_uz_zy*giyz_p-Del_uz_zz*gizz_p;
				Dgam_z_uxy=-Del_ux_zx*gixy_p-Del_ux_zy*giyy_p-Del_ux_zz*gizy_p-Del_uy_zx*gixx_p-Del_uy_zy*giyx_p-Del_uy_zz*gizx_p;
				Dgam_z_uxz=-Del_ux_zx*gixz_p-Del_ux_zy*giyz_p-Del_ux_zz*gizz_p-Del_uz_zx*gixx_p-Del_uz_zy*giyx_p-Del_uz_zz*gizx_p;
				Dgam_z_uyz=-Del_uy_zx*gixz_p-Del_uy_zy*giyz_p-Del_uy_zz*gizz_p-Del_uz_zx*gixy_p-Del_uz_zy*giyy_p-Del_uz_zz*gizy_p;
				Dgam_z_uyx=Dgam_z_uxy;
				Dgam_z_uzx=Dgam_z_uxz;
				Dgam_z_uzy=Dgam_z_uyz;
				
				double zgx_p,zgy_p,zgz_p,zgx_x_p,zgy_y_p,zgz_z_p;
				
				zgx_p=  get_bv(l,k,j,21);
				zgy_p=  get_bv(l,k,j,22);
				zgz_p=  get_bv(l,k,j,23);
				
				//\del_x \Gamma^i
				zgx_x_p=(   -(get_bv(l,k,j+2,21)-get_bv(l,k,j-2,21))
						+8.*(get_bv(l,k,j+1,21)-get_bv(l,k,j-1,21))
						)*dxi12;
				zgy_y_p=(  -(get_bv(l,k+2,j,22)-get_bv(l,k-2,j,22))
						+8.*(get_bv(l,k+1,j,22)-get_bv(l,k-1,j,22))
						)*dyi12;
				zgz_z_p=(  -(get_bv(l+2,k,j,23)-get_bv(l-2,k,j,23))
						+8.*(get_bv(l+1,k,j,23)-get_bv(l-1,k,j,23))
						)*dzi12;
				
				double DivGam=zgx_x_p+zgy_y_p+zgz_z_p+Gam_ux_xx*zgx_p+Gam_uy_yy*zgy_p+Gam_uz_zz*zgz_p;
				
				//Ricci
				double Ricci=DivGam
							+0.5*(gixx_p*(Dgam_x_uxx*Dgam_x_xx+Dgam_x_uyy*Dgam_x_yy+Dgam_x_uzz*Dgam_x_zz
									+2.*(Dgam_x_uxy*Dgam_x_xy+Dgam_x_uxz*Dgam_x_xz+Dgam_x_uyz*Dgam_x_yz))
								+giyy_p*(Dgam_y_uxx*Dgam_y_xx+Dgam_y_uyy*Dgam_y_yy+Dgam_y_uzz*Dgam_y_zz
									+2.*(Dgam_y_uxy*Dgam_y_xy+Dgam_y_uxz*Dgam_y_xz+Dgam_y_uyz*Dgam_y_yz))
								+gizz_p*(Dgam_z_uxx*Dgam_z_xx+Dgam_z_uyy*Dgam_z_yy+Dgam_z_uzz*Dgam_z_zz
									+2.*(Dgam_z_uxy*Dgam_z_xy+Dgam_z_uxz*Dgam_z_xz+Dgam_z_uyz*Dgam_z_yz))
								+gixy_p*(Dgam_x_uxx*Dgam_y_xx+Dgam_x_uyy*Dgam_y_yy+Dgam_x_uzz*Dgam_y_zz
									+2.*(Dgam_x_uxy*Dgam_y_xy+Dgam_x_uxz*Dgam_y_xz+Dgam_x_uyz*Dgam_y_yz))
								+gixz_p*(Dgam_x_uxx*Dgam_z_xx+Dgam_x_uyy*Dgam_z_yy+Dgam_x_uzz*Dgam_z_zz
									+2.*(Dgam_x_uxy*Dgam_z_xy+Dgam_x_uxz*Dgam_z_xz+Dgam_x_uyz*Dgam_z_yz))
								+giyz_p*(Dgam_y_uxx*Dgam_z_xx+Dgam_y_uyy*Dgam_z_yy+Dgam_y_uzz*Dgam_z_zz
									+2.*(Dgam_y_uxy*Dgam_z_xy+Dgam_y_uxz*Dgam_z_xz+Dgam_y_uyz*Dgam_z_yz))
								+gixy_p*(Dgam_y_uxx*Dgam_x_xx+Dgam_y_uyy*Dgam_x_yy+Dgam_y_uzz*Dgam_x_zz
									+2.*(Dgam_y_uxy*Dgam_x_xy+Dgam_y_uxz*Dgam_x_xz+Dgam_y_uyz*Dgam_x_yz))
								+gixz_p*(Dgam_z_uxx*Dgam_x_xx+Dgam_z_uyy*Dgam_x_yy+Dgam_z_uzz*Dgam_x_zz
									+2.*(Dgam_z_uxy*Dgam_x_xy+Dgam_z_uxz*Dgam_x_xz+Dgam_z_uyz*Dgam_x_yz))
								+giyz_p*(Dgam_z_uxx*Dgam_y_xx+Dgam_z_uyy*Dgam_y_yy+Dgam_z_uzz*Dgam_y_zz
									+2.*(Dgam_z_uxy*Dgam_y_xy+Dgam_z_uxz*Dgam_y_xz+Dgam_z_uyz*Dgam_y_yz)))
									
								-gixx_p*(Dgam_x_uxx*Dgam_x_xx+Dgam_x_uyy*Dgam_y_yx+Dgam_x_uzz*Dgam_z_zx
										+Dgam_x_uxy*Dgam_x_yx+Dgam_x_uxz*Dgam_x_zx+Dgam_x_uyz*Dgam_y_zx
										+Dgam_x_uyx*Dgam_y_xx+Dgam_x_uzx*Dgam_z_xx+Dgam_x_uzy*Dgam_z_yx)
								-giyy_p*(Dgam_y_uxx*Dgam_x_xy+Dgam_y_uyy*Dgam_y_yy+Dgam_y_uzz*Dgam_z_zy
										+Dgam_y_uxy*Dgam_x_yy+Dgam_y_uxz*Dgam_x_zy+Dgam_y_uyz*Dgam_y_zy
										+Dgam_y_uyx*Dgam_y_xy+Dgam_y_uzx*Dgam_z_xy+Dgam_y_uzy*Dgam_z_yy)
								-gizz_p*(Dgam_z_uxx*Dgam_x_xz+Dgam_z_uyy*Dgam_y_yz+Dgam_z_uzz*Dgam_z_zz
										+Dgam_z_uxy*Dgam_x_yz+Dgam_z_uxz*Dgam_x_zz+Dgam_z_uyz*Dgam_y_zz
										+Dgam_z_uyx*Dgam_y_xz+Dgam_z_uzx*Dgam_z_xz+Dgam_z_uzy*Dgam_z_yz)
								-gixy_p*(Dgam_x_uxx*Dgam_x_xy+Dgam_x_uyy*Dgam_y_yy+Dgam_x_uzz*Dgam_z_zy
										+Dgam_x_uxy*Dgam_x_yy+Dgam_x_uxz*Dgam_x_zy+Dgam_x_uyz*Dgam_y_zy
										+Dgam_x_uyx*Dgam_y_xy+Dgam_x_uzx*Dgam_z_xy+Dgam_x_uzy*Dgam_z_yy)
								-gixz_p*(Dgam_x_uxx*Dgam_x_xz+Dgam_x_uyy*Dgam_y_yz+Dgam_x_uzz*Dgam_z_zz
										+Dgam_x_uxy*Dgam_x_yz+Dgam_x_uxz*Dgam_x_zz+Dgam_x_uyz*Dgam_y_zz
										+Dgam_x_uyx*Dgam_y_xz+Dgam_x_uzx*Dgam_z_xz+Dgam_x_uzy*Dgam_z_yz)
								-giyz_p*(Dgam_y_uxx*Dgam_x_xz+Dgam_y_uyy*Dgam_y_yz+Dgam_y_uzz*Dgam_z_zz
										+Dgam_y_uxy*Dgam_x_yz+Dgam_y_uxz*Dgam_x_zz+Dgam_y_uyz*Dgam_y_zz
										+Dgam_y_uyx*Dgam_y_xz+Dgam_y_uzx*Dgam_z_xz+Dgam_y_uzy*Dgam_z_yz)
								-giyx_p*(Dgam_y_uxx*Dgam_x_xx+Dgam_y_uyy*Dgam_y_yx+Dgam_y_uzz*Dgam_z_zx
										+Dgam_y_uxy*Dgam_x_yx+Dgam_y_uxz*Dgam_x_zx+Dgam_y_uyz*Dgam_y_zx
										+Dgam_y_uyx*Dgam_y_xx+Dgam_y_uzx*Dgam_z_xx+Dgam_y_uzy*Dgam_z_yx)
								-gizx_p*(Dgam_z_uxx*Dgam_x_xx+Dgam_z_uyy*Dgam_y_yx+Dgam_z_uzz*Dgam_z_zx
										+Dgam_z_uxy*Dgam_x_yx+Dgam_z_uxz*Dgam_x_zx+Dgam_z_uyz*Dgam_y_zx
										+Dgam_z_uyx*Dgam_y_xx+Dgam_z_uzx*Dgam_z_xx+Dgam_z_uzy*Dgam_z_yx)
								-gizy_p*(Dgam_z_uxx*Dgam_x_xy+Dgam_z_uyy*Dgam_y_yy+Dgam_z_uzz*Dgam_z_zy
										+Dgam_z_uxy*Dgam_x_yy+Dgam_z_uxz*Dgam_x_zy+Dgam_z_uyz*Dgam_y_zy
										+Dgam_z_uyx*Dgam_y_xy+Dgam_z_uzx*Dgam_z_xy+Dgam_z_uzy*Dgam_z_yy)
										
								-gixx_p*(Del_ux_xx*Del_ux_xx+Del_uy_yx*Del_uy_yx+Del_uz_zx*Del_uz_zx
										+2.*(Del_ux_yx*Del_uy_xx+Del_ux_zx*Del_uz_xx+Del_uy_zx*Del_uz_yx))
								-giyy_p*(Del_ux_xy*Del_ux_xy+Del_uy_yy*Del_uy_yy+Del_uz_zy*Del_uz_zy
										+2.*(Del_ux_yy*Del_uy_xy+Del_ux_zy*Del_uz_xy+Del_uy_zy*Del_uz_yy))
								-gizz_p*(Del_ux_xz*Del_ux_xz+Del_uy_yz*Del_uy_yz+Del_uz_zz*Del_uz_zz
										+2.*(Del_ux_yz*Del_uy_xz+Del_ux_zz*Del_uz_xz+Del_uy_zz*Del_uz_yz))
								-2.*
								(gixy_p*(Del_ux_xx*Del_ux_xy+Del_uy_yx*Del_uy_yy+Del_uz_zx*Del_uz_zy
										+2.*(Del_ux_yx*Del_uy_xy+Del_ux_zx*Del_uz_xy+Del_uy_zx*Del_uz_yy))
								+gixz_p*(Del_ux_xx*Del_ux_xz+Del_uy_yx*Del_uy_yz+Del_uz_zx*Del_uz_zz
										+2.*(Del_ux_yx*Del_uy_xz+Del_ux_zx*Del_uz_xz+Del_uy_zx*Del_uz_yz))
								+giyz_p*(Del_ux_xy*Del_ux_xz+Del_uy_yy*Del_uy_yz+Del_uz_zy*Del_uz_zz
										+2.*(Del_ux_yy*Del_uy_xz+Del_ux_zy*Del_uz_xz+Del_uy_zy*Del_uz_yz)));
								
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

				//source term for psi
				set_coef(l,k,j,3)=0.125*Ricci;
				set_coef(l,k,j,4)=-0.125*aaaa;
			}
		}
	}
	
	return;
}

void Fmv::set_source_solveham()
{
	double sumnume=0.;
	double sumdeno=0.;

	#pragma omp parallel for \
	reduction(+:sumnume,sumdeno)
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

				//modified by chulmoon on 2016.2.9 for non-Cartesian inhomogeneous grid
				double fac=facl*fack*facj*sqrt(get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l));

				double gxx_p,gyy_p,gzz_p,
				gxy_p,gxz_p,gyz_p,
				gyx_p,gzx_p,gzy_p;
				
				double det_p,det_pi;
				double gixx_p,giyy_p,gizz_p,
				gixy_p,gixz_p,giyz_p,
				giyx_p,gizx_p,gizy_p;
				
				double wa_xx,wa_yy,wa_zz,
				wa_xy,wa_xz,wa_yz,
				wa_p,wa_x_p,wa_y_p,wa_z_p,wawa;
				
				double rho=0.;
				
			
				gxx_p=  get_bv(l,k,j, 7)+get_flat_df2x(j);
				gyy_p=  get_bv(l,k,j, 8)+get_flat_df2y(k);
				gzz_p=  get_bv(l,k,j, 9)+get_flat_df2z(l);
				gxy_p=  get_bv(l,k,j,10);
				gxz_p=  get_bv(l,k,j,11);
				gyz_p=  get_bv(l,k,j,12);
				gyx_p=	gxy_p;
				gzx_p=	gxz_p;
				gzy_p=	gyz_p;
				

				// tilde{gamma}^{ij} -> gi
				det_p=gxx_p*gyy_p*gzz_p+gxy_p*gyz_p*gzx_p+gxz_p*gyx_p*gzy_p
						-gxz_p*gyy_p*gzx_p-gxy_p*gyx_p*gzz_p-gxx_p*gyz_p*gzy_p;
				if(det_p<1.e-16) 
				det_p=1.e-16;

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
				
				wa_p=   get_bv(l,k,j,13);
				
				wa_x_p=get_f_x(l,k,j,13);
				wa_y_p=get_f_y(l,k,j,13);
				wa_z_p=get_f_z(l,k,j,13);
				
				wawa=gixx_p*wa_x_p*wa_x_p +gixy_p*wa_x_p*wa_y_p +gixz_p*wa_x_p*wa_z_p
					+giyx_p*wa_y_p*wa_x_p +giyy_p*wa_y_p*wa_y_p +giyz_p*wa_y_p*wa_z_p
					+gizx_p*wa_z_p*wa_x_p +gizy_p*wa_z_p*wa_y_p +gizz_p*wa_z_p*wa_z_p;

				//\del_i\del_j \psi
				wa_xx=get_f_xx(l,k,j,13);
				wa_yy=get_f_yy(l,k,j,13);
				wa_zz=get_f_zz(l,k,j,13);
				wa_xy=get_f_xy(l,k,j,13);
				wa_xz=get_f_xz(l,k,j,13);
				wa_yz=get_f_yz(l,k,j,13);

				// \tilde Laplacian \psi
				double wa_lap=gixx_p*wa_xx +giyy_p*wa_yy +gizz_p*wa_zz
						+2.*(gixy_p*wa_xy +gixz_p*wa_xz +giyz_p*wa_yz)
						-get_coef(l,k,j,0)*wa_x_p
						-get_coef(l,k,j,1)*wa_y_p
						-get_coef(l,k,j,2)*wa_z_p;
				
				double dfxm2=1./get_flat_df2x(j);
				double dfym2=1./get_flat_df2y(k);
				double dfzm2=1./get_flat_df2z(l);
				double ddfpdfx=get_flat_Gamx(j);
				double ddfpdfy=get_flat_Gamy(k);
				double ddfpdfz=get_flat_Gamz(l);
				
				double lhs=(wa_xx-ddfpdfx*wa_x_p)*dfxm2
							+(wa_yy-ddfpdfy*wa_y_p)*dfym2
							+(wa_zz-ddfpdfz*wa_z_p)*dfzm2;
				
				if(fluidevo)
				{
					double sqgam;
					sqgam=sqrt(det_p)*exp(6.*wa_p);
					rho+=get_bv(l,k,j,24)/sqgam;
				}
				
				if(scalarevo)
				{
					double phi,Pi,phi_x,phi_y,phi_z,dphidphi,Vpot;
					phi=get_bv(l,k,j,nsc);
					Pi=get_bv(l,k,j,nscp);
				
					phi_x=get_f_x(l,k,j,nsc);
					phi_y=get_f_y(l,k,j,nsc);
					phi_z=get_f_z(l,k,j,nsc);
				
					dphidphi=gixx_p*phi_x*phi_x +gixy_p*phi_x*phi_y +gixz_p*phi_x*phi_z
						+giyx_p*phi_y*phi_x +giyy_p*phi_y*phi_y +giyz_p*phi_y*phi_z
						+gizx_p*phi_z*phi_x +gizy_p*phi_z*phi_y +gizz_p*phi_z*phi_z;

					Vpot=funcV(phi);
					rho+=0.5*(pow(Pi,2)+exp(-4.*wa_p)*dphidphi)+Vpot;
				}
				
				set_coef(l,k,j,5)=-wawa-wa_lap+get_coef(l,k,j,3)
								+get_coef(l,k,j,4)*exp(4.*wa_p)
								-pi2*exp(4.*wa_p)*rho;
				
				double k2P5=pow(get_bv(l,k,j,20),2)/12.*exp(5.*wa_p); 
				
				sumnume-=get_coef(l,k,j,5)*exp(wa_p)*fac;
				sumdeno+=k2P5*fac;
				
				set_coef(l,k,j,5)=get_coef(l,k,j,5)
								//+k2P5/exp(wa_p)	will be added after modification of wa_p for integ. cond.
								+lhs;
												
				set_coef(l,k,j,6)=abs(wawa)+abs(wa_lap)+abs(get_coef(l,k,j,3))
								+abs(get_coef(l,k,j,4)
								//+pow(get_bv(l,k,j,20),2)/12.		will be added after modification of wa_p for integ. cond.
								)*exp(4.*wa_p)
								+pi2*exp(4.*wa_p)*abs(rho)
								+abs(lhs);

			}
		}
	}
	double ofac2=sumnume/sumdeno;
	
	if(ofac2<0.)
	{
		cout << "ofac2<0 ofac2=" << ofac2 << endl;
		exit(0);
	}
	
	double ofac=sqrt(ofac2);
	
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				set_bv(l,k,j,20)=get_bv(l,k,j,20)*ofac;
				set_coef(l,k,j,5)+=pow(get_bv(l,k,j,20),2)/12.*exp(4.*get_bv(l,k,j,13));
				set_coef(l,k,j,6)+=abs(pow(get_bv(l,k,j,20),2)/12.*exp(4.*get_bv(l,k,j,13)));
			}
		}
	}
	
	boundary_reflection_even(20);

	return;
}

//iteration with Jacobi method for hamiltonian constraint
void Fmv::solveham(int msstep,double cc,double acc)
//cc: convergence condition
//acc: acceleration/deceleration parameter
{
	//for convergence cond. 
	double errm=0.;						// max error
	double convcond=cc;					// convergence condition
	bool conv=false;					// convergence or not
	
	//set coefficients
	coefficients_solveham();
	//set source term
	set_source_solveham();
	
	//Jacobi method iteration
	for(int n=0;n<msstep;n++)
	{
		boundary_reflection_even(13);
		set_psi();
		errm=0.;
	
		#pragma omp parallel for default(shared) reduction(max:errm)
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int j=jli;j<=jui;j++)
				{
					double wa_x_p=(    -(get_psi(l,k,j+2)-get_psi(l,k,j-2))
							+8.*(get_psi(l,k,j+1)-get_psi(l,k,j-1))
							)/12.;
					double wa_y_p=(   -(get_psi(l,k+2,j)-get_psi(l,k-2,j))
							+8.*(get_psi(l,k+1,j)-get_psi(l,k-1,j))
							)/12.;
					double wa_z_p=(   -(get_psi(l+2,k,j)-get_psi(l-2,k,j))
							+8.*(get_psi(l+1,k,j)-get_psi(l-1,k,j))
							)/12.;

					//\del_i\del_j \psi
					double wa_xx=(    -(    get_psi(l  ,k  ,j+2) + get_psi(l  ,k  ,j-2))
							+16.*( get_psi(l  ,k  ,j+1) + get_psi(l  ,k  ,j-1))
							-30.*  get_psi(l  ,k  ,j  )
							)/12.;
					double wa_yy=(    -(    get_psi(l  ,k+2,j  ) + get_psi(l  ,k-2,j  ))
							+16.*( get_psi(l  ,k+1,j  ) + get_psi(l  ,k-1,j  ))
							-30.*  get_psi(l  ,k  ,j  )
							)/12.;
					double wa_zz=(    -(    get_psi(l+2,k  ,j  ) + get_psi(l-2,k  ,j  ))
							+16.*( get_psi(l+1,k  ,j  ) + get_psi(l-1,k  ,j  ))
							-30.*  get_psi(l  ,k  ,j  )
							)/12.;

					double dfxm2=1./get_flat_df2x(j);
					double dfym2=1./get_flat_df2y(k);
					double dfzm2=1./get_flat_df2z(l);
					double ddfpdfx=get_flat_Gamx(j);
					double ddfpdfy=get_flat_Gamy(k);
					double ddfpdfz=get_flat_Gamz(l);

					double dev=(
						dfxm2*wa_xx
						+dfym2*wa_yy
						+dfzm2*wa_zz
						-dfxm2*ddfpdfx*wa_x_p*dx
						-dfym2*ddfpdfy*wa_y_p*dx
						-dfzm2*ddfpdfz*wa_z_p*dx
						-get_coef(l,k,j,5)*dx*dx
						)/(2.5*(dfxm2+dfym2+dfzm2));
					
					double absv=abs(get_psi(l,k,j));

					double deno=absv+0.5*abs((dx*dx*get_coef(l,k,j,6))/(dfxm2+dfym2+dfzm2));
					double err=abs(dev)/deno;
					
					if(err>1.0e6)
					{
						cout << "large error in solveham" << endl;
						cout << "j=" << j << endl
						<< "k=" << k << endl 
						<< "l=" << l << endl
						<< "n=" << n << endl 
						<< "err=" << err << endl 
						<< "dev=" << dev << endl 
						<< endl;
						
						cout << " get_psi(l,k,j)=" << get_psi(l,k,j) << endl 
							<< " get_psi(l,k,j+1)=" << get_psi(l,k,j+1) << endl 
							<< " get_psi(l,k,j-1)=" << get_psi(l,k,j-1) << endl 
							<< " get_psi(l,k+1,j)=" << get_psi(l,k+1,j) << endl 
							<< " get_psi(l,k-1,j)=" << get_psi(l,k-1,j) << endl 
							<< " get_psi(l+1,k,j)=" << get_psi(l+1,k,j) << endl 
							<< " get_psi(l-1,k,j)=" << get_psi(l-1,k,j) << endl 
							<< " get_psi(l,k-1,j-1)=" << get_psi(l,k-1,j-1) << endl 
							<< " get_psi(l,k+1,j+1)=" << get_psi(l,k+1,j+1) << endl 
							<< " get_coef(l,k,j,0)=" << get_coef(l,k,j,0) << endl 
							<< " get_coef(l,k,j,1)=" << get_coef(l,k,j,1) << endl 
							<< " get_coef(l,k,j,2)=" << get_coef(l,k,j,2) << endl 
							<< " get_coef(l,k,j,3)=" << get_coef(l,k,j,3) << endl 
							<< " get_coef(l,k,j,4)=" << get_coef(l,k,j,4) << endl 
							<< " get_coef(l,k,j,5)=" << get_coef(l,k,j,5) << endl 
							<< " get_coef(l,k,j,6)=" << get_coef(l,k,j,6) << endl 
							<< endl;
					
						exit(1);
					}
					
					
					set_bv(l,k,j,13)=get_bv(l,k,j,13)+acc*dev;
					
					if (errm<err)
					errm=err;
				}
			}
		}

		if(errm<convcond)
		{
			cout << "solveham step=" << n 
				<< "errm solveham=" << errm << endl;
			conv=true;
			break;
			
		}
		else
		{
				cout << "solveham step=" << n 
				<< "errm solveham=" << errm << endl;
		}

		boundary_reflection_even(13);
		set_source_solveham();
	}
	//Jacobi method iteration end
	
	if(!conv)
	{
		cout << "not converge in solveham step" 
			<< " errm=" << errm << endl;
	}
	
	return;
}

//iteration with SOR for hamiltonian constraint no parallelization
void Fmv::solveham_sor_mono(int msstep,double cc,double acc)
//cc: convergence condition
//acc: acceleration/deceleration parameter for SOR. 
{
	//for convergence cond. 
	double errm=0.;						// max error
	double convcond=cc;					// convergence condition
	bool conv=false;					// convergence or not
	
	//set coefficients
	coefficients_solveham();
	//set source term
	set_source_solveham();
		
	for(int n=0;n<msstep;n++)
	{
		errm=0.;
		
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int j=jli;j<=jui;j++)
				{
					double wa_p=   get_bv(l,k,j,13);
					
					double wa_x_p=(    -(get_bv(l,k,j+2,13)-get_bv(l,k,j-2,13))
							+8.*(get_bv(l,k,j+1,13)-get_bv(l,k,j-1,13))
							)/12.;
					double wa_y_p=(   -(get_bv(l,k+2,j,13)-get_bv(l,k-2,j,13))
							+8.*(get_bv(l,k+1,j,13)-get_bv(l,k-1,j,13))
							)/12.;
					double wa_z_p=(   -(get_bv(l+2,k,j,13)-get_bv(l-2,k,j,13))
							+8.*(get_bv(l+1,k,j,13)-get_bv(l-1,k,j,13))
							)/12.;

					//\del_i\del_j \psi
					double wa_xx=(    -(    get_bv(l  ,k  ,j+2,13) + get_bv(l  ,k  ,j-2,13))
							+16.*( get_bv(l  ,k  ,j+1,13) + get_bv(l  ,k  ,j-1,13))
							-30.*  get_bv(l  ,k  ,j  ,13)
							)/12.;
					double wa_yy=(    -(    get_bv(l  ,k+2,j  ,13) + get_bv(l  ,k-2,j  ,13))
							+16.*( get_bv(l  ,k+1,j  ,13) + get_bv(l  ,k-1,j  ,13))
							-30.*  get_bv(l  ,k  ,j  ,13)
							)/12.;
					double wa_zz=(    -(    get_bv(l+2,k  ,j  ,13) + get_bv(l-2,k  ,j  ,13))
							+16.*( get_bv(l+1,k  ,j  ,13) + get_bv(l-1,k  ,j  ,13))
							-30.*  get_bv(l  ,k  ,j  ,13)
							)/12.;
					
					double dfxm2=1./get_flat_df2x(j);
					double dfym2=1./get_flat_df2y(k);
					double dfzm2=1./get_flat_df2z(l);
					double ddfpdfx=get_flat_Gamx(j);
					double ddfpdfy=get_flat_Gamy(k);
					double ddfpdfz=get_flat_Gamz(l);
					
					
					double dev=(
						dfxm2*wa_xx
						+dfym2*wa_yy
						+dfzm2*wa_zz
						-dfxm2*ddfpdfx*wa_x_p*dx
						-dfym2*ddfpdfy*wa_y_p*dx
						-dfzm2*ddfpdfz*wa_z_p*dx
						-get_coef(l,k,j,5)*dx*dx
						)/(2.5*(dfxm2+dfym2+dfzm2));
					
					double absv=abs(wa_p);

					double deno=absv+0.5*abs((dx*dx*get_coef(l,k,j,6))/(dfxm2+dfym2+dfzm2));
					double err=abs(dev)/deno;
					
					if(err>1.0e6)
					{
						cout << "large error in solveham" << endl;
						cout << "j=" << j << endl
						<< "k=" << k << endl 
						<< "l=" << l << endl
						<< "n=" << n << endl 
						<< "err=" << err << endl 
						<< "dev=" << dev << endl 
						//<< "norm=" << norm 
						<< endl;
						
						cout << " get_bv(l,k,j,13)=" << get_bv(l,k,j,13) << endl 
							<< " get_bv(l,k,j+1,13)=" << get_bv(l,k,j+1,13) << endl 
							<< " get_bv(l,k,j-1,13)=" << get_bv(l,k,j-1,13) << endl 
							<< " get_bv(l,k+1,j,13)=" << get_bv(l,k+1,j,13) << endl 
							<< " get_bv(l,k-1,j,13)=" << get_bv(l,k-1,j,13) << endl 
							<< " get_bv(l+1,k,j,13)=" << get_bv(l+1,k,j,13) << endl 
							<< " get_bv(l-1,k,j,13)=" << get_bv(l-1,k,j,13) << endl 
							<< " get_bv(l,k-1,j-1,13)=" << get_bv(l,k-1,j-1,13) << endl 
							<< " get_bv(l,k+1,j+1,13)=" << get_bv(l,k+1,j+1,13) << endl 
							<< " get_coef(l,k,j,0)=" << get_coef(l,k,j,0) << endl 
							<< " get_coef(l,k,j,1)=" << get_coef(l,k,j,1) << endl 
							<< " get_coef(l,k,j,2)=" << get_coef(l,k,j,2) << endl 
							<< " get_coef(l,k,j,3)=" << get_coef(l,k,j,3) << endl 
							<< " get_coef(l,k,j,4)=" << get_coef(l,k,j,4) << endl 
							<< " get_coef(l,k,j,5)=" << get_coef(l,k,j,5) << endl 
							<< " get_coef(l,k,j,6)=" << get_coef(l,k,j,6) << endl 
							<< endl;
					
						exit(1);
					}
					
					
					set_bv(l,k,j,13)=wa_p+acc*dev;

					if (errm<err)
					errm=err;
				}
			}
		}
		
		if(errm<convcond)
		{
			cout << "solveham step=" << n 
				<< "errm solveham=" << errm << endl;
			conv=true;
			break;
			
		}
		else
		{
				cout << "solveham step=" << n 
			   << "errm solveham=" << errm << endl;
		}

		boundary_reflection_even(13);
		set_source_solveham();
		
	}
	
	if(!conv)
	{
		cout << "not converge in solveham step" 
			 << " errm=" << errm << endl;
	}

	return;
}

// setting source terms of all constraints(Hamiltonian and momentum constraints)
void Fmv::set_source()
{
	double sumnume=0.;
	double sumdeno=0.;

	#pragma omp parallel for reduction(+:sumnume,sumdeno)
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
				
				////////// variables /////////////////////////////////////////////////////////////////////////////////////////////
				double hamc,fac;																								//

				double 
				gxx_p,gyy_p,gzz_p,
				gxy_p,gxz_p,gyz_p,
				gyx_p,gzx_p,gzy_p,wa_p,
				akxx_p,akyy_p,akzz_p,
				akxy_p,akxz_p,akyz_p,
				akyx_p,akzx_p,akzy_p;

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
				double rct_xx_p,rct_yy_p,rct_zz_p,
				rct_xy_p,rct_xz_p,rct_yz_p,
				rct_yx_p,rct_zx_p,rct_zy_p;
				double aaaa,ricci;
				
				double rho=0.;
				double p_x=0.;
				double p_y=0.;
				double p_z=0.;

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
				
				double Daxx_x,Dayy_x,Dazz_x,Daxy_x,Daxz_x,Dayz_x;
				double Daxx_y,Dayy_y,Dazz_y,Daxy_y,Daxz_y,Dayz_y;
				double Daxx_z,Dayy_z,Dazz_z,Daxy_z,Daxz_z,Dayz_z;
				
				double Da_x,Da_y,Da_z,M_x,M_y,M_z;
				
				double M_ux,M_uy,M_uz;
				
				double wx_p,wy_p,wz_p;

				double wx_x_p,wy_x_p,wz_x_p;
				double wx_y_p,wy_y_p,wz_y_p;
				double wx_z_p,wy_z_p,wz_z_p;

				double wx_xx,wx_yy,wx_zz,
				wx_xy,wx_xz,wx_yz,
				wx_yx,wx_zx,wx_zy;
				double wy_xx,wy_yy,wy_zz,
				wy_xy,wy_xz,wy_yz,
				wy_yx,wy_zx,wy_zy;
				double wz_xx,wz_yy,wz_zz,
				wz_xy,wz_xz,wz_yz,
				wz_yx,wz_zx,wz_zy;
				
				double Dw_x_ux,Dw_x_uy,Dw_x_uz,Dw_y_ux,Dw_y_uy,Dw_y_uz,Dw_z_ux,Dw_z_uy,Dw_z_uz;
				double Dw_uxx,Dw_uxy,Dw_uxz,Dw_uyx,Dw_uyy,Dw_uyz,Dw_uzx,Dw_uzy,Dw_uzz,divDw;
				double wx_lap,wy_lap,wz_lap;
				double divw_x,divw_y,divw_z;
				
				double Lw_x_ux,Lw_y_ux,Lw_z_ux,Lw_x_uy,Lw_y_uy,Lw_z_uy,Lw_x_uz,Lw_y_uz,Lw_z_uz;
				double Lw_uxx,Lw_uxy,Lw_uxz,Lw_uyx,Lw_uyy,Lw_uyz,Lw_uzx,Lw_uzy,Lw_uzz;
				double Lw_xx,Lw_xy,Lw_xz,Lw_yx,Lw_yy,Lw_yz,Lw_zx,Lw_zy,Lw_zz;
				
				double Del_uxx_x,Del_uxx_y,Del_uxx_z,
					Del_uxy_x,Del_uxy_y,Del_uxy_z,
					Del_uxz_x,Del_uxz_y,Del_uxz_z,
					Del_uyx_x,Del_uyx_y,Del_uyx_z,
					Del_uyy_x,Del_uyy_y,Del_uyy_z,
					Del_uyz_x,Del_uyz_y,Del_uyz_z,
					Del_uzx_x,Del_uzx_y,Del_uzx_z,
					Del_uzy_x,Del_uzy_y,Del_uzy_z,
					Del_uzz_x,Del_uzz_y,Del_uzz_z,
					DGw_xx,DGw_xy,DGw_xz,
					DGw_yx,DGw_yy,DGw_yz,
					DGw_zx,DGw_zy,DGw_zz;																						//
				////////// variables /////////////////////////////////////////////////////////////////////////////////////////////

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

				akxx_p= get_bv(l,k,j,14);
				akyy_p= get_bv(l,k,j,15);
				akzz_p= get_bv(l,k,j,16);
				akxy_p= get_bv(l,k,j,17);
				akxz_p= get_bv(l,k,j,18);
				akyz_p= get_bv(l,k,j,19);
				akyx_p=	akxy_p;
				akzx_p=	akxz_p;
				akzy_p=	akyz_p;
				//ek_p=   get_bv(l,k,j,20);
				
				zgx_p=  get_bv(l,k,j,21);
				zgy_p=  get_bv(l,k,j,22);
				zgz_p=  get_bv(l,k,j,23);
				
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
				
				//\del_x \Gamma^i
				zgx_x_p=get_f_x(l,k,j,21);
				zgy_x_p=get_f_x(l,k,j,22);
				zgz_x_p=get_f_x(l,k,j,23);

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
				
				// \del_y \Gamma^i
				zgx_y_p=get_f_y(l,k,j,21);
				zgy_y_p=get_f_y(l,k,j,22);
				zgz_y_p=get_f_y(l,k,j,23);

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
				
				// \del_z \Gamma^i
				zgx_z_p=get_f_z(l,k,j,21);
				zgy_z_p=get_f_z(l,k,j,22);
				zgz_z_p=get_f_z(l,k,j,23);

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
				
				// wx,wy,wz: variables to be solved for momentum constraints
				wx_p=   get_wvec(l,k,j, 0);
				wy_p=   get_wvec(l,k,j, 1);
				wz_p=   get_wvec(l,k,j, 2);
				
				wx_x_p=(   -(get_wvec(l,k,j+2,0)-get_wvec(l,k,j-2,0))
						+8.*(get_wvec(l,k,j+1,0)-get_wvec(l,k,j-1,0))
						)*dxi12;
				wy_x_p=(-   (get_wvec(l,k,j+2,1)-get_wvec(l,k,j-2,1))
						+8.*(get_wvec(l,k,j+1,1)-get_wvec(l,k,j-1,1))
						)*dxi12;
				wz_x_p=(   -(get_wvec(l,k,j+2,2)-get_wvec(l,k,j-2,2))
						+8.*(get_wvec(l,k,j+1,2)-get_wvec(l,k,j-1,2))
						)*dxi12;

				wx_y_p=(   -(get_wvec(l,k+2,j,0)-get_wvec(l,k-2,j,0))
						+8.*(get_wvec(l,k+1,j,0)-get_wvec(l,k-1,j,0))
						)*dyi12;
				wy_y_p=(   -(get_wvec(l,k+2,j,1)-get_wvec(l,k-2,j,1))
						+8.*(get_wvec(l,k+1,j,1)-get_wvec(l,k-1,j,1))
						)*dyi12;
				wz_y_p=(   -(get_wvec(l,k+2,j,2)-get_wvec(l,k-2,j,2))
						+8.*(get_wvec(l,k+1,j,2)-get_wvec(l,k-1,j,2))
						)*dyi12;

				// \del_z w^i
				wx_z_p=(   -(get_wvec(l+2,k,j,0)-get_wvec(l-2,k,j,0))
						+8.*(get_wvec(l+1,k,j,0)-get_wvec(l-1,k,j,0))
						)*dzi12;
				wy_z_p=(   -(get_wvec(l+2,k,j,1)-get_wvec(l-2,k,j,1))
						+8.*(get_wvec(l+1,k,j,1)-get_wvec(l-1,k,j,1))
						)*dzi12;
				wz_z_p=(   -(get_wvec(l+2,k,j,2)-get_wvec(l-2,k,j,2))
						+8.*(get_wvec(l+1,k,j,2)-get_wvec(l-1,k,j,2))
						)*dzi12;

				// \del_i\del_j w^x
				wx_xx=(    -(     get_wvec(l  ,k  ,j+2,0) + get_wvec(l  ,k  ,j-2,0))
						+16.*( get_wvec(l  ,k  ,j+1,0) + get_wvec(l  ,k  ,j-1,0))
						-30.*  get_wvec(l  ,k  ,j  ,0)
						)*dxi*dxi12;
				wx_yy=(    -(     get_wvec(l  ,k+2,j  ,0) + get_wvec(l  ,k-2,j  ,0))
						+16.*( get_wvec(l  ,k+1,j  ,0) + get_wvec(l  ,k-1,j  ,0))
						-30.*  get_wvec(l  ,k  ,j  ,0)
						)*dyi*dyi12;
				wx_zz=(    -(     get_wvec(l+2,k  ,j  ,0) + get_wvec(l-2,k  ,j  ,0))
						+16.*( get_wvec(l+1,k  ,j  ,0) + get_wvec(l-1,k  ,j  ,0))
						-30.*  get_wvec(l  ,k  ,j  ,0)
						)*dzi*dzi12;

				wx_xy=( -( ( -(   get_wvec(l  ,k+2,j+2,0) - get_wvec(l  ,k+2,j-2,0))
						+8.*(get_wvec(l  ,k+2,j+1,0) - get_wvec(l  ,k+2,j-1,0)))
						-(   -(get_wvec(l  ,k-2,j+2,0) - get_wvec(l  ,k-2,j-2,0))
						+8.*(get_wvec(l  ,k-2,j+1,0) - get_wvec(l  ,k-2,j-1,0))))
						+8.*(  (-(get_wvec(l  ,k+1,j+2,0) - get_wvec(l  ,k+1,j-2,0))
						+8.*(get_wvec(l  ,k+1,j+1,0) - get_wvec(l  ,k+1,j-1,0)))
						-(-(get_wvec(l  ,k-1,j+2,0) - get_wvec(l  ,k-1,j-2,0))
						+8.*(get_wvec(l  ,k-1,j+1,0) - get_wvec(l  ,k-1,j-1,0))))
						)*dxi12*dyi12;
				wx_xz=( -( ( -(   get_wvec(l+2,k  ,j+2,0) - get_wvec(l+2,k  ,j-2,0))
						+8.*(get_wvec(l+2,k  ,j+1,0) - get_wvec(l+2,k  ,j-1,0)))
						-(   -(get_wvec(l-2,k  ,j+2,0) - get_wvec(l-2,k  ,j-2,0))
						+8.*(get_wvec(l-2,k  ,j+1,0) - get_wvec(l-2,k  ,j-1,0))))
						+8.*(  (-(get_wvec(l+1,k  ,j+2,0) - get_wvec(l+1,k  ,j-2,0))
						+8.*(get_wvec(l+1,k  ,j+1,0) - get_wvec(l+1,k  ,j-1,0)))
						-(-(get_wvec(l-1,k  ,j+2,0) - get_wvec(l-1,k  ,j-2,0))
						+8.*(get_wvec(l-1,k  ,j+1,0) - get_wvec(l-1,k  ,j-1,0))))
						)*dxi12*dzi12;
				wx_yz=( -( ( -(   get_wvec(l+2,k+2,j  ,0) - get_wvec(l+2,k-2,j  ,0))
						+8.*(get_wvec(l+2,k+1,j  ,0) - get_wvec(l+2,k-1,j  ,0)))
						-(   -(get_wvec(l-2,k+2,j  ,0) - get_wvec(l-2,k-2,j  ,0))
						+8.*(get_wvec(l-2,k+1,j  ,0) - get_wvec(l-2,k-1,j  ,0))))
						+8.*(  (-(get_wvec(l+1,k+2,j  ,0) - get_wvec(l+1,k-2,j  ,0))
						+8.*(get_wvec(l+1,k+1,j  ,0) - get_wvec(l+1,k-1,j  ,0)))
						-(-(get_wvec(l-1,k+2,j  ,0) - get_wvec(l-1,k-2,j  ,0))
						+8.*(get_wvec(l-1,k+1,j  ,0) - get_wvec(l-1,k-1,j  ,0))))
						)*dyi12*dzi12;
				wx_yx=wx_xy;
				wx_zx=wx_xz;
				wx_zy=wx_yz;

				// \del_i\del_j w^y
				wy_xx=(    -(    get_wvec(l  ,k  ,j+2,1) + get_wvec(l  ,k  ,j-2,1))
						+16.*( get_wvec(l  ,k  ,j+1,1) + get_wvec(l  ,k  ,j-1,1))
						-30.*  get_wvec(l  ,k  ,j  ,1)
						)*dxi*dxi12;
				wy_yy=(    -(    get_wvec(l  ,k+2,j  ,1) + get_wvec(l  ,k-2,j  ,1))
						+16.*( get_wvec(l  ,k+1,j  ,1) + get_wvec(l  ,k-1,j  ,1))
						-30.*  get_wvec(l  ,k  ,j  ,1)
						)*dyi*dyi12;
				wy_zz=(    -(    get_wvec(l+2,k  ,j  ,1) + get_wvec(l-2,k  ,j  ,1))
						+16.*( get_wvec(l+1,k  ,j  ,1) + get_wvec(l-1,k  ,j  ,1))
						-30.*  get_wvec(l  ,k  ,j  ,1)
						)*dzi*dzi12;

				wy_xy=( -( ( -(      get_wvec(l  ,k+2,j+2,1) - get_wvec(l  ,k+2,j-2,1))
						+8.*(get_wvec(l  ,k+2,j+1,1) - get_wvec(l  ,k+2,j-1,1)))
						-(   -(get_wvec(l  ,k-2,j+2,1) - get_wvec(l  ,k-2,j-2,1))
						+8.*(get_wvec(l  ,k-2,j+1,1) - get_wvec(l  ,k-2,j-1,1))))
						+8.*(  (-(get_wvec(l  ,k+1,j+2,1) - get_wvec(l  ,k+1,j-2,1))
						+8.*(get_wvec(l  ,k+1,j+1,1) - get_wvec(l  ,k+1,j-1,1)))
						-(-(get_wvec(l  ,k-1,j+2,1) - get_wvec(l  ,k-1,j-2,1))
						+8.*(get_wvec(l  ,k-1,j+1,1) - get_wvec(l  ,k-1,j-1,1))))
						)*dxi12*dyi12;
				wy_xz=( -( ( -(      get_wvec(l+2,k  ,j+2,1) - get_wvec(l+2,k  ,j-2,1))
						+8.*(get_wvec(l+2,k  ,j+1,1) - get_wvec(l+2,k  ,j-1,1)))
						-(   -(get_wvec(l-2,k  ,j+2,1) - get_wvec(l-2,k  ,j-2,1))
						+8.*(get_wvec(l-2,k  ,j+1,1) - get_wvec(l-2,k  ,j-1,1))))
						+8.*(  (-(get_wvec(l+1,k  ,j+2,1) - get_wvec(l+1,k  ,j-2,1))
						+8.*(get_wvec(l+1,k  ,j+1,1) - get_wvec(l+1,k  ,j-1,1)))
						-(-(get_wvec(l-1,k  ,j+2,1) - get_wvec(l-1,k  ,j-2,1))
						+8.*(get_wvec(l-1,k  ,j+1,1) - get_wvec(l-1,k  ,j-1,1))))
						)*dxi12*dzi12;
				wy_yz=( -( ( -(      get_wvec(l+2,k+2,j  ,1) - get_wvec(l+2,k-2,j  ,1))
						+8.*(get_wvec(l+2,k+1,j  ,1) - get_wvec(l+2,k-1,j  ,1)))
						-(   -(get_wvec(l-2,k+2,j  ,1) - get_wvec(l-2,k-2,j  ,1))
						+8.*(get_wvec(l-2,k+1,j  ,1) - get_wvec(l-2,k-1,j  ,1))))
						+8.*(  (-(get_wvec(l+1,k+2,j  ,1) - get_wvec(l+1,k-2,j  ,1))
						+8.*(get_wvec(l+1,k+1,j  ,1) - get_wvec(l+1,k-1,j  ,1)))
						-(-(get_wvec(l-1,k+2,j  ,1) - get_wvec(l-1,k-2,j  ,1))
						+8.*(get_wvec(l-1,k+1,j  ,1) - get_wvec(l-1,k-1,j  ,1))))
						)*dyi12*dzi12;
				wy_yx=wy_xy;
				wy_zx=wy_xz;
				wy_zy=wy_yz;

				// \del_i\del_j w^z
				wz_xx=(    -(    get_wvec(l  ,k  ,j+2,2) + get_wvec(l  ,k  ,j-2,2))
						+16.*( get_wvec(l  ,k  ,j+1,2) + get_wvec(l  ,k  ,j-1,2))
						-30.*  get_wvec(l  ,k  ,j  ,2)
						)*dxi*dxi12;
				wz_yy=(    -(    get_wvec(l  ,k+2,j  ,2) + get_wvec(l  ,k-2,j  ,2))
						+16.*( get_wvec(l  ,k+1,j  ,2) + get_wvec(l  ,k-1,j  ,2))
						-30.*  get_wvec(l  ,k  ,j  ,2)
						)*dyi*dyi12;
				wz_zz=(    -(    get_wvec(l+2,k  ,j  ,2) + get_wvec(l-2,k  ,j  ,2))
						+16.*( get_wvec(l+1,k  ,j  ,2) + get_wvec(l-1,k  ,j  ,2))
						-30.*  get_wvec(l  ,k  ,j  ,2)
						)*dzi*dzi12;

				wz_xy=( -( ( -(      get_wvec(l  ,k+2,j+2,2) - get_wvec(l  ,k+2,j-2,2))
						+8.*(get_wvec(l  ,k+2,j+1,2) - get_wvec(l  ,k+2,j-1,2)))
						-(   -(get_wvec(l  ,k-2,j+2,2) - get_wvec(l  ,k-2,j-2,2))
						+8.*(get_wvec(l  ,k-2,j+1,2) - get_wvec(l  ,k-2,j-1,2))))
						+8.*(  (-(get_wvec(l  ,k+1,j+2,2) - get_wvec(l  ,k+1,j-2,2))
						+8.*(get_wvec(l  ,k+1,j+1,2) - get_wvec(l  ,k+1,j-1,2)))
						-(-(get_wvec(l  ,k-1,j+2,2) - get_wvec(l  ,k-1,j-2,2))
						+8.*(get_wvec(l  ,k-1,j+1,2) - get_wvec(l  ,k-1,j-1,2))))
						)*dxi12*dyi12;
				wz_xz=( -( ( -(      get_wvec(l+2,k  ,j+2,2) - get_wvec(l+2,k  ,j-2,2))
						+8.*(get_wvec(l+2,k  ,j+1,2) - get_wvec(l+2,k  ,j-1,2)))
						-(   -(get_wvec(l-2,k  ,j+2,2) - get_wvec(l-2,k  ,j-2,2))
						+8.*(get_wvec(l-2,k  ,j+1,2) - get_wvec(l-2,k  ,j-1,2))))
						+8.*(  (-(get_wvec(l+1,k  ,j+2,2) - get_wvec(l+1,k  ,j-2,2))
						+8.*(get_wvec(l+1,k  ,j+1,2) - get_wvec(l+1,k  ,j-1,2)))
						-(-(get_wvec(l-1,k  ,j+2,2) - get_wvec(l-1,k  ,j-2,2))
						+8.*(get_wvec(l-1,k  ,j+1,2) - get_wvec(l-1,k  ,j-1,2))))
						)*dxi12*dzi12;
				wz_yz=( -( ( -(      get_wvec(l+2,k+2,j  ,2) - get_wvec(l+2,k-2,j  ,2))
						+8.*(get_wvec(l+2,k+1,j  ,2) - get_wvec(l+2,k-1,j  ,2)))
						-(   -(get_wvec(l-2,k+2,j  ,2) - get_wvec(l-2,k-2,j  ,2))
						+8.*(get_wvec(l-2,k+1,j  ,2) - get_wvec(l-2,k-1,j  ,2))))
						+8.*(  (-(get_wvec(l+1,k+2,j  ,2) - get_wvec(l+1,k-2,j  ,2))
						+8.*(get_wvec(l+1,k+1,j  ,2) - get_wvec(l+1,k-1,j  ,2)))
						-(-(get_wvec(l-1,k+2,j  ,2) - get_wvec(l-1,k-2,j  ,2))
						+8.*(get_wvec(l-1,k+1,j  ,2) - get_wvec(l-1,k-1,j  ,2))))
						)*dyi12*dzi12;
				wz_yx=wz_xy;
				wz_zx=wz_xz;
				wz_zy=wz_yz;
				
				// tilde{gamma}^{ij} -> gi
				det_p=gxx_p*gyy_p*gzz_p+gxy_p*gyz_p*gzx_p+gxz_p*gyx_p*gzy_p
						-gxz_p*gyy_p*gzx_p-gxy_p*gyx_p*gzz_p-gxx_p*gyz_p*gzy_p;
				if(det_p<1.e-16) 
				 det_p=1.e-16;
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
				
				
				// \tilde D_i D_j \psi = \psi,ij -\tilde{\Gamma}^k_ij \psi_k
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
				// see other files for previous version
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
				
				rct_xx_p=(0.5*(-gDDg_xx
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
						-Del_uy_xx*Del_ux_yx-Del_uz_xx*Del_ux_zx-Del_uz_xy*Del_uy_zx)*exp(-4.*wa_p);
				
				rct_yy_p=(0.5*(-gDDg_yy
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
						)*exp(-4.*wa_p);

				rct_zz_p=(0.5*(-gDDg_zz
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
						-Del_uy_zx*Del_ux_yz-Del_uz_zx*Del_ux_zz-Del_uz_zy*Del_uy_zz)*exp(-4.*wa_p);

				rct_xy_p=(0.5*(-gDDg_xy
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
						-Del_uy_xx*Del_ux_yy-Del_uz_xx*Del_ux_zy-Del_uz_xy*Del_uy_zy)*exp(-4.*wa_p);

				rct_xz_p=(0.5*(-gDDg_xz
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
						-Del_uy_xx*Del_ux_yz-Del_uz_xx*Del_ux_zz-Del_uz_xy*Del_uy_zz)*exp(-4.*wa_p);

				rct_yz_p=(0.5*(-gDDg_yz
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
						-Del_uy_yx*Del_ux_yz-Del_uz_yx*Del_ux_zz-Del_uz_yy*Del_uy_zz)*exp(-4.*wa_p);
				
				rct_yx_p=rct_xy_p;
				rct_zx_p=rct_xz_p;
				rct_zy_p=rct_yz_p;
				
				rc_xx_p=rct_xx_p+(-2.*wa_cdxx -2.*wa_lap*gxx_p -4.*wawa*gxx_p +4.*wa_x_p*wa_x_p)*exp(-4.*wa_p);
				
				rc_yy_p=rct_yy_p+(-2.*wa_cdyy -2.*wa_lap*gyy_p -4.*wawa*gyy_p +4.*wa_y_p*wa_y_p)*exp(-4.*wa_p);

				rc_zz_p=rct_zz_p+(-2.*wa_cdzz -2.*wa_lap*gzz_p -4.*wawa*gzz_p +4.*wa_z_p*wa_z_p)*exp(-4.*wa_p);

				rc_xy_p=rct_yx_p+(-2.*wa_cdxy -2.*wa_lap*gxy_p -4.*wawa*gxy_p +4.*wa_x_p*wa_y_p)*exp(-4.*wa_p);

				rc_xz_p=rct_xz_p+(-2.*wa_cdxz -2.*wa_lap*gxz_p -4.*wawa*gxz_p +4.*wa_x_p*wa_z_p)*exp(-4.*wa_p);

				rc_yz_p=rct_yz_p+(-2.*wa_cdyz -2.*wa_lap*gyz_p -4.*wawa*gyz_p +4.*wa_y_p*wa_z_p)*exp(-4.*wa_p);
				
				rc_yx_p=rc_xy_p;
				rc_zx_p=rc_xz_p;
				rc_zy_p=rc_yz_p;
				
				// R 
				ricci=gixx_p*rc_xx_p +gixy_p*rc_xy_p +gixz_p*rc_xz_p 
						+giyx_p*rc_yx_p +giyy_p*rc_yy_p +giyz_p*rc_yz_p 
						+gizx_p*rc_zx_p +gizy_p*rc_zy_p +gizz_p*rc_zz_p;
				
				if(fluidevo)
				{
					double sqgam=sqrt(det_p)*exp(6.*wa_p);
					
					rho+=get_bv(l,k,j,24)/sqgam;
					
					p_x+=get_bv(l,k,j,25)/sqgam;
					p_y+=get_bv(l,k,j,26)/sqgam;
					p_z+=get_bv(l,k,j,27)/sqgam;
				}
				
				if(scalarevo)				
				{
					double phi=get_bv(l,k,j,nsc);
					double Pi=get_bv(l,k,j,nscp);
					
					double phi_x=get_f_x(l,k,j,nsc);
					double phi_y=get_f_y(l,k,j,nsc);
					double phi_z=get_f_z(l,k,j,nsc);
					
					// (\tilde D_i \phi)(\tilde D^i \phi)
					double dphidphi=gixx_p*phi_x*phi_x +gixy_p*phi_x*phi_y +gixz_p*phi_x*phi_z
						+giyx_p*phi_y*phi_x +giyy_p*phi_y*phi_y +giyz_p*phi_y*phi_z
						+gizx_p*phi_z*phi_x +gizy_p*phi_z*phi_y +gizz_p*phi_z*phi_z;
					
					double Vpot=funcV(phi);
					
					rho+=(0.5*(pow(Pi,2)+exp(-4.*wa_p)*dphidphi)+Vpot);
					p_x+=Pi*phi_x;
					p_y+=Pi*phi_y;
					p_z+=Pi*phi_z;
				}
				
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
				
				M_x=Da_x+6.*(akx_ux_p*wa_x_p+akx_uy_p*wa_y_p+akx_uz_p*wa_z_p)-2./3.*ek_x_p-pi8*p_x;
				M_y=Da_y+6.*(aky_ux_p*wa_x_p+aky_uy_p*wa_y_p+aky_uz_p*wa_z_p)-2./3.*ek_y_p-pi8*p_y;
				M_z=Da_z+6.*(akz_ux_p*wa_x_p+akz_uy_p*wa_y_p+akz_uz_p*wa_z_p)-2./3.*ek_z_p-pi8*p_z;
				
				M_ux=gixx_p*M_x+gixy_p*M_y+gixz_p*M_z;
				M_uy=giyx_p*M_x+giyy_p*M_y+giyz_p*M_z;
				M_uz=gizx_p*M_x+gizy_p*M_y+gizz_p*M_z;
				
				//\cal D_i w^j
				//Dw_1_u2=w2_1_p (+Gam_u2_11*w1_p if 1==2);
				Dw_x_ux=wx_x_p+Gam_ux_xx*wx_p;
				Dw_x_uy=wy_x_p;
				Dw_x_uz=wz_x_p;
				
				Dw_y_ux=wx_y_p;
				Dw_y_uy=wy_y_p+Gam_uy_yy*wy_p;
				Dw_y_uz=wz_y_p;
				
				Dw_z_ux=wx_z_p;
				Dw_z_uy=wy_z_p;
				Dw_z_uz=wz_z_p+Gam_uz_zz*wz_p;
			
				//Delta_uij_k=\tilde \gamma^jl Delta_ui_lk
				//Del_u12_3=gi2x_p*Del_u1_x3+gi2y_p*Del_u1_y3+gi2z_p*Del_u1_z3;
				Del_uxx_x=gixx_p*Del_ux_xx+gixy_p*Del_ux_yx+gixz_p*Del_ux_zx;
				Del_uxx_y=gixx_p*Del_ux_xy+gixy_p*Del_ux_yy+gixz_p*Del_ux_zy;
				Del_uxx_z=gixx_p*Del_ux_xz+gixy_p*Del_ux_yz+gixz_p*Del_ux_zz;
				
				Del_uxy_x=giyx_p*Del_ux_xx+giyy_p*Del_ux_yx+giyz_p*Del_ux_zx;
				Del_uxy_y=giyx_p*Del_ux_xy+giyy_p*Del_ux_yy+giyz_p*Del_ux_zy;
				Del_uxy_z=giyx_p*Del_ux_xz+giyy_p*Del_ux_yz+giyz_p*Del_ux_zz;
				
				Del_uxz_x=gizx_p*Del_ux_xx+gizy_p*Del_ux_yx+gizz_p*Del_ux_zx;
				Del_uxz_y=gizx_p*Del_ux_xy+gizy_p*Del_ux_yy+gizz_p*Del_ux_zy;
				Del_uxz_z=gizx_p*Del_ux_xz+gizy_p*Del_ux_yz+gizz_p*Del_ux_zz;
				
				Del_uyx_x=gixx_p*Del_uy_xx+gixy_p*Del_uy_yx+gixz_p*Del_uy_zx;
				Del_uyx_y=gixx_p*Del_uy_xy+gixy_p*Del_uy_yy+gixz_p*Del_uy_zy;
				Del_uyx_z=gixx_p*Del_uy_xz+gixy_p*Del_uy_yz+gixz_p*Del_uy_zz;
				
				Del_uyy_x=giyx_p*Del_uy_xx+giyy_p*Del_uy_yx+giyz_p*Del_uy_zx;
				Del_uyy_y=giyx_p*Del_uy_xy+giyy_p*Del_uy_yy+giyz_p*Del_uy_zy;
				Del_uyy_z=giyx_p*Del_uy_xz+giyy_p*Del_uy_yz+giyz_p*Del_uy_zz;
				
				Del_uyz_x=gizx_p*Del_uy_xx+gizy_p*Del_uy_yx+gizz_p*Del_uy_zx;
				Del_uyz_y=gizx_p*Del_uy_xy+gizy_p*Del_uy_yy+gizz_p*Del_uy_zy;
				Del_uyz_z=gizx_p*Del_uy_xz+gizy_p*Del_uy_yz+gizz_p*Del_uy_zz;
				
				Del_uzx_x=gixx_p*Del_uz_xx+gixy_p*Del_uz_yx+gixz_p*Del_uz_zx;
				Del_uzx_y=gixx_p*Del_uz_xy+gixy_p*Del_uz_yy+gixz_p*Del_uz_zy;
				Del_uzx_z=gixx_p*Del_uz_xz+gixy_p*Del_uz_yz+gixz_p*Del_uz_zz;
				
				Del_uzy_x=giyx_p*Del_uz_xx+giyy_p*Del_uz_yx+giyz_p*Del_uz_zx;
				Del_uzy_y=giyx_p*Del_uz_xy+giyy_p*Del_uz_yy+giyz_p*Del_uz_zy;
				Del_uzy_z=giyx_p*Del_uz_xz+giyy_p*Del_uz_yz+giyz_p*Del_uz_zz;
				
				Del_uzz_x=gizx_p*Del_uz_xx+gizy_p*Del_uz_yx+gizz_p*Del_uz_zx;
				Del_uzz_y=gizx_p*Del_uz_xy+gizy_p*Del_uz_yy+gizz_p*Del_uz_zy;
				Del_uzz_z=gizx_p*Del_uz_xz+gizy_p*Del_uz_yz+gizz_p*Del_uz_zz;
				
				//\tilde D^i w^j
				//Dw_u12=gi1x_p*Dw_x_u2+gi1y_p*Dw_y_u2+gi1z_p*Dw_z_u2+Del_u21_x*wx_p+Del_u21_y*wy_p+Del_u21_z*wz_p;
				Dw_uxx=gixx_p*Dw_x_ux+gixy_p*Dw_y_ux+gixz_p*Dw_z_ux+Del_uxx_x*wx_p+Del_uxx_y*wy_p+Del_uxx_z*wz_p;
				Dw_uxy=gixx_p*Dw_x_uy+gixy_p*Dw_y_uy+gixz_p*Dw_z_uy+Del_uyx_x*wx_p+Del_uyx_y*wy_p+Del_uyx_z*wz_p;
				Dw_uxz=gixx_p*Dw_x_uz+gixy_p*Dw_y_uz+gixz_p*Dw_z_uz+Del_uzx_x*wx_p+Del_uzx_y*wy_p+Del_uzx_z*wz_p;
				
				Dw_uyx=giyx_p*Dw_x_ux+giyy_p*Dw_y_ux+giyz_p*Dw_z_ux+Del_uxy_x*wx_p+Del_uxy_y*wy_p+Del_uxy_z*wz_p;
				Dw_uyy=giyx_p*Dw_x_uy+giyy_p*Dw_y_uy+giyz_p*Dw_z_uy+Del_uyy_x*wx_p+Del_uyy_y*wy_p+Del_uyy_z*wz_p;
				Dw_uyz=giyx_p*Dw_x_uz+giyy_p*Dw_y_uz+giyz_p*Dw_z_uz+Del_uzy_x*wx_p+Del_uzy_y*wy_p+Del_uzy_z*wz_p;
				
				Dw_uzx=gizx_p*Dw_x_ux+gizy_p*Dw_y_ux+gizz_p*Dw_z_ux+Del_uxz_x*wx_p+Del_uxz_y*wy_p+Del_uxz_z*wz_p;
				Dw_uzy=gizx_p*Dw_x_uy+gizy_p*Dw_y_uy+gizz_p*Dw_z_uy+Del_uyz_x*wx_p+Del_uyz_y*wy_p+Del_uyz_z*wz_p;
				Dw_uzz=gizx_p*Dw_x_uz+gizy_p*Dw_y_uz+gizz_p*Dw_z_uz+Del_uzz_x*wx_p+Del_uzz_y*wy_p+Del_uzz_z*wz_p;
				
				divDw=Dw_x_ux+Dw_y_uy+Dw_z_uz;
				
				//(\tilde Lw)^ij
				//Lw_u12=Dw_u12+Dw_u21-2./3.*gi12*(Dw_x_ux+Dw_x_ux+Dw_x_ux);
				Lw_uxx=Dw_uxx+Dw_uxx-2./3.*gixx_p*divDw;
				Lw_uxy=Dw_uxy+Dw_uyx-2./3.*gixy_p*divDw;
				Lw_uxz=Dw_uxz+Dw_uzx-2./3.*gixz_p*divDw;
				
				Lw_uyx=Dw_uyx+Dw_uxy-2./3.*giyx_p*divDw;
				Lw_uyy=Dw_uyy+Dw_uyy-2./3.*giyy_p*divDw;
				Lw_uyz=Dw_uyz+Dw_uzy-2./3.*giyz_p*divDw;
				
				Lw_uzx=Dw_uzx+Dw_uxz-2./3.*gizx_p*divDw;
				Lw_uzy=Dw_uzy+Dw_uyz-2./3.*gizy_p*divDw;
				Lw_uzz=Dw_uzz+Dw_uzz-2./3.*gizz_p*divDw;
				
				//\tilde gam^jk \cal D_j \calD_k w^i
				wx_lap=gixx_p*wx_xx +gixy_p*wx_xy +gixz_p*wx_xz 
						+giyx_p*wx_yx +giyy_p*wx_yy +giyz_p*wx_yz 
						+gizx_p*wx_zx +gizy_p*wx_zy +gizz_p*wx_zz
						
						+gixx_p*delG_x_ux_xx*wx_p
						+Gam_ux_xx*gixx_p*wx_x_p+Gam_ux_xx*gixy_p*wx_y_p+Gam_ux_xx*gixz_p*wx_z_p
						-Gam_uy_yy*giyy_p*wx_y_p
						-Gam_uz_zz*gizz_p*wx_z_p;

				wy_lap=gixx_p*wy_xx +gixy_p*wy_xy +gixz_p*wy_xz 
						+giyx_p*wy_yx +giyy_p*wy_yy +giyz_p*wy_yz 
						+gizx_p*wy_zx +gizy_p*wy_zy +gizz_p*wy_zz
						
						+giyy_p*delG_y_uy_yy*wy_p
						+Gam_uy_yy*giyx_p*wy_x_p+Gam_uy_yy*giyy_p*wy_y_p+Gam_uy_yy*giyz_p*wy_z_p
						-Gam_ux_xx*gixx_p*wy_x_p-Gam_uz_zz*gizz_p*wy_z_p;
						
				wz_lap=gixx_p*wz_xx +gixy_p*wz_xy +gixz_p*wz_xz 
						+giyx_p*wz_yx +giyy_p*wz_yy +giyz_p*wz_yz 
						+gizx_p*wz_zx +gizy_p*wz_zy +gizz_p*wz_zz
						
						+gizz_p*delG_z_uz_zz*wz_p
						+Gam_uz_zz*gizx_p*wz_x_p+Gam_uz_zz*gizy_p*wz_y_p+Gam_uz_zz*gizz_p*wz_z_p
						-Gam_ux_xx*gixx_p*wz_x_p-Gam_uy_yy*giyy_p*wz_y_p;
				
				//\tilde D_i \tilde D_j w^j=\cal D_i \cal D_j w^j
				divw_x=wx_xx +wy_yx +wz_zx
						+delG_x_ux_xx*wx_p+Gam_ux_xx*wx_x_p+Gam_uy_yy*wy_x_p+Gam_uz_zz*wz_x_p;
				divw_y=wx_xy +wy_yy +wz_zy
						+delG_y_uy_yy*wy_p+Gam_ux_xx*wx_y_p+Gam_uy_yy*wy_y_p+Gam_uz_zz*wz_y_p;
				divw_z=wx_xz +wy_yz +wz_zz
						+delG_z_uz_zz*wz_p+Gam_ux_xx*wx_z_p+Gam_uy_yy*wy_z_p+Gam_uz_zz*wz_z_p;
				
				//DGw_12=\cal D_1 \tilde \gamma_23 w^3
				//DGw_12=Dgam_1_2x*wx_p+Dgam_1_2y*wy_p+Dgam_1_2z*wz_p
				
				DGw_xx=Dgam_x_xx*wx_p+Dgam_x_xy*wy_p+Dgam_x_xz*wz_p;
				DGw_xy=Dgam_x_yx*wx_p+Dgam_x_yy*wy_p+Dgam_x_yz*wz_p;
				DGw_xz=Dgam_x_zx*wx_p+Dgam_x_zy*wy_p+Dgam_x_zz*wz_p;
			
				DGw_yx=Dgam_y_xx*wx_p+Dgam_y_xy*wy_p+Dgam_y_xz*wz_p;
				DGw_yy=Dgam_y_yx*wx_p+Dgam_y_yy*wy_p+Dgam_y_yz*wz_p;
				DGw_yz=Dgam_y_zx*wx_p+Dgam_y_zy*wy_p+Dgam_y_zz*wz_p;
				
				DGw_zx=Dgam_z_xx*wx_p+Dgam_z_xy*wy_p+Dgam_z_xz*wz_p;
				DGw_zy=Dgam_z_yx*wx_p+Dgam_z_yy*wy_p+Dgam_z_yz*wz_p;
				DGw_zz=Dgam_z_zx*wx_p+Dgam_z_zy*wy_p+Dgam_z_zz*wz_p;
				
				//Lw_1_u2==g1x_p*Lw_ux2+g1y_p*Lw_uy2+g1z_p*Lw_uz2;
				Lw_x_ux=gxx_p*Lw_uxx+gxy_p*Lw_uyx+gxz_p*Lw_uzx;
				Lw_x_uy=gxx_p*Lw_uxy+gxy_p*Lw_uyy+gxz_p*Lw_uzy;
				Lw_x_uz=gxx_p*Lw_uxz+gxy_p*Lw_uyz+gxz_p*Lw_uzz;
				
				Lw_y_ux=gyx_p*Lw_uxx+gyy_p*Lw_uyx+gyz_p*Lw_uzx;
				Lw_y_uy=gyx_p*Lw_uxy+gyy_p*Lw_uyy+gyz_p*Lw_uzy;
				Lw_y_uz=gyx_p*Lw_uxz+gyy_p*Lw_uyz+gyz_p*Lw_uzz;
				
				Lw_z_ux=gzx_p*Lw_uxx+gzy_p*Lw_uyx+gzz_p*Lw_uzx;
				Lw_z_uy=gzx_p*Lw_uxy+gzy_p*Lw_uyy+gzz_p*Lw_uzy;
				Lw_z_uz=gzx_p*Lw_uxz+gzy_p*Lw_uyz+gzz_p*Lw_uzz;
				
				//Lw_12=g2x_p*Lw_1_ux+g2y_p*Lw_1_uy+g2z_p*Lw_1_uz;
				Lw_xx=gxx_p*Lw_x_ux+gxy_p*Lw_x_uy+gxz_p*Lw_x_uz;
				Lw_xy=gyx_p*Lw_x_ux+gyy_p*Lw_x_uy+gyz_p*Lw_x_uz;
				Lw_xz=gzx_p*Lw_x_ux+gzy_p*Lw_x_uy+gzz_p*Lw_x_uz;
				
				Lw_yx=Lw_xy;
				Lw_yy=gyx_p*Lw_y_ux+gyy_p*Lw_y_uy+gyz_p*Lw_y_uz;
				Lw_yz=gzx_p*Lw_y_ux+gzy_p*Lw_y_uy+gzz_p*Lw_y_uz;
				
				Lw_zx=Lw_xz;
				Lw_zy=Lw_yz;
				Lw_zz=gzx_p*Lw_z_ux+gzy_p*Lw_z_uy+gzz_p*Lw_z_uz;
				
				M_ux=wx_lap
				
					-zgx_p*Dw_x_ux-zgy_p*Dw_y_ux-zgz_p*Dw_z_ux
					
					+2.*(
					//Delta_ux1_2*Dw_1_u2
					Del_uxx_x*Dw_x_ux
					+Del_uxx_y*Dw_x_uy
					+Del_uxx_z*Dw_x_uz
					
					+Del_uxy_x*Dw_y_ux
					+Del_uxy_y*Dw_y_uy
					+Del_uxy_z*Dw_y_uz
				
					+Del_uxz_x*Dw_z_ux
					+Del_uxz_y*Dw_z_uy
					+Del_uxz_z*Dw_z_uz
					)
					
					-0.5*(
					//gix1_p*DGw_21*zg2_p
					gixx_p*DGw_xx*zgx_p
					+gixx_p*DGw_yx*zgy_p
					+gixx_p*DGw_zx*zgz_p
					
					+gixy_p*DGw_xy*zgx_p
					+gixy_p*DGw_yy*zgy_p
					+gixy_p*DGw_zy*zgz_p
					
					+gixz_p*DGw_xz*zgx_p
					+gixz_p*DGw_yz*zgy_p
					+gixz_p*DGw_zz*zgz_p
					)
					
					+0.5*(
					//DGam_1_ux*w1_p
					DGam_x_ux*wx_p+DGam_y_ux*wy_p+DGam_z_ux*wz_p
					)
					
					-0.5*(
					//gix1_p*g23_p*DGam_1_u3*w2_p
					gixx_p*gxx_p*DGam_x_ux*wx_p
					+gixx_p*gxy_p*DGam_x_uy*wx_p
					+gixx_p*gxz_p*DGam_x_uz*wx_p
					
					+gixx_p*gyx_p*DGam_x_ux*wy_p
					+gixx_p*gyy_p*DGam_x_uy*wy_p
					+gixx_p*gyz_p*DGam_x_uz*wy_p
					
					+gixx_p*gzx_p*DGam_x_ux*wz_p
					+gixx_p*gzy_p*DGam_x_uy*wz_p
					+gixx_p*gzz_p*DGam_x_uz*wz_p
					
					+gixy_p*gxx_p*DGam_y_ux*wx_p
					+gixy_p*gxy_p*DGam_y_uy*wx_p
					+gixy_p*gxz_p*DGam_y_uz*wx_p
					
					+gixy_p*gyx_p*DGam_y_ux*wy_p
					+gixy_p*gyy_p*DGam_y_uy*wy_p
					+gixy_p*gyz_p*DGam_y_uz*wy_p
					
					+gixy_p*gzx_p*DGam_y_ux*wz_p
					+gixy_p*gzy_p*DGam_y_uy*wz_p
					+gixy_p*gzz_p*DGam_y_uz*wz_p
					
					+gixz_p*gxx_p*DGam_z_ux*wx_p
					+gixz_p*gxy_p*DGam_z_uy*wx_p
					+gixz_p*gxz_p*DGam_z_uz*wx_p
					
					+gixz_p*gyx_p*DGam_z_ux*wy_p
					+gixz_p*gyy_p*DGam_z_uy*wy_p
					+gixz_p*gyz_p*DGam_z_uz*wy_p
					
					+gixz_p*gzx_p*DGam_z_ux*wz_p
					+gixz_p*gzy_p*DGam_z_uy*wz_p
					+gixz_p*gzz_p*DGam_z_uz*wz_p
					)
					
					+0.5*(
					//gix1*gDDg_21*w2_p
					gixx_p*gDDg_xx*wx_p
					+gixx_p*gDDg_xy*wy_p
					+gixx_p*gDDg_xz*wz_p
					
					+gixy_p*gDDg_xy*wx_p
					+gixy_p*gDDg_yy*wy_p
					+gixy_p*gDDg_yz*wz_p
					
					+gixz_p*gDDg_xz*wx_p
					+gixz_p*gDDg_yz*wy_p
					+gixz_p*gDDg_zz*wz_p
					)
					
					+
					//Del_u12_3*Del_ux_12*w3_p
					Del_uxx_x*Del_ux_xx*wx_p
					+Del_uxx_y*Del_ux_xx*wy_p
					+Del_uxx_z*Del_ux_xx*wz_p
					
					+Del_uxy_x*Del_ux_xy*wx_p
					+Del_uxy_y*Del_ux_xy*wy_p
					+Del_uxy_z*Del_ux_xy*wz_p
					
					+Del_uxz_x*Del_ux_xz*wx_p
					+Del_uxz_y*Del_ux_xz*wy_p
					+Del_uxz_z*Del_ux_xz*wz_p
					
					+Del_uyx_x*Del_ux_yx*wx_p
					+Del_uyx_y*Del_ux_yx*wy_p
					+Del_uyx_z*Del_ux_yx*wz_p
					
					+Del_uyy_x*Del_ux_yy*wx_p
					+Del_uyy_y*Del_ux_yy*wy_p
					+Del_uyy_z*Del_ux_yy*wz_p
					
					+Del_uyz_x*Del_ux_yz*wx_p
					+Del_uyz_y*Del_ux_yz*wy_p
					+Del_uyz_z*Del_ux_yz*wz_p
					
					+Del_uzx_x*Del_ux_zx*wx_p
					+Del_uzx_y*Del_ux_zx*wy_p
					+Del_uzx_z*Del_ux_zx*wz_p
					
					+Del_uzy_x*Del_ux_zy*wx_p
					+Del_uzy_y*Del_ux_zy*wy_p
					+Del_uzy_z*Del_ux_zy*wz_p
					
					+Del_uzz_x*Del_ux_zz*wx_p
					+Del_uzz_y*Del_ux_zz*wy_p
					+Del_uzz_z*Del_ux_zz*wz_p
					
					-(
					//gix1_p*Del_u23_1*DGw_32
					gixx_p*Del_uxx_x*DGw_xx
					+gixx_p*Del_uxy_x*DGw_yx
					+gixx_p*Del_uxz_x*DGw_zx
					
					+gixx_p*Del_uyx_x*DGw_xy
					+gixx_p*Del_uyy_x*DGw_yy
					+gixx_p*Del_uyz_x*DGw_zy
					
					+gixx_p*Del_uzx_x*DGw_xz
					+gixx_p*Del_uzy_x*DGw_yz
					+gixx_p*Del_uzz_x*DGw_zz
					
					+gixy_p*Del_uxx_y*DGw_xx
					+gixy_p*Del_uxy_y*DGw_yx
					+gixy_p*Del_uxz_y*DGw_zx
					
					+gixy_p*Del_uyx_y*DGw_xy
					+gixy_p*Del_uyy_y*DGw_yy
					+gixy_p*Del_uyz_y*DGw_zy
					
					+gixy_p*Del_uzx_y*DGw_xz
					+gixy_p*Del_uzy_y*DGw_yz
					+gixy_p*Del_uzz_y*DGw_zz
					
					+gixz_p*Del_uxx_z*DGw_xx
					+gixz_p*Del_uxy_z*DGw_yx
					+gixz_p*Del_uxz_z*DGw_zx
					
					+gixz_p*Del_uyx_z*DGw_xy
					+gixz_p*Del_uyy_z*DGw_yy
					+gixz_p*Del_uyz_z*DGw_zy
					
					+gixz_p*Del_uzx_z*DGw_xz
					+gixz_p*Del_uzy_z*DGw_yz
					+gixz_p*Del_uzz_z*DGw_zz
					)
					
					+1./3.*(gixx_p*divw_x+gixy_p*divw_y+gixz_p*divw_z)
					
					+exp(4.*wa_p)*(
					//gix1_p*rc_12_p*w2_p
					gixx_p*rct_xx_p*wx_p
					+gixx_p*rct_xy_p*wy_p
					+gixx_p*rct_xz_p*wz_p
					
					+gixy_p*rct_yx_p*wx_p
					+gixy_p*rct_yy_p*wy_p
					+gixy_p*rct_yz_p*wz_p
					
					+gixz_p*rct_zx_p*wx_p
					+gixz_p*rct_zy_p*wy_p
					+gixz_p*rct_zz_p*wz_p
					)
					
					+6.*(
					//Lw_ux1*wa_1_p
					Lw_uxx*wa_x_p
					+Lw_uxy*wa_y_p
					+Lw_uxz*wa_z_p
					)
					
					+M_ux;
					
				M_uy=wy_lap
				
					-zgx_p*Dw_x_uy-zgy_p*Dw_y_uy-zgz_p*Dw_z_uy
					
					+2.*(
					//Delta_ux1_2*Dw_1_u2
					Del_uyx_x*Dw_x_ux
					+Del_uyx_y*Dw_x_uy
					+Del_uyx_z*Dw_x_uz
					
					+Del_uyy_x*Dw_y_ux
					+Del_uyy_y*Dw_y_uy
					+Del_uyy_z*Dw_y_uz
				
					+Del_uyz_x*Dw_z_ux
					+Del_uyz_y*Dw_z_uy
					+Del_uyz_z*Dw_z_uz
					)
					
					-0.5*(
					//gix1_p*DGw_21*zg2_p
					giyx_p*DGw_xx*zgx_p
					+giyx_p*DGw_yx*zgy_p
					+giyx_p*DGw_zx*zgz_p
					
					+giyy_p*DGw_xy*zgx_p
					+giyy_p*DGw_yy*zgy_p
					+giyy_p*DGw_zy*zgz_p
					
					+giyz_p*DGw_xz*zgx_p
					+giyz_p*DGw_yz*zgy_p
					+giyz_p*DGw_zz*zgz_p
					)
					
					+0.5*(
					//DGam_1_ux*w1_p
					DGam_x_uy*wx_p+DGam_y_uy*wy_p+DGam_z_uy*wz_p
					)
					
					-0.5*(
					//gix1_p*g23_p*DGam_1_u3*w2_p
					giyx_p*gxx_p*DGam_x_ux*wx_p
					+giyx_p*gxy_p*DGam_x_uy*wx_p
					+giyx_p*gxz_p*DGam_x_uz*wx_p
					
					+giyx_p*gyx_p*DGam_x_ux*wy_p
					+giyx_p*gyy_p*DGam_x_uy*wy_p
					+giyx_p*gyz_p*DGam_x_uz*wy_p
					
					+giyx_p*gzx_p*DGam_x_ux*wz_p
					+giyx_p*gzy_p*DGam_x_uy*wz_p
					+giyx_p*gzz_p*DGam_x_uz*wz_p
					
					+giyy_p*gxx_p*DGam_y_ux*wx_p
					+giyy_p*gxy_p*DGam_y_uy*wx_p
					+giyy_p*gxz_p*DGam_y_uz*wx_p
					
					+giyy_p*gyx_p*DGam_y_ux*wy_p
					+giyy_p*gyy_p*DGam_y_uy*wy_p
					+giyy_p*gyz_p*DGam_y_uz*wy_p
					
					+giyy_p*gzx_p*DGam_y_ux*wz_p
					+giyy_p*gzy_p*DGam_y_uy*wz_p
					+giyy_p*gzz_p*DGam_y_uz*wz_p
					
					+giyz_p*gxx_p*DGam_z_ux*wx_p
					+giyz_p*gxy_p*DGam_z_uy*wx_p
					+giyz_p*gxz_p*DGam_z_uz*wx_p
					
					+giyz_p*gyx_p*DGam_z_ux*wy_p
					+giyz_p*gyy_p*DGam_z_uy*wy_p
					+giyz_p*gyz_p*DGam_z_uz*wy_p
					
					+giyz_p*gzx_p*DGam_z_ux*wz_p
					+giyz_p*gzy_p*DGam_z_uy*wz_p
					+giyz_p*gzz_p*DGam_z_uz*wz_p
					)
					
					+0.5*(
					//gix1*gDDg_21*w2_p
					giyx_p*gDDg_xx*wx_p
					+giyx_p*gDDg_xy*wy_p
					+giyx_p*gDDg_xz*wz_p
					
					+giyy_p*gDDg_xy*wx_p
					+giyy_p*gDDg_yy*wy_p
					+giyy_p*gDDg_yz*wz_p
					
					+giyz_p*gDDg_xz*wx_p
					+giyz_p*gDDg_yz*wy_p
					+giyz_p*gDDg_zz*wz_p
					)
					
					+
					//Del_u12_3*Del_ux_12*w3_p
					Del_uxx_x*Del_uy_xx*wx_p
					+Del_uxx_y*Del_uy_xx*wy_p
					+Del_uxx_z*Del_uy_xx*wz_p
					
					+Del_uxy_x*Del_uy_xy*wx_p
					+Del_uxy_y*Del_uy_xy*wy_p
					+Del_uxy_z*Del_uy_xy*wz_p
					
					+Del_uxz_x*Del_uy_xz*wx_p
					+Del_uxz_y*Del_uy_xz*wy_p
					+Del_uxz_z*Del_uy_xz*wz_p
					
					+Del_uyx_x*Del_uy_yx*wx_p
					+Del_uyx_y*Del_uy_yx*wy_p
					+Del_uyx_z*Del_uy_yx*wz_p
					
					+Del_uyy_x*Del_uy_yy*wx_p
					+Del_uyy_y*Del_uy_yy*wy_p
					+Del_uyy_z*Del_uy_yy*wz_p
					
					+Del_uyz_x*Del_uy_yz*wx_p
					+Del_uyz_y*Del_uy_yz*wy_p
					+Del_uyz_z*Del_uy_yz*wz_p
					
					+Del_uzx_x*Del_uy_zx*wx_p
					+Del_uzx_y*Del_uy_zx*wy_p
					+Del_uzx_z*Del_uy_zx*wz_p
					
					+Del_uzy_x*Del_uy_zy*wx_p
					+Del_uzy_y*Del_uy_zy*wy_p
					+Del_uzy_z*Del_uy_zy*wz_p
					
					+Del_uzz_x*Del_uy_zz*wx_p
					+Del_uzz_y*Del_uy_zz*wy_p
					+Del_uzz_z*Del_uy_zz*wz_p
					
					-(
					//gix1_p*Del_u23_1*DGw_32
					giyx_p*Del_uxx_x*DGw_xx
					+giyx_p*Del_uxy_x*DGw_yx
					+giyx_p*Del_uxz_x*DGw_zx
					
					+giyx_p*Del_uyx_x*DGw_xy
					+giyx_p*Del_uyy_x*DGw_yy
					+giyx_p*Del_uyz_x*DGw_zy
					
					+giyx_p*Del_uzx_x*DGw_xz
					+giyx_p*Del_uzy_x*DGw_yz
					+giyx_p*Del_uzz_x*DGw_zz
					
					+giyy_p*Del_uxx_y*DGw_xx
					+giyy_p*Del_uxy_y*DGw_yx
					+giyy_p*Del_uxz_y*DGw_zx
					
					+giyy_p*Del_uyx_y*DGw_xy
					+giyy_p*Del_uyy_y*DGw_yy
					+giyy_p*Del_uyz_y*DGw_zy
					
					+giyy_p*Del_uzx_y*DGw_xz
					+giyy_p*Del_uzy_y*DGw_yz
					+giyy_p*Del_uzz_y*DGw_zz
					
					+giyz_p*Del_uxx_z*DGw_xx
					+giyz_p*Del_uxy_z*DGw_yx
					+giyz_p*Del_uxz_z*DGw_zx
					
					+giyz_p*Del_uyx_z*DGw_xy
					+giyz_p*Del_uyy_z*DGw_yy
					+giyz_p*Del_uyz_z*DGw_zy
					
					+giyz_p*Del_uzx_z*DGw_xz
					+giyz_p*Del_uzy_z*DGw_yz
					+giyz_p*Del_uzz_z*DGw_zz
					)
					
					+1./3.*(giyx_p*divw_x+giyy_p*divw_y+giyz_p*divw_z)
					
					+exp(4.*wa_p)*(
					//gix1_p*rc_12_p*w2_p
					giyx_p*rct_xx_p*wx_p
					+giyx_p*rct_xy_p*wy_p
					+giyx_p*rct_xz_p*wz_p
					
					+giyy_p*rct_yx_p*wx_p
					+giyy_p*rct_yy_p*wy_p
					+giyy_p*rct_yz_p*wz_p
					
					+giyz_p*rct_zx_p*wx_p
					+giyz_p*rct_zy_p*wy_p
					+giyz_p*rct_zz_p*wz_p
					)
					
					+6.*(
					//Lw_ux1*wa_1_p
					Lw_uyx*wa_x_p
					+Lw_uyy*wa_y_p
					+Lw_uyz*wa_z_p
					)

					+M_uy;
					
				M_uz=wz_lap
				
					-zgx_p*Dw_x_uz-zgy_p*Dw_y_uz-zgz_p*Dw_z_uz
					
					+2.*(
					//Delta_ux1_2*Dw_1_u2
					Del_uzx_x*Dw_x_ux
					+Del_uzx_y*Dw_x_uy
					+Del_uzx_z*Dw_x_uz
					
					+Del_uzy_x*Dw_y_ux
					+Del_uzy_y*Dw_y_uy
					+Del_uzy_z*Dw_y_uz
				
					+Del_uzz_x*Dw_z_ux
					+Del_uzz_y*Dw_z_uy
					+Del_uzz_z*Dw_z_uz
					)
					
					-0.5*(
					//gix1_p*DGw_21*zg2_p
					gizx_p*DGw_xx*zgx_p
					+gizx_p*DGw_yx*zgy_p
					+gizx_p*DGw_zx*zgz_p
					
					+gizy_p*DGw_xy*zgx_p
					+gizy_p*DGw_yy*zgy_p
					+gizy_p*DGw_zy*zgz_p
					
					+gizz_p*DGw_xz*zgx_p
					+gizz_p*DGw_yz*zgy_p
					+gizz_p*DGw_zz*zgz_p
					)
					
					+0.5*(
					//DGam_1_ux*w1_p
					DGam_x_uz*wx_p+DGam_y_uz*wy_p+DGam_z_uz*wz_p
					)
					
					-0.5*(
					//gix1_p*g23_p*DGam_1_u3*w2_p
					gizx_p*gxx_p*DGam_x_ux*wx_p
					+gizx_p*gxy_p*DGam_x_uy*wx_p
					+gizx_p*gxz_p*DGam_x_uz*wx_p
					
					+gizx_p*gyx_p*DGam_x_ux*wy_p
					+gizx_p*gyy_p*DGam_x_uy*wy_p
					+gizx_p*gyz_p*DGam_x_uz*wy_p
					
					+gizx_p*gzx_p*DGam_x_ux*wz_p
					+gizx_p*gzy_p*DGam_x_uy*wz_p
					+gizx_p*gzz_p*DGam_x_uz*wz_p
					
					+gizy_p*gxx_p*DGam_y_ux*wx_p
					+gizy_p*gxy_p*DGam_y_uy*wx_p
					+gizy_p*gxz_p*DGam_y_uz*wx_p
					
					+gizy_p*gyx_p*DGam_y_ux*wy_p
					+gizy_p*gyy_p*DGam_y_uy*wy_p
					+gizy_p*gyz_p*DGam_y_uz*wy_p
					
					+gizy_p*gzx_p*DGam_y_ux*wz_p
					+gizy_p*gzy_p*DGam_y_uy*wz_p
					+gizy_p*gzz_p*DGam_y_uz*wz_p
					
					+gizz_p*gxx_p*DGam_z_ux*wx_p
					+gizz_p*gxy_p*DGam_z_uy*wx_p
					+gizz_p*gxz_p*DGam_z_uz*wx_p
					
					+gizz_p*gyx_p*DGam_z_ux*wy_p
					+gizz_p*gyy_p*DGam_z_uy*wy_p
					+gizz_p*gyz_p*DGam_z_uz*wy_p
					
					+gizz_p*gzx_p*DGam_z_ux*wz_p
					+gizz_p*gzy_p*DGam_z_uy*wz_p
					+gizz_p*gzz_p*DGam_z_uz*wz_p
					)
					
					+0.5*(
					//gix1*gDDg_21*w2_p
					gizx_p*gDDg_xx*wx_p
					+gizx_p*gDDg_xy*wy_p
					+gizx_p*gDDg_xz*wz_p
					
					+gizy_p*gDDg_xy*wx_p
					+gizy_p*gDDg_yy*wy_p
					+gizy_p*gDDg_yz*wz_p
					
					+gizz_p*gDDg_xz*wx_p
					+gizz_p*gDDg_yz*wy_p
					+gizz_p*gDDg_zz*wz_p
					)
					
					+
					//Del_u12_3*Del_ux_12*w3_p
					Del_uxx_x*Del_uz_xx*wx_p
					+Del_uxx_y*Del_uz_xx*wy_p
					+Del_uxx_z*Del_uz_xx*wz_p
					
					+Del_uxy_x*Del_uz_xy*wx_p
					+Del_uxy_y*Del_uz_xy*wy_p
					+Del_uxy_z*Del_uz_xy*wz_p
					
					+Del_uxz_x*Del_uz_xz*wx_p
					+Del_uxz_y*Del_uz_xz*wy_p
					+Del_uxz_z*Del_uz_xz*wz_p
					
					+Del_uyx_x*Del_uz_yx*wx_p
					+Del_uyx_y*Del_uz_yx*wy_p
					+Del_uyx_z*Del_uz_yx*wz_p
					
					+Del_uyy_x*Del_uz_yy*wx_p
					+Del_uyy_y*Del_uz_yy*wy_p
					+Del_uyy_z*Del_uz_yy*wz_p
					
					+Del_uyz_x*Del_uz_yz*wx_p
					+Del_uyz_y*Del_uz_yz*wy_p
					+Del_uyz_z*Del_uz_yz*wz_p
					
					+Del_uzx_x*Del_uz_zx*wx_p
					+Del_uzx_y*Del_uz_zx*wy_p
					+Del_uzx_z*Del_uz_zx*wz_p
					
					+Del_uzy_x*Del_uz_zy*wx_p
					+Del_uzy_y*Del_uz_zy*wy_p
					+Del_uzy_z*Del_uz_zy*wz_p
					
					+Del_uzz_x*Del_uz_zz*wx_p
					+Del_uzz_y*Del_uz_zz*wy_p
					+Del_uzz_z*Del_uz_zz*wz_p
					
					-(
					//gix1_p*Del_u23_1*DGw_32
					gizx_p*Del_uxx_x*DGw_xx
					+gizx_p*Del_uxy_x*DGw_yx
					+gizx_p*Del_uxz_x*DGw_zx
					
					+gizx_p*Del_uyx_x*DGw_xy
					+gizx_p*Del_uyy_x*DGw_yy
					+gizx_p*Del_uyz_x*DGw_zy
					
					+gizx_p*Del_uzx_x*DGw_xz
					+gizx_p*Del_uzy_x*DGw_yz
					+gizx_p*Del_uzz_x*DGw_zz
					
					+gizy_p*Del_uxx_y*DGw_xx
					+gizy_p*Del_uxy_y*DGw_yx
					+gizy_p*Del_uxz_y*DGw_zx
					
					+gizy_p*Del_uyx_y*DGw_xy
					+gizy_p*Del_uyy_y*DGw_yy
					+gizy_p*Del_uyz_y*DGw_zy
					
					+gizy_p*Del_uzx_y*DGw_xz
					+gizy_p*Del_uzy_y*DGw_yz
					+gizy_p*Del_uzz_y*DGw_zz
					
					+gizz_p*Del_uxx_z*DGw_xx
					+gizz_p*Del_uxy_z*DGw_yx
					+gizz_p*Del_uxz_z*DGw_zx
					
					+gizz_p*Del_uyx_z*DGw_xy
					+gizz_p*Del_uyy_z*DGw_yy
					+gizz_p*Del_uyz_z*DGw_zy
					
					+gizz_p*Del_uzx_z*DGw_xz
					+gizz_p*Del_uzy_z*DGw_yz
					+gizz_p*Del_uzz_z*DGw_zz
					)
					
					+1./3.*(gizx_p*divw_x+gizy_p*divw_y+gizz_p*divw_z)
					
					+exp(4.*wa_p)*(
					//gix1_p*rc_12_p*w2_p
					gizx_p*rct_xx_p*wx_p
					+gizx_p*rct_xy_p*wy_p
					+gizx_p*rct_xz_p*wz_p
					
					+gizy_p*rct_yx_p*wx_p
					+gizy_p*rct_yy_p*wy_p
					+gizy_p*rct_yz_p*wz_p
					
					+gizz_p*rct_zx_p*wx_p
					+gizz_p*rct_zy_p*wy_p
					+gizz_p*rct_zz_p*wz_p
					)
					
					+6.*(
					//Lw_ux1*wa_1_p
					Lw_uzx*wa_x_p
					+Lw_uzy*wa_y_p
					+Lw_uzz*wa_z_p
					)

					+M_uz;
				
				//////////// updating \tilde A_ij ////////////////////////////////
				//ak12_p= ak12_p+Lw_12;											//
				akxx_p= akxx_p+Lw_xx;
				akxy_p= akxy_p+Lw_xy;
				akxz_p= akxz_p+Lw_xz;
				
				akyx_p= akyx_p+Lw_yx;
				akyy_p= akyy_p+Lw_yy;
				akyz_p= akyz_p+Lw_yz;
				
				akzx_p= akzx_p+Lw_zx;
				akzy_p= akzy_p+Lw_zy;
				akzz_p= akzz_p+Lw_zz;
				
				set_bv(l,k,j,14)=akxx_p;
				set_bv(l,k,j,15)=akyy_p;
				set_bv(l,k,j,16)=akzz_p;
				set_bv(l,k,j,17)=akxy_p;
				set_bv(l,k,j,18)=akxz_p;
				set_bv(l,k,j,19)=akyz_p;										//
				//////////// updating \tilde A_ij ////////////////////////////////
				
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
				
				//hamc = ricci + 2.*pow(ek_p+2.*th_p,2)/3. - aaaa - 2.*lambda - pi16*get_enemom(l,k,j,0);
				hamc = ricci  - aaaa - 2.*lambda-pi16*rho;
				
				double dfxm2=1./get_flat_df2x(j);
				double dfym2=1./get_flat_df2y(k);
				double dfzm2=1./get_flat_df2z(l);
				double ddfpdfx=get_flat_Gamx(j);
				double ddfpdfy=get_flat_Gamy(k);
				double ddfpdfz=get_flat_Gamz(l);
				
				double lhs=(wa_xx-ddfpdfx*wa_x_p)*dfxm2
							+(wa_yy-ddfpdfy*wa_y_p)*dfym2
							+(wa_zz-ddfpdfz*wa_z_p)*dfzm2;
				
				set_coef(l,k,j,5)=0.125*exp(4.*wa_p)*hamc;
				
				double k2P5=pow(get_bv(l,k,j,20),2)/12.*exp(5.*wa_p); 
				
				//modified by chulmoon on 2016.2.9 for non-Cartesian inhomogeneous grid
			 	fac=facl*fack*facj*sqrt(get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l));

				sumnume-=get_coef(l,k,j,5)*exp(wa_p)*fac;
				sumdeno+=k2P5*fac;

				set_coef(l,k,j,5)=get_coef(l,k,j,5)+lhs;
				
				
				lhs=(wx_xx-ddfpdfx*wx_x_p)*dfxm2
							+(wx_yy-ddfpdfy*wx_y_p)*dfym2
							+(wx_zz-ddfpdfz*wx_z_p)*dfzm2;
				
				set_coef(l,k,j,7)=-M_ux+lhs;
				
				lhs=(wy_xx-ddfpdfx*wy_x_p)*dfxm2
							+(wy_yy-ddfpdfy*wy_y_p)*dfym2
							+(wy_zz-ddfpdfz*wy_z_p)*dfzm2;
				
				set_coef(l,k,j,8)=-M_uy+lhs;
				
				lhs=(wz_xx-ddfpdfx*wz_x_p)*dfxm2
							+(wz_yy-ddfpdfy*wz_y_p)*dfym2
							+(wz_zz-ddfpdfz*wz_z_p)*dfzm2;
				
				set_coef(l,k,j,9)=-M_uz+lhs;
				
				
			}
		}
	}
	
	double ofac2=sumnume/sumdeno;
	
	if(ofac2<0.)
	{
		cout << "ofac2<0 ofac2=" << ofac2 << endl;
		exit(0);
	}
	
	double ofac=sqrt(ofac2);
	
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				set_bv(l,k,j,20)=get_bv(l,k,j,20)*ofac;
				set_coef(l,k,j,5)+=pow(get_bv(l,k,j,20),2)/12.*exp(4.*get_bv(l,k,j,13));
				set_coef(l,k,j,6)+=abs(pow(get_bv(l,k,j,20),2)/12.*exp(4.*get_bv(l,k,j,13)));
			}
		}
	}
	
	boundary_reflection_even(20);
	
	return;
}

//solving all constraints(hamiltonian and momentum constraints)
//SOR iteration for all variables
//cc: convergence condition
//acc: acceleration/deceleration parameter for SOR. acc<1 may work 
void Fmv::solveconst(int msstep,double cc,double acc)
{
	/////////////////////////// variables definition /////////////////////////////////
	double errm[4];									// maximum errors
	double convcond=cc;								// convergence condition
	double outputerr=0.1;							// output threshold
	bool conv=false;								// convergence or not
	double errmm=0.;								// maximum error among four componensts
	/////////////////////////// variables definition /////////////////////////////////
	
	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);
	
	//set source terms
	set_source();
	
	//SOR iteration
	for(int n=0;n<msstep;n++)
	{
		#pragma omp parallel sections
		{
		
		#pragma omp section
		{
			int nn=0;
		
			errm[nn]=0.;
		
			for(int l=lli;l<=lui;l++)
			{
				for(int k=kli;k<=kui;k++)
				{
					for(int j=jli;j<=jui;j++)
					{
						double wa_p=   get_bv(l,k,j,13);
						
						double wa_x_p=(    -(get_bv(l,k,j+2,13)-get_bv(l,k,j-2,13))
								+8.*(get_bv(l,k,j+1,13)-get_bv(l,k,j-1,13))
								)/12.;
						double wa_y_p=(   -(get_bv(l,k+2,j,13)-get_bv(l,k-2,j,13))
								+8.*(get_bv(l,k+1,j,13)-get_bv(l,k-1,j,13))
								)/12.;
						double wa_z_p=(   -(get_bv(l+2,k,j,13)-get_bv(l-2,k,j,13))
								+8.*(get_bv(l+1,k,j,13)-get_bv(l-1,k,j,13))
								)/12.;

						double wa_xx=(    -(    get_bv(l  ,k  ,j+2,13) + get_bv(l  ,k  ,j-2,13))
								+16.*( get_bv(l  ,k  ,j+1,13) + get_bv(l  ,k  ,j-1,13))
								-30.*  get_bv(l  ,k  ,j  ,13)
								)/12.;
						double wa_yy=(    -(    get_bv(l  ,k+2,j  ,13) + get_bv(l  ,k-2,j  ,13))
								+16.*( get_bv(l  ,k+1,j  ,13) + get_bv(l  ,k-1,j  ,13))
								-30.*  get_bv(l  ,k  ,j  ,13)
								)/12.;
						double wa_zz=(    -(    get_bv(l+2,k  ,j  ,13) + get_bv(l-2,k  ,j  ,13))
								+16.*( get_bv(l+1,k  ,j  ,13) + get_bv(l-1,k  ,j  ,13))
								-30.*  get_bv(l  ,k  ,j  ,13)
								)/12.;
						
						double dfxm2=1./get_flat_df2x(j);
						double dfym2=1./get_flat_df2y(k);
						double dfzm2=1./get_flat_df2z(l);
						double ddfpdfx=get_flat_Gamx(j);
						double ddfpdfy=get_flat_Gamy(k);
						double ddfpdfz=get_flat_Gamz(l);
						
						
						double dev=(
							dfxm2*wa_xx
							+dfym2*wa_yy
							+dfzm2*wa_zz
							-dfxm2*ddfpdfx*wa_x_p*dx
							-dfym2*ddfpdfy*wa_y_p*dx
							-dfzm2*ddfpdfz*wa_z_p*dx
							-get_coef(l,k,j,5)*dx*dx
							)/(2.5*(dfxm2+dfym2+dfzm2));

						double err=abs(dev);
						
						if(err>1.0e6)
						{
							cout << "large error in solveham" << endl;
							cout << "j=" << j << endl
							<< "k=" << k << endl 
							<< "l=" << l << endl
							<< "n=" << n << endl 
							<< "err=" << err << endl 
							<< "dev=" << dev << endl 
							//<< "norm=" << norm 
							<< endl;
							
							cout << " get_bv(l,k,j,13)=" << get_bv(l,k,j,13) << endl 
								<< " get_bv(l,k,j+1,13)=" << get_bv(l,k,j+1,13) << endl 
								<< " get_bv(l,k,j-1,13)=" << get_bv(l,k,j-1,13) << endl 
								<< " get_bv(l,k+1,j,13)=" << get_bv(l,k+1,j,13) << endl 
								<< " get_bv(l,k-1,j,13)=" << get_bv(l,k-1,j,13) << endl 
								<< " get_bv(l+1,k,j,13)=" << get_bv(l+1,k,j,13) << endl 
								<< " get_bv(l-1,k,j,13)=" << get_bv(l-1,k,j,13) << endl 
								<< " get_bv(l,k-1,j-1,13)=" << get_bv(l,k-1,j-1,13) << endl 
								<< " get_bv(l,k+1,j+1,13)=" << get_bv(l,k+1,j+1,13) << endl 
								<< " get_coef(l,k,j,5)=" << get_coef(l,k,j,5) << endl 
								<< endl;
						
							exit(1);
						}
						
						set_bv(l,k,j,13)=wa_p+acc*dev;

						if (errm[nn]<err)
						errm[nn]=err;
					}
				}
			}
		}
		
		#pragma omp section
		{
			int nn=1;
		
			errm[nn]=0.;
		
			for(int l=lli;l<=lui;l++)
			{
				for(int k=kli;k<=kui;k++)
				{
					for(int j=jli;j<=jui;j++)
					{
						double wa_x_p=(    -(get_wvec(l,k,j+2,0)-get_wvec(l,k,j-2,0))
								+8.*(get_wvec(l,k,j+1,0)-get_wvec(l,k,j-1,0))
								)/12.;
						double wa_y_p=(   -(get_wvec(l,k+2,j,0)-get_wvec(l,k-2,j,0))
								+8.*(get_wvec(l,k+1,j,0)-get_wvec(l,k-1,j,0))
								)/12.;
						double wa_z_p=(   -(get_wvec(l+2,k,j,0)-get_wvec(l-2,k,j,0))
								+8.*(get_wvec(l+1,k,j,0)-get_wvec(l-1,k,j,0))
								)/12.;

						double wa_xx=(    -(    get_wvec(l  ,k  ,j+2,0) + get_wvec(l  ,k  ,j-2,0))
								+16.*( get_wvec(l  ,k  ,j+1,0) + get_wvec(l  ,k  ,j-1,0))
								-30.*  get_wvec(l  ,k  ,j  ,0)
								)/12.;
						double wa_yy=(    -(    get_wvec(l  ,k+2,j  ,0) + get_wvec(l  ,k-2,j  ,0))
								+16.*( get_wvec(l  ,k+1,j  ,0) + get_wvec(l  ,k-1,j  ,0))
								-30.*  get_wvec(l  ,k  ,j  ,0)
								)/12.;
						double wa_zz=(    -(    get_wvec(l+2,k  ,j  ,0) + get_wvec(l-2,k  ,j  ,0))
								+16.*( get_wvec(l+1,k  ,j  ,0) + get_wvec(l-1,k  ,j  ,0))
								-30.*  get_wvec(l  ,k  ,j  ,0)
								)/12.;
						
						double dfxm2=1./get_flat_df2x(j);
						double dfym2=1./get_flat_df2y(k);
						double dfzm2=1./get_flat_df2z(l);
						double ddfpdfx=get_flat_Gamx(j);
						double ddfpdfy=get_flat_Gamy(k);
						double ddfpdfz=get_flat_Gamz(l);
						
						double dev=(
							dfxm2*wa_xx
							+dfym2*wa_yy
							+dfzm2*wa_zz
							-dfxm2*ddfpdfx*wa_x_p*dx
							-dfym2*ddfpdfy*wa_y_p*dx
							-dfzm2*ddfpdfz*wa_z_p*dx
							-get_coef(l,k,j,7)*dx*dx
							)/(2.5*(dfxm2+dfym2+dfzm2));
						
						double err=abs(dev);
						
						if(err>outputerr)
						{
							cout << "large error in solveham" << endl;
							cout << "j=" << j << endl
							<< "k=" << k << endl 
							<< "l=" << l << endl
							<< "n=" << n << endl 
							<< "err=" << err << endl 
							<< "dev=" << dev << endl 
							//<< "norm=" << norm 
							<< endl;
							
							cout << " get_wvec(l,k,j,0)=" << get_wvec(l,k,j,0) << endl 
								<< " get_wvec(l,k,j+1,0)=" << get_wvec(l,k,j+1,0) << endl 
								<< " get_wvec(l,k,j-1,0)=" << get_wvec(l,k,j-1,0) << endl 
								<< " get_wvec(l,k+1,j,0)=" << get_wvec(l,k+1,j,0) << endl 
								<< " get_wvec(l,k-1,j,0)=" << get_wvec(l,k-1,j,0) << endl 
								<< " get_wvec(l+1,k,j,0)=" << get_wvec(l+1,k,j,0) << endl 
								<< " get_wvec(l-1,k,j,0)=" << get_wvec(l-1,k,j,0) << endl 
								<< " get_wvec(l,k-1,j-1,0)=" << get_wvec(l,k-1,j-1,0) << endl 
								<< " get_wvec(l,k+1,j+1,0)=" << get_wvec(l,k+1,j+1,0) << endl 
								<< " get_coef(l,k,j,7)=" << get_coef(l,k,j,7) << endl 
								<< endl;
						
							exit(1);
						}
						
						
						set_wvec(l,k,j,0)=acc*dev;

						if (errm[nn]<err)
						errm[nn]=err;
					}
				}
			}
		}
		
		#pragma omp section
		{
			int nn=2;
		
			errm[nn]=0.;
		
			for(int l=lli;l<=lui;l++)
			{
				for(int k=kli;k<=kui;k++)
				{
					for(int j=jli;j<=jui;j++)
					{
						double wa_x_p=(    -(get_wvec(l,k,j+2,1)-get_wvec(l,k,j-2,1))
								+8.*(get_wvec(l,k,j+1,1)-get_wvec(l,k,j-1,1))
								)/12.;
						double wa_y_p=(   -(get_wvec(l,k+2,j,1)-get_wvec(l,k-2,j,1))
								+8.*(get_wvec(l,k+1,j,1)-get_wvec(l,k-1,j,1))
								)/12.;
						double wa_z_p=(   -(get_wvec(l+2,k,j,1)-get_wvec(l-2,k,j,1))
								+8.*(get_wvec(l+1,k,j,1)-get_wvec(l-1,k,j,1))
								)/12.;

						double wa_xx=(    -(    get_wvec(l  ,k  ,j+2,1) + get_wvec(l  ,k  ,j-2,1))
								+16.*( get_wvec(l  ,k  ,j+1,1) + get_wvec(l  ,k  ,j-1,1))
								-30.*  get_wvec(l  ,k  ,j  ,1)
								)/12.;
						double wa_yy=(    -(    get_wvec(l  ,k+2,j  ,1) + get_wvec(l  ,k-2,j  ,1))
								+16.*( get_wvec(l  ,k+1,j  ,1) + get_wvec(l  ,k-1,j  ,1))
								-30.*  get_wvec(l  ,k  ,j  ,1)
								)/12.;
						double wa_zz=(    -(    get_wvec(l+2,k  ,j  ,1) + get_wvec(l-2,k  ,j  ,1))
								+16.*( get_wvec(l+1,k  ,j  ,1) + get_wvec(l-1,k  ,j  ,1))
								-30.*  get_wvec(l  ,k  ,j  ,1)
								)/12.;
						
						double dfxm2=1./get_flat_df2x(j);
						double dfym2=1./get_flat_df2y(k);
						double dfzm2=1./get_flat_df2z(l);
						double ddfpdfx=get_flat_Gamx(j);
						double ddfpdfy=get_flat_Gamy(k);
						double ddfpdfz=get_flat_Gamz(l);
						
						double dev=(
							dfxm2*wa_xx
							+dfym2*wa_yy
							+dfzm2*wa_zz
							-dfxm2*ddfpdfx*wa_x_p*dx
							-dfym2*ddfpdfy*wa_y_p*dx
							-dfzm2*ddfpdfz*wa_z_p*dx
							-get_coef(l,k,j,8)*dx*dx
							)/(2.5*(dfxm2+dfym2+dfzm2));
						
						double err=abs(dev);
						
						if(err>outputerr)
						{
							cout << "large error in solveham" << endl;
							cout << "j=" << j << endl
							<< "k=" << k << endl 
							<< "l=" << l << endl
							<< "n=" << n << endl 
							<< "err=" << err << endl 
							<< "dev=" << dev << endl 
							//<< "norm=" << norm 
							<< endl;
							
							cout << " get_wvec(l,k,j,1)=" << get_wvec(l,k,j,1) << endl 
								<< " get_wvec(l,k,j+1,1)=" << get_wvec(l,k,j+1,1) << endl 
								<< " get_wvec(l,k,j-1,1)=" << get_wvec(l,k,j-1,1) << endl 
								<< " get_wvec(l,k+1,j,1)=" << get_wvec(l,k+1,j,1) << endl 
								<< " get_wvec(l,k-1,j,1)=" << get_wvec(l,k-1,j,1) << endl 
								<< " get_wvec(l+1,k,j,1)=" << get_wvec(l+1,k,j,1) << endl 
								<< " get_wvec(l-1,k,j,1)=" << get_wvec(l-1,k,j,1) << endl 
								<< " get_wvec(l,k-1,j-1,1)=" << get_wvec(l,k-1,j-1,1) << endl 
								<< " get_wvec(l,k+1,j+1,1)=" << get_wvec(l,k+1,j+1,1) << endl 
								<< " get_coef(l,k,j,8)=" << get_coef(l,k,j,8) << endl 
								<< endl;
						
							exit(1);
						}
						
						set_wvec(l,k,j,1)=acc*dev;

						if (errm[nn]<err)
						errm[nn]=err;
					}
				}
			}
		}
		
		#pragma omp section
		{
			int nn=3;
		
			errm[nn]=0.;
		
			for(int l=lli;l<=lui;l++)
			{
				for(int k=kli;k<=kui;k++)
				{
					for(int j=jli;j<=jui;j++)
					{
						if(get_bflag(l,k,j)!=0)
						continue;

						//double wa_p=   get_wvec(l,k,j,2);
						
						double wa_x_p=(    -(get_wvec(l,k,j+2,2)-get_wvec(l,k,j-2,2))
								+8.*(get_wvec(l,k,j+1,2)-get_wvec(l,k,j-1,2))
								)/12.;
						double wa_y_p=(   -(get_wvec(l,k+2,j,2)-get_wvec(l,k-2,j,2))
								+8.*(get_wvec(l,k+1,j,2)-get_wvec(l,k-1,j,2))
								)/12.;
						double wa_z_p=(   -(get_wvec(l+2,k,j,2)-get_wvec(l-2,k,j,2))
								+8.*(get_wvec(l+1,k,j,2)-get_wvec(l-1,k,j,2))
								)/12.;

						double wa_xx=(    -(    get_wvec(l  ,k  ,j+2,2) + get_wvec(l  ,k  ,j-2,2))
								+16.*( get_wvec(l  ,k  ,j+1,2) + get_wvec(l  ,k  ,j-1,2))
								-30.*  get_wvec(l  ,k  ,j  ,2)
								)/12.;
						double wa_yy=(    -(    get_wvec(l  ,k+2,j  ,2) + get_wvec(l  ,k-2,j  ,2))
								+16.*( get_wvec(l  ,k+1,j  ,2) + get_wvec(l  ,k-1,j  ,2))
								-30.*  get_wvec(l  ,k  ,j  ,2)
								)/12.;
						double wa_zz=(    -(    get_wvec(l+2,k  ,j  ,2) + get_wvec(l-2,k  ,j  ,2))
								+16.*( get_wvec(l+1,k  ,j  ,2) + get_wvec(l-1,k  ,j  ,2))
								-30.*  get_wvec(l  ,k  ,j  ,2)
								)/12.;
						
						double dfxm2=1./get_flat_df2x(j);
						double dfym2=1./get_flat_df2y(k);
						double dfzm2=1./get_flat_df2z(l);
						double ddfpdfx=get_flat_Gamx(j);
						double ddfpdfy=get_flat_Gamy(k);
						double ddfpdfz=get_flat_Gamz(l);
						
						
						double dev=(
							dfxm2*wa_xx
							+dfym2*wa_yy
							+dfzm2*wa_zz
							-dfxm2*ddfpdfx*wa_x_p*dx
							-dfym2*ddfpdfy*wa_y_p*dx
							-dfzm2*ddfpdfz*wa_z_p*dx
							-get_coef(l,k,j,9)*dx*dx
							)/(2.5*(dfxm2+dfym2+dfzm2));
						
						double err=abs(dev);
						
						if(err>outputerr)
						{
							cout << "large error in solveham" << endl;
							cout << "j=" << j << endl
							<< "k=" << k << endl 
							<< "l=" << l << endl
							<< "n=" << n << endl 
							<< "err=" << err << endl 
							<< "dev=" << dev << endl 
							//<< "norm=" << norm 
							<< endl;
							
							cout << " get_wvec(l,k,j,2)=" << get_wvec(l,k,j,2) << endl 
								<< " get_wvec(l,k,j+1,2)=" << get_wvec(l,k,j+1,2) << endl 
								<< " get_wvec(l,k,j-1,2)=" << get_wvec(l,k,j-1,2) << endl 
								<< " get_wvec(l,k+1,j,2)=" << get_wvec(l,k+1,j,2) << endl 
								<< " get_wvec(l,k-1,j,2)=" << get_wvec(l,k-1,j,2) << endl 
								<< " get_wvec(l+1,k,j,2)=" << get_wvec(l+1,k,j,2) << endl 
								<< " get_wvec(l-1,k,j,2)=" << get_wvec(l-1,k,j,2) << endl 
								<< " get_wvec(l,k-1,j-1,2)=" << get_wvec(l,k-1,j-1,2) << endl 
								<< " get_wvec(l,k+1,j+1,2)=" << get_wvec(l,k+1,j+1,2) << endl 
								<< " get_coef(l,k,j,9)=" << get_coef(l,k,j,9) << endl 
								<< endl;
						
							exit(1);
						}
						
						set_wvec(l,k,j,2)=acc*dev;

						if (errm[nn]<err)
						errm[nn]=err;
					}
				}
			}
		}
		}
		#pragma omp barrier
		
		errmm=0.;
		int mm=0;
		
		for(int nn=0;nn<4;nn++)
		{
			//cout << "errm[nn]=" << errm[nn] << endl;
			if(errm[nn]>errmm)
			{
				errmm=errm[nn];
				mm=nn;
			}
		}
		
		if(errmm<convcond)
		{
			cout << "solveconst step=" << n 
				<< " errm solveconst=" << errmm 
				<< " for mm=" << mm 
				<< endl;
			conv=true;
			break;
			
		}
		else
		{
			cout << "not yet solveconst step=" << n 
				<< " errm solveconst=" << errmm 
				<< " for mm=" << mm 
				<< endl;
		}
		
		boundary_reflection_even(13);
		boundary_w();
		set_source();
		
	}
	
	
	if(!conv)
	{
		cout << "not converge in solveham step" 
			<< " errmm=" << errmm << endl;
	}
	return;
}




