/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* BSSN EVOLUTION :: BSSN evolution Class of COSMOS                                                      */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
//#include "omp.h"
#include <algorithm>

void Fmv0::BSSN(int itype)
{
	//time deviation and fraction for Rungekutta start
	double frac;
	double tt;
	
	switch(itype){
	case 1:
		dt=0.5*dt0;
		frac=1./6.;
		tt=get_t();
		break;
	case 2:
		dt=0.5*dt0;
		frac=1./3.;
		tt=get_t()+dt;
		break;
	case 3:
		dt=dt0;
		frac=1./3.;
		tt=get_t()+0.5*dt;
		break;
	case 4:
		dt=dt0;
		frac=1./6.;
		tt=get_t()+dt;
		break;
	default:
		frac=1./6.;
		tt=get_t();
	}
	//time deviation and fraction for Rungekutta end
		
	// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
	// 24:E, 25:S_x, 26:S_y, 27:S_z, 28:N_B(D), 
	// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
	
	set_zero_d();	//bvd=0
	
	//preparation associated with \beta^i
	#pragma omp barrier
	BSSN_adv();

	//setting for fluid flux
	#pragma omp barrier
	if(fluidevo)
	{
		//getting primitive variables from dynamical variables
		dyntoprim();

		//flux calculation
		#pragma omp barrier
		flux_fill();
	}

	#pragma omp parallel for 
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				if(get_bflag(l,k,j)!=0)
				continue;

				//definitions of variables start
				double alpha_p,															//lapse
				gxx_p,gyy_p,gzz_p,gxy_p,gxz_p,gyz_p,gyx_p,gzx_p,gzy_p,					//conformal metric
				wa_p,																	//conformal factor(e^(4 wa_p))
				akxx_p,akyy_p,akzz_p,akxy_p,akxz_p,akyz_p,akyx_p,akzx_p,akzy_p,ek_p;	//extrinsic curvature
				double bx_p,by_p,bz_p;													//shift
				double zgx_p,zgy_p,zgz_p;												//\tilde Gamma^j(see https://arxiv.org/abs/gr-qc/0703035v1)
				
				//x-derivatives
				double alpha_x_p,
				bx_x_p,by_x_p,bz_x_p,
				dxx_x_p,dyy_x_p,dzz_x_p,
				dxy_x_p,dxz_x_p,dyz_x_p,
				dyx_x_p,dzx_x_p,dzy_x_p,wa_x_p,
				ek_x_p,
				zgx_x_p,zgy_x_p,zgz_x_p;

				//y-derivatives
				double alpha_y_p,bx_y_p,by_y_p,bz_y_p,
				dxx_y_p,dyy_y_p,dzz_y_p,
				dxy_y_p,dxz_y_p,dyz_y_p,
				dyx_y_p,dzx_y_p,dzy_y_p,wa_y_p,
				ek_y_p,
				zgx_y_p,zgy_y_p,zgz_y_p;
				
				//z-derivatives
				double alpha_z_p,bx_z_p,by_z_p,bz_z_p,
				dxx_z_p,dyy_z_p,dzz_z_p,
				dxy_z_p,dxz_z_p,dyz_z_p,
				dyx_z_p,dzx_z_p,dzy_z_p,wa_z_p,
				ek_z_p,
				zgx_z_p,zgy_z_p,zgz_z_p;
				
				//second derivatives of lapse
				double alpha_xx,alpha_yy,alpha_zz,
				alpha_xy,alpha_xz,alpha_yz;
				
				//second derivatives of shift
				double bx_xx,bx_yy,bx_zz,
				bx_xy,bx_xz,bx_yz,
				bx_yx,bx_zx,bx_zy;
				double by_xx,by_yy,by_zz,
				by_xy,by_xz,by_yz,
				by_yx,by_zx,by_zy;
				double bz_xx,bz_yy,bz_zz,
				bz_xy,bz_xz,bz_yz,
				bz_yx,bz_zx,bz_zy;
				
				//second derivatives of conformal factor
				double wa_xx,wa_yy,wa_zz,
				wa_xy,wa_xz,wa_yz;
				
				//second derivatives of conformal metric
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
				
				//conformal metric determinant and inverse
				double det_p,det_pi;
				
				//inverse conformal metric
				double gixx_p,giyy_p,gizz_p,
				gixy_p,gixz_p,giyz_p,
				giyx_p,gizx_p,gizy_p;
				
				//\tilde A_i^j
				double akx_ux_p,aky_uy_p,akz_uz_p,
				akx_uy_p,akx_uz_p,aky_uz_p,
				aky_ux_p,akz_ux_p,akz_uy_p;
				
				//\tilde A^{ij}
				double ak_ixx,ak_iyy,ak_izz,
				ak_ixy,ak_ixz,ak_iyz,
				ak_iyx,ak_izx,ak_izy;
				
				//  \tilde{\Gamma}_{ijk}=\tilde{\gamma}_{im}\tilde{\Gamma}^m_{jk} -> crd{i}_{jk}
				double crdx_xx,crdx_yy,crdx_zz,
				crdx_xy,crdx_xz,crdx_yz;
				double crdy_xx,crdy_yy,crdy_zz,
				crdy_xy,crdy_xz,crdy_yz;
				double crdz_xx,crdz_yy,crdz_zz,
				crdz_xy,crdz_xz,crdz_yz;
				
				//  \tilde{\Gamma}^i_{jk} -> cr{i}_{jk}
				double crx_xx,crx_yy,crx_zz,
				crx_xy,crx_xz,crx_yz;
				double cry_xx,cry_yy,cry_zz,
				cry_xy,cry_xz,cry_yz;
				double crz_xx,crz_yy,crz_zz,
				crz_xy,crz_xz,crz_yz;

				//tilde Gamma0^i =- cal D_j tilde gamma ^ji=tilde gamma^jk Delta^i_jk		//(eq(9.80) in https://arxiv.org/abs/gr-qc/0703035v1
				double gamma0_x,gamma0_y,gamma0_z;
				
				// \tilde D_i D_j \psi = \psi,ij -\tilde{\Gamma}^k_ij \psi_k
				// \tilde Laplacian \psi
				// (\tilde D_i \psi)(\tilde D^i \psi)
				double wa_cdxx,wa_cdyy,wa_cdzz,
				wa_cdxy,wa_cdxz,wa_cdyz,
				wa_cdyx,wa_cdzx,wa_cdzy,
				wa_lap,wawa;
				
				// exp(-4\psi)* R_{jk} -> rc_{ij}; 
				double rc_xx_p,rc_yy_p,rc_zz_p,
				rc_xy_p,rc_xz_p,rc_yz_p,
				rc_yx_p,rc_zx_p,rc_zy_p;
				
				//\tilde A_i^j \tilde A_j^i
				// R (not \tilde R)
				double aaaa,ricci;
				
				//D_i beta^i and 2/3 D_i beta^i
				double divbeta,divbeta23;
				
				//time forward values
				double fwa,fgxx,fgyy,fgzz,fgxy,fgxz,fgyz;
				
				// \tilde{\gamma}^{ij} alpha_i \phi_j -> alphawa_ip
				// exp(-4*\phi) D_i D_j \alpha
				// Lap \alpha
				double alphawa_ip,alpha_cdxx4,alpha_cdyy4,alpha_cdzz4,
					alpha_cdxy4,alpha_cdxz4,alpha_cdyz4,alpha_cdyx4,alpha_cdzx4,alpha_cdzy4,
					alpha_cdtr;
				
				//temporary variables for fak
				double ftmp_xx,ftmp_yy,ftmp_zz,
				ftmp_xy,ftmp_xz,ftmp_yz,
				ftmp_yx,ftmp_zx,ftmp_zy,ftmptr,
				trs;
				
				//time forward values for extrinsic curvature
				double fakxx,fakyy,fakzz,
				fakxy,fakxz,fakyz,fek;
				
				//Laplacian for shift and derivative of divergence of shift
				double bx_lap,by_lap,bz_lap,
				divb_x,divb_y,divb_z;
				
				//time forward values for \tilde \Gamma^i
				double fzgx,fzgy,fzgz;
				
				//time forward values for gauge
				double fbx,fby,fbz,fbbx,fbby,fbbz,falpha;
				
				//stress-energy tensor components
				double Ene=0.;
				double px;
				double py;
				double pz;
				double p_x=0.;
				double p_y=0.;
				double p_z=0.;
				double sxx=0.;
				double sxy=0.;
				double sxz=0.;
				double syy=0.;
				double syz=0.;
				double szz=0.;

				//reference metric variables
				double fxx,fyy,fzz,dfxx,dfyy,dfzz;
				
				//bar Gamma_ui_jk for inhomogeneous grid 
				double Gam_ux_xx,Gam_uy_yy,Gam_uz_zz;
				
				//double delG_1_u2_34
				double delG_x_ux_xx,delG_y_uy_yy,delG_z_uz_zz;
				
				//cal D_i tilde gamma_jk
				double Dgam_x_xx,Dgam_x_yy,Dgam_x_zz,Dgam_x_xy,Dgam_x_xz,Dgam_x_yz,Dgam_x_yx,Dgam_x_zx,Dgam_x_zy,
						Dgam_y_xx,Dgam_y_yy,Dgam_y_zz,Dgam_y_xy,Dgam_y_xz,Dgam_y_yz,Dgam_y_yx,Dgam_y_zx,Dgam_y_zy,
						Dgam_z_xx,Dgam_z_yy,Dgam_z_zz,Dgam_z_xy,Dgam_z_xz,Dgam_z_yz,Dgam_z_yx,Dgam_z_zx,Dgam_z_zy;
				
				//Delta_ui_jk=tilde Gamma - bar Gamma
				double Del_ux_xx,Del_ux_yy,Del_ux_zz,Del_ux_xy,Del_ux_xz,Del_ux_yz,Del_ux_yx,Del_ux_zx,Del_ux_zy,
						Del_uy_xx,Del_uy_yy,Del_uy_zz,Del_uy_xy,Del_uy_xz,Del_uy_yz,Del_uy_yx,Del_uy_zx,Del_uy_zy,
						Del_uz_xx,Del_uz_yy,Del_uz_zz,Del_uz_xy,Del_uz_xz,Del_uz_yz,Del_uz_yx,Del_uz_zx,Del_uz_zy;
				
				//cal D_i tilde gamma^jk
				double Dgam_x_uxx,Dgam_x_uyy,Dgam_x_uzz,Dgam_x_uxy,Dgam_x_uxz,Dgam_x_uyz,/*Dgam_x_uyx,Dgam_x_uzx,Dgam_x_uzy,*/
						Dgam_y_uxx,Dgam_y_uyy,Dgam_y_uzz,Dgam_y_uxy,Dgam_y_uxz,Dgam_y_uyz,/*Dgam_y_uyx,Dgam_y_uzx,Dgam_y_uzy,*/
						Dgam_z_uxx,Dgam_z_uyy,Dgam_z_uzz,Dgam_z_uxy,Dgam_z_uxz,Dgam_z_uyz/*,Dgam_z_uyx,Dgam_z_uzx,Dgam_z_uzy*/;
				
				//\tilde gamma^kl \del_k \del_l \tilde \gamma_ij
				double lapgam_xx,lapgam_yy,lapgam_zz,
						lapgam_xy,lapgam_xz,lapgam_yz;
				
				//\tilde gamma^kl \cal D_k \cal D_l \tilde \gamma_ij
				double gDDg_xx,gDDg_yy,gDDg_zz,
						gDDg_xy,gDDg_xz,gDDg_yz;
				
				//\cal D_i \tilde Gamma^j
				double DGam_x_ux,DGam_y_uy,DGam_z_uz,
						DGam_x_uy,DGam_x_uz,DGam_y_uz,
						DGam_y_ux,DGam_z_ux,DGam_z_uy;
				
				//exp(-4*wa), exp(-8*wa), exp(-12*wa)
				double ewa4i,ewa8i,ewa12i;
				
				//fij for inhomogeneous grid 
				fxx=get_flat_df2x(j);
				fyy=get_flat_df2y(k);
				fzz=get_flat_df2z(l);

				//substitution of geometrical variables start
				alpha_p=get_bv(l,k,j, 0);
				bx_p=   get_bv(l,k,j, 1);
				by_p=   get_bv(l,k,j, 2);
				bz_p=   get_bv(l,k,j, 3);
				gxx_p=  get_bv(l,k,j, 7)+fxx;
				gyy_p=  get_bv(l,k,j, 8)+fyy;
				gzz_p=  get_bv(l,k,j, 9)+fzz;
				gxy_p=  get_bv(l,k,j,10);
				gxz_p=  get_bv(l,k,j,11);
				gyz_p=  get_bv(l,k,j,12);
				gyx_p=	gxy_p;
				gzx_p=	gxz_p;
				gzy_p=	gyz_p;
				wa_p=   get_bv(l,k,j,13);
				
				ewa4i=exp(-4.*wa_p);
				ewa8i=exp(-8.*wa_p);
				ewa12i=exp(-12.*wa_p);

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
				//substitution of geometrical variables end
				
				//bar Gamma_ui_jk for inhomogeneous grid 
				Gam_ux_xx=get_flat_Gamx(j);
				Gam_uy_yy=get_flat_Gamy(k);
				Gam_uz_zz=get_flat_Gamz(l);
				
				//calculation of derivatives of flat Christoffel
				//double delG_1_u2_34
				delG_x_ux_xx=get_flat_dGamx(j);
				delG_y_uy_yy=get_flat_dGamy(k);
				delG_z_uz_zz=get_flat_dGamz(l);

				// first derivatives 4-th order
				// \del_x Gauge (\alpha \beta^i)
				alpha_x_p=get_f_x(l,k,j,0);
				bx_x_p=get_f_x(l,k,j,1);
				by_x_p=get_f_x(l,k,j,2);
				bz_x_p=get_f_x(l,k,j,3);
				
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

				// \del_y Gauge (\alpha \beta^i)
				alpha_y_p=get_f_y(l,k,j,0);
				bx_y_p=get_f_y(l,k,j,1);
				by_y_p=get_f_y(l,k,j,2);
				bz_y_p=get_f_y(l,k,j,3);
				
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
				
				// \del_z Gauge (\alpha \beta^i)
				alpha_z_p=get_f_z(l,k,j,0);
				bx_z_p=get_f_z(l,k,j,1);
				by_z_p=get_f_z(l,k,j,2);
				bz_z_p=get_f_z(l,k,j,3);
				
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
				
				alpha_xx=get_f_xx(l,k,j,0);
				alpha_yy=get_f_yy(l,k,j,0);
				alpha_zz=get_f_zz(l,k,j,0);
				
				alpha_xy=get_f_xy(l,k,j,0);
				alpha_xz=get_f_xz(l,k,j,0);
				alpha_yz=get_f_yz(l,k,j,0);
				
				// \del_i\del_j \beta^x
				bx_xx=get_f_xx(l,k,j,1);
				bx_yy=get_f_yy(l,k,j,1);
				bx_zz=get_f_zz(l,k,j,1);
				
				bx_xy=get_f_xy(l,k,j,1);
				bx_xz=get_f_xz(l,k,j,1);
				bx_yz=get_f_yz(l,k,j,1);
				bx_yx=bx_xy;
				bx_zx=bx_xz;
				bx_zy=bx_yz;

				// \del_i\del_j \beta^y
				by_xx=get_f_xx(l,k,j,2);
				by_yy=get_f_yy(l,k,j,2);
				by_zz=get_f_zz(l,k,j,2);
				
				by_xy=get_f_xy(l,k,j,2);
				by_xz=get_f_xz(l,k,j,2);
				by_yz=get_f_yz(l,k,j,2);
				by_yx=by_xy;
				by_zx=by_xz;
				by_zy=by_yz;

				// \del_i\del_j \beta^z
				bz_xx=get_f_xx(l,k,j,3);
				bz_yy=get_f_yy(l,k,j,3);
				bz_zz=get_f_zz(l,k,j,3);
				
				bz_xy=get_f_xy(l,k,j,3);
				bz_xz=get_f_xz(l,k,j,3);
				bz_yz=get_f_yz(l,k,j,3);
				bz_yx=bz_xy;
				bz_zx=bz_xz;
				bz_zy=bz_yz;

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
				det_p=gxx_p*gyy_p*gzz_p+gxy_p*gyz_p*gzx_p+gxz_p*gyx_p*gzy_p
						-gxz_p*gyy_p*gzx_p-gxy_p*gyx_p*gzz_p-gxx_p*gyz_p*gzy_p;
				if(det_p<1.e-16) 
				det_p=1.e-16;
				//	det_p=get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l);
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

				//  \tilde{A}^{ij} -> ak_i{ij}
				ak_ixx=gixx_p*akx_ux_p +gixy_p*aky_ux_p +gixz_p*akz_ux_p;
				ak_ixy=gixx_p*akx_uy_p +gixy_p*aky_uy_p +gixz_p*akz_uy_p;
				ak_ixz=gixx_p*akx_uz_p +gixy_p*aky_uz_p +gixz_p*akz_uz_p;
				ak_iyx=giyx_p*akx_ux_p +giyy_p*aky_ux_p +giyz_p*akz_ux_p;
				ak_iyy=giyx_p*akx_uy_p +giyy_p*aky_uy_p +giyz_p*akz_uy_p;
				ak_iyz=giyx_p*akx_uz_p +giyy_p*aky_uz_p +giyz_p*akz_uz_p;
				ak_izx=gizx_p*akx_ux_p +gizy_p*aky_ux_p +gizz_p*akz_ux_p;
				ak_izy=gizx_p*akx_uy_p +gizy_p*aky_uy_p +gizz_p*akz_uy_p;
				ak_izz=gizx_p*akx_uz_p +gizy_p*aky_uz_p +gizz_p*akz_uz_p;

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
				
				rc_xx_p=(0.5*(-gDDg_xx
							+gxx_p*DGam_x_ux+gxy_p*DGam_x_uy+gxz_p*DGam_x_uz
							+gxx_p*DGam_x_ux+gxy_p*DGam_x_uy+gxz_p*DGam_x_uz)
						-0.5*(Dgam_x_xx*Dgam_x_uxx+Dgam_y_yx*Dgam_x_uyy+Dgam_z_zx*Dgam_x_uzz
							+Dgam_x_yx*Dgam_x_uxy+Dgam_x_zx*Dgam_x_uxz+Dgam_y_zx*Dgam_x_uyz
							+Dgam_y_xx*Dgam_x_uxy+Dgam_z_xx*Dgam_x_uxz+Dgam_z_yx*Dgam_x_uyz
							
							+Dgam_x_xx*Dgam_x_uxx+Dgam_y_yx*Dgam_x_uyy+Dgam_z_zx*Dgam_x_uzz
							+Dgam_x_yx*Dgam_x_uxy+Dgam_x_zx*Dgam_x_uxz+Dgam_y_zx*Dgam_x_uyz
							+Dgam_y_xx*Dgam_x_uxy+Dgam_z_xx*Dgam_x_uxz+Dgam_z_yx*Dgam_x_uyz
							
							-zgx_p*Dgam_x_xx-zgy_p*Dgam_y_xx-zgz_p*Dgam_z_xx)
							//-gamma0_x*Dgam_x_xx-gamma0_y*Dgam_y_xx-gamma0_z*Dgam_z_xx)
						
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
							
							-zgx_p*Dgam_x_yy-zgy_p*Dgam_y_yy-zgz_p*Dgam_z_yy)
							//-gamma0_x*Dgam_x_yy-gamma0_y*Dgam_y_yy-gamma0_z*Dgam_z_yy)
						
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
							
							-zgx_p*Dgam_x_zz-zgy_p*Dgam_y_zz-zgz_p*Dgam_z_zz)
							//-gamma0_x*Dgam_x_zz-gamma0_y*Dgam_y_zz-gamma0_z*Dgam_z_zz)
						
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
							
							-zgx_p*Dgam_x_xy-zgy_p*Dgam_y_xy-zgz_p*Dgam_z_xy)
							//-gamma0_x*Dgam_x_xy-gamma0_y*Dgam_y_xy-gamma0_z*Dgam_z_xy)
						
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
							
							-zgx_p*Dgam_x_xz-zgy_p*Dgam_y_xz-zgz_p*Dgam_z_xz)
							//-gamma0_x*Dgam_x_xz-gamma0_y*Dgam_y_xz-gamma0_z*Dgam_z_xz)
						
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
							
							-zgx_p*Dgam_x_yz-zgy_p*Dgam_y_yz-zgz_p*Dgam_z_yz)
							//-gamma0_x*Dgam_x_yz-gamma0_y*Dgam_y_yz-gamma0_z*Dgam_z_yz)
						
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

				//  \tilde{A}^{ij} -> ak_i{ij}
				ak_ixx=gixx_p*akx_ux_p +gixy_p*aky_ux_p +gixz_p*akz_ux_p;
				ak_ixy=gixx_p*akx_uy_p +gixy_p*aky_uy_p +gixz_p*akz_uy_p;
				ak_ixz=gixx_p*akx_uz_p +gixy_p*aky_uz_p +gixz_p*akz_uz_p;
				ak_iyx=giyx_p*akx_ux_p +giyy_p*aky_ux_p +giyz_p*akz_ux_p;
				ak_iyy=giyx_p*akx_uy_p +giyy_p*aky_uy_p +giyz_p*akz_uy_p;
				ak_iyz=giyx_p*akx_uz_p +giyy_p*aky_uz_p +giyz_p*akz_uz_p;
				ak_izx=gizx_p*akx_ux_p +gizy_p*aky_ux_p +gizz_p*akz_ux_p;
				ak_izy=gizx_p*akx_uy_p +gizy_p*aky_uy_p +gizz_p*akz_uy_p;
				ak_izz=gizx_p*akx_uz_p +gizy_p*aky_uz_p +gizz_p*akz_uz_p;

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
				
				//D_i beta^i
				divbeta=bx_x_p +by_y_p +bz_z_p+Gam_ux_xx*bx_p+Gam_uy_yy*by_p+Gam_uz_zz*bz_p;
				//divbeta=bx_x_p +by_y_p +bz_z_p;
				divbeta23=2./3.*divbeta;

				////////////// for scalar field ///////////////////
				if(scalarevo)
				{
					//time forward values for scalar field and conjugate momentum
					double fphi,fPi;
					
					//scalar field and conjugate momentum, those derivatives
					double phi,Pi,phi_x,phi_y,phi_z,phi_xx,phi_yy,phi_zz,phi_xy,phi_xz,phi_yz;
					
					// \tilde D_i D_j \phi = \phi,ij -\tilde{\Gamma}^k_ij \phi_k
					double phi_cdxx,phi_cdyy,phi_cdzz,
					phi_cdxy,phi_cdxz,phi_cdyz,
					phi_cdyx,phi_cdzx,phi_cdzy;
					
					// \tilde Laplacian \psi
					// (\tilde D_i \phi)(\tilde D^i \phi)
					// (\tilde D_i \psi)(\tilde D^i \phi)
					// (\tilde D_i \alpha)(\tilde D^i \phi)
					double phi_lap,dphidphi,dpsidphi,dalphadphi;
					
					//scalar field potential and derivative
					double Vpot,dVpot;
					
					//substitution start
					phi=get_bv(l,k,j,nsc);
					Pi=get_bv(l,k,j,nscp);
					
					phi_x=get_f_x(l,k,j,nsc);
					phi_y=get_f_y(l,k,j,nsc);
					phi_z=get_f_z(l,k,j,nsc);
					phi_xx=get_f_xx(l,k,j,nsc);
					phi_yy=get_f_yy(l,k,j,nsc);
					phi_zz=get_f_zz(l,k,j,nsc);
					phi_xy=get_f_xy(l,k,j,nsc);
					phi_xz=get_f_xz(l,k,j,nsc);
					phi_yz=get_f_yz(l,k,j,nsc);
				
					// \tilde D_i D_j \phi = \phi,ij -\tilde{\Gamma}^k_ij \phi_k
					// ok even for non-Cartesian
					phi_cdxx=phi_xx-(crx_xx*phi_x +cry_xx*phi_y +crz_xx*phi_z);
					phi_cdyy=phi_yy-(crx_yy*phi_x +cry_yy*phi_y +crz_yy*phi_z);
					phi_cdzz=phi_zz-(crx_zz*phi_x +cry_zz*phi_y +crz_zz*phi_z);
					phi_cdxy=phi_xy-(crx_xy*phi_x +cry_xy*phi_y +crz_xy*phi_z);
					phi_cdxz=phi_xz-(crx_xz*phi_x +cry_xz*phi_y +crz_xz*phi_z);
					phi_cdyz=phi_yz-(crx_yz*phi_x +cry_yz*phi_y +crz_yz*phi_z);
					phi_cdyx=phi_cdxy;
					phi_cdzx=phi_cdxz;
					phi_cdzy=phi_cdyz;
					
					
					// \tilde Laplacian \psi
					phi_lap=gixx_p*phi_cdxx +gixy_p*phi_cdxy +gixz_p*phi_cdxz
							+giyx_p*phi_cdyx +giyy_p*phi_cdyy +giyz_p*phi_cdyz
							+gizx_p*phi_cdzx +gizy_p*phi_cdzy +gizz_p*phi_cdzz;
					
					// (\tilde D_i \phi)(\tilde D^i \phi)
					dphidphi=gixx_p*phi_x*phi_x +gixy_p*phi_x*phi_y +gixz_p*phi_x*phi_z
						+giyx_p*phi_y*phi_x +giyy_p*phi_y*phi_y +giyz_p*phi_y*phi_z
						+gizx_p*phi_z*phi_x +gizy_p*phi_z*phi_y +gizz_p*phi_z*phi_z;
					
					// (\tilde D_i \psi)(\tilde D^i \phi)
					dpsidphi=gixx_p*wa_x_p*phi_x +gixy_p*wa_x_p*phi_y +gixz_p*wa_x_p*phi_z
						+giyx_p*wa_y_p*phi_x +giyy_p*wa_y_p*phi_y +giyz_p*wa_y_p*phi_z
						+gizx_p*wa_z_p*phi_x +gizy_p*wa_z_p*phi_y +gizz_p*wa_z_p*phi_z;
					
					// (\tilde D_i \palpha)(\tilde D^i \phi)
					dalphadphi=gixx_p*alpha_x_p*phi_x +gixy_p*alpha_x_p*phi_y +gixz_p*alpha_x_p*phi_z
						+giyx_p*alpha_y_p*phi_x +giyy_p*alpha_y_p*phi_y +giyz_p*alpha_y_p*phi_z
						+gizx_p*alpha_z_p*phi_x +gizy_p*alpha_z_p*phi_y +gizz_p*alpha_z_p*phi_z;
					
					Vpot=funcV(phi);
					dVpot=funcdV(phi);

					///////////////////////////////
					// Scalar field evolution 
					///////////////////////////////
					
					fphi=-alpha_p*Pi;
					fPi=-ewa4i*(alpha_p*(2.*dpsidphi+phi_lap)+dalphadphi)+alpha_p*(ek_p*Pi+dVpot);
					
					set_dbv(l,k,j,nsc)=get_dbv(l,k,j,nsc)+fphi;
					set_dbv(l,k,j,nscp)=get_dbv(l,k,j,nscp)+fPi;
					
					//////////////////////////////////////////////////////
					// scalar field contribution to stress energy tensor
					//////////////////////////////////////////////////////
					
					Ene=0.5*(pow(Pi,2)+ewa4i*dphidphi)+Vpot;
					
					p_x=Pi*phi_x;
					p_y=Pi*phi_y;
					p_z=Pi*phi_z;
					
					sxx=phi_x*phi_x+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxx_p;
					syy=phi_y*phi_y+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gyy_p;
					szz=phi_z*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gzz_p;
					sxy=phi_x*phi_y+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxy_p;
					sxz=phi_x*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxz_p;
					syz=phi_y*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gyz_p;
				}
				
				if(fluidevo)
				{
					//fluid stress-energy tensor components
					double fEne,fpx,fpy,fpz,fp_x,fp_y,fp_z,fsxx,fsxy,fsxz,fsyy,fsyz,fszz;

					//E+P and P 
					double EpP,press;

					//\sqrt(\gamma)
					double sqgam;
					
					//time forward values for fluid variables
					double fS0,fSx,fSy,fSz,fDen;
					
					//extrinsic curvature componens
					double exKuu_xx,exKuu_xy;
					double exKuu_xz,exKuu_yy,exKuu_yz,exKuu_zz;
					
					//S_ij K^ij
					double SK;
					
					//fluid dynamical variables
					double S0,Sy,Sz,Sx;
					
					//stress tensor componets
					double fsuxx,fsuyy,fsuzz,fsuxy,fsuxz,fsuyz;
					
					//\del_i \gamma^jk;
					double dx_gamxx,dx_gamyy,dx_gamzz,dx_gamxy,dx_gamxz,dx_gamyz;
					double dy_gamxx,dy_gamyy,dy_gamzz,dy_gamxy,dy_gamxz,dy_gamyz;
					double dz_gamxx,dz_gamyy,dz_gamzz,dz_gamxy,dz_gamxz,dz_gamyz;
					
					//substitution start
					sqgam=sqrt(det_p/ewa12i);
					
					//energy momentum tensor of fluid
					fEne=get_bv(l,k,j,24)/sqgam;
					fp_x=get_bv(l,k,j,25)/sqgam;
					fp_y=get_bv(l,k,j,26)/sqgam;
					fp_z=get_bv(l,k,j,27)/sqgam;
					
					//energy momentum tensor
					Ene+=fEne;
					p_x+=fp_x;
					p_y+=fp_y;
					p_z+=fp_z;
				
					fpx=(gixx_p*fp_x+gixy_p*fp_y+gixz_p*fp_z)*ewa4i;
					fpy=(giyx_p*fp_x+giyy_p*fp_y+giyz_p*fp_z)*ewa4i;
					fpz=(gizx_p*fp_x+gizy_p*fp_y+gizz_p*fp_z)*ewa4i;
				
					press=pres(get_primv(l,k,j,0));
					
					EpP=fEne+press;
					
					press/=ewa4i;
					fsxx=fp_x*fp_x/EpP+press*gxx_p;
					fsyy=fp_y*fp_y/EpP+press*gyy_p;
					fszz=fp_z*fp_z/EpP+press*gzz_p;
					fsxy=fp_x*fp_y/EpP+press*gxy_p;
					fsxz=fp_x*fp_z/EpP+press*gxz_p;
					fsyz=fp_y*fp_z/EpP+press*gyz_p;
					
					//stress tensor components
					sxx+=fsxx;
					syy+=fsyy;
					szz+=fszz;
					sxy+=fsxy;
					sxz+=fsxz;
					syz+=fsyz;
					
					///////////////////////////////
					// fluid evolution
					///////////////////////////////
						
					fS0=-dxi*(get_flux_x(l,k,j+1,0)-get_flux_x(l,k,j,0));
					fSx=-dxi*(get_flux_x(l,k,j+1,1)-get_flux_x(l,k,j,1));
					fSy=-dxi*(get_flux_x(l,k,j+1,2)-get_flux_x(l,k,j,2));
					fSz=-dxi*(get_flux_x(l,k,j+1,3)-get_flux_x(l,k,j,3));
					fDen=-dxi*(get_flux_x(l,k,j+1,4)-get_flux_x(l,k,j,4));
					
					fS0=fS0-dyi*(get_flux_y(l,k+1,j,0)-get_flux_y(l,k,j,0));
					fSx=fSx-dyi*(get_flux_y(l,k+1,j,1)-get_flux_y(l,k,j,1));
					fSy=fSy-dyi*(get_flux_y(l,k+1,j,2)-get_flux_y(l,k,j,2));
					fSz=fSz-dyi*(get_flux_y(l,k+1,j,3)-get_flux_y(l,k,j,3));
					fDen=fDen-dyi*(get_flux_y(l,k+1,j,4)-get_flux_y(l,k,j,4));
					
					fS0=fS0-dzi*(get_flux_z(l+1,k,j,0)-get_flux_z(l,k,j,0));
					fSx=fSx-dzi*(get_flux_z(l+1,k,j,1)-get_flux_z(l,k,j,1));
					fSy=fSy-dzi*(get_flux_z(l+1,k,j,2)-get_flux_z(l,k,j,2));
					fSz=fSz-dzi*(get_flux_z(l+1,k,j,3)-get_flux_z(l,k,j,3));
					fDen=fDen-dzi*(get_flux_z(l+1,k,j,4)-get_flux_z(l,k,j,4));
					
					//extrinsic curvature
					exKuu_xx=(ak_ixx+gixx_p*ek_p/3.)*ewa4i;
					exKuu_xy=(ak_ixy+gixy_p*ek_p/3.)*ewa4i;
					exKuu_xz=(ak_ixz+gixz_p*ek_p/3.)*ewa4i;
					//exKuu_yx=(ak_iyx+giyx_p*ek_p/3.)*ewa4i;
					exKuu_yy=(ak_iyy+giyy_p*ek_p/3.)*ewa4i;
					exKuu_yz=(ak_iyz+giyz_p*ek_p/3.)*ewa4i;
					//exKuu_zx=(ak_izx+gizx_p*ek_p/3.)*ewa4i;
					//exKuu_zy=(ak_izy+gizy_p*ek_p/3.)*ewa4i;
					exKuu_zz=(ak_izz+gizz_p*ek_p/3.)*ewa4i;
					
					SK=fsxx*exKuu_xx+fsyy*exKuu_yy+fszz*exKuu_zz
						+2.*(fsxy*exKuu_xy+fsxz*exKuu_xz+fsyz*exKuu_yz);
					
					S0=get_bv(l,k,j,24);
					Sx=get_bv(l,k,j,25);
					Sy=get_bv(l,k,j,26);
					Sz=get_bv(l,k,j,27);
					//Den=get_bv(l,k,j,28);
					
					press*=ewa8i;
					fsuxx=fpx*fpx/EpP+press*gixx_p;
					fsuyy=fpy*fpy/EpP+press*giyy_p;
					fsuzz=fpz*fpz/EpP+press*gizz_p;
					fsuxy=fpx*fpy/EpP+press*gixy_p;
					fsuxz=fpx*fpz/EpP+press*gixz_p;
					fsuyz=fpy*fpz/EpP+press*giyz_p;
					
					//d1_gam23=(4.*g23_p*wa_1_p+2.*d23_1_p)/ewa4i;
					dx_gamxx=(4.*gxx_p*wa_x_p+2.*dxx_x_p)/ewa4i;
					dx_gamyy=(4.*gyy_p*wa_x_p+2.*dyy_x_p)/ewa4i;
					dx_gamzz=(4.*gzz_p*wa_x_p+2.*dzz_x_p)/ewa4i;
					dx_gamxy=(4.*gxy_p*wa_x_p+2.*dxy_x_p)/ewa4i;
					dx_gamxz=(4.*gxz_p*wa_x_p+2.*dxz_x_p)/ewa4i;
					dx_gamyz=(4.*gyz_p*wa_x_p+2.*dyz_x_p)/ewa4i;
					
					dy_gamxx=(4.*gxx_p*wa_y_p+2.*dxx_y_p)/ewa4i;
					dy_gamyy=(4.*gyy_p*wa_y_p+2.*dyy_y_p)/ewa4i;
					dy_gamzz=(4.*gzz_p*wa_y_p+2.*dzz_y_p)/ewa4i;
					dy_gamxy=(4.*gxy_p*wa_y_p+2.*dxy_y_p)/ewa4i;
					dy_gamxz=(4.*gxz_p*wa_y_p+2.*dxz_y_p)/ewa4i;
					dy_gamyz=(4.*gyz_p*wa_y_p+2.*dyz_y_p)/ewa4i;
					
					dz_gamxx=(4.*gxx_p*wa_z_p+2.*dxx_z_p)/ewa4i;
					dz_gamyy=(4.*gyy_p*wa_z_p+2.*dyy_z_p)/ewa4i;
					dz_gamzz=(4.*gzz_p*wa_z_p+2.*dzz_z_p)/ewa4i;
					dz_gamxy=(4.*gxy_p*wa_z_p+2.*dxy_z_p)/ewa4i;
					dz_gamxz=(4.*gxz_p*wa_z_p+2.*dxz_z_p)/ewa4i;
					dz_gamyz=(4.*gyz_p*wa_z_p+2.*dyz_z_p)/ewa4i;
					
					fS0=fS0+sqgam*(-fpx*alpha_x_p-fpy*alpha_y_p-fpz*alpha_z_p
								+alpha_p*SK);
					
					fSx=fSx-S0*alpha_x_p+(Sx*bx_x_p+Sy*by_x_p+Sz*bz_x_p)
							+0.5*alpha_p*sqgam*(fsuxx*dx_gamxx+fsuyy*dx_gamyy+fsuzz*dx_gamzz
											+2.*(fsuxy*dx_gamxy+fsuxz*dx_gamxz+fsuyz*dx_gamyz));
					
					fSy=fSy-S0*alpha_y_p+(Sx*bx_y_p+Sy*by_y_p+Sz*bz_y_p)
							+0.5*alpha_p*sqgam*(fsuxx*dy_gamxx+fsuyy*dy_gamyy+fsuzz*dy_gamzz
											+2.*(fsuxy*dy_gamxy+fsuxz*dy_gamxz+fsuyz*dy_gamyz));
					
					fSz=fSz-S0*alpha_z_p+(Sx*bx_z_p+Sy*by_z_p+Sz*bz_z_p)
							+0.5*alpha_p*sqgam*(fsuxx*dz_gamxx+fsuyy*dz_gamyy+fsuzz*dz_gamzz
											+2.*(fsuxy*dz_gamxy+fsuxz*dz_gamxz+fsuyz*dz_gamyz));
					
					set_dbv(l,k,j,24)=fS0;
					set_dbv(l,k,j,25)=fSx;
					set_dbv(l,k,j,26)=fSy;
					set_dbv(l,k,j,27)=fSz;
					set_dbv(l,k,j,28)=fDen;
				}
				
				px=ewa4i*(gixx_p*p_x+gixy_p*p_y+gixz_p*p_z);
				py=ewa4i*(gixy_p*p_x+giyy_p*p_y+giyz_p*p_z);
				pz=ewa4i*(gixz_p*p_x+giyz_p*p_y+gizz_p*p_z);

				///////////////////////////
				// Wa,\tilde{\gamma}_{ij}
				///////////////////////////
				
				fwa=-1./6.*(alpha_p*ek_p -divbeta);
				set_dbv(l,k,j,13)=get_dbv(l,k,j,13) +fwa;
				
				fgxx=-2.*alpha_p*akxx_p
					+gxx_p*bx_x_p +gxy_p*by_x_p +gxz_p*bz_x_p 
					+gxx_p*bx_x_p +gxy_p*by_x_p +gxz_p*bz_x_p 
					-gxx_p*divbeta23;
					//+bx_p*dfxx;
				fgyy=-2.*alpha_p*akyy_p 
					+gyx_p*bx_y_p +gyy_p*by_y_p +gyz_p*bz_y_p 
					+gyx_p*bx_y_p +gyy_p*by_y_p +gyz_p*bz_y_p 
					-gyy_p*divbeta23;
					//+by_p*dfyy;
				fgzz=-2.*alpha_p*akzz_p 
					+gzx_p*bx_z_p +gzy_p*by_z_p +gzz_p*bz_z_p 
					+gzx_p*bx_z_p +gzy_p*by_z_p +gzz_p*bz_z_p 
					-gzz_p*divbeta23;
					//+bz_p*dfzz;
				fgxy=-2.*alpha_p*akxy_p 
					+gxx_p*bx_y_p +gxy_p*by_y_p +gxz_p*bz_y_p 
					+gyx_p*bx_x_p +gyy_p*by_x_p +gyz_p*bz_x_p 
					-gxy_p*divbeta23;
				fgxz=-2.*alpha_p*akxz_p 
					+gxx_p*bx_z_p +gxy_p*by_z_p +gxz_p*bz_z_p 
					+gzx_p*bx_x_p +gzy_p*by_x_p +gzz_p*bz_x_p 
					-gxz_p*divbeta23;
				fgyz=-2.*alpha_p*akyz_p 
					+gyx_p*bx_z_p +gyy_p*by_z_p +gyz_p*bz_z_p 
					+gzx_p*bx_y_p +gzy_p*by_y_p +gzz_p*bz_y_p 
					-gyz_p*divbeta23;
				set_dbv(l,k,j, 7) =get_dbv(l,k,j, 7) +fgxx;
				set_dbv(l,k,j, 8) =get_dbv(l,k,j, 8) +fgyy;
				set_dbv(l,k,j, 9) =get_dbv(l,k,j, 9) +fgzz;
				set_dbv(l,k,j,10) =get_dbv(l,k,j,10) +fgxy;
				set_dbv(l,k,j,11) =get_dbv(l,k,j,11) +fgxz;
				set_dbv(l,k,j,12) =get_dbv(l,k,j,12) +fgyz;
				
				
				//////////////////////////////////////
				// \tilde{\gamma}^{ij} alpha_i \phi_j
				//////////////////////////////////////
				alphawa_ip=gixx_p*alpha_x_p*wa_x_p
							+gixy_p*alpha_x_p*wa_y_p
							+gixz_p*alpha_x_p*wa_z_p 
							+giyx_p*alpha_y_p*wa_x_p
							+giyy_p*alpha_y_p*wa_y_p 
							+giyz_p*alpha_y_p*wa_z_p 
							+gizx_p*alpha_z_p*wa_x_p
							+gizy_p*alpha_z_p*wa_y_p 
							+gizz_p*alpha_z_p*wa_z_p;

				// exp(-4*\phi) D_i D_j \alpha
				alpha_cdxx4=( alpha_xx
							-( crx_xx*alpha_x_p 
							+cry_xx*alpha_y_p 
							+crz_xx*alpha_z_p)
							-2.*( alpha_x_p*wa_x_p
							+alpha_x_p*wa_x_p
							-gxx_p*alphawa_ip) )*ewa4i;
				alpha_cdyy4=(alpha_yy
							-( crx_yy*alpha_x_p 
							+cry_yy*alpha_y_p 
							+crz_yy*alpha_z_p)
							-2.*( alpha_y_p*wa_y_p
							+alpha_y_p*wa_y_p 
							-gyy_p*alphawa_ip))*ewa4i;
				alpha_cdzz4=(alpha_zz
							-( crx_zz*alpha_x_p 
							+cry_zz*alpha_y_p 
							+crz_zz*alpha_z_p)
							-2.*( alpha_z_p*wa_z_p
							+alpha_z_p*wa_z_p 
							-gzz_p*alphawa_ip))*ewa4i;
				alpha_cdxy4=(alpha_xy
							-( crx_xy*alpha_x_p 
							+cry_xy*alpha_y_p 
							+crz_xy*alpha_z_p)
							-2.*( alpha_x_p*wa_y_p
							+alpha_y_p*wa_x_p 
							-gxy_p*alphawa_ip))*ewa4i;
				alpha_cdxz4=(alpha_xz
							-( crx_xz*alpha_x_p 
							+cry_xz*alpha_y_p 
							+crz_xz*alpha_z_p) 
							-2.*( alpha_x_p*wa_z_p
							+alpha_z_p*wa_x_p 
							-gxz_p*alphawa_ip))*ewa4i;
				alpha_cdyz4=(alpha_yz
							-( crx_yz*alpha_x_p 
							+cry_yz*alpha_y_p 
							+crz_yz*alpha_z_p)
							-2.*( alpha_y_p*wa_z_p
							+alpha_z_p*wa_y_p 
							-gyz_p*alphawa_ip))*ewa4i;
				alpha_cdyx4=alpha_cdxy4;
				alpha_cdzx4=alpha_cdxz4;
				alpha_cdzy4=alpha_cdyz4;
				

				ftmp_xx=-alpha_cdxx4 +alpha_p*(rc_xx_p-pi8*sxx*ewa4i);
				ftmp_yy=-alpha_cdyy4 +alpha_p*(rc_yy_p-pi8*syy*ewa4i);
				ftmp_zz=-alpha_cdzz4 +alpha_p*(rc_zz_p-pi8*szz*ewa4i);
				ftmp_xy=-alpha_cdxy4 +alpha_p*(rc_xy_p-pi8*sxy*ewa4i);
				ftmp_xz=-alpha_cdxz4 +alpha_p*(rc_xz_p-pi8*sxz*ewa4i);
				ftmp_yz=-alpha_cdyz4 +alpha_p*(rc_yz_p-pi8*syz*ewa4i);
				
				
				ftmp_yx=ftmp_xy;
				ftmp_zx=ftmp_xz;
				ftmp_zy=ftmp_yz;


				ftmptr=gixx_p*ftmp_xx +gixy_p*ftmp_xy +gixz_p*ftmp_xz
						+giyx_p*ftmp_yx +giyy_p*ftmp_yy +giyz_p*ftmp_yz
						+gizx_p*ftmp_zx +gizy_p*ftmp_zy +gizz_p*ftmp_zz;


				fakxx=ftmp_xx -1./3.*ftmptr*gxx_p
						+alpha_p*(-2.*( akx_ux_p*akxx_p 
						+akx_uy_p*akyx_p 
						+akx_uz_p*akzx_p) +ek_p*akxx_p )
						+akxx_p*bx_x_p +akxy_p*by_x_p +akxz_p*bz_x_p
						+akxx_p*bx_x_p +akxy_p*by_x_p +akxz_p*bz_x_p
						-akxx_p*divbeta23;
				fakyy=ftmp_yy-1./3.*ftmptr*gyy_p 
						+alpha_p*(-2.*( aky_ux_p*akxy_p 
						+aky_uy_p*akyy_p 
						+aky_uz_p*akzy_p) +ek_p*akyy_p ) 
						+akyx_p*bx_y_p +akyy_p*by_y_p +akyz_p*bz_y_p 
						+akyx_p*bx_y_p +akyy_p*by_y_p +akyz_p*bz_y_p 
						-akyy_p*divbeta23;
				fakzz=ftmp_zz-1./3.*ftmptr*gzz_p 
						+alpha_p*(-2.*( akz_ux_p*akxz_p 
						+akz_uy_p*akyz_p 
						+akz_uz_p*akzz_p) +ek_p*akzz_p ) 
						+akzx_p*bx_z_p +akzy_p*by_z_p +akzz_p*bz_z_p 
						+akzx_p*bx_z_p +akzy_p*by_z_p +akzz_p*bz_z_p 
						-akzz_p*divbeta23;
				fakxy=ftmp_xy-1./3.*ftmptr*gxy_p 
						+alpha_p*(-2.*( akx_ux_p*akxy_p 
						+akx_uy_p*akyy_p 
						+akx_uz_p*akzy_p) +ek_p*akxy_p ) 
						+akxx_p*bx_y_p +akxy_p*by_y_p +akxz_p*bz_y_p 
						+akyx_p*bx_x_p +akyy_p*by_x_p +akyz_p*bz_x_p 
						-akxy_p*divbeta23;
				fakxz=ftmp_xz-1./3.*ftmptr*gxz_p 
						+alpha_p*(-2.*( akx_ux_p*akxz_p 
						+akx_uy_p*akyz_p 
						+akx_uz_p*akzz_p) +ek_p*akxz_p ) 
						+akxx_p*bx_z_p +akxy_p*by_z_p +akxz_p*bz_z_p 
						+akzx_p*bx_x_p +akzy_p*by_x_p +akzz_p*bz_x_p 
						-akxz_p*divbeta23;
				fakyz=ftmp_yz-1./3.*ftmptr*gyz_p 
						+alpha_p*(-2.*( aky_ux_p*akxz_p 
						+aky_uy_p*akyz_p 
						+aky_uz_p*akzz_p) +ek_p*akyz_p ) 
						+akyx_p*bx_z_p +akyy_p*by_z_p +akyz_p*bz_z_p 
						+akzx_p*bx_y_p +akzy_p*by_y_p +akzz_p*bz_y_p 
						-akyz_p*divbeta23;
				set_dbv(l,k,j,14)=get_dbv(l,k,j,14) +fakxx;
				set_dbv(l,k,j,15)=get_dbv(l,k,j,15) +fakyy;
				set_dbv(l,k,j,16)=get_dbv(l,k,j,16) +fakzz;
				set_dbv(l,k,j,17)=get_dbv(l,k,j,17) +fakxy;
				set_dbv(l,k,j,18)=get_dbv(l,k,j,18) +fakxz;
				set_dbv(l,k,j,19)=get_dbv(l,k,j,19) +fakyz;

				// \Delta \alpha
				alpha_cdtr=gixx_p*alpha_cdxx4
							+gixy_p*alpha_cdxy4
							+gixz_p*alpha_cdxz4
							+giyx_p*alpha_cdyx4
							+giyy_p*alpha_cdyy4
							+giyz_p*alpha_cdyz4
							+gizx_p*alpha_cdzx4
							+gizy_p*alpha_cdzy4
							+gizz_p*alpha_cdzz4;

				trs=(gixx_p*sxx +gixy_p*sxy +gixz_p*sxz
						+giyx_p*sxy +giyy_p*syy +giyz_p*syz
						+gizx_p*sxz +gizy_p*syz +gizz_p*szz)*ewa4i;
				
				fek=-alpha_cdtr +alpha_p*(aaaa + pow(ek_p,2)/3.+pi4*(Ene+trs));
				set_dbv(l,k,j,20)=get_dbv(l,k,j,20) +fek;
				
				
				//added gi11_p*delG_1_u1_11*b1_p
				//				+2.*(Gam_u1_11*gi1x_p*b1_x_p+Gam_u1_11*gi1y_p*b1_y_p+Gam_u1_11*gi1z_p*b1_z_p)
				//				-Gam_ux_xx*gixx_p*b1_x_p-Gam_uy_yy*giyy_p*b1_y_p-Gam_uz_zz*gizz_p*b1_z_p;
				bx_lap=gixx_p*bx_xx +gixy_p*bx_xy +gixz_p*bx_xz 
						+giyx_p*bx_yx +giyy_p*bx_yy +giyz_p*bx_yz 
						+gizx_p*bx_zx +gizy_p*bx_zy +gizz_p*bx_zz
						
						+gixx_p*delG_x_ux_xx*bx_p
						+2.*(Gam_ux_xx*gixx_p*bx_x_p+Gam_ux_xx*gixy_p*bx_y_p+Gam_ux_xx*gixz_p*bx_z_p)
						
						-Gam_ux_xx*gixx_p*bx_x_p
						-Gam_uy_yy*giyy_p*bx_y_p
						-Gam_uz_zz*gizz_p*bx_z_p;

				by_lap=gixx_p*by_xx +gixy_p*by_xy +gixz_p*by_xz 
						+giyx_p*by_yx +giyy_p*by_yy +giyz_p*by_yz 
						+gizx_p*by_zx +gizy_p*by_zy +gizz_p*by_zz
						
						+giyy_p*delG_y_uy_yy*by_p
						+2.*(Gam_uy_yy*giyx_p*by_x_p+Gam_uy_yy*giyy_p*by_y_p+Gam_uy_yy*giyz_p*by_z_p)
						
						-Gam_ux_xx*gixx_p*by_x_p
						-Gam_uy_yy*giyy_p*by_y_p
						-Gam_uz_zz*gizz_p*by_z_p;

				bz_lap=gixx_p*bz_xx +gixy_p*bz_xy +gixz_p*bz_xz 
						+giyx_p*bz_yx +giyy_p*bz_yy +giyz_p*bz_yz 
						+gizx_p*bz_zx +gizy_p*bz_zy +gizz_p*bz_zz
						
						+gizz_p*delG_z_uz_zz*bz_p
						+2.*(Gam_uz_zz*gizx_p*bz_x_p+Gam_uz_zz*gizy_p*bz_y_p+Gam_uz_zz*gizz_p*bz_z_p)
						-Gam_ux_xx*gixx_p*bz_x_p
						-Gam_uy_yy*giyy_p*bz_y_p
						-Gam_uz_zz*gizz_p*bz_z_p;
				
				//modified by chulmoon for non-Cartesian inhomogeneous grid
				//added delG_1_u1_11*b1_p+Gam_ux_xx*bx_1_p+Gam_uy_yy*by_1_p+Gam_uz_zz*bz_1_p
				divb_x=bx_xx +by_yx +bz_zx
						+delG_x_ux_xx*bx_p+Gam_ux_xx*bx_x_p+Gam_uy_yy*by_x_p+Gam_uz_zz*bz_x_p;
				divb_y=bx_xy +by_yy +bz_zy
						+delG_y_uy_yy*by_p+Gam_ux_xx*bx_y_p+Gam_uy_yy*by_y_p+Gam_uz_zz*bz_y_p;
				divb_z=bx_xz +by_yz +bz_zz
						+delG_z_uz_zz*bz_p+Gam_ux_xx*bx_z_p+Gam_uy_yy*by_z_p+Gam_uz_zz*bz_z_p;
				
				
				//cr._.. -> Del_u._..
				fzgx=-2.*(ak_ixx*alpha_x_p +ak_ixy*alpha_y_p +ak_ixz*alpha_z_p ) 
					+2.*alpha_p*( Del_ux_xx*ak_ixx +Del_ux_xy*ak_ixy +Del_ux_xz*ak_ixz 
					+Del_ux_yx*ak_iyx +Del_ux_yy*ak_iyy +Del_ux_yz*ak_iyz 
					+Del_ux_zx*ak_izx +Del_ux_zy*ak_izy +Del_ux_zz*ak_izz 
					-2./3.*(gixx_p*ek_x_p 
					+gixy_p*ek_y_p 
					+gixz_p*ek_z_p)
					+6.*( ak_ixx*wa_x_p 
					+ak_ixy*wa_y_p 
					+ak_ixz*wa_z_p)  
					-pi8*px/ewa4i   ) 
					//-(bx_x_p*gamma0_x +bx_y_p*gamma0_y +bx_z_p*gamma0_z) 
					//+divbeta23*gamma0_x 
					-(bx_x_p*zgx_p +bx_y_p*zgy_p +bx_z_p*zgz_p) 
					+divbeta23*zgx_p
					+1./3.*(gixx_p*divb_x +gixy_p*divb_y +gixz_p*divb_z  ) + bx_lap;
				
				fzgy=-2.*(ak_iyx*alpha_x_p +ak_iyy*alpha_y_p +ak_iyz*alpha_z_p ) 
					+2.*alpha_p*( Del_uy_xx*ak_ixx +Del_uy_xy*ak_ixy +Del_uy_xz*ak_ixz 
					+Del_uy_yx*ak_iyx +Del_uy_yy*ak_iyy +Del_uy_yz*ak_iyz 
					+Del_uy_zx*ak_izx +Del_uy_zy*ak_izy +Del_uy_zz*ak_izz 
					-2./3.*(giyx_p*ek_x_p 
					+giyy_p*ek_y_p 
					+giyz_p*ek_z_p) 
					+6.*( ak_iyx*wa_x_p 
					+ak_iyy*wa_y_p 
					+ak_iyz*wa_z_p)    
					-pi8*py/ewa4i  ) 
					//-(by_x_p*gamma0_x +by_y_p*gamma0_y +by_z_p*gamma0_z) 
					//+divbeta23*gamma0_y 
					-(by_x_p*zgx_p +by_y_p*zgy_p +by_z_p*zgz_p) 
					+divbeta23*zgy_p 
					+1./3.*(giyx_p*divb_x +giyy_p*divb_y +giyz_p*divb_z  )  + by_lap;

				fzgz=-2.*(ak_izx*alpha_x_p +ak_izy*alpha_y_p +ak_izz*alpha_z_p ) 
					+2.*alpha_p*( Del_uz_xx*ak_ixx +Del_uz_xy*ak_ixy +Del_uz_xz*ak_ixz 
					+Del_uz_yx*ak_iyx +Del_uz_yy*ak_iyy +Del_uz_yz*ak_iyz 
					+Del_uz_zx*ak_izx +Del_uz_zy*ak_izy +Del_uz_zz*ak_izz 
					-2./3.*(gizx_p*ek_x_p
					+gizy_p*ek_y_p
					+gizz_p*ek_z_p)
					+6.*( ak_izx*wa_x_p
					+ak_izy*wa_y_p
					+ak_izz*wa_z_p)    
					-pi8*pz/ewa4i   )
					//-(bz_x_p*gamma0_x +bz_y_p*gamma0_y +bz_z_p*gamma0_z) 
					//+divbeta23*gamma0_z 
					-(bz_x_p*zgx_p +bz_y_p*zgy_p +bz_z_p*zgz_p) 
					+divbeta23*zgz_p 
					+1./3.*(gizx_p*divb_x +gizy_p*divb_y +gizz_p*divb_z )   + bz_lap;
				
				set_dbv(l,k,j,21)=get_dbv(l,k,j,21) +fzgx;
				set_dbv(l,k,j,22)=get_dbv(l,k,j,22) +fzgy;
				set_dbv(l,k,j,23)=get_dbv(l,k,j,23) +fzgz;
				

				///////////////////////////////
				// Gauge
				///////////////////////////////
				
				//modified 1+log
				falpha=-etaa*(ek_p+2./(1.+fluidw)/tt)*alpha_p;
				
				//synchronous
				//falpha=0.;
					
				//modified harmonic
				//falpha=-(ek_p-get_bv(lui,kui,jui,20))*pow(alpha_p,2);
				
				//Gamma driver
				fbx=etabb*get_bv(l,k,j,4);
				fby=etabb*get_bv(l,k,j,5);
				fbz=etabb*get_bv(l,k,j,6);
				
				//fbbx=fzgx -etab*get_bv(l,k,j,4);
				//fbby=fzgy -etab*get_bv(l,k,j,5);
				//fbbz=fzgz -etab*get_bv(l,k,j,6);
				
				//for Lattice Uni
				fbbx=fzgx -2./(1.+fluidw)/tt*get_bv(l,k,j,4);
				fbby=fzgy -2./(1.+fluidw)/tt*get_bv(l,k,j,5);
				fbbz=fzgz -2./(1.+fluidw)/tt*get_bv(l,k,j,6);
				
				//zero shift gauge
				//fbx=0.;
				//fby=0.;
				//fbz=0.;
				//fbbx=0.;
				//fbby=0.;
				//fbbz=0.;
				
				set_dbv(l,k,j,0)=get_dbv(l,k,j,0) +falpha;
				
				set_dbv(l,k,j,1)=get_dbv(l,k,j,1) +fbx;
				set_dbv(l,k,j,2)=get_dbv(l,k,j,2) +fby;
				set_dbv(l,k,j,3)=get_dbv(l,k,j,3) +fbz;

				set_dbv(l,k,j,4)=get_dbv(l,k,j,4) +fbbx;
				set_dbv(l,k,j,5)=get_dbv(l,k,j,5) +fbby;
				set_dbv(l,k,j,6)=get_dbv(l,k,j,6) +fbbz;
						
				if(itype==1)
				{
					////////////////////////////////////// constraints start //////////////////////////////////////
					double hamc,nham;			//hamiltonian constraint, normalized hamiltonian constraint

					//\tilde A_ij,k=\del_k \tilde A_ij for constraint and curvature invariants
					double daxx_x,dayy_x,dazz_x,daxy_x,daxz_x,dayz_x,
						daxx_y,dayy_y,dazz_y,daxy_y,daxz_y,dayz_y,
						daxx_z,dayy_z,dazz_z,daxy_z,daxz_z,dayz_z;
					
					
					//\cal D_i \tilde A_jk for constraints
					double Daxx_x,Dayy_x,Dazz_x,Daxy_x,Daxz_x,Dayz_x;
					double Daxx_y,Dayy_y,Dazz_y,Daxy_y,Daxz_y,Dayz_y;
					double Daxx_z,Dayy_z,Dazz_z,Daxy_z,Daxz_z,Dayz_z;
					
					//\tilde D^1 \tilde A_x1 for constraints
					double Da_x,Da_y,Da_z;
					
					//momentum constraints and normalized ones
					double M_x,M_y,M_z,normM_x,normM_y,normM_z;
					//double M_ux,M_uy,M_uz;
					//double normM_ux,normM_uy,normM_uz;
					

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
						
						//-gamma0_x*akxx_p-gamma0_y*akxy_p-gamma0_z*akxz_p;
						-zgx_p*akxx_p-zgy_p*akxy_p-zgz_p*akxz_p;
						
					
					Da_y=gixx_p*Daxy_x+gixy_p*Daxy_y+gixz_p*Daxy_z
						+giyx_p*Dayy_x+giyy_p*Dayy_y+giyz_p*Dayy_z
						+gizx_p*Dayz_x+gizy_p*Dayz_y+gizz_p*Dayz_z
						
						-Del_ux_xy*akx_ux_p-Del_ux_yy*akx_uy_p-Del_ux_zy*akx_uz_p
						-Del_uy_xy*aky_ux_p-Del_uy_yy*aky_uy_p-Del_uy_zy*aky_uz_p
						-Del_uz_xy*akz_ux_p-Del_uz_yy*akz_uy_p-Del_uz_zy*akz_uz_p
						
						//-gamma0_x*akyx_p-gamma0_y*akyy_p-gamma0_z*akyz_p;
						-zgx_p*akyx_p-zgy_p*akyy_p-zgz_p*akyz_p;
						
					Da_z=gixx_p*Daxz_x+gixy_p*Daxz_y+gixz_p*Daxz_z
						+giyx_p*Dayz_x+giyy_p*Dayz_y+giyz_p*Dayz_z
						+gizx_p*Dazz_x+gizy_p*Dazz_y+gizz_p*Dazz_z
					
						-Del_ux_xz*akx_ux_p-Del_ux_yz*akx_uy_p-Del_ux_zz*akx_uz_p
						-Del_uy_xz*aky_ux_p-Del_uy_yz*aky_uy_p-Del_uy_zz*aky_uz_p
						-Del_uz_xz*akz_ux_p-Del_uz_yz*akz_uy_p-Del_uz_zz*akz_uz_p
						
						-zgx_p*akzx_p-zgy_p*akzy_p-zgz_p*akzz_p;
					
					M_x=Da_x+6.*(akx_ux_p*wa_x_p+akx_uy_p*wa_y_p+akx_uz_p*wa_z_p)-2./3.*ek_x_p-pi8*p_x;
					M_y=Da_y+6.*(aky_ux_p*wa_x_p+aky_uy_p*wa_y_p+aky_uz_p*wa_z_p)-2./3.*ek_y_p-pi8*p_y;
					M_z=Da_z+6.*(akz_ux_p*wa_x_p+akz_uy_p*wa_y_p+akz_uz_p*wa_z_p)-2./3.*ek_z_p-pi8*p_z;
					
					normM_x=abs(Da_x)+6.*(abs(akx_ux_p*wa_x_p)+abs(akx_uy_p*wa_y_p)+abs(akx_uz_p*wa_z_p))+2./3.*abs(ek_x_p)+pi8*abs(p_x);
					normM_y=abs(Da_y)+6.*(abs(aky_ux_p*wa_x_p)+abs(aky_uy_p*wa_y_p)+abs(aky_uz_p*wa_z_p))+2./3.*abs(ek_y_p)+pi8*abs(p_y);
					normM_z=abs(Da_z)+6.*(abs(akz_ux_p*wa_x_p)+abs(akz_uy_p*wa_y_p)+abs(akz_uz_p*wa_z_p))+2./3.*abs(ek_z_p)+pi8*abs(p_z);

					//M_ux=gixx_p*M_x+gixy_p*M_y+gixz_p*M_z;
					//M_uy=giyx_p*M_x+giyy_p*M_y+giyz_p*M_z;
					//M_uz=gizx_p*M_x+gizy_p*M_y+gizz_p*M_z;
					
					//normM_ux=abs(gixx_p*normM_x)+abs(gixy_p*normM_y)+abs(gixz_p*normM_z);
					//normM_uy=abs(giyx_p*normM_x)+abs(giyy_p*normM_y)+abs(giyz_p*normM_z);
					//normM_uz=abs(gizx_p*normM_x)+abs(gizy_p*normM_y)+abs(gizz_p*normM_z);
					
						
					//hamc = ricci + 2.*pow(ek_p,2)/3. - aaaa - 2.*lambda - pi16*get_enemom(l,k,j,0);
					hamc = ricci + 2.*pow(ek_p,2)/3. - aaaa - 2.*lambda-pi16*Ene;
					
					double ricci_psi=-8.*(wawa+wa_lap)*ewa4i;
					
					//nham=abs(wa_p*dxi*dxi*exp(-4.*wa_p))+abs(pi16*get_enemom(l,k,j,0));
					//nham=abs(wa_p*dxi*dxi*exp(-4.*wa_p))+2.*pow(ek_p,2)/3.+1.+abs(ricci-ricci_psi)+abs(ricci_psi);
					
					//nham=2.*pow(ek_p,2)/3.+1.+abs(ricci-ricci_psi)+abs(ricci_psi)+abs(wa_p*dxi*dxi*exp(-4.*wa_p));
					nham=2.*pow(ek_p,2)/3.+1.+abs(ricci-ricci_psi)+(abs(8.*wawa)+abs(8.*wa_lap))*ewa4i+abs(aaaa)+pi16*abs(Ene);
					
					set_con(l,k,j,0)=hamc/nham;
					
					set_con(l,k,j,1)=hamc;
					set_con(l,k,j,2)=M_x/(normM_x+1.);
					set_con(l,k,j,3)=M_x;
					set_con(l,k,j,4)=M_y/(normM_y+1.);
					set_con(l,k,j,5)=M_y;
					set_con(l,k,j,6)=M_z/(normM_z+1.);
					set_con(l,k,j,7)=M_z;
					set_con(l,k,j,8)=gamma0_x-zgx_p;
					set_con(l,k,j,9)=gamma0_y-zgy_p;
					set_con(l,k,j,10)=gamma0_z-zgz_p;
					
					////////////////////////////////////// constraints end /////////////////////////////////////////
				
					if(curveval)
					{
						
						//////////////////////////////////////////////
						//calculation of Kretschmann invariant start//
						//////////////////////////////////////////////
						
						double Kinv;
						//gamma gamma gamma gamma R4
						double tanRiem_xyxy,tanRiem_xyxz,tanRiem_xyyz,
								tanRiem_xzxz,tanRiem_xzyz,tanRiem_yzyz;
								//tanRiem_xzxy,tanRiem_yzxy,tanRiem_yzxz,
						//upper indeces
						double utanRiem_xyxy,utanRiem_xyxz,utanRiem_xyyz,
								utanRiem_xzxz,utanRiem_xzyz,utanRiem_yzyz,
								utanRiem_xzxy,utanRiem_yzxy,utanRiem_yzxz;
						//gamma gamma gamma n R4
						double nRiem_xy_x,nRiem_xz_x,nRiem_yz_x,
								nRiem_xy_y,nRiem_xz_y,nRiem_yz_y,
								nRiem_xy_z,nRiem_xz_z,nRiem_yz_z;
						//upper indeces
						double unRiem_xy_x,unRiem_xz_x,unRiem_yz_x,
								unRiem_xy_y,unRiem_xz_y,unRiem_yz_y,
								unRiem_xy_z,unRiem_xz_z,unRiem_yz_z;
						//gamma n gamma n gamma R4
						double nnRiem_xx,nnRiem_yy,nnRiem_zz,
								nnRiem_xy,nnRiem_xz,nnRiem_yz;
						//upper indeces
						double unnRiem_xx,unnRiem_yy,unnRiem_zz,
								unnRiem_xy,unnRiem_xz,unnRiem_yz;
						//R3
						double Riem3_xyxy,Riem3_xyxz,Riem3_xyyz,
								Riem3_xzxz,Riem3_xzyz,Riem3_yzyz;
						//upper indeces
						double uRiem3_xyxy,uRiem3_xyxz,uRiem3_xyyz,
								uRiem3_xzxz,uRiem3_xzyz,uRiem3_yzyz;
						//K_ij 
						double exK_xx,exK_xy,exK_xz,exK_yx,exK_yy,exK_yz,exK_zx,exK_zy,exK_zz;
						//K_ij,k=\del_k K_ij
						double /*exK_xx_x,*/exK_xy_x,exK_xz_x,exK_yx_x,exK_yy_x,exK_yz_x,exK_zx_x,exK_zy_x,exK_zz_x,
								exK_xx_y,exK_xy_y,exK_xz_y,/*exK_yy_y,*/exK_yz_y,exK_zx_y,exK_zy_y,exK_zz_y,
								exK_xx_z,exK_xy_z,exK_xz_z,exK_yx_z,exK_yy_z,exK_yz_z/*,exK_zz_z*/;
						//Gamma^i_jk
						double Gamma_x_xx,Gamma_x_yy,Gamma_x_xy,Gamma_x_xz,Gamma_x_yz,Gamma_x_zz,
								Gamma_y_xx,Gamma_y_yy,Gamma_y_xy,Gamma_y_xz,Gamma_y_yz,Gamma_y_zz,
								Gamma_z_xx,Gamma_z_yy,Gamma_z_xy,Gamma_z_xz,Gamma_z_yz,Gamma_z_zz;
						//Gamma^l_ij K_lk
						double /*GammaK_xx_x,*/GammaK_yy_x,GammaK_zz_x,GammaK_xy_x,GammaK_xz_x,GammaK_yz_x,GammaK_yx_x,GammaK_zx_x,GammaK_zy_x,
								GammaK_xx_y,/*GammaK_yy_y,*/GammaK_zz_y,GammaK_xy_y,GammaK_xz_y,GammaK_yz_y,/*GammaK_yx_y,*/GammaK_zx_y,GammaK_zy_y,
								GammaK_xx_z,GammaK_yy_z,/*GammaK_zz_z,*/GammaK_xy_z,GammaK_xz_z,GammaK_yz_z,GammaK_yx_z/*,GammaK_zx_z,GammaK_zy_z*/;
						//D^i wa
						double wa_ux,wa_uy,wa_uz;
						//K^i_j
						double exKu_xx,exKu_xy,exKu_xz,exKu_yx,exKu_yy,exKu_yz,exKu_zx,exKu_zy,exKu_zz;
						//Ki^k K_kj
						double KK_xx,KK_yy,KK_zz,KK_xy,KK_xz,KK_yz;
						//R3^ij
						double uR_xx,uR_yy,uR_zz,uR_xy,uR_xz,uR_yz,uR_yx,uR_zx,uR_zy;
						//K^ij
						double uK_xx,uK_yy,uK_zz,uK_xy,uK_xz,uK_yz,uK_yx,uK_zx,uK_zy;
						//upper RC_ij
						double uRC_xx,uRC_yy,uRC_zz,uRC_xy,uRC_xz,uRC_yz;
						//RC_ij
						double RC_xx,RC_yy,RC_zz,RC_xy,RC_xz,RC_yz;
						//RC
						double RC;
						//RC_i
						double RC_x,RC_y,RC_z;
						double RijRij;
						double RR;
						
						//extrinsic curvature
						exK_xx=(akxx_p+gxx_p*ek_p/3.)/ewa4i;
						exK_xy=(akxy_p+gxy_p*ek_p/3.)/ewa4i;
						exK_xz=(akxz_p+gxz_p*ek_p/3.)/ewa4i;
						exK_yx=(akyx_p+gyx_p*ek_p/3.)/ewa4i;
						exK_yy=(akyy_p+gyy_p*ek_p/3.)/ewa4i;
						exK_yz=(akyz_p+gyz_p*ek_p/3.)/ewa4i;
						exK_zx=(akzx_p+gzx_p*ek_p/3.)/ewa4i;
						exK_zy=(akzy_p+gzy_p*ek_p/3.)/ewa4i;
						exK_zz=(akzz_p+gzz_p*ek_p/3.)/ewa4i;
						
						//3-dim Riemann tensor Riem3_1234=exp(8.*wa_p)*(g13_p*rc_42_p-g14_p*rc_32_p+g24_p*rc_31_p-g23_p*rc_41_p
						//								-0.5*ricci*(g13_p*g42_p-g14_p*g32_p))
						Riem3_xyxy=(gxx_p*rc_yy_p-gxy_p*rc_xy_p+gyy_p*rc_xx_p-gxy_p*rc_xy_p-0.5*ricci*(gxx_p*gyy_p-gxy_p*gxy_p))/ewa8i;
						Riem3_xyxz=(gxx_p*rc_zy_p-gxz_p*rc_xy_p+gyz_p*rc_xx_p-gxy_p*rc_zx_p-0.5*ricci*(gxx_p*gzy_p-gxz_p*gxy_p))/ewa8i;
						Riem3_xyyz=(gxy_p*rc_zy_p-gxz_p*rc_yy_p+gyz_p*rc_yx_p-gxy_p*rc_zx_p-0.5*ricci*(gxy_p*gzy_p-gxz_p*gyy_p))/ewa8i;
						//Riem3_xzxy=Riem3_xyxz;
						Riem3_xzxz=(gxx_p*rc_zz_p-gxz_p*rc_xz_p+gzz_p*rc_xx_p-gzx_p*rc_zx_p-0.5*ricci*(gxx_p*gzz_p-gxz_p*gxz_p))/ewa8i;
						Riem3_xzyz=(gxy_p*rc_zz_p-gxz_p*rc_yz_p+gzz_p*rc_yx_p-gzy_p*rc_zx_p-0.5*ricci*(gxy_p*gzz_p-gxz_p*gyz_p))/ewa8i;
						//Riem3_yzxy=Riem3_xyyz;
						//Riem3_yzxz=Riem3_xzyz;
						Riem3_yzyz=(gyy_p*rc_zz_p-gyz_p*rc_yz_p+gzz_p*rc_yy_p-gzy_p*rc_zy_p-0.5*ricci*(gyy_p*gzz_p-gyz_p*gyz_p))/ewa8i;
						
						//K_ij,k=exK_12_3=4.*wa_3_p*exK_12+exp(4.*wa_p)*(da12_3+d12_3_p*ek_p/3.+g12_p*ek_3_p/3.)
						//exK_xx_x=4.*wa_x_p*exK_xx+exp(4.*wa_p)*(daxx_x+2.*dxx_x_p*ek_p/3.+gxx_p*ek_x_p/3.);
						exK_yy_x=4.*wa_x_p*exK_yy+(dayy_x+2.*dyy_x_p*ek_p/3.+gyy_p*ek_x_p/3.)/ewa4i;
						exK_zz_x=4.*wa_x_p*exK_zz+(dazz_x+2.*dzz_x_p*ek_p/3.+gzz_p*ek_x_p/3.)/ewa4i;
						exK_xy_x=4.*wa_x_p*exK_xy+(daxy_x+2.*dxy_x_p*ek_p/3.+gxy_p*ek_x_p/3.)/ewa4i;
						exK_xz_x=4.*wa_x_p*exK_xz+(daxz_x+2.*dxz_x_p*ek_p/3.+gxz_p*ek_x_p/3.)/ewa4i;
						exK_yz_x=4.*wa_x_p*exK_yz+(dayz_x+2.*dyz_x_p*ek_p/3.+gyz_p*ek_x_p/3.)/ewa4i;
						exK_yx_x=exK_xy_x;
						exK_zx_x=exK_xz_x;
						exK_zy_x=exK_yz_x;
						
						exK_xx_y=4.*wa_y_p*exK_xx+(daxx_y+2.*dxx_y_p*ek_p/3.+gxx_p*ek_y_p/3.)/ewa4i;
						//exK_yy_y=4.*wa_y_p*exK_yy+(dayy_y+2.*dyy_y_p*ek_p/3.+gyy_p*ek_y_p/3.)/ewa4i;
						exK_zz_y=4.*wa_y_p*exK_zz+(dazz_y+2.*dzz_y_p*ek_p/3.+gzz_p*ek_y_p/3.)/ewa4i;
						exK_xy_y=4.*wa_y_p*exK_xy+(daxy_y+2.*dxy_y_p*ek_p/3.+gxy_p*ek_y_p/3.)/ewa4i;
						exK_xz_y=4.*wa_y_p*exK_xz+(daxz_y+2.*dxz_y_p*ek_p/3.+gxz_p*ek_y_p/3.)/ewa4i;
						exK_yz_y=4.*wa_y_p*exK_yz+(dayz_y+2.*dyz_y_p*ek_p/3.+gyz_p*ek_y_p/3.)/ewa4i;
						exK_zx_y=exK_xz_y;
						exK_zy_y=exK_yz_y;
						
						exK_xx_z=4.*wa_z_p*exK_xx+(daxx_z+2.*dxx_z_p*ek_p/3.+gxx_p*ek_z_p/3.)/ewa4i;
						exK_yy_z=4.*wa_z_p*exK_yy+(dayy_z+2.*dyy_z_p*ek_p/3.+gyy_p*ek_z_p/3.)/ewa4i;
						//exK_zz_z=4.*wa_z_p*exK_zz+(dazz_z+2.*dzz_z_p*ek_p/3.+gzz_p*ek_z_p/3.)/ewa4i;
						exK_xy_z=4.*wa_z_p*exK_xy+(daxy_z+2.*dxy_z_p*ek_p/3.+gxy_p*ek_z_p/3.)/ewa4i;
						exK_xz_z=4.*wa_z_p*exK_xz+(daxz_z+2.*dxz_z_p*ek_p/3.+gxz_p*ek_z_p/3.)/ewa4i;
						exK_yz_z=4.*wa_z_p*exK_yz+(dayz_z+2.*dyz_z_p*ek_p/3.+gyz_p*ek_z_p/3.)/ewa4i;
						exK_yx_z=exK_xy_z;
						
						//wa_u1=gi1x_p*wa_x_p+gi1y_p*wa_y_p+gi1z_p*wa_z_p
						wa_ux=gixx_p*wa_x_p+gixy_p*wa_y_p+gixz_p*wa_z_p;
						wa_uy=giyx_p*wa_x_p+giyy_p*wa_y_p+giyz_p*wa_z_p;
						wa_uz=gizx_p*wa_x_p+gizy_p*wa_y_p+gizz_p*wa_z_p;
						
						//Gamma_1_23=cr1_23+2.*(delta13*wa_1_p+delta12*wa_2_p-g23_p*wa_u1)
						Gamma_x_xx=crx_xx+2.*(wa_x_p+wa_x_p-gxx_p*wa_ux);
						Gamma_x_xy=crx_xy+2.*(wa_y_p-gxy_p*wa_ux);
						Gamma_x_xz=crx_xz+2.*(wa_z_p-gxz_p*wa_ux);
						Gamma_x_yz=crx_yz+2.*(-gyz_p*wa_ux);
						Gamma_x_yy=crx_yy+2.*(-gyy_p*wa_ux);
						Gamma_x_zz=crx_zz+2.*(-gzz_p*wa_ux);
						
						Gamma_y_xx=cry_xx+2.*(-gxx_p*wa_uy);
						Gamma_y_yy=cry_yy+2.*(wa_y_p+wa_y_p-gyy_p*wa_uy);
						Gamma_y_zz=cry_zz+2.*(-gzz_p*wa_uy);
						Gamma_y_xy=cry_xy+2.*(wa_x_p-gxy_p*wa_uy);
						Gamma_y_xz=cry_xz+2.*(-gxz_p*wa_uy);
						Gamma_y_yz=cry_yz+2.*(wa_z_p-gyz_p*wa_uy);
						
						Gamma_z_xx=crz_xx+2.*(-gxx_p*wa_uz);
						Gamma_z_yy=crz_yy+2.*(-gyy_p*wa_uz);
						Gamma_z_zz=crz_zz+2.*(2.*wa_z_p-gzz_p*wa_uz);
						Gamma_z_xy=crz_xy+2.*(-gxy_p*wa_uz);
						Gamma_z_xz=crz_xz+2.*(wa_x_p-gxz_p*wa_uz);
						Gamma_z_yz=crz_yz+2.*(wa_y_p-gyz_p*wa_uz);
						
						//GammaK_12_3=Gamma_x_12*exK_x3+Gamma_y_12*exK_y3+Gamma_z_12*exK_z3
						//GammaK_xx_x=Gamma_x_xx*exK_xx+Gamma_y_xx*exK_yx+Gamma_z_xx*exK_zx;
						GammaK_yy_x=Gamma_x_yy*exK_xx+Gamma_y_yy*exK_yx+Gamma_z_yy*exK_zx;
						GammaK_zz_x=Gamma_x_zz*exK_xx+Gamma_y_zz*exK_yx+Gamma_z_zz*exK_zx;
						GammaK_xy_x=Gamma_x_xy*exK_xx+Gamma_y_xy*exK_yx+Gamma_z_xy*exK_zx;
						GammaK_xz_x=Gamma_x_xz*exK_xx+Gamma_y_xz*exK_yx+Gamma_z_xz*exK_zx;
						GammaK_yz_x=Gamma_x_yz*exK_xx+Gamma_y_yz*exK_yx+Gamma_z_yz*exK_zx;
						GammaK_yx_x=GammaK_xy_x;
						GammaK_zx_x=GammaK_xz_x;
						GammaK_zy_x=GammaK_yz_x;
					
						GammaK_xx_y=Gamma_x_xx*exK_xy+Gamma_y_xx*exK_yy+Gamma_z_xx*exK_zy;
						//GammaK_yy_y=Gamma_x_yy*exK_xy+Gamma_y_yy*exK_yy+Gamma_z_yy*exK_zy;
						GammaK_zz_y=Gamma_x_zz*exK_xy+Gamma_y_zz*exK_yy+Gamma_z_zz*exK_zy;
						GammaK_xy_y=Gamma_x_xy*exK_xy+Gamma_y_xy*exK_yy+Gamma_z_xy*exK_zy;
						GammaK_xz_y=Gamma_x_xz*exK_xy+Gamma_y_xz*exK_yy+Gamma_z_xz*exK_zy;
						GammaK_yz_y=Gamma_x_yz*exK_xy+Gamma_y_yz*exK_yy+Gamma_z_yz*exK_zy;
						//GammaK_yx_y=GammaK_xy_y;
						GammaK_zx_y=GammaK_xz_y;
						GammaK_zy_y=GammaK_yz_y;

						//GammaK_12_3=Gamma_x_12*exK_x3+Gamma_y_12*exK_y3+Gamma_z_12*exK_z3;
						GammaK_xx_z=Gamma_x_xx*exK_xz+Gamma_y_xx*exK_yz+Gamma_z_xx*exK_zz;
						GammaK_yy_z=Gamma_x_yy*exK_xz+Gamma_y_yy*exK_yz+Gamma_z_yy*exK_zz;
						//GammaK_zz_z=Gamma_x_zz*exK_xz+Gamma_y_zz*exK_yz+Gamma_z_zz*exK_zz;
						GammaK_xy_z=Gamma_x_xy*exK_xz+Gamma_y_xy*exK_yz+Gamma_z_xy*exK_zz;
						GammaK_xz_z=Gamma_x_xz*exK_xz+Gamma_y_xz*exK_yz+Gamma_z_xz*exK_zz;
						GammaK_yz_z=Gamma_x_yz*exK_xz+Gamma_y_yz*exK_yz+Gamma_z_yz*exK_zz;
						GammaK_yx_z=GammaK_xy_z;
						//GammaK_zx_z=GammaK_xz_z;
						//GammaK_zy_z=GammaK_yz_z;


						//tanRiem_1234=Riem3_1234+exK_13*exK_24-exK_23*exK_14;
						tanRiem_xyxy=Riem3_xyxy+exK_xx*exK_yy-exK_yx*exK_xy;
						tanRiem_xyxz=Riem3_xyxz+exK_xx*exK_yz-exK_yx*exK_xz;
						tanRiem_xyyz=Riem3_xyyz+exK_xy*exK_yz-exK_yy*exK_xz;
						//tanRiem_xzxy=tanRiem_xyxz;
						tanRiem_xzxz=Riem3_xzxz+exK_xx*exK_zz-exK_zx*exK_xz;
						tanRiem_xzyz=Riem3_xzyz+exK_xy*exK_zz-exK_zy*exK_xz;
						//tanRiem_yzxy=tanRiem_xyyz;
						//tanRiem_yzxz=tanRiem_xzyz;
						tanRiem_yzyz=Riem3_yzyz+exK_yy*exK_zz-exK_zy*exK_yz;
						
						
						//nRiem_12_3=-exK_23_1+exK_13_2+GammaK_13_2-GammaK_23_1;
						nRiem_xy_x=-exK_yx_x+exK_xx_y+GammaK_xx_y-GammaK_yx_x;
						nRiem_xz_x=-exK_zx_x+exK_xx_z+GammaK_xx_z-GammaK_zx_x;
						nRiem_yz_x=-exK_zx_y+exK_yx_z+GammaK_yx_z-GammaK_zx_y;
						nRiem_xy_y=-exK_yy_x+exK_xy_y+GammaK_xy_y-GammaK_yy_x;
						nRiem_xz_y=-exK_zy_x+exK_xy_z+GammaK_xy_z-GammaK_zy_x;
						nRiem_yz_y=-exK_zy_y+exK_yy_z+GammaK_yy_z-GammaK_zy_y;
						nRiem_xy_z=-exK_yz_x+exK_xz_y+GammaK_xz_y-GammaK_yz_x;
						nRiem_xz_z=-exK_zz_x+exK_xz_z+GammaK_xz_z-GammaK_zz_x;
						nRiem_yz_z=-exK_zz_y+exK_yz_z+GammaK_yz_z-GammaK_zz_y;
						
						
						//exKu_12=gi1x_p*exK_x2+gi1y_p*exKy2+gi1z_p*exKz2;
						exKu_xx=(gixx_p*exK_xx+gixy_p*exK_yx+gixz_p*exK_zx)*ewa4i;
						exKu_yy=(giyx_p*exK_xy+giyy_p*exK_yy+giyz_p*exK_zy)*ewa4i;
						exKu_zz=(gizx_p*exK_xz+gizy_p*exK_yz+gizz_p*exK_zz)*ewa4i;
						exKu_xy=(gixx_p*exK_xy+gixy_p*exK_yy+gixz_p*exK_zy)*ewa4i;
						exKu_xz=(gixx_p*exK_xz+gixy_p*exK_yz+gixz_p*exK_zz)*ewa4i;
						exKu_yx=(giyx_p*exK_xx+giyy_p*exK_yx+giyz_p*exK_zx)*ewa4i;
						exKu_yz=(giyx_p*exK_xz+giyy_p*exK_yz+giyz_p*exK_zz)*ewa4i;
						exKu_zx=(gizx_p*exK_xx+gizy_p*exK_yx+gizz_p*exK_zx)*ewa4i;
						exKu_zy=(gizx_p*exK_xy+gizy_p*exK_yy+gizz_p*exK_zy)*ewa4i;
						
						
						//KK_12=exK_1x*exKu_x2+exK_1y*exKu_y2+exK_1z*exKu_z2;
						KK_xx=exK_xx*exKu_xx+exK_xy*exKu_yx+exK_xz*exKu_zx;
						KK_yy=exK_yx*exKu_xy+exK_yy*exKu_yy+exK_yz*exKu_zy;
						KK_zz=exK_zx*exKu_xz+exK_zy*exKu_yz+exK_zz*exKu_zz;
						KK_xy=exK_xx*exKu_xy+exK_xy*exKu_yy+exK_xz*exKu_zy;
						KK_xz=exK_xx*exKu_xz+exK_xy*exKu_yz+exK_xz*exKu_zz;
						KK_yz=exK_yx*exKu_xz+exK_yy*exKu_yz+exK_yz*exKu_zz;
						//KK_yx=KK_xy;
						//KK_zx=KK_xz;
						//KK_zy=KK_yz;
						
						
						//nnRiem_12=-KK_12+exp(4.*wa_p)*rc_12_p+ek_p*exK_12+pi4*((trs-Ene)*g12_p-2.*s12);
						nnRiem_xx=-KK_xx+rc_xx_p/ewa4i+ek_p*exK_xx+pi4*((trs-Ene)*gxx_p-2.*sxx);
						nnRiem_yy=-KK_yy+rc_yy_p/ewa4i+ek_p*exK_yy+pi4*((trs-Ene)*gyy_p-2.*syy);
						nnRiem_zz=-KK_zz+rc_zz_p/ewa4i+ek_p*exK_zz+pi4*((trs-Ene)*gzz_p-2.*szz);
						nnRiem_xy=-KK_xy+rc_xy_p/ewa4i+ek_p*exK_xy+pi4*((trs-Ene)*gxy_p-2.*sxy);
						nnRiem_xz=-KK_xz+rc_xz_p/ewa4i+ek_p*exK_xz+pi4*((trs-Ene)*gxz_p-2.*sxz);
						nnRiem_yz=-KK_yz+rc_yz_p/ewa4i+ek_p*exK_yz+pi4*((trs-Ene)*gyz_p-2.*syz);

						
						//unnRiem_12=exp(-8.*wa_p)*(
						//			gi1x_p*gi2x_p*nnRiem_xx+gi1y_p*gi2y_p*nnRiem_yy+gi1z_p*gi2z_p*nnRiem_zz
						//			+gi1x_p*gi2y_p*nnRiem_xy+gi1x_p*gi2z_p*nnRiem_xz+gi1y_p*gi2z_p*nnRiem_yz
						//			+gi1y_p*gi2x_p*nnRiem_xy+gi1z_p*gi2x_p*nnRiem_xz+gi1z_p*gi2y_p*nnRiem_yz);
						unnRiem_xx=ewa8i*(
									gixx_p*gixx_p*nnRiem_xx+gixy_p*gixy_p*nnRiem_yy+gixz_p*gixz_p*nnRiem_zz
									+gixx_p*gixy_p*nnRiem_xy+gixx_p*gixz_p*nnRiem_xz+gixy_p*gixz_p*nnRiem_yz
									+gixy_p*gixx_p*nnRiem_xy+gixz_p*gixx_p*nnRiem_xz+gixz_p*gixy_p*nnRiem_yz);
						
						unnRiem_yy=ewa8i*(
									giyx_p*giyx_p*nnRiem_xx+giyy_p*giyy_p*nnRiem_yy+giyz_p*giyz_p*nnRiem_zz
									+giyx_p*giyy_p*nnRiem_xy+giyx_p*giyz_p*nnRiem_xz+giyy_p*giyz_p*nnRiem_yz
									+giyy_p*giyx_p*nnRiem_xy+giyz_p*giyx_p*nnRiem_xz+giyz_p*giyy_p*nnRiem_yz);
						
						unnRiem_zz=ewa8i*(
									gizx_p*gizx_p*nnRiem_xx+gizy_p*gizy_p*nnRiem_yy+gizz_p*gizz_p*nnRiem_zz
									+gizx_p*gizy_p*nnRiem_xy+gizx_p*gizz_p*nnRiem_xz+gizy_p*gizz_p*nnRiem_yz
									+gizy_p*gizx_p*nnRiem_xy+gizz_p*gizx_p*nnRiem_xz+gizz_p*gizy_p*nnRiem_yz);
						
						unnRiem_xy=ewa8i*(
									gixx_p*giyx_p*nnRiem_xx+gixy_p*giyy_p*nnRiem_yy+gixz_p*giyz_p*nnRiem_zz
									+gixx_p*giyy_p*nnRiem_xy+gixx_p*giyz_p*nnRiem_xz+gixy_p*giyz_p*nnRiem_yz
									+gixy_p*giyx_p*nnRiem_xy+gixz_p*giyx_p*nnRiem_xz+gixz_p*giyy_p*nnRiem_yz);
						
						unnRiem_xz=ewa8i*(
									gixx_p*gizx_p*nnRiem_xx+gixy_p*gizy_p*nnRiem_yy+gixz_p*gizz_p*nnRiem_zz
									+gixx_p*gizy_p*nnRiem_xy+gixx_p*gizz_p*nnRiem_xz+gixy_p*gizz_p*nnRiem_yz
									+gixy_p*gizx_p*nnRiem_xy+gixz_p*gizx_p*nnRiem_xz+gixz_p*gizy_p*nnRiem_yz);
						
						unnRiem_yz=ewa8i*(
									giyx_p*gizx_p*nnRiem_xx+giyy_p*gizy_p*nnRiem_yy+giyz_p*gizz_p*nnRiem_zz
									+giyx_p*gizy_p*nnRiem_xy+giyx_p*gizz_p*nnRiem_xz+giyy_p*gizz_p*nnRiem_yz
									+giyy_p*gizx_p*nnRiem_xy+giyz_p*gizx_p*nnRiem_xz+giyz_p*gizy_p*nnRiem_yz);
						
						
						//unRiem_12_3=exp(-12.*wa_p)*(
						//			+gi1x_p*gi2y_p*gi3x_p*nRiem_xy_x+gi1x_p*gi2z_p*gi3x_p*nRiem_xz_x+gi1y_p*gi2z_p*gi3x_p*nRiem_yz_x
						//			-gi1y_p*gi2x_p*gi3x_p*nRiem_xy_x-gi1z_p*gi2x_p*gi3x_p*nRiem_xz_x-gi1z_p*gi2y_p*gi3x_p*nRiem_yz_x
						//			+gi1x_p*gi2y_p*gi3y_p*nRiem_xy_y+gi1x_p*gi2z_p*gi3y_p*nRiem_xz_y+gi1y_p*gi2z_p*gi3y_p*nRiem_yz_y
						//			-gi1y_p*gi2x_p*gi3y_p*nRiem_xy_y-gi1z_p*gi2x_p*gi3y_p*nRiem_xz_y-gi1z_p*gi2y_p*gi3y_p*nRiem_yz_y
						//			+gi1x_p*gi2y_p*gi3z_p*nRiem_xy_z+gi1x_p*gi2z_p*gi3z_p*nRiem_xz_z+gi1y_p*gi2z_p*gi3z_p*nRiem_yz_z
						//			-gi1y_p*gi2x_p*gi3z_p*nRiem_xy_z-gi1z_p*gi2x_p*gi3z_p*nRiem_xz_z-gi1z_p*gi2y_p*gi3z_p*nRiem_yz_z);
						unRiem_xy_x=ewa12i*(
									+gixx_p*giyy_p*gixx_p*nRiem_xy_x+gixx_p*giyz_p*gixx_p*nRiem_xz_x+gixy_p*giyz_p*gixx_p*nRiem_yz_x
									-gixy_p*giyx_p*gixx_p*nRiem_xy_x-gixz_p*giyx_p*gixx_p*nRiem_xz_x-gixz_p*giyy_p*gixx_p*nRiem_yz_x
									+gixx_p*giyy_p*gixy_p*nRiem_xy_y+gixx_p*giyz_p*gixy_p*nRiem_xz_y+gixy_p*giyz_p*gixy_p*nRiem_yz_y
									-gixy_p*giyx_p*gixy_p*nRiem_xy_y-gixz_p*giyx_p*gixy_p*nRiem_xz_y-gixz_p*giyy_p*gixy_p*nRiem_yz_y
									+gixx_p*giyy_p*gixz_p*nRiem_xy_z+gixx_p*giyz_p*gixz_p*nRiem_xz_z+gixy_p*giyz_p*gixz_p*nRiem_yz_z
									-gixy_p*giyx_p*gixz_p*nRiem_xy_z-gixz_p*giyx_p*gixz_p*nRiem_xz_z-gixz_p*giyy_p*gixz_p*nRiem_yz_z);
						unRiem_xz_x=ewa12i*(
									+gixx_p*gizy_p*gixx_p*nRiem_xy_x+gixx_p*gizz_p*gixx_p*nRiem_xz_x+gixy_p*gizz_p*gixx_p*nRiem_yz_x
									-gixy_p*gizx_p*gixx_p*nRiem_xy_x-gixz_p*gizx_p*gixx_p*nRiem_xz_x-gixz_p*gizy_p*gixx_p*nRiem_yz_x
									+gixx_p*gizy_p*gixy_p*nRiem_xy_y+gixx_p*gizz_p*gixy_p*nRiem_xz_y+gixy_p*gizz_p*gixy_p*nRiem_yz_y
									-gixy_p*gizx_p*gixy_p*nRiem_xy_y-gixz_p*gizx_p*gixy_p*nRiem_xz_y-gixz_p*gizy_p*gixy_p*nRiem_yz_y
									+gixx_p*gizy_p*gixz_p*nRiem_xy_z+gixx_p*gizz_p*gixz_p*nRiem_xz_z+gixy_p*gizz_p*gixz_p*nRiem_yz_z
									-gixy_p*gizx_p*gixz_p*nRiem_xy_z-gixz_p*gizx_p*gixz_p*nRiem_xz_z-gixz_p*gizy_p*gixz_p*nRiem_yz_z);
						unRiem_yz_x=ewa12i*(
									+giyx_p*gizy_p*gixx_p*nRiem_xy_x+giyx_p*gizz_p*gixx_p*nRiem_xz_x+giyy_p*gizz_p*gixx_p*nRiem_yz_x
									-giyy_p*gizx_p*gixx_p*nRiem_xy_x-giyz_p*gizx_p*gixx_p*nRiem_xz_x-giyz_p*gizy_p*gixx_p*nRiem_yz_x
									+giyx_p*gizy_p*gixy_p*nRiem_xy_y+giyx_p*gizz_p*gixy_p*nRiem_xz_y+giyy_p*gizz_p*gixy_p*nRiem_yz_y
									-giyy_p*gizx_p*gixy_p*nRiem_xy_y-giyz_p*gizx_p*gixy_p*nRiem_xz_y-giyz_p*gizy_p*gixy_p*nRiem_yz_y
									+giyx_p*gizy_p*gixz_p*nRiem_xy_z+giyx_p*gizz_p*gixz_p*nRiem_xz_z+giyy_p*gizz_p*gixz_p*nRiem_yz_z
									-giyy_p*gizx_p*gixz_p*nRiem_xy_z-giyz_p*gizx_p*gixz_p*nRiem_xz_z-giyz_p*gizy_p*gixz_p*nRiem_yz_z);
						unRiem_xy_y=ewa12i*(
									+gixx_p*giyy_p*giyx_p*nRiem_xy_x+gixx_p*giyz_p*giyx_p*nRiem_xz_x+gixy_p*giyz_p*giyx_p*nRiem_yz_x
									-gixy_p*giyx_p*giyx_p*nRiem_xy_x-gixz_p*giyx_p*giyx_p*nRiem_xz_x-gixz_p*giyy_p*giyx_p*nRiem_yz_x
									+gixx_p*giyy_p*giyy_p*nRiem_xy_y+gixx_p*giyz_p*giyy_p*nRiem_xz_y+gixy_p*giyz_p*giyy_p*nRiem_yz_y
									-gixy_p*giyx_p*giyy_p*nRiem_xy_y-gixz_p*giyx_p*giyy_p*nRiem_xz_y-gixz_p*giyy_p*giyy_p*nRiem_yz_y
									+gixx_p*giyy_p*giyz_p*nRiem_xy_z+gixx_p*giyz_p*giyz_p*nRiem_xz_z+gixy_p*giyz_p*giyz_p*nRiem_yz_z
									-gixy_p*giyx_p*giyz_p*nRiem_xy_z-gixz_p*giyx_p*giyz_p*nRiem_xz_z-gixz_p*giyy_p*giyz_p*nRiem_yz_z);
						unRiem_xz_y=ewa12i*(
									+gixx_p*gizy_p*giyx_p*nRiem_xy_x+gixx_p*gizz_p*giyx_p*nRiem_xz_x+gixy_p*gizz_p*giyx_p*nRiem_yz_x
									-gixy_p*gizx_p*giyx_p*nRiem_xy_x-gixz_p*gizx_p*giyx_p*nRiem_xz_x-gixz_p*gizy_p*giyx_p*nRiem_yz_x
									+gixx_p*gizy_p*giyy_p*nRiem_xy_y+gixx_p*gizz_p*giyy_p*nRiem_xz_y+gixy_p*gizz_p*giyy_p*nRiem_yz_y
									-gixy_p*gizx_p*giyy_p*nRiem_xy_y-gixz_p*gizx_p*giyy_p*nRiem_xz_y-gixz_p*gizy_p*giyy_p*nRiem_yz_y
									+gixx_p*gizy_p*giyz_p*nRiem_xy_z+gixx_p*gizz_p*giyz_p*nRiem_xz_z+gixy_p*gizz_p*giyz_p*nRiem_yz_z
									-gixy_p*gizx_p*giyz_p*nRiem_xy_z-gixz_p*gizx_p*giyz_p*nRiem_xz_z-gixz_p*gizy_p*giyz_p*nRiem_yz_z);
						unRiem_yz_y=ewa12i*(
									+giyx_p*gizy_p*giyx_p*nRiem_xy_x+giyx_p*gizz_p*giyx_p*nRiem_xz_x+giyy_p*gizz_p*giyx_p*nRiem_yz_x
									-giyy_p*gizx_p*giyx_p*nRiem_xy_x-giyz_p*gizx_p*giyx_p*nRiem_xz_x-giyz_p*gizy_p*giyx_p*nRiem_yz_x
									+giyx_p*gizy_p*giyy_p*nRiem_xy_y+giyx_p*gizz_p*giyy_p*nRiem_xz_y+giyy_p*gizz_p*giyy_p*nRiem_yz_y
									-giyy_p*gizx_p*giyy_p*nRiem_xy_y-giyz_p*gizx_p*giyy_p*nRiem_xz_y-giyz_p*gizy_p*giyy_p*nRiem_yz_y
									+giyx_p*gizy_p*giyz_p*nRiem_xy_z+giyx_p*gizz_p*giyz_p*nRiem_xz_z+giyy_p*gizz_p*giyz_p*nRiem_yz_z
									-giyy_p*gizx_p*giyz_p*nRiem_xy_z-giyz_p*gizx_p*giyz_p*nRiem_xz_z-giyz_p*gizy_p*giyz_p*nRiem_yz_z);
						unRiem_xy_z=ewa12i*(
									+gixx_p*giyy_p*gizx_p*nRiem_xy_x+gixx_p*giyz_p*gizx_p*nRiem_xz_x+gixy_p*giyz_p*gizx_p*nRiem_yz_x
									-gixy_p*giyx_p*gizx_p*nRiem_xy_x-gixz_p*giyx_p*gizx_p*nRiem_xz_x-gixz_p*giyy_p*gizx_p*nRiem_yz_x
									+gixx_p*giyy_p*gizy_p*nRiem_xy_y+gixx_p*giyz_p*gizy_p*nRiem_xz_y+gixy_p*giyz_p*gizy_p*nRiem_yz_y
									-gixy_p*giyx_p*gizy_p*nRiem_xy_y-gixz_p*giyx_p*gizy_p*nRiem_xz_y-gixz_p*giyy_p*gizy_p*nRiem_yz_y
									+gixx_p*giyy_p*gizz_p*nRiem_xy_z+gixx_p*giyz_p*gizz_p*nRiem_xz_z+gixy_p*giyz_p*gizz_p*nRiem_yz_z
									-gixy_p*giyx_p*gizz_p*nRiem_xy_z-gixz_p*giyx_p*gizz_p*nRiem_xz_z-gixz_p*giyy_p*gizz_p*nRiem_yz_z);
						unRiem_xz_z=ewa12i*(
									+gixx_p*gizy_p*gizx_p*nRiem_xy_x+gixx_p*gizz_p*gizx_p*nRiem_xz_x+gixy_p*gizz_p*gizx_p*nRiem_yz_x
									-gixy_p*gizx_p*gizx_p*nRiem_xy_x-gixz_p*gizx_p*gizx_p*nRiem_xz_x-gixz_p*gizy_p*gizx_p*nRiem_yz_x
									+gixx_p*gizy_p*gizy_p*nRiem_xy_y+gixx_p*gizz_p*gizy_p*nRiem_xz_y+gixy_p*gizz_p*gizy_p*nRiem_yz_y
									-gixy_p*gizx_p*gizy_p*nRiem_xy_y-gixz_p*gizx_p*gizy_p*nRiem_xz_y-gixz_p*gizy_p*gizy_p*nRiem_yz_y
									+gixx_p*gizy_p*gizz_p*nRiem_xy_z+gixx_p*gizz_p*gizz_p*nRiem_xz_z+gixy_p*gizz_p*gizz_p*nRiem_yz_z
									-gixy_p*gizx_p*gizz_p*nRiem_xy_z-gixz_p*gizx_p*gizz_p*nRiem_xz_z-gixz_p*gizy_p*gizz_p*nRiem_yz_z);
						unRiem_yz_z=ewa12i*(
									+giyx_p*gizy_p*gizx_p*nRiem_xy_x+giyx_p*gizz_p*gizx_p*nRiem_xz_x+giyy_p*gizz_p*gizx_p*nRiem_yz_x
									-giyy_p*gizx_p*gizx_p*nRiem_xy_x-giyz_p*gizx_p*gizx_p*nRiem_xz_x-giyz_p*gizy_p*gizx_p*nRiem_yz_x
									+giyx_p*gizy_p*gizy_p*nRiem_xy_y+giyx_p*gizz_p*gizy_p*nRiem_xz_y+giyy_p*gizz_p*gizy_p*nRiem_yz_y
									-giyy_p*gizx_p*gizy_p*nRiem_xy_y-giyz_p*gizx_p*gizy_p*nRiem_xz_y-giyz_p*gizy_p*gizy_p*nRiem_yz_y
									+giyx_p*gizy_p*gizz_p*nRiem_xy_z+giyx_p*gizz_p*gizz_p*nRiem_xz_z+giyy_p*gizz_p*gizz_p*nRiem_yz_z
									-giyy_p*gizx_p*gizz_p*nRiem_xy_z-giyz_p*gizx_p*gizz_p*nRiem_xz_z-giyz_p*gizy_p*gizz_p*nRiem_yz_z);
						
						
						uR_xx=ewa4i*(
									gixx_p*gixx_p*rc_xx_p+gixy_p*gixy_p*rc_yy_p+gixz_p*gixz_p*rc_zz_p
									+gixx_p*gixy_p*rc_xy_p+gixx_p*gixz_p*rc_xz_p+gixy_p*gixz_p*rc_yz_p
									+gixy_p*gixx_p*rc_xy_p+gixz_p*gixx_p*rc_xz_p+gixz_p*gixy_p*rc_yz_p);
						uR_yy=ewa4i*(
									giyx_p*giyx_p*rc_xx_p+giyy_p*giyy_p*rc_yy_p+giyz_p*giyz_p*rc_zz_p
									+giyx_p*giyy_p*rc_xy_p+giyx_p*giyz_p*rc_xz_p+giyy_p*giyz_p*rc_yz_p
									+giyy_p*giyx_p*rc_xy_p+giyz_p*giyx_p*rc_xz_p+giyz_p*giyy_p*rc_yz_p);
						uR_zz=ewa4i*(
									gizx_p*gizx_p*rc_xx_p+gizy_p*gizy_p*rc_yy_p+gizz_p*gizz_p*rc_zz_p
									+gizx_p*gizy_p*rc_xy_p+gizx_p*gizz_p*rc_xz_p+gizy_p*gizz_p*rc_yz_p
									+gizy_p*gizx_p*rc_xy_p+gizz_p*gizx_p*rc_xz_p+gizz_p*gizy_p*rc_yz_p);
						uR_xy=ewa4i*(
									gixx_p*giyx_p*rc_xx_p+gixy_p*giyy_p*rc_yy_p+gixz_p*giyz_p*rc_zz_p
									+gixx_p*giyy_p*rc_xy_p+gixx_p*giyz_p*rc_xz_p+gixy_p*giyz_p*rc_yz_p
									+gixy_p*giyx_p*rc_xy_p+gixz_p*giyx_p*rc_xz_p+gixz_p*giyy_p*rc_yz_p);
						uR_xz=ewa4i*(
									gixx_p*gizx_p*rc_xx_p+gixy_p*gizy_p*rc_yy_p+gixz_p*gizz_p*rc_zz_p
									+gixx_p*gizy_p*rc_xy_p+gixx_p*gizz_p*rc_xz_p+gixy_p*gizz_p*rc_yz_p
									+gixy_p*gizx_p*rc_xy_p+gixz_p*gizx_p*rc_xz_p+gixz_p*gizy_p*rc_yz_p);
						uR_yz=ewa4i*(
									giyx_p*gizx_p*rc_xx_p+giyy_p*gizy_p*rc_yy_p+giyz_p*gizz_p*rc_zz_p
									+giyx_p*gizy_p*rc_xy_p+giyx_p*gizz_p*rc_xz_p+giyy_p*gizz_p*rc_yz_p
									+giyy_p*gizx_p*rc_xy_p+giyz_p*gizx_p*rc_xz_p+giyz_p*gizy_p*rc_yz_p);
						uR_yx=uR_xy;
						uR_zx=uR_xz;
						uR_zy=uR_yz;

						uK_xx=ewa8i*(
									gixx_p*gixx_p*exK_xx+gixy_p*gixy_p*exK_yy+gixz_p*gixz_p*exK_zz
									+gixx_p*gixy_p*exK_xy+gixx_p*gixz_p*exK_xz+gixy_p*gixz_p*exK_yz
									+gixy_p*gixx_p*exK_xy+gixz_p*gixx_p*exK_xz+gixz_p*gixy_p*exK_yz);
						uK_yy=ewa8i*(
									giyx_p*giyx_p*exK_xx+giyy_p*giyy_p*exK_yy+giyz_p*giyz_p*exK_zz
									+giyx_p*giyy_p*exK_xy+giyx_p*giyz_p*exK_xz+giyy_p*giyz_p*exK_yz
									+giyy_p*giyx_p*exK_xy+giyz_p*giyx_p*exK_xz+giyz_p*giyy_p*exK_yz);
						uK_zz=ewa8i*(
									gizx_p*gizx_p*exK_xx+gizy_p*gizy_p*exK_yy+gizz_p*gizz_p*exK_zz
									+gizx_p*gizy_p*exK_xy+gizx_p*gizz_p*exK_xz+gizy_p*gizz_p*exK_yz
									+gizy_p*gizx_p*exK_xy+gizz_p*gizx_p*exK_xz+gizz_p*gizy_p*exK_yz);
						uK_xy=ewa8i*(
									gixx_p*giyx_p*exK_xx+gixy_p*giyy_p*exK_yy+gixz_p*giyz_p*exK_zz
									+gixx_p*giyy_p*exK_xy+gixx_p*giyz_p*exK_xz+gixy_p*giyz_p*exK_yz
									+gixy_p*giyx_p*exK_xy+gixz_p*giyx_p*exK_xz+gixz_p*giyy_p*exK_yz);
						uK_xz=ewa8i*(
									gixx_p*gizx_p*exK_xx+gixy_p*gizy_p*exK_yy+gixz_p*gizz_p*exK_zz
									+gixx_p*gizy_p*exK_xy+gixx_p*gizz_p*exK_xz+gixy_p*gizz_p*exK_yz
									+gixy_p*gizx_p*exK_xy+gixz_p*gizx_p*exK_xz+gixz_p*gizy_p*exK_yz);
						uK_yz=ewa8i*(
									giyx_p*gizx_p*exK_xx+giyy_p*gizy_p*exK_yy+giyz_p*gizz_p*exK_zz
									+giyx_p*gizy_p*exK_xy+giyx_p*gizz_p*exK_xz+giyy_p*gizz_p*exK_yz
									+giyy_p*gizx_p*exK_xy+giyz_p*gizx_p*exK_xz+giyz_p*gizy_p*exK_yz);
						uK_yx=uK_xy;
						uK_zx=uK_xz;
						uK_zy=uK_yz;
						
						
						//uRiem3_xyxy=exp(-4.*wa_p)*(gi13_p*uR_42-gi14_p*uR_32+gi24_p*uR_31-gi23_p*uR_41)
						//								-0.5*ricci*exp(-8.*wa_p)*(gi13_p*gi42_p-gi14_p*gi32_p);
						uRiem3_xyxy=ewa4i*(gixx_p*uR_yy-gixy_p*uR_xy+giyy_p*uR_xx-gixy_p*uR_xy)-0.5*ricci*ewa8i*(gixx_p*giyy_p-gixy_p*gixy_p);
						uRiem3_xyxz=ewa4i*(gixx_p*uR_zy-gixz_p*uR_xy+giyz_p*uR_xx-gixy_p*uR_zx)-0.5*ricci*ewa8i*(gixx_p*gizy_p-gixz_p*gixy_p);
						uRiem3_xyyz=ewa4i*(gixy_p*uR_zy-gixz_p*uR_yy+giyz_p*uR_yx-gixy_p*uR_zx)-0.5*ricci*ewa8i*(gixy_p*gizy_p-gixz_p*giyy_p);
						//uRiem3_xzxy=uRiem3_xyxz;
						uRiem3_xzxz=ewa4i*(gixx_p*uR_zz-gixz_p*uR_xz+gizz_p*uR_xx-gizx_p*uR_zx)-0.5*ricci*ewa8i*(gixx_p*gizz_p-gixz_p*gixz_p);
						uRiem3_xzyz=ewa4i*(gixy_p*uR_zz-gixz_p*uR_yz+gizz_p*uR_yx-gizy_p*uR_zx)-0.5*ricci*ewa8i*(gixy_p*gizz_p-gixz_p*giyz_p);
						//uRiem3_yzxy=uRiem3_xyyz;
						//uRiem3_yzxz=uRiem3_xzyz;
						uRiem3_yzyz=ewa4i*(giyy_p*uR_zz-giyz_p*uR_yz+gizz_p*uR_yy-gizy_p*uR_zy)-0.5*ricci*ewa8i*(giyy_p*gizz_p-giyz_p*giyz_p);
						
						
						//utanRiem_1234=uRiem3_1234+uK_13*uK_24-uK_23*uK_14;
						utanRiem_xyxy=uRiem3_xyxy+uK_xx*uK_yy-uK_yx*uK_xy;
						utanRiem_xyxz=uRiem3_xyxz+uK_xx*uK_yz-uK_yx*uK_xz;
						utanRiem_xyyz=uRiem3_xyyz+uK_xy*uK_yz-uK_yy*uK_xz;
						utanRiem_xzxy=utanRiem_xyxz;
						utanRiem_xzxz=uRiem3_xzxz+uK_xx*uK_zz-uK_zx*uK_xz;
						utanRiem_xzyz=uRiem3_xzyz+uK_xy*uK_zz-uK_zy*uK_xz;
						utanRiem_yzxy=utanRiem_xyyz;
						utanRiem_yzxz=utanRiem_xzyz;
						utanRiem_yzyz=uRiem3_yzyz+uK_yy*uK_zz-uK_zy*uK_yz;
						
						//
						Kinv=4.*(
								unnRiem_xx*nnRiem_xx+unnRiem_yy*nnRiem_yy+unnRiem_zz*nnRiem_zz
								+2.*(unnRiem_xy*nnRiem_xy+unnRiem_xz*nnRiem_xz+unnRiem_yz*nnRiem_yz))
							-4.*(
								2.*(unRiem_xy_x*nRiem_xy_x+unRiem_xz_x*nRiem_xz_x+unRiem_yz_x*nRiem_yz_x
									+unRiem_xy_y*nRiem_xy_y+unRiem_xz_y*nRiem_xz_y+unRiem_yz_y*nRiem_yz_y
									+unRiem_xy_z*nRiem_xy_z+unRiem_xz_z*nRiem_xz_z+unRiem_yz_z*nRiem_yz_z))
							+4.*(
								utanRiem_xyxy*tanRiem_xyxy+utanRiem_xzxz*tanRiem_xzxz+utanRiem_yzyz*tanRiem_yzyz
									+2.*(utanRiem_xyxz*tanRiem_xyxz+utanRiem_xyyz*tanRiem_xyyz+utanRiem_xzyz*tanRiem_xzyz));
						
						set_outv(l,k,j,0)=Kinv;
						
						//upper RC_ij
						//uRC_12=exp(4.*wa_p)*(gxx_p*utanRiem_x1x2+gyy_p*utanRiem_y1y2+gzz_p*utanRiem_z1z2
						//						+gxy_p*utanRiem_x1y2+gyz_p*utanRiem_y1z2+gzx_p*utanRiem_z1x2
						//						+gxy_p*utanRiem_y1x2+gyz_p*utanRiem_z1y2+gzx_p*utanRiem_x1z2);
						
						uRC_xx=(gyy_p*utanRiem_xyxy+gzz_p*utanRiem_xzxz
												+gyz_p*utanRiem_xyxz+gyz_p*utanRiem_xzxy)/ewa4i;
						uRC_yy=(gxx_p*utanRiem_xyxy+gzz_p*utanRiem_yzyz
												-gzx_p*utanRiem_yzxy-gzx_p*utanRiem_xyyz)/ewa4i;
						uRC_zz=(gxx_p*utanRiem_xzxz+gyy_p*utanRiem_yzyz
												+gxy_p*utanRiem_xzyz+gxy_p*utanRiem_yzxz)/ewa4i;
						uRC_xy=(gzz_p*utanRiem_xzyz+gyz_p*utanRiem_xyyz
												-gzx_p*utanRiem_xzxy-gxy_p*utanRiem_xyxy)/ewa4i;
						uRC_xz=(-gyy_p*utanRiem_xyyz-gzx_p*utanRiem_xzxz
												-gxy_p*utanRiem_xyxz-gyz_p*utanRiem_xzyz)/ewa4i;
						uRC_yz=(gxx_p*utanRiem_xyxz+gxy_p*utanRiem_xyyz
												-gzx_p*utanRiem_yzxz-gyz_p*utanRiem_yzyz)/ewa4i;
						
						//RC_ij
						
						RC_xx=(
									gxx_p*gxx_p*uRC_xx+gxy_p*gxy_p*uRC_yy+gxz_p*gxz_p*uRC_zz
									+gxx_p*gxy_p*uRC_xy+gxx_p*gxz_p*uRC_xz+gxy_p*gxz_p*uRC_yz
									+gxy_p*gxx_p*uRC_xy+gxz_p*gxx_p*uRC_xz+gxz_p*gxy_p*uRC_yz)/ewa8i;
						RC_yy=(
									gyx_p*gyx_p*uRC_xx+gyy_p*gyy_p*uRC_yy+gyz_p*gyz_p*uRC_zz
									+gyx_p*gyy_p*uRC_xy+gyx_p*gyz_p*uRC_xz+gyy_p*gyz_p*uRC_yz
									+gyy_p*gyx_p*uRC_xy+gyz_p*gyx_p*uRC_xz+gyz_p*gyy_p*uRC_yz)/ewa8i;
						RC_zz=(
									gzx_p*gzx_p*uRC_xx+gzy_p*gzy_p*uRC_yy+gzz_p*gzz_p*uRC_zz
									+gzx_p*gzy_p*uRC_xy+gzx_p*gzz_p*uRC_xz+gzy_p*gzz_p*uRC_yz
									+gzy_p*gzx_p*uRC_xy+gzz_p*gzx_p*uRC_xz+gzz_p*gzy_p*uRC_yz)/ewa8i;
						RC_xy=(
									gxx_p*gyx_p*uRC_xx+gxy_p*gyy_p*uRC_yy+gxz_p*gyz_p*uRC_zz
									+gxx_p*gyy_p*uRC_xy+gxx_p*gyz_p*uRC_xz+gxy_p*gyz_p*uRC_yz
									+gxy_p*gyx_p*uRC_xy+gxz_p*gyx_p*uRC_xz+gxz_p*gyy_p*uRC_yz)/ewa8i;
						RC_xz=(
									gxx_p*gzx_p*uRC_xx+gxy_p*gzy_p*uRC_yy+gxz_p*gzz_p*uRC_zz
									+gxx_p*gzy_p*uRC_xy+gxx_p*gzz_p*uRC_xz+gxy_p*gzz_p*uRC_yz
									+gxy_p*gzx_p*uRC_xy+gxz_p*gzx_p*uRC_xz+gxz_p*gzy_p*uRC_yz)/ewa8i;
						RC_yz=(
									gyx_p*gzx_p*uRC_xx+gyy_p*gzy_p*uRC_yy+gyz_p*gzz_p*uRC_zz
									+gyx_p*gzy_p*uRC_xy+gyx_p*gzz_p*uRC_xz+gyy_p*gzz_p*uRC_yz
									+gyy_p*gzx_p*uRC_xy+gyz_p*gzx_p*uRC_xz+gyz_p*gzy_p*uRC_yz)/ewa8i;
						//RC_yx=RC_xy;
						//RC_zx=RC_xz;
						//RC_zy=RC_yz;

						//RC
						
						RC=(gxx_p*unnRiem_xx+gyy_p*unnRiem_yy+gzz_p*unnRiem_zz
										+2.0*(gxy_p*unnRiem_xy+gxz_p*unnRiem_xz+gyz_p*unnRiem_yz))/ewa4i;
						
						//RC_i
						
						//RC_1=exp(-4.*wa_p)*(gixx_p*nRiem_x1_x+giyy_p*nRiem_y1_y+gizz_p*nRiem_z1_z
						//					+gixy_p*nRiem_x1_y+gixz_p*nRiem_x1_z+giyz_p*nRiem_y1_z
						//					+giyx_p*nRiem_y1_x+gizx_p*nRiem_z1_x+gizy_p*nRiem_z1_y);
						RC_x=ewa4i*(-giyy_p*nRiem_xy_y-gizz_p*nRiem_xz_z
											-giyz_p*nRiem_xy_z-giyx_p*nRiem_xy_x
											-gizx_p*nRiem_xz_x-gizy_p*nRiem_xz_y);
						RC_y=ewa4i*(gixx_p*nRiem_xy_x-gizz_p*nRiem_yz_z
											+gixy_p*nRiem_xy_y+gixz_p*nRiem_xy_z
											-gizx_p*nRiem_yz_x-gizy_p*nRiem_yz_y);
						RC_z=ewa4i*(gixx_p*nRiem_xz_x+giyy_p*nRiem_yz_y
											+gixy_p*nRiem_xz_y+gixz_p*nRiem_xz_z
											+giyz_p*nRiem_yz_z+giyx_p*nRiem_yz_x);
						
						RijRij=(unnRiem_xx-uRC_xx)*(nnRiem_xx-RC_xx)
										+(unnRiem_yy-uRC_yy)*(nnRiem_yy-RC_yy)
										+(unnRiem_zz-uRC_zz)*(nnRiem_zz-RC_zz)
										+2.*(unnRiem_xy-uRC_xy)*(nnRiem_xy-RC_xy)
										+2.*(unnRiem_xz-uRC_xz)*(nnRiem_xz-RC_xz)
										+2.*(unnRiem_yz-uRC_yz)*(nnRiem_yz-RC_yz)
										-2.*ewa4i
											*(gixx_p*RC_x*RC_x+giyy_p*RC_y*RC_y+gizz_p*RC_z*RC_z
												+2.*(gixy_p*RC_x*RC_y+gixz_p*RC_x*RC_z+giyz_p*RC_y*RC_z))
										+RC*RC;
										
						RR=-2.*RC
									+(gxx_p*uRC_xx+gyy_p*uRC_yy+gzz_p*uRC_zz
										+2.0*(gxy_p*uRC_xy+gxz_p*uRC_xz+gyz_p*uRC_yz))/ewa4i;
						
						set_outv(l,k,j,1)=RijRij;
						set_outv(l,k,j,2)=RR*RR;
						set_outv(l,k,j,3)=Kinv-2.*RijRij+RR/3.;
					}
				}
			}
		}
	}
	
	//Kreiss-Oliger dissipation
	if(abs(KOep)>1.0e-10)
	KOdiss();
	
	//excision
	if(exc) 
	excision();
	
	//  4th order Runge-Kutta method
	//  add flux
	runge_kutta(dt0*frac);		//bvr+=dbv*dt
	
	#pragma omp barrier

	if(itype!=4)
	{
		new_bv(dt);				//bv=bv0+dbv*dt
	}
	else
	{
		new_bv4();				//bv=bv0+bvr
	}
	#pragma omp barrier

	//determinant of gamma=1 and tr A_ij=0
	enforce_const();
	return;
}

//guarantee the determinant of gamma=1 and tracelessness of A_ij
//  -> psi(wa) and trK(ek) also change
void Fmv0::enforce_const()
{
	#pragma omp parallel for 
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				if(get_bflag(l,k,j)==-1)
				continue;

				enforce_const_gp(l,k,j);
			}
		}
	}
	return;
}

//enforce const for each grid point
void Fmv0::enforce_const_gp(int l,int k,int j)
{
	double gxx_p,gyy_p,gzz_p,gxy_p,gxz_p,gyz_p,gyx_p,gzx_p,gzy_p;
	double akxx_p,akyy_p,akzz_p,akxy_p,akxz_p,akyz_p;
	double det_p,det_pi, det_p13;
	double gixx_p,gixy_p,gixz_p,giyy_p,giyz_p,gizz_p;
	double tmp_ek, tmp_ek3;
		
	//  calculate determinant of tilde{\gamma_{ij}}
	gxx_p=get_bv(l,k,j, 7)+get_flat_df2x(j);
	gyy_p=get_bv(l,k,j, 8)+get_flat_df2y(k);
	gzz_p=get_bv(l,k,j, 9)+get_flat_df2z(l);
	gxy_p=get_bv(l,k,j,10);
	gxz_p=get_bv(l,k,j,11);
	gyz_p=get_bv(l,k,j,12);
	gyx_p=gxy_p;
	gzx_p=gxz_p;
	gzy_p=gyz_p;
	det_p= gxx_p*gyy_p*gzz_p +gxy_p*gyz_p*gzx_p +gxz_p*gyx_p*gzy_p
			-gxz_p*gyy_p*gzx_p -gxy_p*gyx_p*gzz_p -gxx_p*gyz_p*gzy_p;

	if(det_p<1.e-16) 
		det_p=1.e-16;
	det_pi=1./det_p;
	
	//det_p13=pow(det_p,(-1./3.));

	gixx_p= (gyy_p*gzz_p-gyz_p*gzy_p)*det_pi;
	giyy_p= (gxx_p*gzz_p-gxz_p*gzx_p)*det_pi;
	gizz_p= (gxx_p*gyy_p-gxy_p*gyx_p)*det_pi;
	gixy_p=-(gyx_p*gzz_p-gyz_p*gzx_p)*det_pi;
	gixz_p= (gyx_p*gzy_p-gyy_p*gzx_p)*det_pi;
	giyz_p=-(gxx_p*gzy_p-gxy_p*gzx_p)*det_pi;
	
	double fac=det_p/(get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l));
	det_p13=pow(fac,(-1./3.));
	
	//set_bv(l,k,j,13)=get_bv(l,k,j,13)+log(det_p)/12.;
	set_bv(l,k,j,13)=get_bv(l,k,j,13)+log(fac)/12.;
	
	akxx_p=get_bv(l,k,j,14);
	akyy_p=get_bv(l,k,j,15);
	akzz_p=get_bv(l,k,j,16);
	akxy_p=get_bv(l,k,j,17);
	akxz_p=get_bv(l,k,j,18);
	akyz_p=get_bv(l,k,j,19);

	tmp_ek=gixx_p*akxx_p
			+giyy_p*akyy_p
			+gizz_p*akzz_p
			+2.*( gixy_p*akxy_p
			+gixz_p*akxz_p
			+giyz_p*akyz_p );
	set_bv(l,k,j,20)=get_bv(l,k,j,20) +tmp_ek;
	
	set_bv(l,k,j, 7) = gxx_p*det_p13-get_flat_df2x(j);
	set_bv(l,k,j, 8) = gyy_p*det_p13-get_flat_df2y(k);
	set_bv(l,k,j, 9) = gzz_p*det_p13-get_flat_df2z(l);
	set_bv(l,k,j,10) = gxy_p*det_p13;
	set_bv(l,k,j,11) = gxz_p*det_p13;
	set_bv(l,k,j,12) = gyz_p*det_p13;

	tmp_ek3=tmp_ek/3.;
	set_bv(l,k,j,14)=akxx_p*det_p13 -tmp_ek3*gxx_p*det_p13;
	set_bv(l,k,j,15)=akyy_p*det_p13 -tmp_ek3*gyy_p*det_p13;
	set_bv(l,k,j,16)=akzz_p*det_p13 -tmp_ek3*gzz_p*det_p13;
	set_bv(l,k,j,17)=akxy_p*det_p13 -tmp_ek3*gxy_p*det_p13;
	set_bv(l,k,j,18)=akxz_p*det_p13 -tmp_ek3*gxz_p*det_p13;
	set_bv(l,k,j,19)=akyz_p*det_p13 -tmp_ek3*gyz_p*det_p13;
}

//preparation associated with \beta^i
void Fmv0::BSSN_adv()
{
	//int s=0,e=nn;
	//int s=7,e=nn;
	//nn=24
	int s=0,e=24;

	#pragma omp parallel for 
	for(int j=jli;j<=jui;j++)
	{
		double Gam_ux_xx=get_flat_Gamx(j);
		double fxx=get_flat_df2x(j);
		double dfxx=2.*fxx*Gam_ux_xx;
		for(int l=lli;l<=lui;l++)
		{
			double fzz=get_flat_df2z(l);
			double Gam_uz_zz=get_flat_Gamz(l);
			double dfzz=2.*fzz*Gam_uz_zz;
			for(int k=kli;k<=kui;k++)
			{
				if(get_bflag(l,k,j)!=0)
				 continue;
				double Gam_uy_yy=get_flat_Gamy(k);
				double fyy=get_flat_df2y(k);
				double dfyy=2.*fyy*Gam_uy_yy;

				double bx_p,by_p,bz_p;
			
				// beta^i
				bx_p=get_bv(l,k,j,1);
				by_p=get_bv(l,k,j,2);
				bz_p=get_bv(l,k,j,3);

				for(int i=s;i<e;i++)
				{
					double adv_x,adv_y,adv_z;

					if(bx_p<0){
						adv_x=( -get_bv(l,k,j-3,i)
								+ 6.*get_bv(l,k,j-2,i)
								-18.*get_bv(l,k,j-1,i)
								+10.*get_bv(l,k,j  ,i)
								+ 3.*get_bv(l,k,j+1,i)
								)*dxi12*bx_p;
					}else{
						adv_x=-( -get_bv(l,k,j+3,i)
								+ 6.*get_bv(l,k,j+2,i)
								-18.*get_bv(l,k,j+1,i)
								+10.*get_bv(l,k,j  ,i)
								+ 3.*get_bv(l,k,j-1,i)
								)*dxi12*bx_p;
					}


					if(by_p<0){
						adv_y=( -get_bv(l,k-3,j,i)
								+ 6.*get_bv(l,k-2,j,i)
								-18.*get_bv(l,k-1,j,i)
								+10.*get_bv(l,k  ,j,i)
								+ 3.*get_bv(l,k+1,j,i)
								)*dyi12*by_p;
					}else{
						adv_y=-( -get_bv(l,k+3,j,i)
								+ 6.*get_bv(l,k+2,j,i)
								-18.*get_bv(l,k+1,j,i)
								+10.*get_bv(l,k  ,j,i)
								+ 3.*get_bv(l,k-1,j,i)
								)*dyi12*by_p;
					}


					if(bz_p<0){
						adv_z=( -get_bv(l-3,k,j,i)
								+ 6.*get_bv(l-2,k,j,i)
								-18.*get_bv(l-1,k,j,i)
								+10.*get_bv(l  ,k,j,i)
								+ 3.*get_bv(l+1,k,j,i)
								)*dzi12*bz_p;
					}else{
						adv_z=-( -get_bv(l+3,k,j,i)
								+ 6.*get_bv(l+2,k,j,i)
								-18.*get_bv(l+1,k,j,i)
								+10.*get_bv(l  ,k,j,i)
								+ 3.*get_bv(l-1,k,j,i)
								)*dzi12*bz_p;
					}
					//set_bv(l,k,j,i)= get_bv(l,k,j,i) +adv_x+adv_y+adv_z;
					set_dbv(l,k,j,i)= adv_x+adv_y+adv_z;
				}
				
				set_dbv(l,k,j,7)+=bx_p*dfxx;
				set_dbv(l,k,j,8)+=by_p*dfyy;
				set_dbv(l,k,j,9)+=bz_p*dfzz;
				
			}
		}
	}// for-loop end
	
	
	return;
}

//Kreiss-Oliger dissipation term
void Fmv0::KOdiss()
{
	//int s=0,e=nn;
	//int s=0,e=nn;
	int s=0,e=23;

	#pragma omp parallel for 
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				if(get_bflag(l,k,j)!=0)
				 continue;
				
				double diss_x,diss_y,diss_z;
				
				for(int i=s;i<=e;i++)
				{
					diss_x=( get_bv(l,k,j-3,i)
							 	-6.*get_bv(l,k,j-2,i)
									+ 15.*get_bv(l,k,j-1,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l,k,j+1,i)
								-6.*get_bv(l,k,j+2,i)
							+get_bv(l,k,j+3,i)
							)*dxi;
					diss_y=( get_bv(l,k-3,j,i)
							 	-6.*get_bv(l,k-2,j,i)
									+ 15.*get_bv(l,k-1,j,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l,k+1,j,i)
								-6.*get_bv(l,k+2,j,i)
							+get_bv(l,k+3,j,i)
							)*dyi;
					diss_z=( get_bv(l-3,k,j,i)
							 	-6.*get_bv(l-2,k,j,i)
									+ 15.*get_bv(l-1,k,j,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l+1,k,j,i)
								-6.*get_bv(l+2,k,j,i)
							+get_bv(l+3,k,j,i)
							)*dzi;
				        
					set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + KOep*(diss_x+diss_y+diss_z)/64.;
				}
				
				//for scalar field
				if(scalarevo)
				{
					int i=nsc;
					diss_x=( get_bv(l,k,j-3,i)
							 	-6.*get_bv(l,k,j-2,i)
									+ 15.*get_bv(l,k,j-1,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l,k,j+1,i)
								-6.*get_bv(l,k,j+2,i)
							+get_bv(l,k,j+3,i)
							)*dxi;
					diss_y=( get_bv(l,k-3,j,i)
							 	-6.*get_bv(l,k-2,j,i)
									+ 15.*get_bv(l,k-1,j,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l,k+1,j,i)
								-6.*get_bv(l,k+2,j,i)
							+get_bv(l,k+3,j,i)
							)*dyi;
					diss_z=( get_bv(l-3,k,j,i)
							 	-6.*get_bv(l-2,k,j,i)
									+ 15.*get_bv(l-1,k,j,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l+1,k,j,i)
								-6.*get_bv(l+2,k,j,i)
							+get_bv(l+3,k,j,i)
							)*dzi;
				        
					set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + KOep*(diss_x+diss_y+diss_z)/64.;
					
					
					i=nscp;
					diss_x=( get_bv(l,k,j-3,i)
							 	-6.*get_bv(l,k,j-2,i)
									+ 15.*get_bv(l,k,j-1,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l,k,j+1,i)
								-6.*get_bv(l,k,j+2,i)
							+get_bv(l,k,j+3,i)
							)*dxi;
					diss_y=( get_bv(l,k-3,j,i)
							 	-6.*get_bv(l,k-2,j,i)
									+ 15.*get_bv(l,k-1,j,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l,k+1,j,i)
								-6.*get_bv(l,k+2,j,i)
							+get_bv(l,k+3,j,i)
							)*dyi;
					diss_z=( get_bv(l-3,k,j,i)
							 	-6.*get_bv(l-2,k,j,i)
									+ 15.*get_bv(l-1,k,j,i)
										-20.*get_bv(l,k,j,i)
									+ 15.*get_bv(l+1,k,j,i)
								-6.*get_bv(l+2,k,j,i)
							+get_bv(l+3,k,j,i)
							)*dzi;
				        
					set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + KOep*(diss_x+diss_y+diss_z)/64.;
				
				}
			}
		}
	}// for-loop end

	return;
}

//setting dbv for boundary grids of the excised region
void Fmv0::excision()
{
	//int s=0,e=nn;
	int s=0,e=nn;
	//nn=24
	
	boundary_d_quarter();
	
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				if(get_bflag(l,k,j)!=3)
				continue;
				
				int jj=j;
				int kk=k;
				int ll=l;
				
				if(get_bflag(l,k,j+1)==2)
				jj+=1;
				if(get_bflag(l,k,j-1)==2)
				jj-=1;
				
				if(get_bflag(l,k+1,j)==2)
				kk+=1;
				if(get_bflag(l,k-1,j)==2)
				kk-=1;
				
				if(get_bflag(l+1,k,j)==2)
				ll+=1;
				if(get_bflag(l-1,k,j)==2)
				ll-=1;
		
				for(int i=s;i<e;i++)
				 set_dbv(l,k,j,i)=get_dbv(ll,kk,jj,i);
			}
		}
	}// for-loop end
	
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				if(get_bflag(l,k,j)!=2)
				continue;
				
				int jj=j;
				int kk=k;
				int ll=l;
				
				if(get_bflag(l,k,j+1)==1)
				jj+=1;
				if(get_bflag(l,k,j-1)==1)
				jj-=1;
				
				if(get_bflag(l,k+1,j)==1)
				kk+=1;
				if(get_bflag(l,k-1,j)==1)
				kk-=1;
				
				if(get_bflag(l+1,k,j)==1)
				ll+=1;
				if(get_bflag(l-1,k,j)==1)
				ll-=1;
				
				for(int i=s;i<e;i++)
				 set_dbv(l,k,j,i)=get_dbv(ll,kk,jj,i);
			}
		}
	}// for-loop end
	
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				if(get_bflag(l,k,j)!=1)
				continue;
				
				int jj=j;
				int kk=k;
				int ll=l;
				
				if(get_bflag(l,k,j+1)==0)
				jj+=1;
				if(get_bflag(l,k,j-1)==0)
				jj-=1;
				
				if(get_bflag(l,k+1,j)==0)
				kk+=1;
				if(get_bflag(l,k-1,j)==0)
				kk-=1;
				
				if(get_bflag(l+1,k,j)==0)
				ll+=1;
				if(get_bflag(l-1,k,j)==0)
				ll-=1;
				
				for(int i=s;i<e;i++)
				 set_dbv(l,k,j,i)=get_dbv(ll,kk,jj,i);
			}
		}
	}// for-loop end

	boundary_d_quarter();
	
	return;
}

//flux setting for fluid
void Fmv0::flux_fill()
{
	#pragma omp parallel for 
	for(int j=jli;j<=jui+1;j++)
	{
		double xx=get_x(j)-0.5*dx;
		//double yy,zz;
		
		//metric variables
		double alp_ph,bx_ph,by_ph,bz_ph,wa_ph,gxx_ph,gyy_ph,gzz_ph,gxy_ph,gxz_ph,gyz_ph;
		double det,deti,sqdet;
		double gixx_ph;
		//double giyy_ph,gizz_ph;
		double fxx,fyy,fzz;
		
		//deviations of variable
		double Del2,Del1,Del0,Delpbar,Delmbar;
		
		//right and left dynamical variables 
		double S0L,S0R,SxL,SxR,SyL,SyR,SzL,SzR,DenL,DenR,rhoL,rhoR;
		double VxL,VxR,VyL,VyR,VzL,VzR,UxL,UxR,UyL,UyR,UzL,UzR,U2L,U2R;
		double cs2L,cs2R,ovfacL,ovfacR,lam1tL,lam2tL,lam1tR,lam2tR,lampL,lampR,lammL,lammR,as,prsL,prsR;

		for(int l=lli;l<=lui;l++)
		{
			//zz=get_z(l);
			for(int k=kli;k<=kui;k++)
			{
				//yy=get_y(k);
				
				//substitution for geometrical variables start
				alp_ph=get_ipol_x_lower_mid(l,k,j,0);
				
				bx_ph=get_ipol_x_lower_mid(l,k,j,1);
				by_ph=get_ipol_x_lower_mid(l,k,j,2);
				bz_ph=get_ipol_x_lower_mid(l,k,j,3);
				
				wa_ph=get_ipol_x_lower_mid(l,k,j,13);

				fxx=pow(df(xx),2);
				fyy=get_flat_df2y(k);
				fzz=get_flat_df2z(l);
				
				gxx_ph=get_ipol_x_lower_mid(l,k,j,7)+fxx;
				gyy_ph=get_ipol_x_lower_mid(l,k,j,8)+fyy;
				gzz_ph=get_ipol_x_lower_mid(l,k,j,9)+fzz;
				
				gxy_ph=get_ipol_x_lower_mid(l,k,j,10);
				gxz_ph=get_ipol_x_lower_mid(l,k,j,11);
				gyz_ph=get_ipol_x_lower_mid(l,k,j,12);
				
			//	det=gxx_ph*gyy_ph*gzz_ph+gxy_ph*gyz_ph*gxz_ph+gxz_ph*gxy_ph*gyz_ph
			//			-gxz_ph*gyy_ph*gxz_ph-gxy_ph*gxy_ph*gzz_ph-gxx_ph*gyz_ph*gyz_ph;
				det=fxx*fyy*fzz;

				deti=1./det;
				gixx_ph= (gyy_ph*gzz_ph-gyz_ph*gyz_ph)*deti*exp(-4.*wa_ph);
				
				sqdet= sqrt(det)*exp(6.*wa_ph);
				//substitution for geometrical variables end
				
				//calculation of left and right S0 start
				Del2=get_bv(l,k,j+1,24)-get_bv(l,k,j,24);
				Del1=get_bv(l,k,j,24)-get_bv(l,k,j-1,24);
				Del0=get_bv(l,k,j-1,24)-get_bv(l,k,j-2,24);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				S0L=get_bv(l,k,j-1,24)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				S0R=get_bv(l,k,j,24)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right S0 end
				
				//calculation of left and right Sx start
				Del2=get_bv(l,k,j+1,25)-get_bv(l,k,j,25);
				Del1=get_bv(l,k,j,25)-get_bv(l,k,j-1,25);
				Del0=get_bv(l,k,j-1,25)-get_bv(l,k,j-2,25);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SxL=get_bv(l,k,j-1,25)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SxR=get_bv(l,k,j,25)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sx end

				//calculation of left and right Sy start
				Del2=get_bv(l,k,j+1,26)-get_bv(l,k,j,26);
				Del1=get_bv(l,k,j,26)-get_bv(l,k,j-1,26);
				Del0=get_bv(l,k,j-1,26)-get_bv(l,k,j-2,26);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SyL=get_bv(l,k,j-1,26)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SyR=get_bv(l,k,j,26)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sy end
				
				//calculation of left and right Sz start
				Del2=get_bv(l,k,j+1,27)-get_bv(l,k,j,27);
				Del1=get_bv(l,k,j,27)-get_bv(l,k,j-1,27);
				Del0=get_bv(l,k,j-1,27)-get_bv(l,k,j-2,27);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SzL=get_bv(l,k,j-1,27)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SzR=get_bv(l,k,j,27)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sz end

				//calculation of left and right Density start (Den:Density=\rho_* in my note)
				Del2=get_bv(l,k,j+1,28)-get_bv(l,k,j,28);
				Del1=get_bv(l,k,j,28)-get_bv(l,k,j-1,28);
				Del0=get_bv(l,k,j-1,28)-get_bv(l,k,j-2,28);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				DenL=get_bv(l,k,j-1,28)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				DenR=get_bv(l,k,j,28)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Density end
				
				//calculation of left and right rho(primitive) start
				Del2=get_primv(l,k,j+1,0)-get_primv(l,k,j,0);
				Del1=get_primv(l,k,j,0)-get_primv(l,k,j-1,0);
				Del0=get_primv(l,k,j-1,0)-get_primv(l,k,j-2,0);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				rhoL=get_primv(l,k,j-1,0)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				rhoR=get_primv(l,k,j,0)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right rho(primitive) end
				
				
				//calculation of left and right Vx(primitive) start
				Del2=get_primv(l,k,j+1,1)-get_primv(l,k,j,1);
				Del1=get_primv(l,k,j,1)-get_primv(l,k,j-1,1);
				Del0=get_primv(l,k,j-1,1)-get_primv(l,k,j-2,1);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VxL=get_primv(l,k,j-1,1)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VxR=get_primv(l,k,j,1)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vx(primitive) end
				
				
				//calculation of left and right Vy(primitive) start
				Del2=get_primv(l,k,j+1,2)-get_primv(l,k,j,2);
				Del1=get_primv(l,k,j,2)-get_primv(l,k,j-1,2);
				Del0=get_primv(l,k,j-1,2)-get_primv(l,k,j-2,2);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VyL=get_primv(l,k,j-1,2)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VyR=get_primv(l,k,j,2)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vx(primitive) end
				
				
				//calculation of left and right Vz(primitive) start
				Del2=get_primv(l,k,j+1,3)-get_primv(l,k,j,3);
				Del1=get_primv(l,k,j,3)-get_primv(l,k,j-1,3);
				Del0=get_primv(l,k,j-1,3)-get_primv(l,k,j-2,3);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VzL=get_primv(l,k,j-1,3)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VzR=get_primv(l,k,j,3)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vz(primitive) end

				//calculation of left and right epsilon(primitive) start
				//Del2=get_primv(l,k,j+2,4)-get_primv(l,k,j+1,4);
				//Del1=get_primv(l,k,j+1,4)-get_primv(l,k,j,4);
				//Del0=get_primv(l,k,j,4)-get_primv(l,k,j-1,4);
				
				//Delpbar=minmod(Del1,b_minmod*Del0);
				//Delmbar=minmod(Del0,b_minmod*Del1);
				
				//epsLx=get_primv(l,k,j,4)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				//Delpbar=minmod(Del2,b_minmod*Del1);
				//Delmbar=minmod(Del1,b_minmod*Del2);
				
				//epsRx=get_primv(l,k,j+1,4)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right epsilon(primitive) end
				
				//calculation of left and right U_i start
				UxL=(VxL+bx_ph)/alp_ph;
				UxR=(VxR+bx_ph)/alp_ph;
				UyL=(VyL+by_ph)/alp_ph;
				UyR=(VyR+by_ph)/alp_ph;
				UzL=(VzL+bz_ph)/alp_ph;
				UzR=(VzR+bz_ph)/alp_ph;
				
				U2L=(UxL*UxL*gxx_ph+UyL*UyL*gyy_ph+UzL*UzL*gzz_ph
					+2.*(UxL*UyL*gxy_ph+UxL*UzL*gxz_ph+UyL*UzL*gyz_ph))*exp(4.*wa_ph);
				U2R=(UxR*UxR*gxx_ph+UyR*UyR*gyy_ph+UzR*UzR*gzz_ph
					+2.*(UxR*UyR*gxy_ph+UxR*UzR*gxz_ph+UyR*UzR*gyz_ph))*exp(4.*wa_ph);
				
				//the speed should not exceed 1
				U2L=min(U2L,1.-1.0e-10);
				U2R=min(U2R,1.-1.0e-10);
				//calculation of left and right U_i end
				
				//sound velocity
				cs2L=dpres(rhoL);
				cs2R=dpres(rhoR);
				
				//calculation of characteristic speeds start
				ovfacL=alp_ph/(1.-U2L*cs2L);
				ovfacR=alp_ph/(1.-U2R*cs2R);
				
				lam1tL=UxL*(1-cs2L);
				lam2tL=sqrt(max(cs2L*(1-U2L)*(gixx_ph*(1-U2L*cs2L)-(1-cs2L)*UxL*UxL),1.0e-16));

				lam1tR=UxR*(1-cs2R);
				lam2tR=sqrt(max(cs2R*(1-U2R)*(gixx_ph*(1-U2R*cs2R)-(1-cs2R)*UxR*UxR),1.0e-16));
				
				lampL=ovfacL*(lam1tL+lam2tL)-bx_ph;
				lampR=ovfacR*(lam1tR+lam2tR)-bx_ph;
				
				lammL=ovfacL*(lam1tL-lam2tL)-bx_ph;
				lammR=ovfacR*(lam1tR-lam2tR)-bx_ph;
				
				as=max({abs(VxL),abs(VxR),abs(lampL),abs(lampR),abs(lammL),abs(lammR)});
				//calculation of characteristic speeds end
				
				//pressure
				prsL=pres(rhoL);
				prsR=pres(rhoR);
				
				//fluxes
				set_flux_x(l,k,j,0)=0.5*(S0L*VxL+S0R*VxR+sqdet*(prsL*(VxL+bx_ph)+prsR*(VxR+bx_ph))-as*(S0R-S0L));
				set_flux_x(l,k,j,1)=0.5*(SxL*VxL+SxR*VxR+sqdet*alp_ph*(prsL+prsR)-as*(SxR-SxL));
				set_flux_x(l,k,j,2)=0.5*(SyL*VxL+SyR*VxR-as*(SyR-SyL));
				set_flux_x(l,k,j,3)=0.5*(SzL*VxL+SzR*VxR-as*(SzR-SzL));
				set_flux_x(l,k,j,4)=0.5*(DenL*VxL+DenR*VxR-as*(DenR-DenL));
				
			}
		}
	}

	#pragma omp parallel for 
	for(int i=kli;i<=kui+1;i++)
	{
		int k=i;
		int j,l;
		double yy=get_y(k)-0.5*dy;
		double zz;
		
		//metric variables
		double alp_ph,bx_ph,by_ph,bz_ph,wa_ph,gxx_ph,gyy_ph,gzz_ph,gxy_ph,gxz_ph,gyz_ph;
		double det,deti,sqdet;
		//double gixx_ph;
		double giyy_ph,gizz_ph;
		double fxx,fyy,fzz;
		
		//deviations of variable
		double Del2,Del1,Del0,Delpbar,Delmbar;
		
		//right and left dynamical variables 
		double S0L,S0R,SxL,SxR,SyL,SyR,SzL,SzR,DenL,DenR,rhoL,rhoR;
		double VxL,VxR,VyL,VyR,VzL,VzR,UxL,UxR,UyL,UyR,UzL,UzR,U2L,U2R;
		double cs2L,cs2R,ovfacL,ovfacR,lam1tL,lam2tL,lam1tR,lam2tR,lampL,lampR,lammL,lammR,as,prsL,prsR;

		for(j=jli;j<=jui;j++)
		{
			//xx=get_x(j);
		
			for(l=lli;l<=lui;l++)
			{
				//zz=get_z(l);

				//substitution for geometrical variables start
				alp_ph=get_ipol_y_lower_mid(l,k,j,0);
				
				bx_ph=get_ipol_y_lower_mid(l,k,j,1);
				by_ph=get_ipol_y_lower_mid(l,k,j,2);
				bz_ph=get_ipol_y_lower_mid(l,k,j,3);
				
				wa_ph=get_ipol_y_lower_mid(l,k,j,13);
				
				fxx=get_flat_df2x(j);
				fyy=pow(df(yy),2);
				fzz=get_flat_df2z(l);
				
				gxx_ph=get_ipol_y_lower_mid(l,k,j,7)+fxx;
				gyy_ph=get_ipol_y_lower_mid(l,k,j,8)+fyy;
				gzz_ph=get_ipol_y_lower_mid(l,k,j,9)+fzz;
				gxy_ph=get_ipol_y_lower_mid(l,k,j,10);
				gxz_ph=get_ipol_y_lower_mid(l,k,j,11);
				gyz_ph=get_ipol_y_lower_mid(l,k,j,12);
				
			//	det=gxx_ph*gyy_ph*gzz_ph+gxy_ph*gyz_ph*gxz_ph+gxz_ph*gxy_ph*gyz_ph
			//			-gxz_ph*gyy_ph*gxz_ph-gxy_ph*gxy_ph*gzz_ph-gxx_ph*gyz_ph*gyz_ph;
				det=fxx*fyy*fzz;
				
				if(det<1.e-16) 
				 det=1.e-16;
				deti=1./det;
				giyy_ph= (gxx_ph*gzz_ph-gxz_ph*gxz_ph)*deti*exp(-4.*wa_ph);
				
				sqdet= sqrt(det)*exp(6.*wa_ph);
				//substitution for geometrical variables end
				
				//calculation of left and right S0 start
				Del2=get_bv(l,k+1,j,24)-get_bv(l,k,j,24);
				Del1=get_bv(l,k,j,24)-get_bv(l,k-1,j,24);
				Del0=get_bv(l,k-1,j,24)-get_bv(l,k-2,j,24);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				S0L=get_bv(l,k-1,j,24)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				S0R=get_bv(l,k,j,24)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right S0 end
				
				//calculation of left and right Sx start
				Del2=get_bv(l,k+1,j,25)-get_bv(l,k,j,25);
				Del1=get_bv(l,k,j,25)-get_bv(l,k-1,j,25);
				Del0=get_bv(l,k-1,j,25)-get_bv(l,k-2,j,25);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SxL=get_bv(l,k-1,j,25)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SxR=get_bv(l,k,j,25)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sx end
				
				//calculation of left and right Sy start
				Del2=get_bv(l,k+1,j,26)-get_bv(l,k,j,26);
				Del1=get_bv(l,k,j,26)-get_bv(l,k-1,j,26);
				Del0=get_bv(l,k-1,j,26)-get_bv(l,k-2,j,26);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SyL=get_bv(l,k-1,j,26)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SyR=get_bv(l,k,j,26)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sy end
				
				//calculation of left and right Sz start
				Del2=get_bv(l,k+1,j,27)-get_bv(l,k,j,27);
				Del1=get_bv(l,k,j,27)-get_bv(l,k-1,j,27);
				Del0=get_bv(l,k-1,j,27)-get_bv(l,k-2,j,27);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SzL=get_bv(l,k-1,j,27)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SzR=get_bv(l,k,j,27)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sz end
				
				//calculation of left and right Density start (Den:Density=\rho_* in my note)
				Del2=get_bv(l,k+1,j,28)-get_bv(l,k,j,28);
				Del1=get_bv(l,k,j,28)-get_bv(l,k-1,j,28);
				Del0=get_bv(l,k-1,j,28)-get_bv(l,k-2,j,28);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				DenL=get_bv(l,k-1,j,28)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				DenR=get_bv(l,k,j,28)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Density end
				
				//calculation of left and right rho(primitive) start
				Del2=get_primv(l,k+1,j,0)-get_primv(l,k,j,0);
				Del1=get_primv(l,k,j,0)-get_primv(l,k-1,j,0);
				Del0=get_primv(l,k-1,j,0)-get_primv(l,k-2,j,0);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				rhoL=get_primv(l,k-1,j,0)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				rhoR=get_primv(l,k,j,0)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right rho(primitive) end
				
				
				//calculation of left and right Vx(primitive) start
				Del2=get_primv(l,k+1,j,1)-get_primv(l,k,j,1);
				Del1=get_primv(l,k,j,1)-get_primv(l,k-1,j,1);
				Del0=get_primv(l,k-1,j,1)-get_primv(l,k-2,j,1);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VxL=get_primv(l,k-1,j,1)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VxR=get_primv(l,k,j,1)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vx(primitive) end
				
				
				//calculation of left and right Vy(primitive) start
				Del2=get_primv(l,k+1,j,2)-get_primv(l,k,j,2);
				Del1=get_primv(l,k,j,2)-get_primv(l,k-1,j,2);
				Del0=get_primv(l,k-1,j,2)-get_primv(l,k-2,j,2);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VyL=get_primv(l,k-1,j,2)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VyR=get_primv(l,k,j,2)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vy(primitive) end
				
				
				//calculation of left and right Vz(primitive) start
				Del2=get_primv(l,k+1,j,3)-get_primv(l,k,j,3);
				Del1=get_primv(l,k,j,3)-get_primv(l,k-1,j,3);
				Del0=get_primv(l,k-1,j,3)-get_primv(l,k-2,j,3);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VzL=get_primv(l,k-1,j,3)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VzR=get_primv(l,k,j,3)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vz(primitive) end

				//calculation of left and right epsilon(primitive) start
				//Del2=get_primv(l,k,j+2,4)-get_primv(l,k,j+1,4);
				//Del1=get_primv(l,k,j+1,4)-get_primv(l,k,j,4);
				//Del0=get_primv(l,k,j,4)-get_primv(l,k,j-1,4);
				
				//Delpbar=minmod(Del1,b_minmod*Del0);
				//Delmbar=minmod(Del0,b_minmod*Del1);
				
				//epsLx=get_primv(l,k,j,4)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				//Delpbar=minmod(Del2,b_minmod*Del1);
				//Delmbar=minmod(Del1,b_minmod*Del2);
				
				//epsRx=get_primv(l,k,j+1,4)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right epsilon(primitive) end
				
				//calculation of left and right U_i start
				UxL=(VxL+bx_ph)/alp_ph;
				UxR=(VxR+bx_ph)/alp_ph;
				UyL=(VyL+by_ph)/alp_ph;
				UyR=(VyR+by_ph)/alp_ph;
				UzL=(VzL+bz_ph)/alp_ph;
				UzR=(VzR+bz_ph)/alp_ph;
				
				U2L=(UxL*UxL*gxx_ph+UyL*UyL*gyy_ph+UzL*UzL*gzz_ph
					+2.*(UxL*UyL*gxy_ph+UxL*UzL*gxz_ph+UyL*UzL*gyz_ph))*exp(4.*wa_ph);
				U2R=(UxR*UxR*gxx_ph+UyR*UyR*gyy_ph+UzR*UzR*gzz_ph
					+2.*(UxR*UyR*gxy_ph+UxR*UzR*gxz_ph+UyR*UzR*gyz_ph))*exp(4.*wa_ph);
				
				//the speed should not exceed 1
				U2L=min(U2L,1.-1.0e-10);
				U2R=min(U2R,1.-1.0e-10);
				//calculation of left and right U_i end

				//the speed should not exceed 1
				cs2L=dpres(rhoL);
				cs2R=dpres(rhoR);
				
				//calculation of characteristic speeds start
				ovfacL=alp_ph/(1.-U2L*cs2L);
				ovfacR=alp_ph/(1.-U2R*cs2R);
				
				lam1tL=UyL*(1-cs2L);
				lam2tL=sqrt(max(cs2L*(1-U2L)*(giyy_ph*(1-U2L*cs2L)-(1-cs2L)*UyL*UyL),1.0e-16));

				lam1tR=UyR*(1-cs2R);
				lam2tR=sqrt(max(cs2R*(1-U2R)*(giyy_ph*(1-U2R*cs2R)-(1-cs2R)*UyR*UyR),1.0e-16));
				
				lampL=ovfacL*(lam1tL+lam2tL)-by_ph;
				lampR=ovfacR*(lam1tR+lam2tR)-by_ph;
				
				lammL=ovfacL*(lam1tL-lam2tL)-by_ph;
				lammR=ovfacR*(lam1tR-lam2tR)-by_ph;
				
				as=max({abs(VyL),abs(VyR),abs(lampL),abs(lampR),abs(lammL),abs(lammR)});
				//calculation of characteristic speeds end
				
				//pressure
				prsL=pres(rhoL);
				prsR=pres(rhoR);
				
				//fluxes
				set_flux_y(l,k,j,0)=0.5*(S0L*VyL+S0R*VyR+sqdet*(prsL*(VyL+by_ph)+prsR*(VyR+by_ph))-as*(S0R-S0L));
				set_flux_y(l,k,j,1)=0.5*(SxL*VyL+SxR*VyR-as*(SxR-SxL));
				set_flux_y(l,k,j,2)=0.5*(SyL*VyL+SyR*VyR+sqdet*alp_ph*(prsL+prsR)-as*(SyR-SyL));
				set_flux_y(l,k,j,3)=0.5*(SzL*VyL+SzR*VyR-as*(SzR-SzL));
				set_flux_y(l,k,j,4)=0.5*(DenL*VyL+DenR*VyR-as*(DenR-DenL));
			}
		}

		l=i;
		zz=get_z(l)-0.5*dz;

		for(k=kli;k<=kui;k++)
		{
			//yy=get_y(k);
			for(j=jli;j<=jui;j++)
			{
				//xx=get_x(j);
			
				//substitution for geometrical variables start
				alp_ph=get_ipol_z_lower_mid(l,k,j,0);
				
				bx_ph=get_ipol_z_lower_mid(l,k,j,1);
				by_ph=get_ipol_z_lower_mid(l,k,j,2);
				bz_ph=get_ipol_z_lower_mid(l,k,j,3);
				
				wa_ph=get_ipol_z_lower_mid(l,k,j,13);
				
				fxx=get_flat_df2x(j);
				fyy=get_flat_df2y(k);
				fzz=pow(df(zz),2);

				gxx_ph=get_ipol_z_lower_mid(l,k,j,7)+fxx;
				gyy_ph=get_ipol_z_lower_mid(l,k,j,8)+fyy;
				gzz_ph=get_ipol_z_lower_mid(l,k,j,9)+fzz;
				gxy_ph=get_ipol_z_lower_mid(l,k,j,10);
				gxz_ph=get_ipol_z_lower_mid(l,k,j,11);
				gyz_ph=get_ipol_z_lower_mid(l,k,j,12);
				
			//	det=gxx_ph*gyy_ph*gzz_ph+gxy_ph*gyz_ph*gxz_ph+gxz_ph*gxy_ph*gyz_ph
			//			-gxz_ph*gyy_ph*gxz_ph-gxy_ph*gxy_ph*gzz_ph-gxx_ph*gyz_ph*gyz_ph;
				det=fxx*fyy*fzz;
				
				if(det<1.e-16) 
				 det=1.e-16;
				deti=1./det;
				gizz_ph= (gxx_ph*gyy_ph-gxy_ph*gxy_ph)*deti*exp(-4.*wa_ph);
				
				sqdet= sqrt(det)*exp(6.*wa_ph);
				//substitution for geometrical variables end

				//calculation of left and right S0 start
				Del2=get_bv(l+1,k,j,24)-get_bv(l,k,j,24);
				Del1=get_bv(l,k,j,24)-get_bv(l-1,k,j,24);
				Del0=get_bv(l-1,k,j,24)-get_bv(l-2,k,j,24);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				S0L=get_bv(l-1,k,j,24)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				S0R=get_bv(l,k,j,24)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right S0 end
				
				
				//calculation of left and right Sx start
				Del2=get_bv(l+1,k,j,25)-get_bv(l,k,j,25);
				Del1=get_bv(l,k,j,25)-get_bv(l-1,k,j,25);
				Del0=get_bv(l-1,k,j,25)-get_bv(l-2,k,j,25);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SxL=get_bv(l-1,k,j,25)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SxR=get_bv(l,k,j,25)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sx end
				
				
				//calculation of left and right Sy start
				Del2=get_bv(l+1,k,j,26)-get_bv(l,k,j,26);
				Del1=get_bv(l,k,j,26)-get_bv(l-1,k,j,26);
				Del0=get_bv(l-1,k,j,26)-get_bv(l-2,k,j,26);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SyL=get_bv(l-1,k,j,26)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SyR=get_bv(l,k,j,26)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sy end
				
				
				//calculation of left and right Sz start
				Del2=get_bv(l+1,k,j,27)-get_bv(l,k,j,27);
				Del1=get_bv(l,k,j,27)-get_bv(l-1,k,j,27);
				Del0=get_bv(l-1,k,j,27)-get_bv(l-2,k,j,27);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				SzL=get_bv(l-1,k,j,27)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				SzR=get_bv(l,k,j,27)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Sz end
				
				
				//calculation of left and right Density start (Den:Density=\rho_* in my note)
				Del2=get_bv(l+1,k,j,28)-get_bv(l,k,j,28);
				Del1=get_bv(l,k,j,28)-get_bv(l-1,k,j,28);
				Del0=get_bv(l-1,k,j,28)-get_bv(l-2,k,j,28);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				DenL=get_bv(l-1,k,j,28)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				DenR=get_bv(l,k,j,28)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Density end
				
				
				//calculation of left and right rho(primitive) start
				Del2=get_primv(l+1,k,j,0)-get_primv(l,k,j,0);
				Del1=get_primv(l,k,j,0)-get_primv(l-1,k,j,0);
				Del0=get_primv(l-1,k,j,0)-get_primv(l-2,k,j,0);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				rhoL=get_primv(l-1,k,j,0)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				rhoR=get_primv(l,k,j,0)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right rho(primitive) end
				
				
				//calculation of left and right Vx(primitive) start
				Del2=get_primv(l+1,k,j,1)-get_primv(l,k,j,1);
				Del1=get_primv(l,k,j,1)-get_primv(l-1,k,j,1);
				Del0=get_primv(l-1,k,j,1)-get_primv(l-2,k,j,1);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VxL=get_primv(l-1,k,j,1)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VxR=get_primv(l,k,j,1)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vx(primitive) end
				
				
				//calculation of left and right Vy(primitive) start
				Del2=get_primv(l+1,k,j,2)-get_primv(l,k,j,2);
				Del1=get_primv(l,k,j,2)-get_primv(l-1,k,j,2);
				Del0=get_primv(l-1,k,j,2)-get_primv(l-2,k,j,2);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VyL=get_primv(l-1,k,j,2)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VyR=get_primv(l,k,j,2)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vy(primitive) end
				
				
				//calculation of left and right Vz(primitive) start
				Del2=get_primv(l+1,k,j,3)-get_primv(l,k,j,3);
				Del1=get_primv(l,k,j,3)-get_primv(l-1,k,j,3);
				Del0=get_primv(l-1,k,j,3)-get_primv(l-2,k,j,3);
				
				Delpbar=minmod(Del1,b_minmod*Del0);
				Delmbar=minmod(Del0,b_minmod*Del1);
				
				VzL=get_primv(l-1,k,j,3)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				Delpbar=minmod(Del2,b_minmod*Del1);
				Delmbar=minmod(Del1,b_minmod*Del2);
				
				VzR=get_primv(l,k,j,3)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right Vz(primitive) end

				//calculation of left and right epsilon(primitive) start
				//Del2=get_primv(l,k,j+2,4)-get_primv(l,k,j+1,4);
				//Del1=get_primv(l,k,j+1,4)-get_primv(l,k,j,4);
				//Del0=get_primv(l,k,j,4)-get_primv(l,k,j-1,4);
				
				//Delpbar=minmod(Del1,b_minmod*Del0);
				//Delmbar=minmod(Del0,b_minmod*Del1);
				
				//epsLx=get_primv(l,k,j,4)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
				
				//Delpbar=minmod(Del2,b_minmod*Del1);
				//Delmbar=minmod(Del1,b_minmod*Del2);
				
				//epsRx=get_primv(l,k,j+1,4)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
				//calculation of left and right epsilon(primitive) end
				
				//calculation of left and right U_i start
				UxL=(VxL+bx_ph)/alp_ph;
				UxR=(VxR+bx_ph)/alp_ph;
				UyL=(VyL+by_ph)/alp_ph;
				UyR=(VyR+by_ph)/alp_ph;
				UzL=(VzL+bz_ph)/alp_ph;
				UzR=(VzR+bz_ph)/alp_ph;
				
				U2L=(UxL*UxL*gxx_ph+UyL*UyL*gyy_ph+UzL*UzL*gzz_ph
					+2.*(UxL*UyL*gxy_ph+UxL*UzL*gxz_ph+UyL*UzL*gyz_ph))*exp(4.*wa_ph);
				U2R=(UxR*UxR*gxx_ph+UyR*UyR*gyy_ph+UzR*UzR*gzz_ph
					+2.*(UxR*UyR*gxy_ph+UxR*UzR*gxz_ph+UyR*UzR*gyz_ph))*exp(4.*wa_ph);

				//the speed should not exceed 1
				U2L=min(U2L,1.-1.0e-10);
				U2R=min(U2R,1.-1.0e-10);
				//calculation of left and right U_i end
				
				//sound velocity
				cs2L=dpres(rhoL);
				cs2R=dpres(rhoR);
				
				//calculation of characteristic speeds start
				ovfacL=alp_ph/(1.-U2L*cs2L);
				ovfacR=alp_ph/(1.-U2R*cs2R);
				
				lam1tL=UzL*(1-cs2L);
				lam2tL=sqrt(max(cs2L*(1-U2L)*(gizz_ph*(1-U2L*cs2L)-(1-cs2L)*UzL*UzL),1.0e-16));
				
				lam1tR=UzR*(1-cs2R);
				lam2tR=sqrt(max(cs2R*(1-U2R)*(gizz_ph*(1-U2R*cs2R)-(1-cs2R)*UzR*UzR),1.0e-16));
				
				
				lampL=ovfacL*(lam1tL+lam2tL)-bz_ph;
				lampR=ovfacR*(lam1tR+lam2tR)-bz_ph;
				
				lammL=ovfacL*(lam1tL-lam2tL)-bz_ph;
				lammR=ovfacR*(lam1tR-lam2tR)-bz_ph;
				
				
				as=max({abs(VzL),abs(VzR),abs(lampL),abs(lampR),abs(lammL),abs(lammR)});
				//calculation of characteristic speeds end
				
				//pressure
				prsL=pres(rhoL);
				prsR=pres(rhoR);
				
				//fluxes
				set_flux_z(l,k,j,0)=0.5*(S0L*VzL+S0R*VzR+sqdet*(prsL*(VzL+bz_ph)+prsR*(VzR+bz_ph))-as*(S0R-S0L));
				set_flux_z(l,k,j,1)=0.5*(SxL*VzL+SxR*VzR-as*(SxR-SxL));
				set_flux_z(l,k,j,2)=0.5*(SyL*VzL+SyR*VzR-as*(SyR-SyL));
				set_flux_z(l,k,j,3)=0.5*(SzL*VzL+SzR*VzR+sqdet*alp_ph*(prsL+prsR)-as*(SzR-SzL));
				set_flux_z(l,k,j,4)=0.5*(DenL*VzL+DenR*VzR-as*(DenR-DenL));
				
				
			}
		}
	}
}

