/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* APARENT HORISON SOLVER :: BSSN evolution Class of COSMOS                                              */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "ahf2d.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

void Ahf2d::find_ah(Fmv0 *fmv,int loopmax,double err_p,double err_eps,ofstream& fout,short int hsign)
{
	int loop=0;
	
	//apparent horizon or cosmological horizon
	//ahorch=1 for apparent horizon (outer trapped)
	//ahorch=3 for cosmological horizon (inner trapped)
	//the value of var is empirically obtained optimized value(AH(CH) cannot be found with var=3(1))
	//Ref.arXiv:1404.1435
	int ahorch=int(var+0.001);
	
	//hsign is direction of the normal null vector
	//Ref.arXiv:1404.1435
	string horizonkind;
	if(ahorch==1 && hsign==1)
	 horizonkind= "AH+";
	else if(ahorch==1 && hsign==-1)
	 horizonkind= "AH-";
	else if(ahorch==3 && hsign==1)
	 horizonkind= "CH+";
	else
	 horizonkind= "CH-";
	
	//reset the initial radius if horizon has not be found yet
	if(error || hb[0][0]<1.0e-4) 
	{
		cout << "reset " << horizonkind << " finder" << endl;
		//set_hini(hini);
		set_hini(0.5*fmv->funcf(fmv->get_xu()));
	}
	
	//value of maximum and minimum expansion
	//to see future- or past-
	//Ref.arXiv:1811.00762
	double expmax=0.;
	double expmin=0.;
	
	for(loop=0;loop<loopmax;loop++)
	{
		expmax=-9999999999999999999.;
		expmin=99999999999999999999.;
		
		bool level_mismatch=false;
		bool too_small=false;

		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			//sin(phi) and cos(phi)
			double sp=sin(phi[k]);
			double cp=cos(phi[k]);

			for(int j=1;j<(mt-1);j++)
			{
				double st=sin(theta[j]);
				double ct=cos(theta[j]);
				double ta=tan(theta[j]);
				
				double hc=hb[k][j];							//previous value as an initial guess
				//for non-Cartesian inhomogeneous grid
				double xc=fmv->ifuncf(hc*st*cp);			//x=r*sin(theta)*cos(phi)
				double yc=fmv->ifuncf(hc*st*sp);			//y=r*sint(theta)*sin(phi)
				double zc=fmv->ifuncf(hc*ct);				//z=r*cos(theta)
				
				//reference metric quantities
				double dfx=fmv->df(xc);
				double dfy=fmv->df(yc);
				double dfz=fmv->df(zc);
				double dfix=1./dfx;
				double dfiy=1./dfy;
				double dfiz=1./dfz;
				double ddfx=fmv->ddf(xc);
				double ddfy=fmv->ddf(yc);
				double ddfz=fmv->ddf(zc);
				
				////////// prep. for coordinate transformation to spherical coordinates //////////
				double dxdr=st*cp*dfix;															//
				double dydr=st*sp*dfiy;
				double dzdr=ct*dfiz;
				double dxdt=hc*ct*cp*dfix;
				double dydt=hc*ct*sp*dfiy;
				double dzdt=-hc*st*dfiz;
				double dxdp=-hc*st*sp*dfix;
				double dydp=hc*st*cp*dfiy;
				double dzdp=0.;
				
				double ddfdfix=ddfx*dfix;
				double ddfdfiy=ddfy*dfiy;
				double ddfdfiz=ddfz*dfiz;
				
				// d2xdrt = \frac{\partial^2 x}{\partial\theta\partial r}
				double d2xdrr=0.-pow(dxdr,2)*ddfdfix;
				double d2xdrt= ct*cp*dfix-dxdr*dxdt*ddfdfix;
				double d2xdrp=-st*sp*dfix-dxdr*dxdp*ddfdfix;
				double d2ydrr=0.-dydr*dydr*ddfdfiy;
				double d2ydrt= ct*sp*dfiy-dydr*dydt*ddfdfiy;
				double d2ydrp= st*cp*dfiy-dydr*dydp*ddfdfiy;
				double d2zdrr=0.-dzdr*dzdr*ddfdfiz;
				double d2zdrt=-st*dfiz-dzdr*dzdt*ddfdfiz;
				double d2zdrp=0.-dzdr*dzdp*ddfdfiz;

				double d2xdtr=    ct*cp*dfix-dxdr*dxdr*ddfdfix;
				double d2xdtt=-hc*st*cp*dfix-dxdt*dxdt*ddfdfix;
				double d2xdtp=-hc*ct*sp*dfix-dxdt*dxdp*ddfdfix;
				double d2ydtr=    ct*sp*dfiy-dydt*dydr*ddfdfiy;
				double d2ydtt=-hc*st*sp*dfiy-dydt*dydt*ddfdfiy;
				double d2ydtp= hc*ct*cp*dfiy-dydt*dydp*ddfdfiy;
				double d2zdtr=   -st*dfiz-dzdt*dzdr*ddfdfiz;
				double d2zdtt=-hc*ct*dfiz-dzdt*dzdt*ddfdfiz;
				double d2zdtp=0.-dzdt*dzdp*ddfdfiz;

				double d2xdpr=   -st*sp*dfix-dxdp*dxdr*ddfdfix;
				double d2xdpt=-hc*ct*sp*dfix-dxdp*dxdt*ddfdfix;
				double d2xdpp=-hc*st*cp*dfix-dxdp*dxdp*ddfdfix;
				double d2ydpr=    st*cp*dfiy-dydp*dydr*ddfdfiy;
				double d2ydpt= hc*ct*cp*dfiy-dydp*dydt*ddfdfiy;
				double d2ydpp=-hc*st*sp*dfiy-dydp*dydp*ddfdfiy;
				double d2zdpr=0.-dzdp*dzdr*ddfdfiz;
				double d2zdpt=0.-dzdp*dzdt*ddfdfiz;
				double d2zdpp=0.-dzdp*dzdp*ddfdfiz;

				double idrdx=st*cp*dfx;
				double idrdy=st*sp*dfy;
				double idrdz=ct*dfz;
				double idtdx=ct*cp/hc*dfx;
				double idtdy=ct*sp/hc*dfy;
				double idtdz=-st/hc*dfz;
				double idpdx=-sp/(st*hc)*dfx;
				double idpdy=cp/(st*hc)*dfy;
				double idpdz=0.;																//
				////////// prep. for coordinate transformation to spherical coordinates //////////
				
				//coordinate deviation inverse
				double dxi=fmv->get_dxi();
				double dyi=fmv->get_dyi();
				double dzi=fmv->get_dzi();

				int jmin=fmv->get_jmin();
				int kmin=fmv->get_kmin();
				int lmin=fmv->get_lmin();
				int jmax=fmv->get_jmax();
				int kmax=fmv->get_kmax();
				int lmax=fmv->get_lmax();
				int jc,kc,lc;
				
				//nearest lower grid 
				jc=int(xc*dxi);
				kc=int(yc*dyi);
				lc=int(zc*dzi);
				
				///////// run off error //////////////////////////////////////////////////////////
				if( (jc>(jmax-3)) || (jc<(jmin+3)) ||											//
				(kc>(kmax-3)) || (kc<(kmin+3)) ||
				(lc>(lmax-3)) || (lc<(lmin+3)) )
				{
					#pragma omp critical (error_output)
					{
						//cout << horizonkind << " finder Error :: Level mismatch   "
						// << jc << "," << kc << "," << lc << ", h=" << hc
						// << endl;
						level_mismatch=true;
					}
					
					continue;
					
				}
				if(fmv->ifuncf(hb[0][0])<fmv->get_dx())
				{
				    #pragma omp critical (error_output)
					{
						//cout << "too small" << endl;
				    	too_small=true;
					}
					
					continue;
				}																				//
				///////// run off error //////////////////////////////////////////////////////////

				///////// metric variables ///////////////////////////////////////////////////////
				double mgxx,mgyy,mgzz,mgxy,mgxz,mgyz,mgyx,mgzx,mgzy,mphi,mpsi;					//
				double makxx,makyy,makzz,makxy,makxz,makyz,makyx,makzx,makzy,mek;
				double mgixx,mgiyy,mgizz,mgixy,mgixz,mgiyz,mgiyx,mgizx,mgizy,mdet;
				double mgrr,mgtt,mgpp,mgrt,mgrp,mgtp;
				double mgirr,mgitt,mgipp,mgirt,mgirp,mgitp,mgitr,mgipr,mgipt;
				double mzgx,mzgy,mzgz;
				double hxx,hyy,hzz;
				
				
				//substitution <- interpolation
				//for non-Cartesian
				hxx  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4, 7);
				hyy  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4, 8);
				hzz  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4, 9);
				mgxx  = hxx+pow(dfx,2);
				mgyy  = hyy+pow(dfy,2);
				mgzz  = hzz+pow(dfz,2);
				mgxy  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,10);
				mgxz  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,11);
				mgyz  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,12);
				mgyx=mgxy;
				mgzx=mgxz;
				mgzy=mgyz;
				mphi  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,13);
				// \psi:=exp{\phi} "phi" is "w" in Fmv 
				mpsi  = exp(mphi);

				mdet=mgxx*mgyy*mgzz +mgxy*mgyz*mgzx +mgxz*mgyx*mgzy
				-mgxz*mgyy*mgzx -mgxy*mgyx*mgzz -mgxx*mgyz*mgzy;
				mgixx= (mgyy*mgzz-mgyz*mgzy)/mdet;
				mgiyy= (mgxx*mgzz-mgxz*mgzx)/mdet;
				mgizz= (mgxx*mgyy-mgxy*mgyx)/mdet;
				mgixy=-(mgyx*mgzz-mgyz*mgzx)/mdet;
				mgixz= (mgyx*mgzy-mgyy*mgzx)/mdet;
				mgiyz=-(mgxx*mgzy-mgxy*mgzx)/mdet;
				mgiyx=mgixy;
				mgizx=mgixz;
				mgizy=mgiyz;
				
				//coordinate transformation
				mgrr=dxdr*dxdr*mgxx +dxdr*dydr*mgxy +dxdr*dzdr*mgxz
					+dydr*dxdr*mgyx +dydr*dydr*mgyy +dydr*dzdr*mgyz
					+dzdr*dxdr*mgzx +dzdr*dydr*mgzy +dzdr*dzdr*mgzz;
				mgtt=dxdt*dxdt*mgxx +dxdt*dydt*mgxy +dxdt*dzdt*mgxz
					+dydt*dxdt*mgyx +dydt*dydt*mgyy +dydt*dzdt*mgyz
					+dzdt*dxdt*mgzx +dzdt*dydt*mgzy +dzdt*dzdt*mgzz;
				mgpp=dxdp*dxdp*mgxx +dxdp*dydp*mgxy +dxdp*dzdp*mgxz
					+dydp*dxdp*mgyx +dydp*dydp*mgyy +dydp*dzdp*mgyz
					+dzdp*dxdp*mgzx +dzdp*dydp*mgzy +dzdp*dzdp*mgzz;
				mgrt=dxdr*dxdt*mgxx +dxdr*dydt*mgxy +dxdr*dzdt*mgxz
					+dydr*dxdt*mgyx +dydr*dydt*mgyy +dydr*dzdt*mgyz
					+dzdr*dxdt*mgzx +dzdr*dydt*mgzy +dzdr*dzdt*mgzz;
				mgrp=dxdr*dxdp*mgxx +dxdr*dydp*mgxy +dxdr*dzdp*mgxz
					+dydr*dxdp*mgyx +dydr*dydp*mgyy +dydr*dzdp*mgyz
					+dzdr*dxdp*mgzx +dzdr*dydp*mgzy +dzdr*dzdp*mgzz;
				mgtp=dxdt*dxdp*mgxx +dxdt*dydp*mgxy +dxdt*dzdp*mgxz
					+dydt*dxdp*mgyx +dydt*dydp*mgyy +dydt*dzdp*mgyz
					+dzdt*dxdp*mgzx +dzdt*dydp*mgzy +dzdt*dzdp*mgzz;
				
				mgirr=idrdx*idrdx*mgixx +idrdx*idrdy*mgixy +idrdx*idrdz*mgixz
					+idrdy*idrdx*mgiyx +idrdy*idrdy*mgiyy +idrdy*idrdz*mgiyz
					+idrdz*idrdx*mgizx +idrdz*idrdy*mgizy +idrdz*idrdz*mgizz;
				mgitt=idtdx*idtdx*mgixx +idtdx*idtdy*mgixy +idtdx*idtdz*mgixz
					+idtdy*idtdx*mgiyx +idtdy*idtdy*mgiyy +idtdy*idtdz*mgiyz
					+idtdz*idtdx*mgizx +idtdz*idtdy*mgizy +idtdz*idtdz*mgizz;
				mgipp=idpdx*idpdx*mgixx +idpdx*idpdy*mgixy +idpdx*idpdz*mgixz
					+idpdy*idpdx*mgiyx +idpdy*idpdy*mgiyy +idpdy*idpdz*mgiyz
					+idpdz*idpdx*mgizx +idpdz*idpdy*mgizy +idpdz*idpdz*mgizz;
				mgirt=idrdx*idtdx*mgixx +idrdx*idtdy*mgixy +idrdx*idtdz*mgixz
					+idrdy*idtdx*mgiyx +idrdy*idtdy*mgiyy +idrdy*idtdz*mgiyz
					+idrdz*idtdx*mgizx +idrdz*idtdy*mgizy +idrdz*idtdz*mgizz;
				mgirp=idrdx*idpdx*mgixx +idrdx*idpdy*mgixy +idrdx*idpdz*mgixz
					+idrdy*idpdx*mgiyx +idrdy*idpdy*mgiyy +idrdy*idpdz*mgiyz
					+idrdz*idpdx*mgizx +idrdz*idpdy*mgizy +idrdz*idpdz*mgizz;
				mgitp=idtdx*idpdx*mgixx +idtdx*idpdy*mgixy +idtdx*idpdz*mgixz
					+idtdy*idpdx*mgiyx +idtdy*idpdy*mgiyy +idtdy*idpdz*mgiyz
					+idtdz*idpdx*mgizx +idtdz*idpdy*mgizy +idtdz*idpdz*mgizz;
				mgipr=mgirp;
				mgitr=mgirt;
				mgipt=mgitp;																	//
				///////// metric variables ///////////////////////////////////////////////////////
				
				//second order finite differences of h(theta,phi)
				double ht = ( hb[k][j+1] -hb[k][j-1] )*dti*0.5;
				double hp = ( hb[k+1][j] -hb[k-1][j] )*dpi*0.5;

				double htt = ( hb[k][j+1] +hb[k][j-1] -2.*hc )*dti2;
				double hpp = ( hb[k+1][j] +hb[k-1][j] -2.*hc )*dpi2;
				double htp =  ( (hb[k+1][j+1] -hb[k+1][j-1])
						-(hb[k-1][j+1] -hb[k-1][j-1]))*dti*dpi*0.25;

				///////// vectors ////////////////////////////////////////////////////////////////
				// \tilde{s}_i :: spherical coordinates											//
				// f(r,theta,phi):=r-h(theta,phi)=0
				double sdr = 1.;//del_r f
				double sdt = -ht;//-del_theta f
				double sdp = -hp;//-del_phi f
				double cn = 1./sqrt( mgirr*sdr*sdr +mgirt*sdr*sdt +mgirp*sdr*sdp
							+mgitr*sdt*sdr +mgitt*sdt*sdt +mgitp*sdt*sdp
							+mgipr*sdp*sdr +mgipt*sdp*sdt +mgipp*sdp*sdp);
							//normalziation

				// \tilde{s}^i = \tilde{\gamma}^{ij}\tilde{s}_{j}
				double sur = mgirr*sdr +mgirt*sdt +mgirp*sdp;
				double sut = mgitr*sdr +mgitt*sdt +mgitp*sdp;
				double sup = mgipr*sdr +mgipt*sdt +mgipp*sdp;

				// // \tilde{s}_I = dx^i/dx^I \tilde{s}_{I} :: Cartesian coordinates
				double sdx = idrdx*sdr +idtdx*sdt +idpdx*sdp;
				double sdy = idrdy*sdr +idtdy*sdt +idpdy*sdp;
				double sdz = idrdz*sdr +idtdz*sdt +idpdz*sdp;
				// // \tilde{s}^I = dx^I/dx^i \tilde{s}^{j} :: Cartesian coordinates
				double sx = mgixx*sdx +mgixy*sdy +mgixz*sdz;
				double sy = mgiyx*sdx +mgiyy*sdy +mgiyz*sdz;
				double sz = mgizx*sdx +mgizy*sdy +mgizz*sdz;									//
				///////// vectors ////////////////////////////////////////////////////////////////

				///////// extrinsic curvature ////////////////////////////////////////////////////
				makxx = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,14);									//
				makyy = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,15);
				makzz = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,16);
				makxy = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,17);
				makxz = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,18);
				makyz = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,19);
				makyx=makxy;
				makzx=makxz;
				makzy=makyz;
				mek   = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,20);									//
				///////// extrinsic curvature ////////////////////////////////////////////////////

				// kssmk :=(K_{ij}\tilde s^i \tilde s^j - K)*(some factor) 
				//        = (C^2\tilde{A}_{ij}\tilde{s}^i\tilde{s}^j -2/3*K)*h^2*\psi^2/C^3
				double kssmk = (makxx*sx*sx +makxy*sx*sy +makxz*sx*sz
								+makyx*sy*sx +makyy*sy*sy +makyz*sy*sz
								+makzx*sz*sx +makzy*sz*sy +makzz*sz*sz)
								*hc*hc*mpsi*mpsi/cn -2.*mek*mpsi*mpsi*hc*hc/(3.*cn*cn*cn);
				
				///////// calculation of D_i s^i /////////////////////////////////////////////////
				double mdxx_x,mdxy_x,mdxz_x,mdyx_x,mdyy_x,mdyz_x,mdzx_x,mdzy_x,mdzz_x;			//
				double mdxx_y,mdxy_y,mdxz_y,mdyx_y,mdyy_y,mdyz_y,mdzx_y,mdzy_y,mdzz_y;
				double mdxx_z,mdxy_z,mdxz_z,mdyx_z,mdyy_z,mdyz_z,mdzx_z,mdzy_z,mdzz_z;
				double mphix,mphiy,mphiz;

				double facx=(xc-fmv->get_x(jc))*dxi;
				double facy=(yc-fmv->get_y(kc))*dyi;
				double facz=(zc-fmv->get_z(lc))*dzi;
				
				//\tilde \gamma_{ij,k}/2  !!half of the derivative!!
				//modified for non-Cartesian 
				mdxx_x = 0.5*fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,hxx,4, 7)+dfx*ddfx;
				mdyy_x = 0.5*fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,hyy,4, 8);
				mdzz_x = 0.5*fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,hzz,4, 9);
				mdxy_x = 0.5*fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,mgxy,4,10);
				mdxz_x = 0.5*fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,mgxz,4,11);
				mdyz_x = 0.5*fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,mgyz,4,12);
				mdyx_x = mdxy_x;
				mdzx_x = mdxz_x;
				mdzy_x = mdyz_x;

				mdxx_y = 0.5*fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,hxx,4, 7);
				mdyy_y = 0.5*fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,hyy,4, 8)+dfy*ddfy;
				mdzz_y = 0.5*fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,hzz,4, 9);
				mdxy_y = 0.5*fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,mgxy,4,10);
				mdxz_y = 0.5*fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,mgxz,4,11);
				mdyz_y = 0.5*fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,mgyz,4,12);
				mdyx_y = mdxy_y;
				mdzx_y = mdxz_y;
				mdzy_y = mdyz_y;

				mdxx_z = 0.5*fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,hxx,4, 7);
				mdyy_z = 0.5*fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,hyy,4, 8);
				mdzz_z = 0.5*fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,hzz,4, 9)+dfz*ddfz;
				mdxy_z = 0.5*fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,mgxy,4,10);
				mdxz_z = 0.5*fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,mgxz,4,11);
				mdyz_z = 0.5*fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,mgyz,4,12);
				mdyx_z = mdxy_z;
				mdzx_z = mdxz_z;
				mdzy_z = mdyz_z;
				
				//\phi_{,i}  !! not half !!
				mphix = fmv->bv_ipol_diff_x(jc,kc,lc,xc,yc,zc,facx,mphi,4,13);
				mphiy = fmv->bv_ipol_diff_y(jc,kc,lc,xc,yc,zc,facy,mphi,4,13);
				mphiz = fmv->bv_ipol_diff_z(jc,kc,lc,xc,yc,zc,facz,mphi,4,13);

				double mpsix=mpsi*mphix;
				double mpsiy=mpsi*mphiy;
				double mpsiz=mpsi*mphiz;

				double mpsir=mpsix*dxdr +mpsiy*dydr +mpsiz*dzdr;
				double mpsit=mpsix*dxdt +mpsiy*dydt +mpsiz*dzdt;
				double mpsip=mpsix*dxdp +mpsiy*dydp +mpsiz*dzdp;

				//-0.5*\del_r C^{-2}=-0.5*\del_r(\tilde s^i \tilde s_i)=C^{-3} \del_r C
				double gss_r= dxdr*(mdxx_x*sx*sx +mdxy_x*sx*sy +mdxz_x*sx*sz
							+mdyx_x*sy*sx +mdyy_x*sy*sy +mdyz_x*sy*sz
							+mdzx_x*sz*sx +mdzy_x*sz*sy +mdzz_x*sz*sz)
						+dydr*(mdxx_y*sx*sx +mdxy_y*sx*sy +mdxz_y*sx*sz
							+mdyx_y*sy*sx +mdyy_y*sy*sy +mdyz_y*sy*sz
							+mdzx_y*sz*sx +mdzy_y*sz*sy +mdzz_y*sz*sz)
						+dzdr*(mdxx_z*sx*sx +mdxy_z*sx*sy +mdxz_z*sx*sz
							+mdyx_z*sy*sx +mdyy_z*sy*sy +mdyz_z*sy*sz
							+mdzx_z*sz*sx +mdzy_z*sz*sy +mdzz_z*sz*sz)
						+d2xdrr*sur*sdx +d2xdtr*sut*sdx +d2xdpr*sup*sdx
						+d2ydrr*sur*sdy +d2ydtr*sut*sdy +d2ydpr*sup*sdy
						+d2zdrr*sur*sdz +d2zdtr*sut*sdz +d2zdpr*sup*sdz;
							
				//-0.5*\del_t C^{-2}+\tilde s^i \del_t \tilde s_i
				//=-0.5*\del_t(\tilde s^i \tilde s_i)+\tilde s^i \del_t \tilde s_i
				//=C^{-3} \del_t C+\tilde s^i \del_t \tilde s_i
				double gss_t= dxdt*(mdxx_x*sx*sx +mdxy_x*sx*sy +mdxz_x*sx*sz
							+mdyx_x*sy*sx +mdyy_x*sy*sy +mdyz_x*sy*sz
							+mdzx_x*sz*sx +mdzy_x*sz*sy +mdzz_x*sz*sz)
						+dydt*(mdxx_y*sx*sx +mdxy_y*sx*sy +mdxz_y*sx*sz
							+mdyx_y*sy*sx +mdyy_y*sy*sy +mdyz_y*sy*sz
							+mdzx_y*sz*sx +mdzy_y*sz*sy +mdzz_y*sz*sz)
						+dzdt*(mdxx_z*sx*sx +mdxy_z*sx*sy +mdxz_z*sx*sz
							+mdyx_z*sy*sx +mdyy_z*sy*sy +mdyz_z*sy*sz
							+mdzx_z*sz*sx +mdzy_z*sz*sy +mdzz_z*sz*sz)
						+d2xdrt*sur*sdx +d2xdtt*sut*sdx +d2xdpt*sup*sdx
						+d2ydrt*sur*sdy +d2ydtt*sut*sdy +d2ydpt*sup*sdy
						+d2zdrt*sur*sdz +d2zdtt*sut*sdz +d2zdpt*sup*sdz;
				
				//-0.5*\del_p C^{-2}+\tilde s^i \del_p \tilde s_i
				//=-0.5*\del_p(\tilde s^i \tilde s_i)+\tilde s^i \del_p \tilde s_i
				//=C^{-3} \del_p C+\tilde s^i \del_p \tilde s_i
				double gss_p= dxdp*(mdxx_x*sx*sx +mdxy_x*sx*sy +mdxz_x*sx*sz
							+mdyx_x*sy*sx +mdyy_x*sy*sy +mdyz_x*sy*sz
							+mdzx_x*sz*sx +mdzy_x*sz*sy +mdzz_x*sz*sz)
						+dydp*(mdxx_y*sx*sx +mdxy_y*sx*sy +mdxz_y*sx*sz
							+mdyx_y*sy*sx +mdyy_y*sy*sy +mdyz_y*sy*sz
							+mdzx_y*sz*sx +mdzy_y*sz*sy +mdzz_y*sz*sz)
						+dzdp*(mdxx_z*sx*sx +mdxy_z*sx*sy +mdxz_z*sx*sz
							+mdyx_z*sy*sx +mdyy_z*sy*sy +mdyz_z*sy*sz
							+mdzx_z*sz*sx +mdzy_z*sz*sy +mdzz_z*sz*sz)
						+d2xdrp*sur*sdx +d2xdtp*sut*sdx +d2xdpp*sup*sdx
						+d2ydrp*sur*sdy +d2ydtp*sut*sdy +d2ydpp*sup*sdy
						+d2zdrp*sur*sdz +d2zdtp*sut*sdz +d2zdpp*sup*sdz;

				// c_i = C_{,i}/C^3
				double c_r=gss_r;
				double c_t=gss_t +htt*sut +htp*sup;
				double c_p=gss_p +htp*sut +hpp*sup;
				
				
				//modified by chulmoon on 2016.3.27
				mzgx  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,21)
						+2.*ddfx*dfix*mgixx+ddfy*dfiy*mgixy+ddfz*dfiz*mgixz;
				mzgy  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,22)
						+2.*ddfy*dfiy*mgiyy+ddfx*dfix*mgixy+ddfz*dfiz*mgiyz;
				mzgz  = fmv->bv_ipol(jc,kc,lc,xc,yc,zc,4,23)
						+2.*ddfz*dfiz*mgizz+ddfy*dfiy*mgiyz+ddfx*dfix*mgixz;
				
				double gammas=sdx*(d2xdrr*mgirr +d2xdrt*mgirt +d2xdrp*mgirp
							+d2xdtr*mgitr +d2xdtt*mgitt +d2xdtp*mgitp
							+d2xdpr*mgipr +d2xdpt*mgipt +d2xdpp*mgipp)
							+sdy*(d2ydrr*mgirr +d2ydrt*mgirt +d2ydrp*mgirp
							+d2ydtr*mgitr +d2ydtt*mgitt +d2ydtp*mgitp
							+d2ydpr*mgipr +d2ydpt*mgipt +d2ydpp*mgipp)
							+sdz*(d2zdrr*mgirr +d2zdrt*mgirt +d2zdrp*mgirp
							+d2zdtr*mgitr +d2zdtt*mgitt +d2zdtp*mgitp
							+d2zdpr*mgipr +d2zdpt*mgipt +d2zdpp*mgipp)
							
							+sur*(d2xdrr*idrdx +d2xdrt*idtdx +d2xdrp*idpdx
							+d2ydrr*idrdy +d2ydrt*idtdy +d2ydrp*idpdy
							+d2zdrr*idrdz +d2zdrt*idtdz +d2zdrp*idpdz)
							+sut*(d2xdtr*idrdx +d2xdtt*idtdx +d2xdtp*idpdx
							+d2ydtr*idrdy +d2ydtt*idtdy +d2ydtp*idpdy
							+d2zdtr*idrdz +d2zdtt*idtdz +d2zdtp*idpdz)
							+sup*(d2xdpr*idrdx +d2xdpt*idtdx +d2xdpp*idpdx
							+d2ydpr*idrdy +d2ydpt*idtdy +d2ydpp*idpdy
							+d2zdpr*idrdz +d2zdpt*idtdz +d2zdpp*idpdz)
							
							+mzgx*sdx +mzgy*sdy +mzgz*sdz; 

				//double xiirr = mgirr-1.;
				//double xiitt = mgitt-1./(hc*hc);
				//double xiipp = mgipp-1./(hc*hc*st*st);

				// disi = D_i s^i * h^2 * \psi^2 / C^3 +linear\Delta h -(2-\eta) h
				// D_i s^i=1/\sqrt{\gamma}\del_i(\sqrt{\gamma}s^i)
				// in Cartesian \sqrt{\gamma}=\psi^6  !!not in spherical coord.!!
				// for non-Cartesian inhomogeneous grid 
				double disi= var*hc
							+(2.*hc*sur +hc*hc*sut/ta)/(cn*cn)
							-(2.*hc -ht/ta)				// sent to the left-hand-side of the equation
							+hc*hc*(-mgitt*htt -mgipp*hpp-2.*mgitp*htp)/(cn*cn)
							-(-htt -hpp/(st*st))			// sent to the left-hand-side of the equation
							+(c_r*sur +c_t*sut +c_p*sup)*hc*hc 
							+4.*(mpsir*sur +mpsit*sut +mpsip*sup)*hc*hc/(mpsi*cn*cn) 
							-gammas*hc*hc/(cn*cn);
				
				src[k][j] = disi +hsign*kssmk;													//
				///////// calculation of D_i s^i /////////////////////////////////////////////////
				
				///////// FOTH or POTH ///////////////////////////////////////////////////////////
				double expv=2.*kssmk*pow(hc*mpsi,-2)*pow(cn,3);									//
				if(expv>expmax)
				expmax=expv;
				if(expv<expmin)
				expmin=expv;																	//
				///////// FOTH or POTH ///////////////////////////////////////////////////////////
				
				///////// calculation for horizon geometry ////////////////////////////////////////
				double da = sqrt( (mgrr*ht*ht +mgtt +2.*mgrt*ht)*(mgrr*hp*hp +mgpp +2.*mgrp*hp)	//
							-pow((mgrr*ht*hp +mgrt*hp +mgrp*ht +mgtp),2) );
				ds[k][j] = pow(mpsi,4)*da;
				
				//for the length arround equator
				if( j==mt-3 )
				{
					dce[k][0] = sqrt(mgrr*hp*hp +mgpp +2.*mgrp*hp)*mpsi*mpsi;
				}
				if( j==mt-2 )
				{
					dce[k][1] = sqrt(mgrr*hp*hp +mgpp +2.*mgrp*hp)*mpsi*mpsi;
					dce[k][2] = dce[k][1];
				}
				//if( j==(mt+2)/2 )
				//{
				//	dce[k][2] = sqrt(mgrr*hp*hp +mgpp +2.*mgrp*hp)*mpsi*mpsi;
				//}
				
				//for length of two meridians
				//phi=0
				if(k==1)
				{
					dcp1[j][0] = sqrt(mgrr*ht*ht +mgtt +2.*mgrt*ht)*mpsi*mpsi;
					dcp1[j][1] = dcp1[j][0];
				}
				if(k==2)
				{
					dcp1[j][2] =sqrt(mgrr*ht*ht +mgtt +2.*mgrt*ht)*mpsi*mpsi;
				}
				
				//phi=pi/2
				if(k==mp/2-1)
				{
					dcp2[j][0] = sqrt(mgrr*ht*ht +mgtt +2.*mgrt*ht)*mpsi*mpsi;
				}	
				if(k==mp/2)
				{
					dcp2[j][1] = sqrt(mgrr*ht*ht +mgtt +2.*mgrt*ht)*mpsi*mpsi;
				}
				if(k==mp/2+1)
				{
					dcp2[j][2] = sqrt(mgrr*ht*ht +mgtt +2.*mgrt*ht)*mpsi*mpsi;
				}																				//
				///////// calculation for horizon geometry ///////////////////////////////////////
			}
		}
		
		if(too_small)
		cout << "too small" << endl;
		if(level_mismatch)
		cout << "level_mismatch" << endl;

		if(too_small || level_mismatch)
		return;

		///// Boundary Condition for the quantities related to AH ////////////////////////////////////////
		boundaryset(ds);																				//
		
		//for(int k=0;k<mp;k++)
		//{
		//	ds[k][0] = ds[k][1];
		//	ds[k][mt-1] = ds[k][mt-2];
		//}
		//for(int j=0;j<mt;j++)
		//{
			//ds[0][j] = ds[mp-2][j];
		//	ds[0][j] = ds[1][j];
			//ds[mp-1][j] = ds[1][j];
		//	ds[mp-1][j] = ds[mp-2][j];
		//}

		//dce[0][0] = dce[mp-2][0];
		//dce[mp-1][0] = dce[1][0];
		//dce[0][1] = dce[mp-2][1];
		//dce[mp-1][1] = dce[1][1];
		//dce[0][2] = dce[mp-2][2];
		//dce[mp-1][2] = dce[1][2];
		
		dce[0][0] = dce[1][0];
		dce[mp-1][0] = dce[mp-2][0];
		dce[0][1] = dce[1][1];
		dce[mp-1][1] = dce[mp-2][1];
		dce[0][2] = dce[1][2];
		dce[mp-1][2] = dce[mp-2][2];
		
		dcp1[0][0] = dcp1[1][0];
		dcp1[mt-1][0] = dcp1[mt-2][0];
		dcp1[0][1] = dcp1[1][1];
		dcp1[mt-1][1] = dcp1[mt-2][1];
		dcp1[0][2] = dcp1[1][2];
		dcp1[mt-1][2] = dcp1[mt-2][2];

		dcp2[0][0] = dcp2[1][0];
		dcp2[mt-1][0] = dcp2[mt-2][0];
		dcp2[0][1] = dcp2[1][1];
		dcp2[mt-1][1] = dcp2[mt-2][1];
		dcp2[0][2] = dcp2[1][2];
		dcp2[mt-1][2] = dcp2[mt-2][2];
		
		
		boundaryset(src);																				//
		///// Boundary Condition for the quantities related to AH ////////////////////////////////////////
		
		/////////////////////////////////////////
		// Poisson Solver
		/////////////////////////////////////////
		poisson(ha, src, err_p);

		double errv=0.;
		double sum=0.;
		#pragma omp parallel for reduction(+:errv,sum)
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				double st=sin(theta[j]);

				double dtmp = dt*dp*st;
				sum += dtmp;
				errv+= (ha[k][j]-hb[k][j])*(ha[k][j]-hb[k][j])*dtmp;
			}
		}
		#pragma omp barrier
		double errtmp=sqrt(errv/sum);

		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				hb[k][j] = fac*ha[k][j] +(1.-fac)*hb[k][j];
			}
		}
		#pragma omp barrier
		
		boundaryset(hb);
		
		cout << "      AHF step=" << loop << " : err=" << errtmp 
		<< " : R(0,0)=" << hb[1][1] 
		<< " : R(Pi/2,0)=" << hb[1][mt-1] 
		<< " : R(Pi/2,Pi/2)=" << hb[mp/2][mt-1] 
		<< " : r(0,0)/zmax=" << fmv->ifuncf(hb[1][1])/ fmv->get_zu()
		<< " : r(Pi/2,0)/xmax=" << fmv->ifuncf(hb[1][mt-1]) / fmv->get_xu()
		<< " : r(Pi/2,Pi/2)/ymax=" << fmv->ifuncf(hb[mp/2][mt-1]) / fmv->get_yu()
		//<< " : area=" << area*8. 
		<< " : expmax=" << expmax
		<< " : expmin=" << expmin
		<<endl;

		if(errtmp <= err_eps)
		{
		  error=false;
		  break;
		}	  
	}
	
	if(loop==loopmax)
	{
		cout <<  horizonkind << " finder Error :: No convergence"<< endl;
		error=true;
		return;
	}
	
	////// Numerical integration for AH area, length of equator and meridian /////////////////////////
	area=0.;																						//
	circe=0.;
	circp1=0.;
	circp2=0.;
	
	#pragma omp parallel for reduction(+:area,circe,circp1,circp2)
	for(int k=0;k<mp;k++)
	{
		for(int j=0;j<mt;j++)
		{
			double facj=1., fack=1.;
			//////////////////////////////////////////////////////////
			// Trapezoidal Rule with 2nd order Lagrange interpolation
			//////////////////////////////////////////////////////////
			if( (j==0) || (j==(mt-1)) )
			{
				facj=0.125*dt;
			}
			else if( (j==1) || (j==(mt-2)) )
			{
				facj=0.875*dt;
			}
			else
			{
				facj=1.*dt;
			}
			
			if( (k==0) || (k==(mp-1)) )
			{
				fack=0.125*dp;
			}
			else if( (k==1) || (k==(mp-2)) )
			{
				fack=0.875*dp;
			}
			else
			{
				fack=1.*dp;
			}

			area += ds[k][j]*facj*fack;
			
			//if(k==1)
			//{
			//	circp1 += dcp1[j]*facj;
			//}
			
			//if(k==mp-1)
			//{
			//	circp2 += dcp2[j]*facj;
			//}
			
			
			if( (k==0) || (k==1) || (k==2) )
			{
				int index=0;
				double facktmp=0.;
				if(k==0)
				{
					facktmp = 0.375;
					index = 0;
				}
				if(k==1)
				{
					facktmp = 0.75;
					index = 1;
				}
				if(k==2)
				{
					facktmp = -0.125;
					index = 2;
				}
				circp1 += dcp1[j][index]*facj*facktmp;
			}
			
			if( (k==mp/2-1) || (k==mp/2) || (k==mp/2+1) )
			{
				int index=0;
				double facktmp=0.;
				if(k==mp/2-1)
				{
					facktmp = 0.375;
					index = 0;
				}
				if(k==mp/2)
				{
					facktmp = 0.75;
					index = 1;
				}
				if(k==mp/2+1)
				{
					facktmp = -0.125;
					index = 2;
				}

				circp2 += dcp2[j][index]*facj*facktmp;
			}
			
			if( (j==mt-3) || (j==mt-2) || (j==mt-1))
			{
				int index=0;
				if(j==mt-3)
				{
					facj = -0.125;
					index = 0;
				}
				if(j==mt-2)
				{
					facj = 0.75;
					index = 1;
				}
				if(j==mt-1)
				{
					facj = 0.375;
					index = 2;
				}

				circe += dce[k][index]*facj*fack;
			}
		}
	}
	#pragma omp barrier																								//
	////// Numerical integration for AH area, length of equator and meridian /////////////////////////
	
	////// output for AH area, length of equator and meridian ////////////////////////////////////////
	area*=4.;																						//
	circe*=2.;
	circp1*=4.;
	circp2*=4.;
	
	mass = circe*0.25/M_PI;
	double imass = sqrt(area/16./M_PI/mass/mass);
	double tmp = pow(1. -area*0.125/(M_PI*mass*mass),2);
	spin = sqrt(abs(1.-tmp));
	if(1.-tmp<0.) spin *= -1.;
	double rp = 0.5*area/circe;
	double tt = fmv->get_t();
	double zh=0.375*ha[1][0]+0.75*ha[1][1]-0.125*ha[1][2];
	double eh=0.375*ha[1][mt]+0.75*ha[1][mt-1]-0.125*ha[1][mt-2];

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << setw(20) << tt														//1
	<< " " << setw(20) << mass													//2
	<< " " << setw(20) << imass													//3
	<< " " << setw(20) << spin//mass											//4
	<< " " << setw(20) << area//(mass*mass)										//5
	<< " " << setw(20) << circe													//6
	<< " " << setw(20) << circp1												//7
	<< " " << setw(20) << circp2												//8
	<< " " << setw(20) << rp													//9
	<< " " << setw(20) << zh													//10
	<< " " << setw(20) << eh													//11
	<< " " << setw(20) << fmv->ifuncf(ha[1][mt-1])/ fmv->get_xu()				//12
	<< " " << setw(20) << fmv->ifuncf(ha[1][1])/ fmv->get_zu()					//13
	<< " " << setw(20) << expmax												//14
	<< " " << setw(20) << expmin												//15
	<< endl;																						//
	////// output for AH area, length of equator and meridian ////////////////////////////////////////
	
	return;
}

// output for horizon figure 
void Ahf2d::print_ah(Fmv0 *fmv,ofstream& fout,double time)
{
	fout << "##print_num=" << fignum << "  t=" << time << endl;
	fignum++;
	
	int plotk=int(mp/20);
	int plotj=int(mt/20);

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);

	for(int k=0;k<mp;k++)
	{	
		if(k%plotk==0 || k==mp-1)
		{
			for(int j=1;j<mt;j++)
			{
				if(j%plotj==0 || j==mt-1)
				{
					double sp = sin(phi[k]);
					double cp = cos(phi[k]);
					double st = sin(theta[j]);
					double ct = cos(theta[j]);
					double har=ha[k][j];
					
					if(k==0 && j!=mt-1)
					{
						sp = 0.;
						cp = 1.;
						
						har=ha[k+1][j]*1.125-ha[k+2][j]*0.125;
					}
					if(k==mp-1 && j!=mt-1)
					{
						sp = 0.;
						cp = -1.;
						har=ha[k-1][j]*1.125-ha[k-2][j]*0.125;
					}
					if(j==mt-1 && k!= 0 && k!= mp-1)
					{
						st = 1.;
						ct = 0.;
						har=ha[k][j-1]*1.125-ha[k][j-2]*0.125;
					}
					if(j==mt-1 && k==0)
					{
						sp = 0.;
						cp = 1.;
						st = 1.;
						ct = 0.;
					
						double fr1=ha[k+2][j-2]*1.125-ha[k+2][j-1]*0.125;
						double fr2=ha[k+1][j-2]*1.125-ha[k+1][j-1]*0.125;
						har=fr2*1.125-fr1*0.125;
					}
					if(j==mt-1 && k==mp-1)
					{
						sp = 0.;
						cp = -1.;
						st = 1.;
						ct = 0.;
					
						double fr1=ha[k-2][j-2]*1.125-ha[k-2][j-1]*0.125;
						double fr2=ha[k-1][j-2]*1.125-ha[k-1][j-1]*0.125;
						har=fr2*1.125-fr1*0.125;
					}
					
					//har*=1.0014;
					
					fout << setw(19) << fmv->ifuncf(har*st*cp) << " " 
						<< setw(19) << fmv->ifuncf(har*st*sp) << " "
					<< setw(19) << fmv->ifuncf(har*ct) << " "
					<< setw(19) << har*st*cp << " " 
					<< setw(19) << har*st*sp << " "
					<< setw(20) << har*ct
					<< endl;
				}
			}
			fout << endl;
		}
		
	}
	fout << endl << endl;
	return;
}

