#ifndef _AHF2D_H_
#define _AHF2D_H_
/**************************************************************************************************************/
/*------------------------------------------------------------------------------------------------------------*/
/* HEADER FILE :: Apparent Horizon Finder Class for COSMOS                                                    */
/*                                                          ver. 1.00  coded by Chulmoon Yoo                  */
/*                                                          based on original code provided by Hirotada Okawa */
/* scheme is IMLUCGS(incomplete modified lower-upper decomposition and conjugate gradient squared) method     */
/* PRD55 2002 (1997)                                                                                          */
/*------------------------------------------------------------------------------------------------------------*/
/**************************************************************************************************************/
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "cosmos.h"
using namespace std;

class Ahf2d
{
	bool error;								//true if not found
	int mt;									//number of theta grid
	int mp;									//number of phi grid
	int fignum;								//number for figure plot
	
	double dt;								//d theta
	double dti;								//1/ dtheta
	double dti2;							//1/ dtheta^2
	double dp;								//d phi
	double dpi;								//1/ dphi
	double dpi2;							//1/ dphi^2
	double var;								//acceleration parameter in PRD55 2002
	double fac;								//factor of iteration update
	double *theta, *phi;					//theta and phi coordinates
	double hini;							//initial radius

	double **ha;							//previous horizon radius
	double **hb;							//updated horizon radius
	double **src;							//source term
	
	double **ds;							//area element
	double **dce;							//for equator length
	double **dcp1;							//for meridian length-1
	double **dcp2;							//for merician length-2
	
	double area;							//area
	double mass;							//mass
	double spin;							//spin
	double circe;							//equator length
	double circp1;							//meridian length-1
	double circp2;							//meridian length-2


	//for IMLUCGS poisson solver (see ``Basics of numerical computation(note in Japanese)" p.90 for names of vectors)
	double **vtmp;							//temporary vector in Poisson solver
	double **vp;							//p-vector in Poisson solver
	double **vr;							//r-vector in Poisson solver
	double **vr0;							//r0-vector in Poisson solver
	double **ve;							//e-vector in Poisson solver
	double **vh;							//h-vector in Poisson solver
	double **vcon;							//temporary used in precondition

	double **cb;							//vector for matrix non-zero components
	double **cc;							//vector for matrix non-zero components
	double **pcd;							//vector for matrix non-zero components
	double **cd;							//vector for the matrix D in IMLU decomposition
	double **ce;							//vector for matrix non-zero components
	double **cf;							//vector for matrix non-zero components

public:
	Ahf2d(int ntheta, int nphi, double va, double fa, double hc)
	{
		// Dphi: whole region of phi
		// Dthe: whole region of theta
		// nphi: inside grid number of phi
		// nthe: inside grid number of theta
		//-----------------------------------------------------
		// b*1**2**...**n*b 
		// where b:boundary, number:grid number, *:half space
		//-----------------------------------------------------
		// Therefore dphi=Dphi/nphi
		//           dthe=Dthe/nthe
		error=true;
		mt=ntheta+2; 					//ntheta +buffer 
		mp=nphi+2; 						//nphi + buffer
										//but buffer is not used
		dt=M_PI*0.5/double(ntheta);		//dtheta
		dti=1./dt;						//1/dtheta
		dti2=1./(dt*dt);				//1/dtheta^2
		dp=M_PI/double(nphi);			//dphi
		dpi=1./dp;						//1/dphi
		dpi2=1./(dp*dp);				//1/dphi^2
		var=va;							//acceleration parameter \eta in PRD55 2002 (1997)
		fac=fa;							//factor of iteration update
		hini=hc;						//initial radius
		fignum=0;						//number for figure plot
		
		///////////////  coordinate setting //////////////////////////////////////
		theta = new double[mt];							//theta coordinate
		phi   = new double[mp];							//phi coordinate
		for(int j=0;j<mt;j++){
			theta[j] = (double(j)-0.5)*dt;
		}
		for(int k=0;k<mp;k++){
			phi[k] = (double(k)-0.5)*dp;
		}																		//
		///////////////  coordinate setting //////////////////////////////////////
		
		
		///////////////  strage allocation and initial values ////////////////////
		ha = new double*[mp];													//
		hb = new double*[mp];
		src= new double*[mp];
		
		ds = new double*[mp];
		dce = new double*[mp];
		dcp1= new double*[mt];
		dcp2= new double*[mt];

		vtmp = new double*[mp];
		vp   = new double*[mp];
		vr   = new double*[mp];
		vr0  = new double*[mp];
		ve   = new double*[mp];
		vh   = new double*[mp];
		vcon = new double*[mp];

		cb   = new double*[mp];
		cc   = new double*[mp];
		pcd  = new double*[mp];
		cd   = new double*[mp];
		ce   = new double*[mp];
		cf   = new double*[mp];
		
		for(int j=0;j<mt;j++)
		{
			dcp1[j]=new double[3];
			dcp2[j]=new double[3];
			for(int i=0;i<3;i++)
			{
				dcp1[j][i] =0.;
				dcp2[j][i] =0.;
			}
		}
		
		for(int k=0;k<mp;k++)
		{
			ha[k] = new double[mt];
			hb[k] = new double[mt];
			src[k]= new double[mt];
			ds[k] = new double[mt];
			dce[k]= new double[3];
			for(int i=0;i<3;i++)
			{
				dce[k][i] =0.;
			}

			vtmp[k] = new double[mt];
			vp[k]   = new double[mt];
			vr[k]   = new double[mt];
			vr0[k]  = new double[mt];
			ve[k]   = new double[mt];
			vh[k]   = new double[mt];
			vcon[k] = new double[mt];

			cb[k] = new double[mt];
			cc[k] = new double[mt];
			pcd[k]= new double[mt];
			cd[k] = new double[mt];
			ce[k] = new double[mt];
			cf[k] = new double[mt];

			for(int j=0;j<mt;j++)
			{
				double st = sin(theta[j]);
				double ct = cos(theta[j]);
				
				double a=hc;
				double b=hc;
				
				ha[k][j] = a*b/sqrt(b*b*ct*ct+a*a*st*st);
				hb[k][j] = a*b/sqrt(b*b*ct*ct+a*a*st*st);

				//ha[k][j] = hc;
				//hb[k][j] = hc;
				src[k][j]= 0.;

				vtmp[k][j] = 0.;
				vp[k][j]   = 0.;
				vr[k][j]   = 0.;
				vr0[k][j]  = 0.;
				ve[k][j]   = 0.;
				vh[k][j]   = 0.;
				vcon[k][j] = 0.;

				cb[k][j] = 0.;
				cc[k][j] = 0.;
				pcd[k][j]= 0.;
				cd[k][j] = 0.;
				ce[k][j] = 0.;
				cf[k][j] = 0.;
			}
		}

		area=0.;
		mass=0.;
		spin=0.;
		circe=0.;
		circp1=0.;
		circp2=0.;																//
		///////////////  strage allocation and initial values ////////////////////

		cout << "         **  Initialize matrix  " << endl;

		////////////////////////////////////////////////////////////////////////////
		//    Construct matrix components of incomplete modified LU decomposition
		//     M=L D U - R (see PRD55 2002 (1997)) with R = 0
		//     M:original matrix
		//     L:lower left triangle with diagonal section = diag(1/d_i)
		//     U:upper right triangle with diagonal section = diag(1/d_i)
		//     D:diagonal given by diag(d_i) <- diag(LUD)=diag(M) 
		//     R:given by LDU-M if needed but R=0 for reflection boundary!!
		/////////////////////////////////////////////////////////////////////////////
		//
		//components needed in M for each line
		//typical organization of a line without boundary consideration
		//
		//		k ->
		//		.....
		//	j	0..0 cc 0*(np-2) cb pcd cf 0*(np-2) ce 0..0
		//	|	0..0 0 cc 0*(np-2) cb pcd cf 0*(np-2) ce 0..0
		//	v	.....
		/////////////////////////////////////////////////////////////////////////////
		
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				double st = sin(theta[j]);
				double ta = tan(theta[j]);

				cb[k][j] = dpi2/(st*st);
				cc[k][j] = dti2 -0.5*dti/ta;

				pcd[k][j] =-2.*dti2 -2.*dpi2/(st*st) -(2.-var);

				ce[k][j] = dti2 +0.5*dti/ta;
				cf[k][j] = dpi2/(st*st);
			}
		}
		
		//reflection boundary for phi (since matrix R=0)
		//buffer values for initial trial must be zero
		//should be commented out for other boundary cond.
		//and also should change the boundary cond. setting just before multiple() below
		//for(int j=1;j<(mt-1);j++)
		//{
		//	pcd[1][j]+=cb[1][j];
		//	cb[1][j]=0.;
		//	pcd[mp-2][j]+=cf[mp-2][j];
		//	cf[mp-2][j]=0.;
		//}
		
		//reflection boundary for theta (since matrix R=0)
		//buffer values for initial trial must be zero
		//should be commented out for other boundary cond.
		//and also should change the boundary cond. setting just before multiple() below
		//for(int k=1;k<(mp-1);k++)
		//{
		//	pcd[k][1]+=cc[k][1];
		//	cc[k][1]=0.;
		//	pcd[k][mt-2]+=ce[k][mt-2];
		//	ce[k][mt-2]=0.;
		//}
		
		//keep pcd(pre cd) and change cd following the IMUL decomposition/////////////////
		for(int k=1;k<(mp-1);k++)														//
		{
			for(int j=1;j<(mt-1);j++)
			{
				cd[k][j] = pcd[k][j];
			}
		}
		
		//calculation of d_i
		for(int k=1;k<(mp-1);k++)
		{
			int km = k-1;
			for(int j=1;j<(mt-1);j++)
			{
				int jm = j-1;
				double tmp = cd[k][j]
						-cc[k][j]*ce[k ][jm]*cd[k ][jm]
						-cb[k][j]*cf[km][j ]*cd[km][j ];
				cd[k][j] = 1./tmp;
			}
		}																				//
		//keep pcd(pre cd) and change cd following the IMUL decomposition/////////////////
		
		
		cout << "         **  IMUL components constructed  " << endl;
	}

	~Ahf2d()
	{
		for(int j=0;j<mt;j++)
		{
			delete[] dcp1[j];
			delete[] dcp2[j];
		}
		
		for(int k=0;k<mp;k++)
		{
			delete[] cb[k];
			delete[] cc[k];
			delete[] pcd[k];
			delete[] cd[k];
			delete[] ce[k];
			delete[] cf[k];

			delete[] vtmp[k];
			delete[] vp[k];
			delete[] vr[k];
			delete[] vr0[k];
			delete[] ve[k];
			delete[] vh[k];
			delete[] vcon[k];

			delete[] ha[k];
			delete[] hb[k];
			delete[] src[k];
			delete[] ds[k];
			delete[] dce[k];
		}

		delete[] cb;
		delete[] cc;
		delete[] pcd;
		delete[] cd;
		delete[] ce;
		delete[] cf;

		delete[] vtmp;
		delete[] vp;
		delete[] vr;
		delete[] vr0;
		delete[] ve;
		delete[] vh;
		delete[] vcon;

		delete[] ha;
		delete[] hb;
		delete[] src;
		delete[] ds;

		delete[] dce;
		delete[] dcp1;
		delete[] dcp2;
		
		delete[] theta;
		delete[] phi;
	}

	///////////////////////////////////////////////
	//   Matrix multiplier
	//      input  ::  x     ( source vector )
	//             ::  ax    ( multiplied vector )
	///////////////////////////////////////////////
	void multiple(double **x,double **ax)
	{
		#pragma omp parallel for
		for(int k=0;k<mp;k++)
		{
			for(int j=0;j<mt;j++)
			{
				ax[k][j]=0.;
			}
		}
		#pragma omp barrier
		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			int kp = k+1;
			int km = k-1;

			for(int j=1;j<(mt-1);j++)
			{
				int jp = j+1;
				int jm = j-1;

				ax[k][j]=x[k][jm]*cc[k][j]+x[km][j]*cb[k][j]
						+x[k][j]*pcd[k][j]+x[kp][j]*cf[k][j]
						+x[k][jp]*ce[k][j];
				
				//the above equation corresponds to the equation below
				//double st=sin(theta[j]);
				//double ta=tan(theta[j]);
				
				//ax[k][j] = ( x[k ][jp] +x[k ][jm] -2.*x[k][j] )*dti2
				//+( x[k ][jp] -x[k ][jm] )*0.5*dti/ta
				//+( x[kp][j ] +x[km][j ] -2.*x[k][j] )*dpi2/(st*st)
				//-(2.-var)*x[k][j];
			}
		}
		#pragma omp barrier
	
		return;
	}

	//////////////////////////////////////////////
	//  Poisson solver
	//     input   :: vsrc   ( source vector )
	//             :: vf     ( solution vector )
	//             :: merr   ( max error )
	//////////////////////////////////////////////
	void poisson(double **vf,double **src,double errmax)
	{
		//iteration parameters
		int it=0,itmax=100000;
		
		//construction of the source term
		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				vr[k][j] = src[k][j];
			}
		}
		#pragma omp barrier

		//U^{-1}D^{-1}L^{-1}S
		precondition(vr);
		
		///////////// calculate the norm of source for reference /////////////////////////
		double norm_src=0.;																//
		#pragma omp parallel for reduction(+:norm_src)
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				norm_src += vr[k][j]*vr[k][j];
			}
		}
		#pragma omp barrier
		
		if(norm_src==0.){
			cout << "             ** Error :: No source term!!" << endl;
		}										//
		///////////// calculate the norm of source for reference /////////////////////////
		
		///////////// boundary setting just before multiple //////////////////////////////
		// buffer value must be 0 for multiple in the case of the matrix R=0            //
		//boundary0(vf);
		// boundary condition setting
		boundary_periodicphi_nosympoletheta(vf);										//
		///////////// boundary setting just before multiple //////////////////////////////

		//M h
		multiple(vf,vtmp);
		
		//U^{-1}D^{-1}L^{-1} (M h)
		precondition(vtmp);
		
		//construction of r0
		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				vr0[k][j] = vr[k][j]-vtmp[k][j];
			}
		}
		#pragma omp barrier

		//construction of r,p,e
		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				vr[k][j] = vr0[k][j];
				vp[k][j] = vr0[k][j];
				ve[k][j] = vr0[k][j];
			}
		}
		#pragma omp barrier
		
		//iteration start for A x = b 
		do
		{
			it++;
			
			//boundary setting if R\neq0
			boundary_periodicphi_nosympoletheta(vp);
			
			//M p
			multiple(vp,vtmp);

			//U^{-1}D^{-1}L^{-1} (M p)
			precondition(vtmp);
			
			////////calculation of alpha//////////////////////////////
			double mu1=0.;											//
			#pragma omp parallel for reduction(+:mu1)
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					mu1 += vr0[k][j]*vr[k][j];
				}
			}
			#pragma omp barrier

			double gamma=0.;
			#pragma omp parallel for reduction(+:gamma)
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					gamma += vr0[k][j]*vtmp[k][j];
				}
			}
			#pragma omp barrier
			double alpha = mu1/gamma;								//
			////////calculation of alpha//////////////////////////////

			
			//h_{k+1}=e_k-alpha_k* A p_k
			//e_k <- e_k + h_{k+1}
			#pragma omp parallel for
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					vh[k][j] = ve[k][j] -alpha*vtmp[k][j];
					ve[k][j] = ve[k][j] + vh[k][j];
				}
			}
			#pragma omp barrier
			
			///////////////// vtmp <- A e_k //////////////////////////
			//boundary setting if R\neq0							//
			boundary_periodicphi_nosympoletheta(ve);				//

			//M p
			multiple(ve,vtmp);

			//U^{-1}D^{-1}L^{-1} (M p)
			precondition(vtmp);										//
			///////////////// vtmp <- A e_k //////////////////////////
			
			
			//x_{k+1}=x_k + alpha_k e_k  !! e_k <- e_k + h_{k+1} before
			//r_{k+1}=r_k - alpha_k vtmp !! vtmp <- A (e_k + h_{k+1})
			#pragma omp parallel for
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					vf[k][j] = vf[k][j] +alpha*ve[k][j];
					vr[k][j]  = vr[k][j]  -alpha*vtmp[k][j];
				}
			}
			#pragma omp barrier
			
			///////// calculation of beta ////////////////////////////
			double mu2=mu1;											//
			mu1=0.;
			#pragma omp parallel for reduction(+:mu1)
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					mu1 += vr0[k][j]*vr[k][j];
				}
			}
			#pragma omp barrier
			double beta = mu1/mu2;									//
			///////// calculation of beta ////////////////////////////
			
			//e_{k+1}=r_{k+1}+beta_{k+1} h_{k+1}
			//p_{k+1}=e_{k+1}+beta_{k+1} (h_{k+1} + beta_{k+1} p_k)
			#pragma omp parallel for
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					ve[k][j] = vr[k][j] +beta*vh[k][j];
					vp[k][j] = ve[k][j] +beta*(vh[k][j]+beta*vp[k][j]);
				}
			}
			#pragma omp barrier
			
			//calculation of norm
			double norm_r=0.;
			#pragma omp parallel for reduction(+:norm_r)
			for(int k=1;k<(mp-1);k++)
			{
				for(int j=1;j<(mt-1);j++)
				{
					norm_r += vr[k][j]*vr[k][j];
				}
			}
			#pragma omp barrier

			/* cout << "                   Poisson error :: " << norm_src*errmax -norm_r */
			/* 	   << " alpha=" << alpha << " beta=" << beta << endl; */

			if(norm_src*errmax-norm_r>0.) break;
		}
		while(it<itmax);
		
		//final boundary setting
		//boundary_half(vf);
		//boundaryset(vf);
        boundary_periodicphi_nosympoletheta(vf);
		return;
	}
	
	//this is operated to vtpm and vtmp=0 in buffer region
	void precondition(double **vs)
	{
		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				vcon[k][j] = vs[k][j];
			}
		}
		#pragma omp barrier
	
		//------------------------
		//       s => C^{-1} s
		//------------------------
		//   C x = s
		// LDU x = s
		//------------------------
		// L z = s
		//   z = L^{-1} s
		for(int k=1;k<(mp-1);k++)
		{
			int km = k-1;
			for(int j=1;j<(mt-1);j++)
			{
				int jm = j-1;
				vcon[k][j]=cd[k][j]*(vs[k][j]
							-cc[k][j]*vcon[k ][jm]
							-cb[k][j]*vcon[km][j ]);
			}
		}

		//------------------------
		// DU x = z 
		//    x = U^{-1}D^{-1}
		for(int k=mp-2;k>=1;k--)
		{
			int kp = k+1;
			for(int j=mt-2;j>=1;j--)
			{
				int jp = j+1;
				vcon[k][j]=vcon[k][j]
							-cd[k][j]*( ce[k][j]*vcon[k ][jp]
							+cf[k][j]*vcon[kp][j ]);
			}
		}
		
		#pragma omp parallel for
		for(int k=1;k<(mp-1);k++)
		{
			for(int j=1;j<(mt-1);j++)
			{
				vs[k][j] = vcon[k][j];
			}
		}
		#pragma omp barrier

		return;
	}
	
	//reflection boundary cond.
	void boundaryset(double **vs)
	{
		#pragma omp parallel for
		for(int k=0;k<mp;k++)
		{
			vs[k][0] = vs[k][1];
			vs[k][mt-1] = vs[k][mt-2];
		}
		#pragma omp barrier

		#pragma omp parallel for
		for(int j=0;j<mt;j++)
		{
			vs[0][j] = vs[1][j];
			vs[mp-1][j] = vs[mp-2][j];
		}
		#pragma omp barrier

		return;
	}

	//reflection boundary cond. half
	void boundary_half(double **vs)
	{
		for(int k=0;k<mp;k++)
		{
			vs[k][0] = vs[k][1];
			vs[k][mt-1] = vs[k][mt-2];
		}
		for(int j=0;j<mt;j++)
		{
			vs[0][j] = vs[mp-2][j];
			vs[mp-1][j] = vs[1][j];
		}
		
		return;
	}

	//periodic boundary condition for phi
	//the poles can be non-symmetric points
	//mp is supposed to be an even integer
	void boundary_periodicphi_nosympoletheta(double **vs)
	{
		#pragma omp parallel for
		for(int k=1;k<mp/2;k++)
		{
			vs[k][0] = vs[k+mp/2-1][1];
			vs[k][mt-1] = vs[k+mp/2-1][mt-2];
		}
		#pragma omp barrier
		#pragma omp parallel for
		for(int k=mp/2;k<(mp-1);k++)
		{
			vs[k][0] = vs[k-mp/2+1][1];
			vs[k][mt-1] = vs[k-mp/2+1][mt-2];
		}
		#pragma omp barrier
		#pragma omp parallel for
		for(int j=0;j<mt;j++)
		{
			vs[0][j] = vs[mp-2][j];
			vs[mp-1][j] = vs[0][j];
		}
		
		return;
	}
	
	//substituting 0 for boundary
	void boundary0(double **vs)
	{
		for(int k=0;k<mp;k++)
		{
			vs[k][0] = 0.;
			vs[k][mt-1] = 0.;
		}
		for(int j=0;j<mt;j++)
		{
			vs[0][j] = 0.;
			vs[mp-1][j] = 0.;
		}
		return;
	}
	
	//retun true if not found
	bool get_error() const
	{
		return error;
	}

	//set true if not found
	void set_error(bool e){
		error=e;
		return;
	}
	
	//initial horizon radius
	void set_hini(double h)
	{
		for(int k=0;k<mp;k++)
		{
			for(int j=0;j<mt;j++)
			{
				double st = sin(theta[j]);
				double ct = cos(theta[j]);
				
				double a=h*0.7;
				double b=h;
				
				ha[k][j] = a*b/sqrt(b*b*ct*ct+a*a*st*st);
				hb[k][j] = a*b/sqrt(b*b*ct*ct+a*a*st*st);
			}
		}
		return;
	}
	
	void find_ah(Fmv0 *fmv,int loopmax,double err_p,double err_eps,ofstream& fout,short int hsign);
	void print_ah(Fmv0 *fmv,ofstream& fout,double time);
	
	//setting horizon flag
	// 0:outside horizon
	//-1:inside  horizon
	void set_hflag(Fmv0 *fmv)
	{
		int jli=fmv->get_jli();
		int kli=fmv->get_kli();
		int lli=fmv->get_lli();
		int jui=fmv->get_jui();
		int kui=fmv->get_kui();
		int lui=fmv->get_lui();
		
		for(int j=jli;j<=jui;j++)
		{
		  double XX=fmv->funcf(fmv->get_x(j));
			for(int k=kli;k<=kui;k++)
			{
			  double YY=fmv->funcf(fmv->get_y(k));
				double rho=sqrt(XX*XX+YY*YY);
				double ph;
				
				if(abs(XX)<1.0e-5)
				 ph=M_PI*0.5;
				else
				 ph=atan2(YY,XX);
				
				for(int l=lli;l<=lui;l++)
				{

				  double ZZ=fmv->funcf(fmv->get_z(l));
					double th;
					if(abs(ZZ)<1.0e-5)
					 th=M_PI*0.5;
					else
					 th=atan2(rho,abs(ZZ));
					
					double rr=sqrt(rho*rho+ZZ*ZZ);

					
					int tc=(int ((th/dt+0.5)));
					int pc=(int ((ph/dp+0.5)));
					
					if(tc>=mt-2)
					 tc=mt-3;
					if(pc>=mp-2)
					 pc=mp-3;
					
					if(rr>hrad(tc,pc,th,ph,3))
					 fmv->set_hflag(l,k,j)=0;
					else
					 fmv->set_hflag(l,k,j)=-1;
				}
			}
		}
		
		fmv->boundary_quarter_hflags();
		
		return;
	}
	
	//return horizon radius with interpolation
	double hrad(int tc,int pc,double th,double ph,int order)
	{
		double ans=0.;
		double *xx,**kk,*ktc;
		xx=new double [order];
		kk=new double *[order];
		ktc=new double [order];
		int buf=order/2-1;
		
		for(int p=0;p<order;p++)
		{
			kk[p]=new double[order];
			xx[p]=0.;
			ktc[p]=0.;
			int pt=pc-buf+p;
			for(int t=0;t<order;t++)
			{
				int tt=tc-buf+t;
				kk[p][t]=ha[pt][tt];
			}
		}
		
		for(int t=0;t<order;t++)
		{
			xx[t]=theta[tc-buf+t];
		}
		for(int p=0;p<order;p++)
		{
			ktc[p]=ipol(th,xx,kk[p],order);
		}
		
		for(int p=0;p<order;p++)
		{
			xx[p]=phi[pc-buf+p];
		}
		ans =ipol(ph,xx,ktc,order);
		
		for(int p=0;p<order;p++)
		{
			delete[] kk[p];
		}
		delete[] xx;
		delete[] kk;
		delete[] ktc;
		
		return ans;
	}
	
	//interpolation by using Lagrange formula
	double ipol( double rr,double *xx,double *yy,int order )
	{
		//rr:position, xx:grid, yy:variable, order: order of interpolation
		double *kk;
		kk = new double[order];
		double ans=0.;
		for(int i=0;i<order;i++)
		{
			kk[i]=1.;
			for(int n=0;n<order;n++)
			{
				if(n==i) continue;
				kk[i] *= (rr-xx[n])/(xx[i]-xx[n]);
			}
			ans += kk[i]*yy[i];
		}
		delete[] kk;
		return ans;
	}
};


#endif
