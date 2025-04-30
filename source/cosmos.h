#ifndef _COSMOS_H_
#define _COSMOS_H_
/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* HEADER FILE :: BSSN evolution Class of COSMOS                                                         */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E 
#define M_E 2.71828182845904523536
#endif

using namespace std;

class Fmv0{

protected:
	///////////////  general parameters and variables ////////////////////////
	bool fluidevo,scalarevo;									// on/off for fluid, scalar evolution
	bool curveval;												// on/off for curvature evaluation
	bool hform;													// horizon formation
	bool exc;													// excision
	bool mrf;													// mesh refinement flag

	short int *rotsign;											// sign for y=0 boundary rotation
	short int *refsign;											// sign for reflection boundary

	int tab;													// tab for refinement boundary
	int jmax,jmin,kmax,kmin,lmax,lmin;							// maximum/minimum number of grid points
	int jui,jli,kui,kli,lui,lli;								// range to evolve
	int nx,ny,nz,nn,nc,nv;										// number of grid points, variables and outputs 
	int nt;														// number of energy momentum tensor
	int npr;													// number of primitive variables of fluids
	int nsc,nscp;												// number for scalar field and its momentum	
	int jhm,khm,lhm;											// grid point of max constraint violation
	int jmm,kmm,lmm;											// grid point of max constraint violation
	int jkm,kkm,lkm;											// grid point of max Kre. inv.
	int jwm,kwm,lwm;											// grid point of max Weyl inv.
	int ipolo;													// order of interpolation
	int exg;													// excision grid number

	int layn;													// layer number
	int jjmin;													// minimum j for constraint check
	int kkmin;													// minimum k for constraint check
	int llmin;													// minimum l for constraint check
	int jjmax;													// maximum j for constraint check
	int kkmax;													// maximum k for constraint check
	int llmax;													// maximum l for constraint check

	double xu,xl, yu,yl, zu,zl;									// maximum/minimum of coordinate length
	double dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,dvol;			// grid intervals and inverse of those
	double dxi4,dyi4,dzi4;										// dx/12,dx/24 for convenience
	double dxi12,dyi12,dzi12,dxi24,dyi24,dzi24;					// dx/12,dx/24 for convenience
	double t,dt,dt0,dtp,dtpp;									// time intervals dtp:previous dtpp:one before previous
	double tmax,cfl,etaa,etab,etabb;							// CFL condition and gauge parameters
	double KOep;												// Kreiss-Oliger dissipation dissipation term
	double amp;													// amplitude of inhomogeneous grid
	double lambda;												// cosmological constant
	double pi2,pi4,pi8,pi16;									// 4*pi,8*pi,16*pi
	///////////////  general parameters and variables ////////////////////////
	
	///////////////  constraints and curvature invariants ////////////////////
	double ham,hammax,mom,mommax;								// values of constraint violation
	double dGam;												// deviation of Gamma: constraint
	double dGammax;												// deviation of Gamma: constraint
	double Kremax;												// value of max Kre. inv.
	double Weylmax;												// value of max Weyl inv.
	///////////////  constraints and curvature invariants ////////////////////
	
	
	///////////////  fluid parameters ////////////////////////////////////////
	double kap_MUSCL,b_minmod;									// parameters for MUSCL interpolation
	double fluidw;												// parameter for fluid EOS
	///////////////  fluid parameters ////////////////////////////////////////
	
	
	///////////////  scalar field parameters /////////////////////////////////
	double scalarm;												// scalar field mass parameter
	///////////////  scalar field parameters /////////////////////////////////
	
	
	///////////////  parameters for initial data and settings ////////////////
	double Hb;													// initial Hubble
	double tini;												// initial time
	double tk2;													// initial trK^2
	double dtk;													// dtrk/dt for maximal slice
	double boxL;												// boxsize
	///////////////  parameters for initial data and settings ////////////////
	
	///////////////  coordinates and buffers /////////////////////////////////
	int *ju,*jl,*ku,*kl,*lu,*ll;								// from innermost grid to outermost buffer grid
	double *x;													// coordinate x
	double *y;													// coordinate y
	double *z;													// coordinate z
	///////////////  coordinates and buffers /////////////////////////////////
	
	///////////////  boundary flags //////////////////////////////////////////
	int ***bflag;												// flag used for excsion
	int ***hflag;												// horizon flag
	int ***fmrflag;												// flag for the existence of higher layer
	///////////////  boundary flags //////////////////////////////////////////

	///////////////  dynamical variables /////////////////////////////////////
	double ****bv;												// variables (BSSN,GAUGE,fluid,[scalar])
	double ****dbv;												// variables (BSSN,GAUGE,fluid,[scalar])
	double ****bv0;												// variables (previous step)
	double ****bv1;												// variables (two steps before)
	double ****bvr;												// variables (sum for Rungekutta)
	///////////////  dynamical variables /////////////////////////////////////
	
	///////////////  fluid fluxes and primitive vars /////////////////////////
	double ****flux_x;											// deposit for x-flux at i-1/2
	double ****flux_y;											// deposit for y-flux at i-1/2
	double ****flux_z;											// deposit for z-flux at i-1/2
	double ****primv;											// primitive variables of fluid
	///////////////  fluid fluxes and primitive vars /////////////////////////
	
	///////////////  temporary storages //////////////////////////////////////
	double ***psi;												// for psi in initial data setting
	///////////////  temporary storages //////////////////////////////////////
	
	///////////////  geometrical variables for inhomogeneous coordinate //////
//	double ***flat_det;											// flat metric determinant
	double *flat_df2x;											// flat metric variables
	double *flat_df2y;											// flat metric variables
	double *flat_df2z;											// flat metric variables
	double *flat_Gamx;											// bar Gamma_ux_xx
	double *flat_Gamy;											// bar Gamma_uy_yy
	double *flat_Gamz;											// bar Gamma_uz_zz
	double *flat_dGamx;											// bar delG_x_ux_xx
	double *flat_dGamy;											// bar delG_y_uy_yy
	double *flat_dGamz;											// bar delG_z_uz_zz
	///////////////  geometrical variables for inhomogeneous coordinate //////
	
	///////////////  constraints and output storages /////////////////////////
	double ****con;												// constraints, etc...
	double ****outv;											// constraints, etc...
	///////////////  constraints and output storages /////////////////////////
	
	///////////////  supplements /////////////////////////////////////////////
	///////////////  supplements /////////////////////////////////////////////
	

public:
	Fmv0(int tabs,int jupper,int jlower,int kupper,int klower,int lupper,int llower,
	double xupper,double xlower,double yupper,double ylower,double zupper,double zlower,double am,bool fld, bool scl, bool cuev)
	{
		
		///////////////  coordinate settings /////////////////////////////////////////////
		tab=tabs;
		ju = new int[tab+1];
		jl = new int[tab+1];
		ku = new int[tab+1];
		kl = new int[tab+1];
		lu = new int[tab+1];
		ll = new int[tab+1];
		for(int n=0;n<=tab;n++)
		{
			ju[n] = jupper + n;
			jl[n] = jlower - n;
			ku[n] = kupper + n;
			kl[n] = klower - n;
			lu[n] = lupper + n;
			ll[n] = llower - n;
		}
		jmax=ju[tab];
		jmin=jl[tab];
		kmax=ku[tab];
		kmin=kl[tab];
		lmax=lu[tab];
		lmin=ll[tab];

		jui=ju[0];
		jli=jl[0];
		kui=ku[0];
		kli=kl[0];
		lui=lu[0];
		lli=ll[0];

		nx=jmax-jmin+1;
		ny=kmax-kmin+1;
		nz=lmax-lmin+1;

		amp=am;

		xu=xupper;
		xl=xlower;
		yu=yupper;
		yl=ylower;
		zu=zupper;
		zl=zlower;
		///////////////  coordinate settings /////////////////////////////////////////////

		///////////////  equations on/off ////////////////////////////////////////////////
		fluidevo=fld;									// for fluid
		scalarevo=scl;									// for sclar
		curveval=cuev;									// for curvareu evaluation
		///////////////  equations on/off ////////////////////////////////////////////////

		///////////////  number of variables /////////////////////////////////////////////
		nn=24;										// for geometry
		
		if(fld)
		nn+=5;										// for fluid
		
		if(scl)
		{
			nsc=nn;
			nscp=nsc+1;
			nn+=2;									// for scalar
		}
		else
		{
			nsc=nscp=0;
		}

		// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
		//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
		//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
		//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
		//  13:wa , -> psi 
		//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
		//  20:ek  , -> trK
		// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
		// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
		// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25

		npr=5;										// for primitive variables of fluid
		// @@@@@@@@@@@@   primitive variables of fluid   @@@@@@@@@@@
		//  0:rho , 1:V^x  ,  2:V^y  ,  3:V^z   ,  4:epsilon

		nc=11;										// for constraints
		// @@@@@@@@@@@@   constraints   @@@@@@@@@@@
		//  0:normalized ham, 1:ham, 2:normalized momx, 3:momx 
		//  4:normalized momy,5:momy,6:normalized momz, 7:momz
		//  8:Gamma^x-D_gamma^x, 9:Gamma^y-D_gamma^y, 10:Gamma^z-D_gamma^z
		
		nv=5;										// for output variables
		///////////////  number of variables /////////////////////////////////////////////

		///////////////  initial sign setting for boundary conditions  ///////////////////
		refsign = new short int[nn];					//sign flips for reflection
		rotsign = new short int[nn];					//sign flips for rotation compared to reflection

		// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
		//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
		//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
		//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
		//  13:wa , -> psi 
		//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
		//  20:ek  , -> trK
		// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
		// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
		// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25

		for(int i=0;i<nn;i++)
		{
			refsign[i]=1; 
			rotsign[i]=1;
    	}

		//Note: the sign flips are not purely due to rotation, but compared to the reflection sign flips on y=const. surface
		//Namely, only one x index => sign flip
		rotsign[1]=-1;
		rotsign[4]=-1;
		rotsign[10]=-1;
		rotsign[11]=-1;
		rotsign[17]=-1;
		rotsign[18]=-1;
		rotsign[21]=-1;
		
		if(fld)
		{
			rotsign[25]=-1;
		}
		///////////////  sign setting for boundary conditions  ///////////////////////////

		///////////////  parameters and temporary variables //////////////////////////////
		ipolo=6;									// order of interpolation
		
		pi2=2.*M_PI;									// 2pi
		pi4=4.*M_PI;									// 4pi
		pi8=8.*M_PI;									// 8pi
		pi16=16.*M_PI;									//16pi
		
		tk2=1.;											// temporary initial value of trK^2

		t=dt=dt0=dtp=dtpp=tmax=0.;						// temporary initial values wrt time
		cfl=0.;											// temporary initial values for cfl
		
		boxL=1.;										// box size is normalized to unity
		jjmin=jli;
		kkmin=kli;
		llmin=lli;
		jjmax=jui;
		kkmax=kui;
		llmax=lui;
		///////////////  parameters and temporary variables //////////////////////////////

		///////////////  spatial coordinate intervals ////////////////////////////////////
		dx=(xu-xl)/(double(ju[0]-jl[0]));
		dy=(yu-yl)/(double(ku[0]-kl[0]));
		dz=(zu-zl)/(double(lu[0]-ll[0]));

		//////////////////////////////////////////////////
		//boundary positions
		//     0 1 2 3 4 5 6 7 8 9 
		//     B B B I I I I B B B --> inner and buffer
		//           |     |       --> boundary position
		//////////////////////////////////////////////////

		dxi=1./dx;
		dyi=1./dy;
		dzi=1./dz;

		dxi2=0.5*dxi;
		dyi2=0.5*dyi;
		dzi2=0.5*dzi;

		dxi4=0.5*dxi2;
		dyi4=0.5*dyi2;
		dzi4=0.5*dzi2;
		dxi12=dxi/12.;
		dyi12=dyi/12.;
		dzi12=dzi/12.;
		dxi24=dxi/24.;
		dyi24=dyi/24.;
		dzi24=dzi/24.;

		dvol=dx*dy*dz;
		///////////////  spatial coordinate intervals ////////////////////////////////////

		///////////////  coordinate strages //////////////////////////////////////////////
		x = new double[nx];																//
		for(int j=jmin;j<=jmax;j++){
			x[j-jmin]=xl +dx*(double(j-jl[0]));
		}
		y = new double[ny];
		for(int k=kmin;k<=kmax;k++){
			y[k-kmin]=yl +dy*(double(k-kl[0]));
		}
		z = new double[nz];
		for(int l=lmin;l<=lmax;l++){
			z[l-lmin]=zl +dz*(double(l-ll[0]));
		}																				//
		///////////////  coordinates strages /////////////////////////////////////////////
		
		///////////////  variable strages ////////////////////////////////////////////////
		bv = new double***[nn];															//
		dbv = new double***[nn];
		bv0 = new double***[nn];
		bv1 = new double***[nn];
		bvr = new double***[nn];
		
		for(int l=0;l<nn;l++){
			bv[l] = new double**[nz];
			dbv[l] = new double**[nz];
			bv0[l] = new double**[nz];
			bv1[l] = new double**[nz];
			bvr[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				bv[l][k] = new double*[ny];
				dbv[l][k] = new double*[ny];
				bv0[l][k] = new double*[ny];
				bv1[l][k] = new double*[ny];
				bvr[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					bv[l][k][j]= new double[nx];
					dbv[l][k][j]= new double[nx];
					bv0[l][k][j]= new double[nx];
					bv1[l][k][j]= new double[nx];
					bvr[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  variable strages ////////////////////////////////////////////////

		///////////////  boundary flags //////////////////////////////////////////////////
		bflag = new int**[nz];															//
		hflag = new int**[nz];
		fmrflag = new int**[nz];
		for(int l=0;l<nz;l++){
			bflag[l] = new int*[ny];
			hflag[l] = new int*[ny];
			fmrflag[l] = new int*[ny];
			for(int k=0;k<ny;k++){
				bflag[l][k] = new int[nx];
				hflag[l][k] = new int[nx];
				fmrflag[l][k] = new int[nx];
			}
		}																				//
		///////////////  boundary flags //////////////////////////////////////////////////

		///////////////  temporary variables for elliptic solvers ////////////////////////
		psi = new double**[nz];															//
		for(int l=0;l<nz;l++)
		{
			psi[l] = new double*[ny];
			for(int k=0;k<ny;k++)
			{
				psi[l][k]= new double[nx];
			}
		}
		///////////////  temporary variables for elliptic solvers ////////////////////////

		///////////////  fluid primitive and fluxes //////////////////////////////////////
		primv = new double***[npr];														//
		for(int l=0;l<npr;l++){
			primv[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				primv[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					primv[l][k][j]= new double[nx];
				}
			}
		}
		
		flux_x = new double***[npr];
		flux_y = new double***[npr];
		flux_z = new double***[npr];
		for(int l=0;l<npr;l++){
			flux_x[l] = new double**[nz];
			flux_y[l] = new double**[nz];
			flux_z[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				flux_x[l][k] = new double*[ny];
				flux_y[l][k] = new double*[ny];
				flux_z[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					flux_x[l][k][j]= new double[nx];
					flux_y[l][k][j]= new double[nx];
					flux_z[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  fluid primitive and fluxes //////////////////////////////////////

		///////////////  constraints and output variables ////////////////////////////////
		con = new double***[nc];														//
		for(int l=0;l<nc;l++){
			con[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				con[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					con[l][k][j]= new double[nx];
				}
			}
		}

		outv = new double***[nv];
		for(int l=0;l<nv;l++){
			outv[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				outv[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					outv[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  constraints and output variables ////////////////////////////////


		///////////////  geometrical variables for inhomogeneous coordinate //////////////
		flat_df2x = new double[nx];														//
		flat_df2y = new double[ny];
		flat_df2z = new double[nz];
		
		flat_Gamx = new double[nx];
		flat_Gamy = new double[ny];
		flat_Gamz = new double[nz];
		
		flat_dGamx = new double[nx];
		flat_dGamy = new double[ny];
		flat_dGamz = new double[nz];													//
		///////////////  geometrical variables for inhomogeneous coordinate //////////////
		
		///////////////  supplements /////////////////////////////////////////////////////
																						//
		///////////////  supplements /////////////////////////////////////////////////////
		
		set_zero_all();

		for(int l=lmin;l<=lmax;l++){
			for(int k=kmin;k<=kmax;k++){
				for(int j=jmin;j<=jmax;j++){
						set_bflag(l,k,j)=0;
						set_hflag(l,k,j)=0;
						set_fmrflag(l,k,j)=0;
				}
			}
		}

		cout.setf(ios_base::fixed, ios_base::floatfield);
		cout.precision(3);
		cout << "          xl : xu  / jli : jui      " << setw(9) << right << xl << " : " << setw(9) << right << xu
		<< "   /   " << setw(4) << right << jli << " : " << setw(4) << right << jui << endl
		<< "          yl : yu  / kli : kui      " << setw(9) << right << yl << " : " << setw(9) << right << yu
		<< "   /   " << setw(4) << right << kli << " : " << setw(4) << right << kui << endl
		<< "          zl : zu  / lli : lui      " << setw(9) << right << zl << " : " << setw(9) << right << zu
		<< "   /   " << setw(4) << right << lli << " : " << setw(4) << right << lui << endl
		<< endl;
		cout.setf(ios_base::scientific, ios_base::floatfield);
		cout.precision(10);
	}


	virtual ~Fmv0(){

		delete[] ju;
		delete[] jl;
		delete[] ku;
		delete[] kl;
		delete[] lu;
		delete[] ll;

		delete[] refsign;
		delete[] rotsign;

		delete[] x;
		delete[] y;
		delete[] z;
		
		delete[] flat_df2x;
		delete[] flat_df2y;
		delete[] flat_df2z;
		
		delete[] flat_Gamx;
		delete[] flat_Gamy;
		delete[] flat_Gamz;
		
		delete[] flat_dGamx;
		delete[] flat_dGamy;
		delete[] flat_dGamz;
		
		for(int l=0;l<nn;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] bv[l][k][j];
					bv[l][k][j]=NULL;
					delete[] dbv[l][k][j];
					dbv[l][k][j]=NULL;
					delete[] bv0[l][k][j];
					bv0[l][k][j]=NULL;
					delete[] bv1[l][k][j];
					bv1[l][k][j]=NULL;
					delete[] bvr[l][k][j];
					bvr[l][k][j]=NULL;
				}
				delete[] bv[l][k];
				bv[l][k]=NULL;
				delete[] dbv[l][k];
				dbv[l][k]=NULL;
				delete[] bv0[l][k];
				bv0[l][k]=NULL;
				delete[] bv1[l][k];
				bv1[l][k]=NULL;
				delete[] bvr[l][k];
				bvr[l][k]=NULL;
			}
			delete[] bv[l];
			bv[l]=NULL;
			delete[] dbv[l];
			dbv[l]=NULL;
			delete[] bv0[l];
			bv0[l]=NULL;
			delete[] bv1[l];
			bv1[l]=NULL;
			delete[] bvr[l];
			bvr[l]=NULL;
		}
		delete[] bv;
		bv=NULL;
		delete[] dbv;
		dbv=NULL;
		delete[] bv0;
		bv0=NULL;
		delete[] bv1;
		bv1=NULL;
		delete[] bvr;
		bvr=NULL;
		
		for(int l=0;l<nc;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] con[l][k][j];
					con[l][k][j]=NULL;
				}
				delete[] con[l][k];
				con[l][k]=NULL;
			}
			delete[] con[l];
			con[l]=NULL;
		}
		delete[] con;
		con=NULL;
		
		for(int l=0;l<nv;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] outv[l][k][j];
					outv[l][k][j]=NULL;
				}
				delete[] outv[l][k];
				outv[l][k]=NULL;
			}
			delete[] outv[l];
			outv[l]=NULL;
		}
		delete[] outv;
		outv=NULL;
		
		//maxslice mod
		for(int k=0;k<nz;k++){
			for(int j=0;j<ny;j++){
				delete[] psi[k][j];
				psi[k][j]=NULL;
			}
			delete[] psi[k];
			psi[k]=NULL;
		}
		delete[] psi;
		psi=NULL;
		
		for(int l=0;l<npr;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] primv[l][k][j];
					primv[l][k][j]=NULL;
				}
				delete[] primv[l][k];
				primv[l][k]=NULL;
			}
			delete[] primv[l];
			primv[l]=NULL;
		}
		delete[] primv;
		primv=NULL;

		for(int l=0;l<npr;l++){
			for(int k=0;k<nz;k++){
				delete[] flux_x[l][k];
				flux_x[l][k]=NULL;
			}
			delete[] flux_x[l];
			flux_x[l]=NULL;
		}
		delete[] flux_x;
		flux_x=NULL;
		
		for(int l=0;l<npr;l++){
			for(int k=0;k<nz;k++){
				delete[] flux_y[l][k];
				flux_y[l][k]=NULL;
			}
			delete[] flux_y[l];
			flux_y[l]=NULL;
		}
		delete[] flux_y;
		flux_y=NULL;

		for(int l=0;l<npr;l++){
			for(int k=0;k<ny;k++){
				delete[] flux_z[l][k];
				flux_z[l][k]=NULL;
			}
			delete[] flux_z[l];
			flux_z[l]=NULL;
		}
		delete[] flux_z;
		flux_z=NULL;
		
	}

	public:

	////////////////////////////////////////
	//  GET func.
	////////////////////////////////////////
	bool get_fluidevo() const{
		return fluidevo;
	}
	bool get_scalarevo() const{
		return scalarevo;
	}
	bool get_hform() const{
		return hform;
	}
	bool get_exc() const{
		return exc;
	}
	bool get_mrf() const{
		return mrf;
	}
	int get_jmax() const{
		return jmax;
	}
	int get_jmin() const{
		return jmin;
	}
	int get_kmax() const{
		return kmax;
	}
	int get_kmin() const{
		return kmin;
	}
	int get_lmax() const{
		return lmax;
	}
	int get_lmin() const{
		return lmin;
	}
	int get_jui() const{
		return jui;
	}
	int get_jli() const{
		return jli;
	}
	int get_kui() const{
		return kui;
	}
	int get_kli() const{
		return kli;
	}
	int get_lui() const{
		return lui;
	}
	int get_lli() const{
		return lli;
	}
	int get_nx() const{
		return nx;
	}
	int get_ny() const{
		return ny;
	}
	int get_nz() const{
		return nz;
	}
	int get_jhm() const{
		return jhm;
	}
	int get_khm() const{
		return khm;
	}
	int get_lhm() const{
		return lhm;
	}
	int get_jmm() const{
		return jmm;
	}
	int get_kmm() const{
		return kmm;
	}
	int get_lmm() const{
		return lmm;
	}
	int get_jkm() const{
		return jkm;
	}
	int get_kkm() const{
		return kkm;
	}
	int get_lkm() const{
		return lkm;
	}
	int get_jwm() const{
		return jwm;
	}
	int get_kwm() const{
		return kwm;
	}
	int get_lwm() const{
		return lwm;
	}
	int get_exg() const{
		return exg;
	}
	int get_layn() const{
		return layn;
	}
	int get_jjmin() const{
		return jjmin;
	}
	int get_kkmin() const{
		return kkmin;
	}
	int get_llmin() const{
		return llmin;
	}
	
	double get_t() const{
		return t;
	}
	double get_dt() const{
		return dt;
	}
	double get_dt0() const{
		return dt0;
	}
	double get_dtp() const{
		return dtp;
	}
	double get_dtpp() const{
		return dtpp;
	}
	double get_tmax() const{
		return tmax;
	}
	double get_cfl() const{
		return cfl;
	}
	double get_dx() const{
		return dx;
	}
	double get_dy() const{
		return dy;
	}
	double get_dz() const{
		return dz;
	}
	double get_dvol() const{
		return dvol;
	}
	double get_dxi() const{
		return dxi;
	}
	double get_dyi() const{
		return dyi;
	}
	double get_dzi() const{
		return dzi;
	}
	double get_dxi2() const{
		return dxi2;
	}
	double get_dyi2() const{
		return dyi2;
	}
	double get_dzi2() const{
		return dzi2;
	}

	double get_dxi4() const{
		return dxi4;
	}
	double get_dyi4() const{
		return dyi4;
	}
	double get_dzi4() const{
		return dzi4;
	}

	double get_xu() const{
		return xu;
	}
	double get_xl() const{
		return xl;
	}
	double get_yu() const{
		return yu;
	}
	double get_yl() const{
		return yl;
	}
	double get_zu() const{
		return zu;
	}
	double get_zl() const{
		return zl;
	}
	double get_ham() const{
		return ham;
	}
	double get_hammax() const{
		return hammax;
	}
	double get_Kremax() const{
		return Kremax;
	}
	double get_Weylmax() const{
		return Weylmax;
	}
	double get_mom() const{
		return mom;
	}
	double get_mommax() const{
		return mommax;
	}
	double get_etaa() const{
		return etaa;
	}
	double get_etab() const{
		return etab;
	}
	double get_etabb() const{
		return etabb;
	}
	double get_lambda() const{
		return lambda;
	}
	double get_tini() const{
		return tini;
	}
	double get_KOep() const{
		return KOep;
	}
	double get_Mkap() const{
		return kap_MUSCL;
	}
	double get_b() const{
		return b_minmod;
	}
	double get_Hb() const{
		return Hb;
	}
	double get_fluidw() const{
		return fluidw;
	}
	double get_scalarm() const{
		return scalarm;
	}

	int get_bflag(int l,int k,int j) const{
		return bflag[l-lmin][k-kmin][j-jmin];
	}
	int get_hflag(int l,int k,int j) const{
		return hflag[l-lmin][k-kmin][j-jmin];
	}
	int get_fmrflag(int l,int k,int j) const{
		return fmrflag[l-lmin][k-kmin][j-jmin];
	}

	double get_x(int j)const{
		return x[j-jmin];
	}
	double get_y(int k)const{
		return y[k-kmin];
	}
	double get_z(int l)const{
		return z[l-lmin];
	}
	double get_ext_x(int j)const{
		return(xl +dx*(double(j-jl[0])));
	}
	double get_ext_y(int k)const{
		return(yl +dy*(double(k-kl[0])));
	}
	double get_ext_z(int l)const{
		return(zl +dz*(double(l-ll[0])));
	}	
	double get_bv(int l,int k,int j,int i) const{
		return bv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_dbv(int l,int k,int j,int i) const{
		return dbv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_bv0(int l,int k,int j,int i) const{
		return bv0[i][l-lmin][k-kmin][j-jmin];
	}
	double get_bv1(int l,int k,int j,int i) const{
		return bv1[i][l-lmin][k-kmin][j-jmin];
	}
	double get_bvr(int l,int k,int j,int i) const{
		return bvr[i][l-lmin][k-kmin][j-jmin];
	}
	double get_con(int l,int k,int j,int i) const{
		return con[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flat_df2x(int j)const{
		return flat_df2x[j-jmin];
	}
	double get_flat_df2y(int k)const{
		return flat_df2y[k-kmin];
	}
	double get_flat_df2z(int l)const{
		return flat_df2z[l-lmin];
	}
	double get_flat_Gamx(int j)const{
		return flat_Gamx[j-jmin];
	}
	double get_flat_Gamy(int k)const{
		return flat_Gamy[k-kmin];
	}
	double get_flat_Gamz(int l)const{
		return flat_Gamz[l-lmin];
	}
	double get_flat_dGamx(int j)const{
		return flat_dGamx[j-jmin];
	}
	double get_flat_dGamy(int k)const{
		return flat_dGamy[k-kmin];
	}
	double get_flat_dGamz(int l)const{
		return flat_dGamz[l-lmin];
	}
	double get_outv(int l,int k,int j,int i) const{
		return outv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_psi(int l,int k,int j) const{
		return psi[l-lmin][k-kmin][j-jmin];
	}
	double get_primv(int l,int k,int j,int i) const{
		return primv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flux_x(int l,int k,int j,int i) const{
		return flux_x[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flux_y(int l,int k,int j,int i) const{
		return flux_y[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flux_z(int l,int k,int j,int i) const{
		return flux_z[i][l-lmin][k-kmin][j-jmin];
	}
	double*** get_bv(int i) const{
		return bv[i];
	}
	double*** get_dbv(int i) const{
		return dbv[i];
	}
	double*** get_bv0(int i) const{
		return bv0[i];
	}
	double*** get_bv1(int i) const{
		return bv1[i];
	}
	double*** get_bvr(int i) const{
		return bvr[i];
	}

	////////////////////////////////////////
	//  GET derivative func.
	////////////////////////////////////////
	double get_f_x(int l,int k,int j,int i) const{
		double w=(   -(get_bv(l,k,j+2,i)-get_bv(l,k,j-2,i))
			+8.*(get_bv(l,k,j+1,i)-get_bv(l,k,j-1,i))
			)*dxi12;
		return w;
	}
	double get_f_y(int l,int k,int j,int i) const{
		double w=(   -(get_bv(l,k+2,j,i)-get_bv(l,k-2,j,i))
			+8.*(get_bv(l,k+1,j,i)-get_bv(l,k-1,j,i))
			)*dyi12;
		return w;
	}
	double get_f_z(int l,int k,int j,int i) const{
		double w=(   -(get_bv(l+2,k,j,i)-get_bv(l-2,k,j,i))
			+8.*(get_bv(l+1,k,j,i)-get_bv(l-1,k,j,i))
			)*dzi12;
		return w;
	}
	double get_f_xx(int l,int k,int j,int i) const{
		double w=(    -(     get_bv(l  ,k  ,j+2,i) + get_bv(l  ,k  ,j-2,i))
			+16.*( get_bv(l  ,k  ,j+1,i) + get_bv(l  ,k  ,j-1,i))
			-30.*  get_bv(l  ,k  ,j  ,i)
			)*dxi*dxi12;
		return w;
	}
	double get_f_yy(int l,int k,int j,int i) const{
		double w=(    -(     get_bv(l  ,k+2,j  ,i) + get_bv(l  ,k-2,j  ,i))
			+16.*( get_bv(l  ,k+1,j  ,i) + get_bv(l  ,k-1,j  ,i))
			-30.*  get_bv(l  ,k  ,j  ,i)
			)*dyi*dyi12;
		return w;
	}
	double get_f_zz(int l,int k,int j,int i) const{
		double w=(    -(     get_bv(l+2,k  ,j  ,i) + get_bv(l-2,k  ,j  ,i))
			+16.*( get_bv(l+1,k  ,j  ,i) + get_bv(l-1,k  ,j  ,i))
			-30.*  get_bv(l  ,k  ,j  ,i)
			)*dzi*dzi12;
		return w;
	}

	double get_f_xy(int l,int k,int j,int i) const{
		double w=( -( ( -(   get_bv(l  ,k+2,j+2,i) - get_bv(l  ,k+2,j-2,i))
						+8.*(get_bv(l  ,k+2,j+1,i) - get_bv(l  ,k+2,j-1,i)))
						-(   -(get_bv(l  ,k-2,j+2,i) - get_bv(l  ,k-2,j-2,i))
						+8.*(get_bv(l  ,k-2,j+1,i) - get_bv(l  ,k-2,j-1,i))))
						+8.*(  (-(get_bv(l  ,k+1,j+2,i) - get_bv(l  ,k+1,j-2,i))
						+8.*(get_bv(l  ,k+1,j+1,i) - get_bv(l  ,k+1,j-1,i)))
						-(-(get_bv(l  ,k-1,j+2,i) - get_bv(l  ,k-1,j-2,i))
						+8.*(get_bv(l  ,k-1,j+1,i) - get_bv(l  ,k-1,j-1,i))))
						)*dxi12*dyi12;
		return w;
	}
	double get_f_xz(int l,int k,int j,int i) const{
		double w=( -( ( -(   get_bv(l+2,k  ,j+2,i) - get_bv(l+2,k  ,j-2,i))
						+8.*(get_bv(l+2,k  ,j+1,i) - get_bv(l+2,k  ,j-1,i)))
						-(   -(get_bv(l-2,k  ,j+2,i) - get_bv(l-2,k  ,j-2,i))
						+8.*(get_bv(l-2,k  ,j+1,i) - get_bv(l-2,k  ,j-1,i))))
						+8.*(  (-(get_bv(l+1,k  ,j+2,i) - get_bv(l+1,k  ,j-2,i))
						+8.*(get_bv(l+1,k  ,j+1,i) - get_bv(l+1,k  ,j-1,i)))
						-(-(get_bv(l-1,k  ,j+2,i) - get_bv(l-1,k  ,j-2,i))
						+8.*(get_bv(l-1,k  ,j+1,i) - get_bv(l-1,k  ,j-1,i))))
						)*dxi12*dzi12;
		return w;
	}
	double get_f_yz(int l,int k,int j,int i) const{
		double w=( -( ( -(   get_bv(l+2,k+2,j  ,i) - get_bv(l+2,k-2,j  ,i))
						+8.*(get_bv(l+2,k+1,j  ,i) - get_bv(l+2,k-1,j  ,i)))
						-(   -(get_bv(l-2,k+2,j  ,i) - get_bv(l-2,k-2,j  ,i))
						+8.*(get_bv(l-2,k+1,j  ,i) - get_bv(l-2,k-1,j  ,i))))
						+8.*(  (-(get_bv(l+1,k+2,j  ,i) - get_bv(l+1,k-2,j  ,i))
						+8.*(get_bv(l+1,k+1,j  ,i) - get_bv(l+1,k-1,j  ,i)))
						-(-(get_bv(l-1,k+2,j  ,i) - get_bv(l-1,k-2,j  ,i))
						+8.*(get_bv(l-1,k+1,j  ,i) - get_bv(l-1,k-1,j  ,i))))
						)*dyi12*dzi12;
		return w;
	}
	double get_ipol_x_lower_mid(int l,int k,int j,int i) const{
		double w=0.0625*(-get_bv(l,k,j-2,i)+9.*get_bv(l,k,j-1,i)+9.*get_bv(l,k,j,i)-get_bv(l,k,j+1,i));
		return w;
	}
	double get_ipol_y_lower_mid(int l,int k,int j,int i) const{
		double w=0.0625*(-get_bv(l,k-2,j,i)+9.*get_bv(l,k-1,j,i)+9.*get_bv(l,k,j,i)-get_bv(l,k+1,j,i));
		return w;
	}
	double get_ipol_z_lower_mid(int l,int k,int j,int i) const{
		double w=0.0625*(-get_bv(l-2,k,j,i)+9.*get_bv(l-1,k,j,i)+9.*get_bv(l,k,j,i)-get_bv(l+1,k,j,i));
		return w;
	}

	double get_2ndf_x(int l,int k,int j,int i) const{
		double w=(   
			get_bv(l,k,j+1,i)-get_bv(l,k,j-1,i)
			)*dxi2;
		return w;
	}
	double get_2ndf_y(int l,int k,int j,int i) const{
		double w=(
			get_bv(l,k+1,j,i)-get_bv(l,k-1,j,i)
			)*dyi2;
		return w;
	}
	double get_2ndf_z(int l,int k,int j,int i) const{
		double w=(get_bv(l+1,k,j,i)-get_bv(l-1,k,j,i)
			)*dzi2;
		return w;
	}
	double get_2ndf_xx(int l,int k,int j,int i) const{
		double w=(
			get_bv(l  ,k  ,j+1,i) + get_bv(l  ,k  ,j-1,i)
			-2.*  get_bv(l  ,k  ,j  ,i)
			)*dxi*dxi;
		return w;
	}
	double get_2ndf_yy(int l,int k,int j,int i) const{
		double w=(
			get_bv(l  ,k+1,j  ,i) + get_bv(l  ,k-1,j  ,i)
			-2.*  get_bv(l  ,k  ,j  ,i)
			)*dyi*dyi;
		return w;
	}
	double get_2ndf_zz(int l,int k,int j,int i) const{
		double w=(
			get_bv(l+1,k  ,j  ,i) + get_bv(l-1,k  ,j  ,i)
			-2.*  get_bv(l  ,k  ,j  ,i)
			)*dzi*dzi;
		return w;
	}

	double get_2ndf_xy(int l,int k,int j,int i) const{
		double w=(get_bv(l  ,k+1,j+1,i) - get_bv(l  ,k+1,j-1,i)
					-get_bv(l  ,k-1,j+1,i) + get_bv(l  ,k-1,j-1,i)
						)*dxi2*dyi2;
		return w;
	}
	double get_2ndf_xz(int l,int k,int j,int i) const{
		double w=(get_bv(l+1,k  ,j+1,i) - get_bv(l+1,k  ,j-1,i)
					-get_bv(l-1,k  ,j+1,i) + get_bv(l-1,k  ,j-1,i)
						)*dxi2*dzi2;
		return w;
	}
	double get_2ndf_yz(int l,int k,int j,int i) const{
		double w=(get_bv(l+1,k+1,j  ,i) - get_bv(l+1,k-1,j  ,i)
					-get_bv(l-1,k+1,j  ,i) + get_bv(l-1,k-1,j  ,i)
						)*dyi2*dzi2;
		return w;
	}

		
	////////////////////////////////////////
	//  SET func.
	////////////////////////////////////////
	void set_fluidevo(bool f){
		fluidevo=f;
		return;
	}
	void set_scalarevo(bool s){
		scalarevo=s;
		return;
	}
	void set_t(double time){
		t=time;
		return;
	}
	void set_dt(double time){
		dt=time;
		return;
	}
	void set_dt0(double time){
		dt0=time;
		return;
	}
	void set_dtp(double time){
		dtp=time;
		return;
	}
	void set_dtpp(double time){
		dtpp=time;
		return;
	}
	void set_tmax(double time){
		tmax=time;
		return;
	}
	void set_cfl(double c){
		cfl=c;
		return;
	}
	void set_etaa(double e){
		etaa=e;
		return;
	}
	void set_etab(double e){
		etab=e;
		return;
	}
	void set_etabb(double e){
		etabb=e;
		return;
	}
	void set_lambda(double l){
		lambda=l;
		return;
	}
	void set_Hb(double hb){
		Hb=hb;
		return;
	}
	void set_tini(double t){
		tini=t;
		return;
	}
	void set_KOep(double l){
		KOep=l;
		return;
	}
	void set_exg(int eg){
		exg=eg;
		return;
	}
	void set_amp(double a){
		amp=a;
		return;
	}
	void set_fluidw(double fw){
		fluidw=fw;
		return;
	}
	void set_scalarm(double sm){
		scalarm=sm;
		return;
	}
	void set_Mkap(double k){
		kap_MUSCL=k;
		return;
	}
	void set_b(double b){
		b_minmod=b;
		return;
	}
	void set_exc(bool e){
		exc=e;
		return;
	}
	void set_mrf(bool s){
		mrf=s;
		return;
	}
	void set_jjmin(int jjj){
		jjmin=jjj;
		return;
	}
	void set_kkmin(int kkk){
		kkmin=kkk;
		return;
	}
	void set_llmin(int lll){
		llmin=lll;
		return;
	}
	int& set_bflag(int l,int k,int j){
		return bflag[l-lmin][k-kmin][j-jmin];
	}
	int& set_hflag(int l,int k,int j){
		return hflag[l-lmin][k-kmin][j-jmin];
	}
	int& set_fmrflag(int l,int k,int j){
		return fmrflag[l-lmin][k-kmin][j-jmin];
	}
	void set_fmrregion(int lmi,int lma,int kmi,int kma,int jmi,int jma)
	{
		for(int l=lmi;l<=lma;l+=1){
			for(int k=kmi;k<=kma;k+=1){
				for(int j=jmi;j<=jma;j+=1){
					set_fmrflag(l,k,j)=1;
				}
			}
		}
	}

	double& set_bv(int l,int k,int j,int i){
		return bv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_dbv(int l,int k,int j,int i){
		return dbv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_bv0(int l,int k,int j,int i){
		return bv0[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_bv1(int l,int k,int j,int i){
		return bv1[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_bvr(int l,int k,int j,int i){
		return bvr[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_con(int l,int k,int j,int i){
		return con[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flat_df2x(int j){
		return flat_df2x[j-jmin];
	}
	double& set_flat_df2y(int k){
		return flat_df2y[k-kmin];
	}
	double& set_flat_df2z(int l){
		return flat_df2z[l-lmin];
	}
	double& set_flat_Gamx(int j){
		return flat_Gamx[j-jmin];
	}
	double& set_flat_Gamy(int k){
		return flat_Gamy[k-kmin];
	}
	double& set_flat_Gamz(int l){
		return flat_Gamz[l-lmin];
	}
	double& set_flat_dGamx(int j){
		return flat_dGamx[j-jmin];
	}
	double& set_flat_dGamy(int k){
		return flat_dGamy[k-kmin];
	}
	double& set_flat_dGamz(int l){
		return flat_dGamz[l-lmin];
	}
	double& set_outv(int l,int k,int j,int i){
		return outv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_primv(int l,int k,int j,int i){
		return primv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flux_x(int l,int k,int j,int i){
		return flux_x[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flux_y(int l,int k,int j,int i){
		return flux_y[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flux_z(int l,int k,int j,int i){
		return flux_z[i][l-lmin][k-kmin][j-jmin];
	}
	double*** set_bv(int i) const{
		return bv[i];
	}
	double*** set_dbv(int i) const{
		return dbv[i];
	}
	double*** set_bv0(int i) const{
		return bv0[i];
	}
	double*** set_bv1(int i) const{
		return bv1[i];
	}
	double*** set_bvr(int i) const{
		return bvr[i];
	}
	double& set_psi(int l,int k,int j){
		return psi[l-lmin][k-kmin][j-jmin];
	}

	////////////////////////////////////////
	//  SET zero or unity func.
	////////////////////////////////////////
	void set_zero_all()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bv[i][l][k][j] = 0.;
						dbv[i][l][k][j] = 0.;
						bv0[i][l][k][j] = 0.;
						bv1[i][l][k][j] = 0.;
					}
				}
			}
		}
	
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nc;i++){
					for(int k=0;k<ny;k++){
						con[i][l][k][j] = 0.;
					}
				}
			}
		}

		#pragma omp parallel for
		for(int j=0;j<nx;j++){
			flat_df2x[j] = 0.;
		}
		#pragma omp parallel for 
		for(int k=0;k<ny;k++){
			flat_df2y[k] = 0.;
		}
		#pragma omp parallel for
		for(int l=0;l<nz;l++){
			flat_df2z[l] = 0.;
		}

		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			flat_Gamx[j] = 0.;
		}
		#pragma omp parallel for 
		for(int k=0;k<ny;k++){
			flat_Gamy[k] = 0.;
		}
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			flat_Gamz[l] = 0.;
		}

		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			flat_dGamx[j] = 0.;
		}
		#pragma omp parallel for 
		for(int k=0;k<ny;k++){
			flat_dGamy[k] = 0.;
		}
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			flat_dGamz[l] = 0.;
		}
	
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nv;i++){
					for(int k=0;k<ny;k++){
						outv[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					psi[l][k][j] = 0.;
				}
			}
		}

		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<npr;i++){
					for(int k=0;k<ny;k++){
						primv[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<npr;i++){
					for(int k=0;k<ny;k++){
						flux_x[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<npr;i++){
					for(int k=0;k<ny;k++){
						flux_y[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<npr;i++){
					for(int k=0;k<ny;k++){
						flux_z[i][l][k][j] = 0.;
					}
				}
			}
		}

		return;
	}
	

	void set_zero_all_exc()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						if(bflag[l][k][j]!=-1)
						continue;
						
						//bv[i][l][k][j] = 0.;
						dbv[i][l][k][j] = 0.;
						//bv0[i][l][k][j] = 0.;
						//bv1[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nc;i++){
					for(int k=0;k<ny;k++){
						if(bflag[l][k][j]!=-1)
						continue;
						
						con[i][l][k][j] = 0.;
					}
				}
			}
		}

		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nv;i++){
					for(int k=0;k<ny;k++){
						if(bflag[l][k][j]!=-1)
						continue;
						
						outv[i][l][k][j] = 0.;
					}
				}
			}
		}

		return;
	}

	void set_zero_primv()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<npr;i++){
					for(int k=0;k<ny;k++){
						primv[i][l][k][j] = 0.;
					}
				}
			}
		}
	}
	
	void set_zero()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bv[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}
	
	void set_zero_0()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bv0[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}
	
	void set_zero_1()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int i=0;i<nn;i++){
				for(int l=0;l<nz;l++){
					for(int k=0;k<ny;k++){
						bv1[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}

	void set_zero_d()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						dbv[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}

	void set_zero_r()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bvr[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}
	
	void set_bflag_zero(){
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					bflag[l][k][j] = 0;
				}
			}
		}
		return;
	}

	void set_hflag_zero()
	{	
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					hflag[l][k][j] = 0;
				}
			}
		}
		return;
	}
	void set_fmrflag_zero()
	{	
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					fmrflag[l][k][j] = 0;
				}
			}
		}
		return;
	}

	////////////////////////////////////////
	//  UPDATE func.
	////////////////////////////////////////

	void setv0()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bv0[i][l][k][j] = bv[i][l][k][j];
					}
				}
			}
		}

		return;
	}
	
	void set01()
	{
		#pragma omp parallel for
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bv1[i][l][k][j] = bv0[i][l][k][j];
					}
				}
			}
		}
	
		return;
	}
	
	void set0v(){
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int i=0;i<nn;i++){
					for(int k=0;k<ny;k++){
						bv[i][l][k][j] = bv0[i][l][k][j];
					}
				}
			}
		}
		return;
	}

	void runge_kutta(double dt)
	{
		#pragma omp parallel for 
		for(int j=jli;j<=jui;j++){
			for(int l=lli;l<=lui;l++){
				for(int i=0;i<nn;i++){
					for(int k=kli;k<=kui;k++){
						set_bvr(l,k,j,i)=get_bvr(l,k,j,i)+get_dbv(l,k,j,i)*dt;
					}
				}
			}
		}
		return;
	}
	
	void new_bv(double dt)
	{
		#pragma omp parallel for 
		for(int j=jli;j<=jui;j++){
			for(int l=lli;l<=lui;l++){
				for(int i=0;i<nn;i++){
					for(int k=kli;k<=kui;k++){
						set_bv(l,k,j,i)=get_bv0(l,k,j,i)+get_dbv(l,k,j,i)*dt;
					}
				}
			}
		}
		return;
	}

	void new_bv4()
	{
		#pragma omp parallel for 
		for(int j=jli;j<=jui;j++){
			for(int l=lli;l<=lui;l++){
				for(int i=0;i<nn;i++){
					for(int k=kli;k<=kui;k++){
						set_bv(l,k,j,i)=get_bv0(l,k,j,i)+get_bvr(l,k,j,i);
					}
				}
			}
		}
		return;
	}
	
	void set_psi()
	{
		#pragma omp parallel for 
		for(int j=0;j<nx;j++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					psi[l][k][j] = bv[13][l][k][j];
				}
			}
		}

		return;
	}

	////////////////////////////////////////
	//  ERROR evaluation
	////////////////////////////////////////

	void check_const()
	//checking constraint start
	{
		ham=0.;
		double hamtmp=0.;
		double momtmp=0.;
		hammax=0.;
		mommax=0.;
		mom=0.;
		lhm=-1;
		khm=-1;
		jhm=-1;
		lmm=-1;
		kmm=-1;
		jmm=-1;
		
		dGam=0.;
		double dGamtmp=0.;
		dGammax=0.;
		
		int excnum=0;

		#pragma omp parallel 
		{
			double hammaxp=0.;
			double mommaxp=0.;
			double lhmp=-1;
			double khmp=-1;
			double jhmp=-1;
			double lmmp=-1;
			double kmmp=-1;
			double jmmp=-1;
			double dGammaxp=0.;
			
			// @@@@@@@@@@@@   constraints   @@@@@@@@@@@
			//  0:normalized ham, 1:ham, 2:normalized momx, 3:momx 
			//  4:normalized momy,5:momy,6:normalized momz, 7:momz
			//  8:Gamma^x-D_gamma^x, 9:Gamma^y-D_gamma^y, 10:Gamma^z-D_gamma^z
		
			#pragma omp for reduction(+:hamtmp,momtmp,dGamtmp,excnum) 
			for(int j=jjmin;j<=jjmax;j++)
			{
				for(int l=lli;l<=llmax;l++)
				{
					for(int k=kli;k<=kkmax;k++)
					{
						if(get_hflag(l,k,j)!=0 || get_fmrflag(l,k,j)!=0)
						{
							excnum++;
							continue;
						}
						
						if(!isfinite(abs(get_con(l,k,j,0))))
						{
							ofstream fileerr("out_err.dat");
							print_all(fileerr);
							fileerr.close();
							exit(0);
						}
						
						hamtmp+= abs(get_con(l,k,j,0));
						momtmp+= (abs(get_con(l,k,j,2))+abs(get_con(l,k,j,4))+abs(get_con(l,k,j,6)));
						dGamtmp+= (abs(get_con(l,k,j,8))+abs(get_con(l,k,j,9))+abs(get_con(l,k,j,10)));
						
						if(hammaxp<abs(get_con(l,k,j,0)))
						{
							hammaxp=abs(get_con(l,k,j,0));
							lhmp=l;
							khmp=k;
							jhmp=j;
						}
						
						
						double momcompmax=max(max(abs(get_con(l,k,j,2)),abs(get_con(l,k,j,4))),abs(get_con(l,k,j,6)));
						
						if(mommaxp<momcompmax)
						{
							mommaxp=momcompmax;
							lmmp=l;
							kmmp=k;
							jmmp=j;
						}

						double dGamcompmax=max(max(abs(get_con(l,k,j,8)),abs(get_con(l,k,j,9))),abs(get_con(l,k,j,10)));
						
						if(dGammaxp<dGamcompmax)
						{
							dGammaxp=dGamcompmax;
						}
					}
				}
			}
			
			if(hammaxp > hammax)
			{
				#pragma omp critical (hamcheck)
				{
					if(hammaxp > hammax) 
					{
						hammax=hammaxp;
						lhm=lhmp;
						khm=khmp;
						jhm=jhmp;
					}
				}
			}

			if(mommaxp > mommax)
			{
				#pragma omp critical (momcheck)
				{
					if(mommaxp > mommax) 
					{
						mommax=mommaxp;
						lmm=lmmp;
						kmm=kmmp;
						jmm=jmmp;
					}
				}
			}
			
			if(dGammaxp > dGammax)
			{
				#pragma omp critical (dGamcheck)
				{
					if(dGammaxp > dGammax) 
					{
						dGammax=dGammaxp;
					}
				}
			}
		}

		#pragma omp barrier

		double fac=((llmax-lli+1)*(kkmax-kli+1)*(jjmax-jjmin+1))-excnum;
		ham=hamtmp/fac;
		mom=momtmp/(3.*fac);
		dGam=dGamtmp/(3.*fac);
		
	}
	
	void check_Kremax()
	//checking Kretschmann start
	{
		Kremax=0.;
		lkm=-1;
		kkm=-1;
		jkm=-1;
		
		#pragma omp parallel 
		{
			double Kremaxp=0.;
			double lkmp=-1;
			double kkmp=-1;
			double jkmp=-1;
		
			#pragma omp for  
			for(int j=jjmin;j<=jjmax;j++)
			{
				for(int l=lli;l<=llmax;l++)
				{
					for(int k=kli;k<=kkmax;k++)
					{
						if(get_hflag(l,k,j)!=0 || get_fmrflag(l,k,j)!=0)
						continue;
						
						if(Kremaxp<abs(get_outv(l,k,j,0)))
						{
							Kremaxp=abs(get_outv(l,k,j,0));
							lkmp=l;
							kkmp=k;
							jkmp=j;
						}
					}
				}
			}
			
			if(Kremaxp > Kremax)
			{
				#pragma omp critical (Krecheck)
				{
					if(Kremaxp > Kremax) 
					{
						Kremax=Kremaxp;
						lkm=lkmp;
						kkm=kkmp;
						jkm=jkmp;
					}
				}
			}
		}
	}

	void check_Weylmax()
	//checking Weyl-square start
	{
		Weylmax=0.;
		lwm=-1;
		kwm=-1;
		jwm=-1;
		
		#pragma omp parallel 
		{
			double Weylmaxp=0.;
			double lwmp=-1;
			double kwmp=-1;
			double jwmp=-1;
			
			#pragma omp for  
			for(int j=jjmin;j<=jjmax;j++)
			{
				for(int l=lli;l<=llmax;l++)
				{
					for(int k=kli;k<=kkmax;k++)
					{
						if(get_hflag(l,k,j)!=0 || get_fmrflag(l,k,j)!=0)
						continue;
						
						if(Weylmaxp<abs(get_outv(l,k,j,3)))
						{
							Weylmaxp=abs(get_outv(l,k,j,3));
							lwmp=l;
							kwmp=k;
							jwmp=j;
						}
					}
				}
			}
			
			if(Weylmaxp > Weylmax)
			{
				#pragma omp critical (Weylcheck)
				{
					if(Weylmaxp > Weylmax) 
					{
						Weylmax=Weylmaxp;
						lwm=lwmp;
						kwm=kwmp;
						jwm=jwmp;
					}
				}
			}
		}
	}
	
	void set_excflags_square()
	{
		exc=true;
		#pragma omp parallel for
		for(int j=-exg+1;j<exg;j++)
		{
			for(int l=lli;l<lli+exg;l++)
			{
				for(int k=kli;k<kli+exg;k++)
				{
					if((j-(exg-1))*(j-(-exg+1))*(k-(kli+exg-1))*(l-(lli+exg-1))==0)
					set_bflag(l,k,j)=1;
					else if((j-(exg-2))*(j-(-exg+2))*(k-(kli+exg-2))*(l-(lli+exg-2))==0)
					set_bflag(l,k,j)=2;
					else if((j-(exg-3))*(j-(-exg+3))*(k-(kli+exg-3))*(l-(lli+exg-3))==0)
					set_bflag(l,k,j)=3;
					else
					set_bflag(l,k,j)=-1;
				}
			}
		}
		boundary_quarter_excflags();
		#pragma omp barrier
		set_zero_all_exc();
	}
	
	////////////////////////////////////////
	//  fluid func.
	////////////////////////////////////////
	
	double pres(double rho);
	double dpres(double rho);
	void get_rhoGam(double Ene, double S,double& rho,double& Gam);
	double minmod(double a,double b);
	double sign(double A);
	void dyntoprim();
	
	
	////////////////////////////////////////
	//  Interpolation func.
	////////////////////////////////////////
	
	double ipol( double rr,double *xx,double *yy,int order );
	double bv0_ipol(int jc,int kc,int lc,double xc,double yc,double zc,int order,int number);
	double bv1_ipol(int jc,int kc,int lc,double xc,double yc,double zc,int order,int number);
	double bv0_ipol_diff_x(int jc,int kc,int lc,double xc,double yc,double zc,double facx,double qh,int order,int number);
	double bv0_ipol_diff_y(int jc,int kc,int lc,double xc,double yc,double zc,double facy,double qh,int order,int number);
	double bv0_ipol_diff_z(int jc,int kc,int lc,double xc,double yc,double zc,double facz,double qh,int order,int number);
	double bv_ipol(int jc,int kc,int lc,double xc,double yc,double zc,int order,int number);
	double bv_ipol_diff_x(int jc,int kc,int lc,	double xc,double yc,double zc,double facx,double qh,int order,int number);
	double bv_ipol_diff_y(int jc,int kc,int lc,double xc,double yc,double zc,double facy,double qh,int order,int number);
	double bv_ipol_diff_z(int jc,int kc,int lc,double xc,double yc,double zc,double facz,double qh,int order,int number);

	////////////////////////////////////////
	//  BSSN func.
	////////////////////////////////////////
	
	void BSSN_adv();
	void BSSN(int itype);
	void enforce_const();
	void enforce_const_gp(int l,int k,int j);
	void KOdiss();								//Kreiss-Oliger dissipation
	void excision();
	void flux_fill();

	////////////////////////////////////////
	//  BOUNDARY func.
	////////////////////////////////////////

	void boundary_quarter();
	void boundary_d_quarter();
	void boundary_quarter_even(int i);
	void boundary_prim_quarter();
	
	void boundary_quarter_excflags();
	void boundary_quarter_hflags();
	void boundary_psi_initial_quarter();
	void boundary_beta_quarter();
	void boundary_quarter_fluid();
	void boundary_quarter_scalar();
	
	////////////////////////////////////////
	//  Initial func.
	////////////////////////////////////////

	void initial_continue(ifstream& fcontinue);
	double funcf(double X);
	double ifuncf(double X);
	double df(double X);
	double ddf(double X);
	double dddf(double X);
	void set_flat();
	void set_Gam();
	void set_enemomini();
	void initial_params(double cfli,double etaai,double etabi,double etabbi,double lambdai,double dt0i,double dtpi,double dtppi,double ti,double tinii,double Hbi,double KOepi,int exgi,double fluidwi,double scalarmi,double kap_MUSCLi,double b_minmodi);
	
	////////////////////////////////////////
	//  OUTPUT func.
	////////////////////////////////////////

	void print_x(ofstream& fn, int k,int l);
	void print_y(ofstream& fn, int j,int l);
	void print_z(ofstream& fn, int j,int k);
	void print_xy(ofstream& fn, int l);
	void print_xz(ofstream& fn, int k);
	void print_yz(ofstream& fn, int j);
	void print_Kremax(ofstream& fout);
	void print_const(ofstream& fout);
	
	void print_all(ofstream& fout);
	void print_3d(ofstream& fout);

	////////////////////////////////////////
	//  Potential func.
	////////////////////////////////////////
	
	double funcV(double p)
	{
		double w=0.5*pow(scalarm*p,2);
		
		return(w);
	}
	
	double funcdV(double p)
	{
		double w=pow(scalarm,2)*p;
		
		return(w);
	}


};

class Fmv : public Fmv0{
private:

public:
	Fmv(int tabs,int jupper,int jlower,int kupper,int klower,int lupper,int llower,
	double xupper,double xlower,double yupper,double ylower,double zupper,double zlower,double am,bool fld, bool scl, bool cuev) : 
	Fmv0(tabs,jupper,jlower,kupper,klower,lupper,llower,
	xupper,xlower,yupper,ylower,zupper,zlower,am,fld, scl, cuev){
		layn=0;
	}

public:
		
	////////////////////////////////////////
	//  BOUNDARY func.
	////////////////////////////////////////

	void boundary_periodic();
	void boundary_reflection();
	void boundary_d_reflection();
	void boundary_d_gs_reflection();
	void boundary_reflection_even(int i);
	void boundary_prim_reflection();
	
	void boundary_reflection_hflags();
	void boundary_w();
	void boundary_psi_initial();
	void boundary_beta();
	void boundary_reflection_fluid();
	void boundary_reflection_scalar();

	void asymcond(int l,int k,int j,int i,double bgv);
	void boundary_asym0();
	void asymcond(int l,int k,int j,int i,double bgv1,double bgv,double dt,int itype);
	void boundary_asym(int itype);

	////////////////////////////////////////
	//  func. for specific initial conditions
	////////////////////////////////////////

	void set_Psi_nonsph(double mu,double k,double xi2,double xi3);
	void set_Psi_nonsph(double mu,double k,double xi2,double xi3,double xit2,double xit3,double w);
	void initial_nonsph(double mu,double k,double xi2,double xi3,double xit2,double xit3,double w);
	void initial_nonsph(double mu,double k,double xi2,double xi3);
	void set_initial_scalar(double mu,double k,double xi2,double xi3);
	void set_initial_fluid(double mu,double k,double xi2,double xi3);
	void initial(double mu);
	void initial(double mu,double k);
	void set_ini_from_Psi();

};


class Fmv1 : public Fmv0{
private: 
	int ljli,ljui,lkli,lkui,llui,llli,intp,llfmrl,lkfmrl,ljfmrl,llfmru,lkfmru,ljfmru;
	Fmv0* llay;
	Fmv1* ulay;

public:
	Fmv1(int tabs,int jupper,int jlower,int kupper,int klower,int lupper,int llower,
	double xupper,double xlower,double yupper,double ylower,double zupper,double zlower,double am, bool fld, bool scl, bool cuev, Fmv0* lolay)
	: Fmv0(tabs,jupper,jlower,kupper,klower,lupper,llower,xupper,xlower,yupper,ylower,zupper,zlower,am,fld, scl, cuev)
	{
		layn=lolay->get_layn()+1;
		llay=lolay;
		intp=3;
		mrf=false;

		jjmin=jli+3;

		llmax=lui-3;
		kkmax=kui-3;
		jjmax=jui-3;

		llfmrl=llay->get_lli();
		lkfmrl=llay->get_kli();
		ljfmrl=int(get_jli()/2)+2;
		llfmru=int(get_lui()/2)-2;
		lkfmru=int(get_kui()/2)-2;
		ljfmru=int(get_jui()/2)-2;

		llay->set_fmrregion(llfmrl-tab,llfmru,lkfmrl-tab,lkfmru,ljfmrl,ljfmru);

		cout << "Fmv1 for layer #" << layn << " constructer done" << endl;
	}

 ~Fmv1()
	{
	}

public:
	int get_ljli() const{
		return ljli;
	}
	
	void set_ljli(int p){
		ljli=p;
		return;
	}
	int get_lkli() const{
		return lkli;
	}
	
	void set_lkli(int p){
		lkli=p;
		return;
	}
	int get_llli() const{
		return llli;
	}
	
	void set_llli(int p){
		llli=p;
		return;
	}
	void set_ulay(Fmv1* uplay){
		ulay=uplay;
		return;
	}
	
	void set_boundary(int btype,int mm);
	void set_fmr_initial();
	void evolve();
	void refine_llay();
	void onestep(int btype);
	
	void tstep_ipol(int l,int k,int j,int ll, int kk, int jj, int i,double aa,double bb,double cc);

};


#endif
