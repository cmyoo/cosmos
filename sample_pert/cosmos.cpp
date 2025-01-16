//**********************************************************************************************//
//----------------------------------------------------------------------------------------------//
//                     @@@@@   @@@@    @@@@    @@@ @@    @@@@    @@@@                           //
//                   @@@      @@   @  @@   @  @@@ @@ @  @@   @  @@   @                          //
//                   @@       @@   @   @@@    @@  @  @  @@   @   @@@                            //
//                   @@       @@   @ @    @@  @@     @  @@   @ @    @@                          //
//                     @@@@@   @@@@   @@@@@   @@     @   @@@@   @@@@@                           //
//----------------------------------------------------------------------------------------------//
//**********************************************************************************************//
//----------------------------------------------------------------------------------------------//
//                 cosmos: quarter region, mesh-refinement                                      //
//                                             ver. 1.00           coded by Chulmoon Yoo        //
//                                                                                              //
//----------------------------------------------------------------------------------------------//
// Main file   :: cosmos.cpp                                                                    //
// header      :: cosmos.h ahf2d.h                                                              //
// definition  :: cosmos_bssn.cpp cosmos_initial.cpp cosmos_output.cpp                          //
//                cosmos_boundary.cpp cosmos_ahf.cpp cosmos_ipol.cpp                            //
//                cosmos_fluid.cpp cosmos_fmr.cpp                                               //
// how to make :: makefile                                                                      //
//----------------------------------------------------------------------------------------------//
//    Compile :: $ make                                                                         //
//               Choose compiler by changing makefile.                                          //
//    OpenMP  :: $ export OMP_NUM_THREADS=#                                                     //
//               Set the number of CPU cores you use.                                           //
//    Execute :: $ taskset -c 0-# ./cosmos                                                      //
//               Set CPU cores directly on the machine.                                         //
//    Clean   :: $ make clean                                                                   //
//               Delete .o files.                                                               //
//**********************************************************************************************//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <cstring>
#include "../cosmos.h"		//header for Einstein solver
#include "../ahf2d.h"			//header for apparent horizon finder

using namespace std;

void printpack(Fmv0 **fmv0,int ln,int pk,int pl,
ofstream& filex,
ofstream& filey,
ofstream& filez,
ofstream& filex0z,
ofstream& filexy0);

void check_continue_file(					//checking parameter consistency with the continue file
ifstream& fcontinue, 								//initial parameter file
bool& fld,									//for fluid evolution 
bool& scl,									//for scalar evolution
int& tab, 									//tab num for buf grids
int& nxmin, 								//min grid num of x
int& nxmax, 								//max grid num of x
int& nymin,									//min grid num of y
int& nymax, 								//max grid num of y
int& nzmin, 								//min grid num of z
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *jbs,									//grid number for fmr region on x-axis
int *kbs,									//grid number for fmr region on y-axis
int *lbs									//grid number for fmr region on z-axis
);

void output_params(							//output all needed parameters for continue
ofstream& fcontinue, 						//initial parameter file
bool& fld,									//for fluid evolution 
bool& scl,									//for scalar evolution
int& tab, 									//tab num for buf grids
int& nxmin, 								//min grid num of x
int& nxmax, 								//max grid num of x
int& nymin,									//min grid num of y
int& nymax, 								//max grid num of y
int& nzmin, 								//min grid num of z
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *jbs,									//grid number for fmr region on x-axis
int *kbs,									//grid number for fmr region on y-axis
int *lbs									//grid number for fmr region on z-axis
);

void initial_read(							//reading inital para
ifstream& fin, 								//initial parameter file
long int& mstep, 							//max step num for time evo
double& tmax, 								//max time to evolve
int& tab, 									//tab num for buf grids
double& amp, 								//amplitude for inhomogeneous grid
int& nxmin, 								//min grid num of x
int& nxmax, 								//max grid num of x
int& nymin,									//min grid num of y
int& nymax, 								//max grid num of y
int& nzmin, 								//min grid num of z
int& nzmax, 								//max grid num of z
double& xmin, 								//min coord val of x
double& xmax, 								//max coord val of x
double& ymin, 								//min coord val of y
double& ymax, 								//max coord val of y
double& zmin, 								//min coord val of z
double& zmax, 								//max coord val of z
double& cfl, 								//CFL para
double& cdt, 								//factor for cosmological time scale(dt=cfl*cdt*1/H)
double& etaa, 								//etaa for gauge
double& etab, 								//etab for gauge
double& etabb, 								//etabb for gauge
double& KOep, 								//factor for Kreiss-Oliger disspation term
int& exg,									//excision grid
int& contn,									//1 for continue
char *file_continue,						//continue file
double& mu, 								//amplitude of the perturbation
double& kk, 								//scale of the perturbation
double& xi2, 								//nonsph parameter 1
double& xi3, 								//nonsph parameter 2
double& w3, 								//nonsph parameter 3
double& mus, 								//amplitude of the perturbation of scalar field
double& kks, 								//scale of the perturbation of scalar field
double& xi2s, 								//nonsph parameter 1 of laplacial
double& xi3s, 								//nonsph parameter 2 of laplacial
double& Hb, 								//initial Hubble
double& fluidw, 							//fluidw
double& Mkap, 								//kappa in MUSCL
double& bminmod, 							//b in minmod function
double& ptintval1,							//1st part print interval boundary time
double& ptintval2,							//2nd
double& changept							//printing interval change time
);

void initial_fmr(							//reading initial parameter for fmr
ifstream& fin,								//initial parameter file for fmr
int& laymax,								//max fmr layer number
int *jbs,									//grid number for fmr region on x-axis
int *kbs,									//grid number for fmr region on y-axis
int *lbs,									//grid number for fmr region on z-axis
double *alp_fmr								//values of the lapse for starting fmr
);

void ahf_read(								//reading initial para for AHF
ifstream& fin, 								//AHF para file
double& AHFstart,							// time to start AHF
int& ahcint, 								//interval step to call AHF
int& ahpint, 								//interval step to print AH
int& ntheta, 								//grid num of theta
int& nphi, 									//grid num of phi
int& ahloop, 								//max step of AHF
double& var, 								//parameter \eta_apparent
double& facn, 								//factor of mix with the next
double& err_p, 								//error for Poisson solver
double& err_e, 								//error for AH equation
double& ahc 								//initial guess of radius
);


int main(int argc,char* argv[])
{
	bool fld,								//for fluid evolution 
		scl,								//for scalar evolution
		cuev;								//for curvature evaluation
		
	long int mstep; 						//max step num for time evo
	int ahloop, 							//max step of AHF
		ahcint, 							//interval step to call AHF
		ahpint, 							//interval step to print AH
		tab, 								//tab num for buf grids
		nxmin, 								//min grid num of x
		nxmax, 								//max grid num of x
		nymin,								//min grid num of y
		nymax, 								//max grid num of y
		nzmin, 								//min grid num of z
		nzmax, 								//max grid num of z
		ntheta, 							//grid num of theta
		nphi, 								//grid num of phi
		pk,									//printing grid number for y axis
		pl;									//printing grid number for z axis
	
	int exg;								//excision grid

	int ln=0;								//current number of higher layers
	
	double amp,								//for inhom grid parameter
		xmin, 								//min coord val of x
		xmax, 								//max coord val of x
		ymin, 								//min coord val of y
		ymax, 								//max coord val of y
		zmin, 								//min coord val of z
		zmax, 								//max coord val of z
		t, 									//time
		tmax, 								//max time to evolve
		var, 								//parameter \eta_apparent
		facn, 								//factor of mix with the next
		err_p, 								//error for Poisson solver
		err_e, 								//error for AH equation
		ahc, 								//initial guess of radius
		cfl, 								//CFL para
		cdt, 								//factor for cosmological time scale(dt=cfl*cdt*1/H)
		etaa, 								//etaa for gauge
		etab, 								//etab for gauge
		etabb, 								//etabb for gauge	
		AHFstart, 							//time to start AHF
		KOep,								//factor for Kreiss-Oliger dissipation term
		mu, 								//amplitude of the perturbation
		kk, 								//scale of the perturbation
		xi2, 								//nonsph parameter 1
		xi3, 								//nonsph parameter 2
		w3, 								//nonsph parameter 3
		mus, 								//amplitude of the perturbation
		kks, 								//scale of the perturbation
		xi2s, 								//nonsph parameter 1
		xi3s, 								//nonsph parameter 2
		Hb,									//initial Hubble
		fluidw, 							//w for fluid EOS
		Mkap,								//kappa for MUSCL
		bminmod,							//b for minmod 
		hammax,								//maximum value of hamiltonian constraint violation
		mommax;								//maximum value of momentum constraint violation

		int contn;							//1(true) for continue
		
		char file_continue[100]; 			//data file of geometry initial

		
	double ptintval1;						//1st part print interval boundary time
	double ptintval2;						//2nd part print interval boundary time
	double changept;						//printing interval change time
	double nexttimeprint;                   //next time of printing
	
	int laymax;								//maximum number of higher layers
	int jbs[10],kbs[10],lbs[10];			//minimum grid numbers for lower layers
	double alp_fmr[10];						//values of lapse at the center for adding a layer
	
	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);
	
	//reading initial parameter file start
	ifstream fin("par_ini.d");				//open par file
	initial_read(fin,
		mstep,tmax, tab, amp, 
		nxmin, nxmax, nymin,
		nymax, nzmin, nzmax,
		xmin, xmax, ymin,
		ymax, zmin, zmax,
		cfl, cdt, etaa, etab, etabb,
		KOep, exg, contn, file_continue, 
        mu,kk,xi2,xi3, w3, 
        mus,kks,xi2s,xi3s, 
        Hb, fluidw, Mkap, bminmod,
		ptintval1,ptintval2,changept);
	fin.close();
	//reading initial parameter file end
	
	//reading fmr parameter file start
	fin.open("par_fmr.d");				//open par file
	initial_fmr(fin,
		laymax,jbs,kbs,lbs,alp_fmr);
	fin.close();
	//reading fmr parameter file end
	
	//reading ahf setting start 
	fin.open("par_ahf.d");
	ahf_read(fin, AHFstart, ahcint, ahpint, 
	ntheta,nphi, ahloop, var, facn, err_p,err_e, ahc);
	fin.close();
	//reading ahf setting end 
	
	//setting for bools start
	fld=true;						// fluid evolution -> true/false
	scl=false;						// scalar evolution -> true/false
	cuev=false;						// curvature evaluation -> true/false
	//setting for bools end
	
	//main class for bssn
	Fmv *fmv=new Fmv(tab,nxmax,nxmin,nymax,nymin,nzmax,nzmin,
					xmax,xmin,ymax,ymin,zmax,zmin,amp,fld,scl,cuev);
	Fmv1 **fmv1;
	Fmv0 **fmv0;

	fmv0 = new Fmv0 *[laymax+1];
	fmv1 = new Fmv1 *[laymax];
	
	fmv0[0]=fmv;

	//apparent horizon finder class
	Ahf2d *ahf=new Ahf2d(ntheta, nphi, var, facn, ahc);
	
	//initial parameter setting start
	double dr=fmv->get_dx();
	double dtini=cfl*dr;
	double lambda=0.;
	double scalarm=0.;
	double tini=2./(3.*(1+fluidw)*Hb);
	t=tini;

	fmv->initial_params(cfl,etaa,etab,etabb,lambda,dtini,dtini,dtini,t,t,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);
	//initial parameter setting end
	
	//check print line number start
	pk=fmv->get_kli();
	pl=fmv->get_lli();
	//pk=20;
	//pl=109;
	//check print line number end
	
	//output files start
	ofstream filex("out_xkl.dat");						//as functions of x on (k,l)
	filex << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	ofstream filey("out_jyl.dat");						//as functions of y on (j,k)
	filey << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	ofstream filez("out_jkz.dat");						//as functions of z on (j,k)
	filez << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	ofstream filex0z("out_xkz.dat");					//as functions of x,z on lower boundary
	filex0z << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;

	ofstream filexy0("out_xyl.dat");					//as functions of x,y on lower boundary
	filexy0 << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;

	ofstream fileah("out_AH.dat");						//time evolution of apparent horizon variables
	fileah << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	ofstream fileahf("out_AHfig.dat");					//for apparent horizon figure
	fileahf << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	ofstream fileconst("out_const.dat");				//constraint evolution

	ofstream file3d;									//for 3d figure
	
	ofstream fileall;									//for all variables to continue
	//output files end
	
	//reading continue or setting initial date start
	if(contn)
	{
		//file preparateion start
		ifstream fcontinue(file_continue);
		if(fcontinue.fail()){
			cout << "Initial data file to continue is not found." << endl;
			exit(1);
		}
	
		//parameter consistency check
		check_continue_file(fcontinue,fld,scl,tab,nxmin,nxmax,nymin,nymax,nzmin,nzmax,laymax,ln,jbs,kbs,lbs);

		//initial date loading 
		fmv->initial_continue(fcontinue);
		cout << "continue" << endl;

		//boundary setting
		#pragma omp barrier
		fmv->boundary_reflection();

		//line spaces
		string buf;
		getline(fcontinue, buf);
		getline(fcontinue, buf);

		//initial data loading for higher layers start
		for(int i=0;i<ln;i++)
		{
			//class preparation start
			fmv1[i]=new Fmv1(tab,2*jbs[i],-2*jbs[i],2*kbs[i],nymin,2*lbs[i],nzmin,fmv0[i]->get_x(jbs[i]),fmv0[i]->get_x(-jbs[i]),fmv0[i]->get_y(kbs[i]),ymin,fmv0[i]->get_z(lbs[i]),zmin,amp,fld, scl, cuev, fmv0[i]);
			fmv0[i+1]=fmv1[i];
			//class preparation end

			//initial setting for the upper layer start
			double deltat=0.5*fmv0[i]->get_dt0();
			fmv1[i]->initial_params(cfl,etaa,etab,etabb,lambda,deltat,deltat,deltat,t,tini,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);
			fmv1[i]->set_fmr_initial();					
			cout << "initial setting done" << endl;
			//initial setting for the upper layer end
			
			//initial data loading
			fmv1[i]->initial_continue(fcontinue);

			//boundary setting 
			fmv1[i]->boundary_quarter();

			//line spaces
			getline(fcontinue, buf);
			getline(fcontinue, buf);

			if(i>0)
			{
				fmv1[i-1]->set_mrf(true);						//mesh refinement flag
				fmv1[i-1]->set_ulay(fmv1[i]);					//upper layer designation
			}
		}
		//initial data loading for higher layers end
		
		fcontinue.close();
		t=fmv->get_t();
		fmv->setv0();
		printpack(fmv0,ln,pk,pl,filex,filey,filez,filex0z,filexy0);
	}
	else
	{
		cout << "no continue" << endl;
		
		//initial data setting start
		//fmv->set_initial_scalar(mus,kks,xi2s,xi3s);
		//#pragma omp barrier
		//fmv->initial_nonsph(mu,kk,xi2,xi3,xi2s,xi3s,w3);
		fmv->initial(mu);
		#pragma omp barrier
		printpack(fmv0,ln,pk,pl,filex,filey,filez,filex0z,filexy0);
		//initial data setting end

		//printpack(fmv0,ln,pk,pl,filex,filey,filez,filex0z,filexy0);
	}
	//reading continue or setting initial date end
	
	//initial output start
	//file3d.open("out_3d.dat",std::ios::out);
	//fmv->print_3d(file3d);
	//file3d.close();
	cout << "ln=" << ln << endl;
	cout << "##Hini=" << abs(fmv->get_bv(nxmax,nymax,nzmax,20)/3) << endl;
	//initial output end
	
	cout << endl 
	<< "          +---------------------------------------+" << endl
	<< "          |      Set up initial condition         |" << endl
	<< "          +---------------------------------------+" << endl << endl;
	
	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);
	
	//time step settings start
	double dtl=cfl*dr;
	double dtt=cdt*cfl/abs(fmv->get_bv(nxmax,nymax,nzmax,20));

	//setting time interval if not continued
	if(!contn)
	{
		fmv->set_dtpp(fmv->get_dt0());			//one before previous time step
		fmv->set_dtp(fmv->get_dt0());			//previous time step
		fmv->set_dt0(min(dtl,dtt));
	}
	//time step settings end

	//other settings for main loop start
	nexttimeprint=t+ptintval1;
	bool AHFflag=false;							//true for actuation
	bool prehorizon=false;						//previous apparent horizon search?
	//other settings for main loop end
	
	//seach for AH if continued
	if(contn)
	{
		cout << "aparent horizon finder" << endl;
		ahf->find_ah( fmv0[ln], ahloop, err_p, err_e, fileah,1 );
	}

	//if AH is found
	if(!prehorizon && !ahf->get_error())
	{
		cout << "horizon found" << endl;
		fmv0[ln]->set_excflags_square();		//setting excision flags
		//fmv->set_excflags_sphere();
		cout << "excflags set" << endl;
		prehorizon=true;
		ahf->set_hflag(fmv0[ln]);				//setting horizon flags
	}

	//main loop start
	cout << "enter the main loop" << endl;
	
	bool changedt=false;
	for(int step=1;step<mstep+1;step++)
	{ 
		fmv->set01();							//xxx1=xxx0
		#pragma omp barrier
		fmv->setv0();							//xxx0=xxx 
		#pragma omp barrier
		fmv->set_zero_r();						//xxxr=0 
		#pragma omp barrier
		
		//each step start
		fmv->BSSN(1);
		#pragma omp barrier

		fmv->boundary_reflection();
		#pragma omp barrier

		//checking constraint start /////////////////////////////////
		#pragma omp barrier
		fmv->check_const();
		#pragma omp barrier
		//fmv->print_const(fileconst);
		hammax=fmv->get_hammax();
		mommax=fmv->get_mommax();
		#pragma omp barrier
		cout << " time=" << t 
			<< " step=" << step+1 << endl
			<< " alp=" << fmv->get_bv(fmv->get_lli(),fmv->get_kli(),0,0) 
			<< " -ek=" << abs(fmv->get_bv(nxmax,nymax,nzmax,20)) 
			<< endl 
			<< " ham=" << fmv->get_ham() 
			<< " hammax=" << fmv->get_hammax() 
			<< "  (j,k,l)=(" 
			<< fmv->get_jhm() << ","  
			<< fmv->get_khm() << "," 
			<< fmv->get_lhm() << ")" 
			<< " r=" << sqrt(pow(fmv->get_x(fmv->get_jhm()),2)+pow(fmv->get_y(fmv->get_khm()),2)+pow(fmv->get_x(fmv->get_lhm()),2))
			<< endl 
			<< " mom=" << fmv->get_mom() 
			<< " mommax=" << fmv->get_mommax() 
			<< "  (j,k,l)=(" 
			<< fmv->get_jmm() << ","  
			<< fmv->get_kmm() << "," 
			<< fmv->get_lmm() << ")" 
			<< " r=" << sqrt(pow(fmv->get_x(fmv->get_jmm()),2)+pow(fmv->get_y(fmv->get_kmm()),2)+pow(fmv->get_x(fmv->get_lmm()),2))
			<< endl;
		//checking constraint end /////////////////////////////////////
		
		//escape judge //////////////////////////////////////////
		if(t>tmax)
		{
			step=mstep+1;
			// printpack(fmv0,ln,pk,pl,filex,filey,filez,filex0z,filexy0);
			break;
		}	
		/////////////////////////////////////////////////////////
		
		fmv->BSSN(2);
		#pragma omp barrier
		
		fmv->boundary_reflection();
		#pragma omp barrier
		
		fmv->BSSN(3);
		#pragma omp barrier
		
		fmv->boundary_reflection();
		#pragma omp barrier
		
		fmv->BSSN(4);
		#pragma omp barrier
		
		fmv->boundary_reflection();
		#pragma omp barrier
		//each step end
		
		fmv->set_dtpp(fmv->get_dtp());			//one before previous time step
		fmv->set_dtp(fmv->get_dt0());			//previous time step

		//evolution of higher layers start /////////////////////////////////////////////
		if(ln!=0)
		{
			fmv1[0]->evolve();
			#pragma omp barrier
			fmv->boundary_quarter();
		}
		#pragma omp barrier
		//evolution of higher layers end /////////////////////////////////////////////

		//constraint check for higher layers start //////////////////////////////////////
		if(ln!=0)
		{
			//int lni=ln;
			for(int lni=1;lni<=ln;lni++)
			{
				hammax=max(hammax,fmv0[lni]->get_hammax());
				mommax=max(mommax,fmv0[lni]->get_mommax());

				cout << " layer number=" << lni << endl
				<< " ham=" << fmv0[lni]->get_ham() 
				<< " hammax=" << fmv0[lni]->get_hammax() 
				<< "  (j,k,l)=(" 
				<< fmv0[lni]->get_jhm() << ","  
				<< fmv0[lni]->get_khm() << "," 
				<< fmv0[lni]->get_lhm() << ")" 
				<< " r=" << sqrt(pow(fmv0[lni]->get_x(fmv0[lni]->get_jhm()),2)+pow(fmv0[lni]->get_y(fmv0[lni]->get_khm()),2)+pow(fmv0[lni]->get_x(fmv0[lni]->get_lhm()),2))
				<< " mom=" << fmv0[lni]->get_mom() 
				<< " mommax=" << fmv0[lni]->get_mommax() 
				<< "  (j,k,l)=(" 
				<< fmv0[lni]->get_jmm() << ","  
				<< fmv0[lni]->get_kmm() << "," 
				<< fmv0[lni]->get_lmm() << ")" 
				<< " r=" << sqrt(pow(fmv0[lni]->get_x(fmv0[lni]->get_jmm()),2)+pow(fmv0[lni]->get_y(fmv0[lni]->get_kmm()),2)+pow(fmv0[lni]->get_x(fmv0[lni]->get_lmm()),2))
				<< endl << endl;
			}
		}
				
		fileconst << setw(20) << t					//1
		<< " "  << setw(20) << hammax				//2
		<< " "  << setw(20) << mommax				//3
		<< endl;
		//constraint check for higher layers end //////////////////////////////////////

		//time forward start ////////////////////////////////////////
		t=fmv->get_t()+fmv->get_dt0();
		fmv->set_t(t);
		dtt=cdt*cfl/abs(fmv->get_bv(nxmax,nymax,nzmax,20));
		fmv->set_dt0(min(dtl,dtt));

		if(changedt==false)
		{
			if(dtl<dtt)
			{
				cout << "dtl<dtt" << " t="<< t <<  endl;
				changedt=true;
			}
		}
		//time forward end /////////////////////////////////////////////

		//output judge start //////////////////////////////////////
		bool printflag=false;
		
		if(t>nexttimeprint-1.0e-10)
		{
			cout << "t=" << t << " nexttimeprint=" << nexttimeprint << endl;
			printflag=true;
			
			if(t+ptintval1>changept)
			nexttimeprint+=ptintval2;
			else
			nexttimeprint+=ptintval1;
		}
		//output judge end //////////////////////////////////////////
		
		//AHF judge 
		if(!AHFflag && (t>AHFstart || ln==laymax))
		{
			AHFflag=true;
		}

		//horizon finder start
		int quot=step/ahcint;
		int rem=step%ahcint;

		if(rem == 0 && AHFflag)
		{
			cout << "aparent horizon finder" << endl;
			ahf->find_ah( fmv0[ln], ahloop, err_p, err_e, fileah,1 );
			
			//print for AH 
			if(quot%ahpint == 0)
			{
				if(!ahf->get_error())
				{
					ahf->print_ah(fmv0[ln],fileahf,t);			//for AH figure
				}
			}
		}

		//if found
		if(!prehorizon && !ahf->get_error())
		{
			cout << "horizon found" << endl;
			printflag=true;
			fmv0[ln]->set_excflags_square();					//excision flag setting 
			cout << "excflags set" << endl;
			prehorizon=true;
			ahf->set_hflag(fmv0[ln]);							//setting horizon flag
			//fmv->set_excflags_sphere();
		}
		//horizon finder end

		//upper layer setting start
		if(ln<laymax && ahf->get_error())						//if AH is already found not going to upper layer
		{	
			if(fmv->get_bv(0,0,0,0)<alp_fmr[ln])				//condition for the next upper layer
			{
				//upper layer class setting
				fmv1[ln]=new Fmv1(tab,2*jbs[ln],-2*jbs[ln],2*kbs[ln],nymin,2*lbs[ln],nzmin,fmv0[ln]->get_x(jbs[ln]),fmv0[ln]->get_x(-jbs[ln]),fmv0[ln]->get_y(kbs[ln]),ymin,fmv0[ln]->get_z(lbs[ln]),zmin,amp,fld, scl, cuev, fmv0[ln]);
				fmv0[ln+1]=fmv1[ln];

				//initial setting for the upper layer start
				double deltat=0.5*fmv0[ln]->get_dt0();
				fmv1[ln]->initial_params(cfl,etaa,etab,etabb,lambda,deltat,deltat,deltat,t,tini,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);
				fmv1[ln]->set_fmr_initial();					
				cout << "initial setting done" << endl;
				//initial setting for the upper layer end

				if(ln>0)
				{
					fmv1[ln-1]->set_mrf(true);					//mesh refinement flag
					fmv1[ln-1]->set_ulay(fmv1[ln]);				//upper layer designation
				}
				ln++;
				cout << "fmr layer #" << ln << " start" << endl;
				printflag=true;
			}
		}
		//upper layer setting end

		//printing start
		if(printflag)
		{
			cout << "printing" << endl;
			printpack(fmv0,ln,pk,pl,filex,filey,filez,filex0z,filexy0);		//output sections
	
			// file3d.open("out_3d.dat",std::ios::app);
			// fmv->print_3d(file3d);
			// file3d.close();

			//output all data start
			fileall.open("out_all.dat", ios::out );
			
			output_params(fileall,fld,scl,tab,nxmin,nxmax,nymin,nymax,nzmin,nzmax,laymax,ln,jbs,kbs,lbs);
			
			for(int i=0;i<=ln;i++)
			{	
				fmv0[i]->print_all(fileall);
			}
			fileall.close();
			//output all data end
		}
		//printing end
	}
	//main loop end
	
	//final print start
	if(mstep!=0)
	{
		printpack(fmv0,ln,pk,pl,filex,filey,filez,filex0z,filexy0);		//output sections

		if(AHFflag)
		{
			ahf->find_ah( fmv0[ln], ahloop, err_p, err_e, fileah,1 );
			
			if(!ahf->get_error())
			{
				ahf->print_ah(fmv0[ln],fileahf,t);
			}
		}
	
		//output all data start
		fileall.open("out_all.dat", ios::out );

		output_params(fileall,fld,scl,tab,nxmin,nxmax,nymin,nymax,nzmin,nzmax,laymax,ln,jbs,kbs,lbs);
		
		for(int i=0;i<=ln;i++)
		{	
			fmv0[i]->print_all(fileall);
		}
		fileall.close();
		//output all data end
	}
	//final print end
		
	//finalize start
	delete ahf;

	for(int i=0;i<=ln;i++)
	{	
		delete fmv0[i];	
		cout << "fmv0-"<< i << "delete" << endl;
	}

	delete[] fmv0;
	cout << "fmv0 delete" << endl;

	delete[] fmv1;
	cout << "fmv1 delete" << endl;
	//finalize end

	return 0;
}

//printing function
void printpack(Fmv0 **fmv0,int ln,int pk,int pl,
ofstream& filex,
ofstream& filey,
ofstream& filez,
ofstream& filex0z,
ofstream& filexy0)
{
	for(int i=0;i<=ln;i++)
	{
		fmv0[i]->dyntoprim();
		fmv0[i]->print_x(filex,pk,pl);
		fmv0[i]->print_y(filey,pk,pl);
		fmv0[i]->print_z(filez,pk,pl);
		fmv0[i]->print_xz(filex0z, pk);
		fmv0[i]->print_xy(filexy0, pl);
	}
}

//initial parameter reading function
void initial_read(ifstream& fin, 
long int& mstep,
double& tmax, 
int& tab, 
double& amp,
int& nxmin, 
int& nxmax, 
int& nymin,
int& nymax, 
int& nzmin, 
int& nzmax,
double& xmin, 
double& xmax, 
double& ymin,
double& ymax, 
double& zmin, 
double& zmax,
double& cfl, 
double& cdt, 
double& etaa, 
double& etab, 
double& etabb,
double& KOep, 
int& exg, 
int& contn, 
char *file_continue, 
double& mu, 
double& kk, 
double& xi2, 
double& xi3, 
double& w3, 
double& mus, 
double& kks, 
double& xi2s, 
double& xi3s, 
double& Hb, 
double& fluidw,
double& Mkap, 
double& bminmod,
double& ptintval1,
double& ptintval2,
double& changept
){
	string buf;
	char cp[100];
	//cp = NULL;

	if(fin.fail()){
		cout << "Parameter File does not exist." << endl;
		abort();
	}

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%ld %s",&mstep,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&tmax,cp);
	}  
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&tab,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&amp,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nxmin,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nxmax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nymin,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nymax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nzmin,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nzmax,cp);
	}

	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xmin,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xmax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ymin,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ymax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&zmin,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&zmax,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&cfl,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&cdt,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etaa,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etab,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etabb,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&KOep,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&exg,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&contn,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%s %s",file_continue,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&mu,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&kk,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xi2,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xi3,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&w3,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&mus,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&kks,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xi2s,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xi3s,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&Hb,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&fluidw,cp);
	}
	getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&Mkap,cp);
	}
	getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&bminmod,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ptintval1,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ptintval2,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&changept,cp);
	}
	return;
}

//AHF parameter reading function
void ahf_read(ifstream& fin, double& AHFstart,int& ahcint, int& ahpint, 
int& ntheta, int& nphi, int& ahloop, 
double& etaah,double& facn,
double& err_p, double& err_e, double& ahc
){
	string buf;
	char cp[100];
	//  cp = NULL;
	if(fin.fail()){
		cout << "Parameter File does not exist." << endl;
		abort();
	}

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&AHFstart,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&ahcint,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&ahpint,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&ntheta,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nphi,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&ahloop,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etaah,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&facn,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%le %s",&err_p,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%le %s",&err_e,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ahc,cp);
	}
}

//mesh refinement parameter reading function
void initial_fmr(							//reading initial parameter for fmr
ifstream& fin,								//initial parameter file for fmr
int& laymax,								//max fmr layer number
int *jbs,									//grid number for fmr region on x-axis
int *kbs,									//grid number for fmr region on y-axis
int *lbs,									//grid number for fmr region on z-axis
double *alp_fmr	){
	string buf;
	char cp[100];
	//cp = NULL;

	if(fin.fail()){
		cout << "Parameter File does not exist." << endl;
		abort();
	}

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&laymax,cp);
	}
	getline(fin, buf);
	
	for(int n=0;n<=laymax;n++)
	{
		getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
			sscanf(buf.c_str(),"%d %s",&jbs[n],cp);
		}
	}
	getline(fin, buf);
	for(int n=0;n<=laymax;n++)
	{
		getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
			sscanf(buf.c_str(),"%d %s",&kbs[n],cp);
		}
	}
	getline(fin, buf);
	for(int n=0;n<=laymax;n++)
	{
		getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
			sscanf(buf.c_str(),"%d %s",&lbs[n],cp);
		}
	}
	getline(fin, buf);
	
	for(int n=0;n<=laymax;n++)
	{
		getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
			sscanf(buf.c_str(),"%lf %s",&alp_fmr[n],cp);
		}
	}
	return;
}


void check_continue_file(					//checking parameter consistency with the continue file
ifstream& fcontinue, 						//initial parameter file
bool& fld,									//for fluid evolution 
bool& scl,									//for scalar evolution
int& tab, 									//tab num for buf grids
int& nxmin, 								//min grid num of x
int& nxmax, 								//max grid num of x
int& nymin,									//min grid num of y
int& nymax, 								//max grid num of y
int& nzmin, 								//min grid num of z
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *jbs,									//grid number for fmr region on x-axis
int *kbs,									//grid number for fmr region on y-axis
int *lbs									//grid number for fmr region on z-axis
){
		string buf;
		//file preparateion end
		int cpar;
		//get the parameters of the continue file
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##fld=%d",&cpar);
		if(cpar!=fld)
		{
			cout << "fld is different from the initial data file" << endl
				<< "fld in data file=" << cpar << endl
				<< "fld in the setting=" << fld << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##scl=%d",&cpar);
		if(cpar!=scl)
		{
			cout << "scl is different from the initial data file" << endl
				<< "scl in data file=" << cpar << endl
				<< "scl in the setting=" << scl << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##tab=%d",&cpar);
		if(cpar!=tab)
		{
			cout << "tab is different from the initial data file" << endl
				<< "tab in data file=" << cpar << endl
				<< "tab in the setting=" << tab << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nxmin=%d",&cpar);
		if(cpar!=nxmin)
		{
			cout << "nxmin is different from the initial data file" << endl
				<< "nxmin in data file=" << cpar << endl
				<< "nxmin in the setting=" << nxmin << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nxmax=%d",&cpar);
		if(cpar!=nxmax)
		{
			cout << "nxmax is different from the initial data file" << endl
				<< "nxmax in data file=" << cpar << endl
				<< "nxmax in the setting=" << nxmax << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nymin=%d",&cpar);
		if(cpar!=nymin)
		{
			cout << "nymin is different from the initial data file" << endl
				<< "nymin in data file=" << cpar << endl
				<< "nymin in the setting=" << nymin << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nymax=%d",&cpar);
		if(cpar!=nymax)
		{
			cout << "nymax is different from the initial data file" << endl
				<< "nymax in data file=" << cpar << endl
				<< "nymax in the setting=" << nymax << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nzmin=%d",&cpar);
		if(cpar!=nzmin)
		{
			cout << "nzmin is different from the initial data file" << endl
				<< "nzmin in data file=" << cpar << endl
				<< "nzmin in the setting=" << nzmin << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nzmax=%d",&cpar);
		if(cpar!=nzmax)
		{
			cout << "nzmax is different from the initial data file" << endl
				<< "nzmax in data file=" << cpar << endl
				<< "nzmax in the setting=" << nzmax << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##ln=%d",&cpar);
		if(cpar>laymax)
		{
			cout << "ln is larger than laymax" << endl
				<< "ln in data file=" << cpar << endl
				<< "laymax in par_fmr.d=" << laymax << endl;
			exit(1);
		}
		ln=cpar;
		for(int i=0;i<ln;i++)
		{
			getline(fcontinue, buf);
			cout << buf << endl;
			sscanf(buf.data(),"##fmrxgnum=%d",&cpar);
			if(cpar!=jbs[i])
			{
				cout << i+1 << "-th fmrxgnum is different from the initial data file" << endl 
					<< "fmrxgnum in data file=" << cpar << endl
					<< "fmrxgnum in the setting=" << jbs[i] << endl;
			exit(1);
			}
		}
		for(int i=0;i<ln;i++)
		{
			getline(fcontinue, buf);
			cout << buf << endl;
			sscanf(buf.data(),"##fmrygnum=%d",&cpar);
			if(cpar!=kbs[i])
			{
				cout << i+1 << "-th fmrygnum is different from the initial data file" << endl 
					<< "fmrygnum in data file=" << cpar << endl
					<< "fmrygnum in the setting=" << kbs[i] << endl;
			exit(1);
			}
		}
		for(int i=0;i<ln;i++)
		{
			getline(fcontinue, buf);
			cout << buf << endl;
			sscanf(buf.data(),"##fmrzgnum=%d",&cpar);
			if(cpar!=lbs[i])
			{
				cout << i+1 << "-th fmrzgnum is different from the initial data file" << endl 
					<< "fmrzgnum in data file=" << cpar << endl
					<< "fmrzgnum in the setting=" << lbs[i] << endl;
			exit(1);
			}
		}
	return;
}

void output_params(							//output all needed parameters for continue
ofstream& fileall, 							//initial parameter file
bool& fld,									//for fluid evolution 
bool& scl,									//for scalar evolution
int& tab, 									//tab num for buf grids
int& nxmin, 								//min grid num of x
int& nxmax, 								//max grid num of x
int& nymin,									//min grid num of y
int& nymax, 								//max grid num of y
int& nzmin, 								//min grid num of z
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *jbs,									//grid number for fmr region on x-axis
int *kbs,									//grid number for fmr region on y-axis
int *lbs									//grid number for fmr region on z-axis
){
	fileall.setf(ios_base::scientific, ios_base::floatfield);
	fileall.precision(10);
	fileall 
	<< "##fld="<< fld << endl
	<< "##scl="<< scl << endl
	<< "##tab="<< tab << endl
	<< "##nxmin="<< nxmin << endl
	<< "##nxmax="<< nxmax << endl
	<< "##nymin="<< nymin << endl
	<< "##nymax="<< nymax << endl
	<< "##nzmin="<< nzmin << endl
	<< "##nzmax="<< nzmax << endl
	<< "##ln="<< ln << endl;
	for(int i=0;i<ln;i++)
	fileall << "##fmrxgnum="<< jbs[i] << endl;
	for(int i=0;i<ln;i++)
	fileall << "##fmrygnum="<< kbs[i] << endl;
	for(int i=0;i<ln;i++)
	fileall << "##fmrzgnum="<< lbs[i] << endl;

}
