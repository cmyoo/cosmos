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
#include "cosmos.h"			//header for Einstein solver
// #include "ahf2d.h"			//header for apparent horizon finder

using namespace std;

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


int main(int argc,char* argv[])
{
	bool fld,								//for fluid evolution 
		scl,								//for scalar evolution
		cuev;								//for curvature evaluation
		
	long int mstep; 						//max step num for time evo

	int	tab, 								//tab num for buf grids
		nxmin, 								//min grid num of x
		nxmax, 								//max grid num of x
		nymin,								//min grid num of y
		nymax, 								//max grid num of y
		nzmin, 								//min grid num of z
	        nzmax; 								//max grid num of z
		// pk,									//printing grid number for y axis
		// pl;									//printing grid number for z axis
	
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
		cfl, 								//CFL para
		cdt, 								//factor for cosmological time scale(dt=cfl*cdt*1/H)
		etaa, 								//etaa for gauge
		etab, 								//etab for gauge
		etabb, 								//etabb for gauge	
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
	
	//setting for bools start
	fld=false;						// fluid evolution -> true/false
	scl=false;						// scalar evolution -> true/false
	cuev=false;						// curvature evaluation -> true/false
	//setting for bools end
	
	//main class for bssn
	Fmv *fmv=new Fmv(tab,nxmax,nxmin,nymax,nymin,nzmax,nzmin,
					xmax,xmin,ymax,ymin,zmax,zmin,amp,fld,scl,cuev);
	
	//initial parameter setting start
	double dr=fmv->get_dx();
	double dtini=cfl*dr;
	double lambda=0.;
	double scalarm=0.;
	double tini=0.;
	t=tini;

	//initial parameter setting 
	fmv->initial_params(cfl,etaa,etab,etabb,lambda,dtini,dtini,dtini,t,t,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);

	//check print line number start
	// pk=fmv->get_kli();
	// pl=fmv->get_lli();
	//pk=20;
	//pl=109;
	//check print line number end
	
	//output files start
	// ofstream filex("out_xkl.dat");						//as functions of x on (k,l)
	// filex << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	// ofstream filey("out_jyl.dat");						//as functions of y on (j,k)
	// filey << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	// ofstream filez("out_jkz.dat");						//as functions of z on (j,k)
	// filez << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	// ofstream filex0z("out_xkz.dat");					//as functions of x,z on lower boundary
	// filex0z << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;

	// ofstream filexy0("out_xyl.dat");					//as functions of x,y on lower boundary
	// filexy0 << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;

	// ofstream fileah("out_AH.dat");						//time evolution of apparent horizon variables
	// fileah << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	// ofstream fileahf("out_AHfig.dat");					//for apparent horizon figure
	// fileahf << "## nmax=" << nxmax << " xmax=" << xmax <<  endl;
	
	ofstream fileconst("out_const.dat");				//constraint evolution

	ofstream file3d("out_3d.dat",std::ios::out);									//for 3d figure
	
	ofstream fileall;									//for all variables to continue
	//output files end
	
	//initial data setting start
	fmv->initial(mu);
        #pragma omp barrier
	//initial data setting end

	//initial output start
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
	double dtt=cfl*dr;

	fmv->set_dtpp(fmv->get_dt0());			//one before previous time step
	fmv->set_dtp(fmv->get_dt0());			//previous time step
	fmv->set_dt0(min(dtl,dtt));
	// //time step settings end

	//other settings for main loop start
	nexttimeprint=t+ptintval1;
	//other settings for main loop end
	

	//main loop start
	cout << "enter the main loop" << endl;
	
	// bool changedt=false;
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
			<< " step=" << step << endl
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

		#pragma omp barrier


		fileconst << setw(20) << t					//1
		<< " "  << setw(20) << hammax				//2
		<< " "  << setw(20) << mommax				//3
		<< endl;

		//time forward start ////////////////////////////////////////
		t=fmv->get_t()+fmv->get_dt0();
		fmv->set_t(t);
		dtt=cfl*dr;
		fmv->set_dt0(min(dtl,dtt));

		// if(changedt==false)
		// {
		// 	if(dtl<dtt)
		// 	{
		// 		cout << "dtl<dtt" << " t="<< t <<  endl;
		// 		changedt=true;
		// 	}
		// }
		//time forward end /////////////////////////////////////////////

		//output judge start //////////////////////////////////////
		// bool printflag=true;
		
		if(t>nexttimeprint-1.0e-10)
		{
			cout << "t=" << t << " nexttimeprint=" << nexttimeprint << endl;
			
			if(t+ptintval1>changept)
			nexttimeprint+=ptintval2;
			else
			nexttimeprint+=ptintval1;
		}
		//output judge end //////////////////////////////////////////
		
		//escape judge //////////////////////////////////////////
		if(t>tmax)
		 break;
		/////////////////////////////////////////////////////////


		// //printing start
		// if(printflag)
		// {
		// 	cout << "printing" << endl;
		// 	fmv->print_3d(file3d);
		// 	// fmv->print_y(filey,pk,pl);
		// 	// fmv->print_z(filez,pk,pl);
		// 	// fmv->print_xz(filex0z, pk);
		// 	// fmv->print_xy(filexy0, pl);
		// }
		// // //printing end
	}
	//main loop end
	
	//final print start
	if(mstep!=0)
	{
	  fmv->print_3d(file3d);
	  // fmv->print_x(filex,pk,pl);
	  // fmv->print_y(filey,pk,pl);
	  // fmv->print_z(filez,pk,pl);
	  // fmv->print_xz(filex0z, pk);
	  // fmv->print_xy(filexy0, pl);
	}
	//final print end
		
	//finalize start
	file3d.close();
	// filex.close();
	// filey.close();
	// filez.close();
	// filex0z.close();
	// filexy0.close();

	//finalize end

	return 0;
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

