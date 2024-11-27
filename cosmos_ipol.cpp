/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* INITERPOLATION FUNCTIONS :: BSSN evolution Class of COSMOS                                            */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

//interpolation for bv0
double Fmv0::bv0_ipol(int jc,int kc,int lc,double xc,double yc,double zc,
				int order,int number)
//jc,kc,lc:gridnumber of nearest smaller point
//xc,yc,zc:coordinate value of the point
//order:order of the interpolation
//number:specifies the variable for interpolation
{
	////////// variables definition and memory allocation ////////////////////////////
	double ans=0.;										// final answer 
	double *xx, ***kkk, **kxc, *kxyc;
	xx  = new double[order];							// coordinate values of grid points
	kkk  = new double**[order];							// values at points in the cubic region
	kxc = new double*[order];							// interpolated for x
	kxyc= new double[order];							// interpolated for x and y
	int buf=order/2-1;									// number of shift to smaller grids
	////////// variables definition and memory allocation ////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////
	//in the interpolation															//
	for(int l=0;l<order;l++)
	{
		kkk[l] = new double*[order];					// memory allocation
		kxc[l]= new double[order];						// memory allocation
		xx[l]  =0.;										// initialize
		kxyc[l]=0.;										// initialize
		int lt=lc-buf+l;
		for(int k=0;k<order;k++)
		{
			kkk[l][k] = new double[order];				// memory allocation
			kxc[l][k] = 0.;								// initialize
			int kt=kc-buf+k;
			for(int j=0;j<order;j++)
			{
				int jt=jc-buf+j;
				kkk[l][k][j]=get_bv0(lt,kt,jt,number);	// get variables to be used
			}
		}
	}																				//
	////////// preparation of values of variables to be used /////////////////////////

	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++)
	{
		xx[j] = get_x(jc-buf+j);
	}
	
	// x-interpolation
	for(int l=0;l<order;l++)
	{
		for(int k=0;k<order;k++)
		{
			kxc[l][k] = ipol( xc, xx, kkk[l][k], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++)
	{
		xx[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int l=0;l<order;l++)
	{
		kxyc[l] = ipol( yc, xx, kxc[l], order );
		
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++)
	{
		xx[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	ans = ipol( zc, xx, kxyc, order );												//
	////////// z-coordinate interpolation ////////////////////////////////////////////
	
	////////// memory release ////////////////////////////////////////////////////////
	for(int l=0;l<order;l++)														//
	{
		for(int k=0;k<order;k++)
		{
			delete[] kkk[l][k];
		}
		delete[] kxc[l];
		delete[] kkk[l];
	}
	delete[] xx;
	delete[] kkk;
	delete[] kxc;
	delete[] kxyc;																	//
	////////// memory release ////////////////////////////////////////////////////////
	
	return ans;
}

//interpolation for bv0
double Fmv0::bv1_ipol(int jc,int kc,int lc,double xc,double yc,double zc,
				int order,int number)
//jc,kc,lc:gridnumber of nearest smaller point
//xc,yc,zc:coordinate value of the point
//order:order of the interpolation
//number:specifies the variable for interpolation
{
	////////// variables definition and memory allocation ////////////////////////////
	double ans=0.;										// final answer 
	double *xx, ***kkk, **kxc, *kxyc;
	xx  = new double[order];							// coordinate values of grid points
	kkk  = new double**[order];							// values at points in the cubic region
	kxc = new double*[order];							// interpolated for x
	kxyc= new double[order];							// interpolated for x and y
	int buf=order/2-1;									// number of shift to smaller grids
	////////// variables definition and memory allocation ////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////
	//in the interpolation															//
	for(int l=0;l<order;l++)
	{
		kkk[l] = new double*[order];					// memory allocation
		kxc[l]= new double[order];						// memory allocation
		xx[l]  =0.;										// initialize
		kxyc[l]=0.;										// initialize
		int lt=lc-buf+l;
		for(int k=0;k<order;k++)
		{
			kkk[l][k] = new double[order];				// memory allocation
			kxc[l][k] = 0.;								// initialize
			int kt=kc-buf+k;
			for(int j=0;j<order;j++)
			{
				int jt=jc-buf+j;
				kkk[l][k][j]=get_bv1(lt,kt,jt,number);	// get variables to be used
			}
		}
	}																				//
	////////// preparation of values of variables to be used /////////////////////////

	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++)
	{
		xx[j] = get_x(jc-buf+j);
	}
	
	// x-interpolation
	for(int l=0;l<order;l++)
	{
		for(int k=0;k<order;k++)
		{
			kxc[l][k] = ipol( xc, xx, kkk[l][k], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++)
	{
		xx[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int l=0;l<order;l++)
	{
		kxyc[l] = ipol( yc, xx, kxc[l], order );
		
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++)
	{
		xx[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	ans = ipol( zc, xx, kxyc, order );												//
	////////// z-coordinate interpolation ////////////////////////////////////////////
	
	////////// memory release ////////////////////////////////////////////////////////
	for(int l=0;l<order;l++)														//
	{
		for(int k=0;k<order;k++)
		{
			delete[] kkk[l][k];
		}
		delete[] kxc[l];
		delete[] kkk[l];
	}
	delete[] xx;
	delete[] kkk;
	delete[] kxc;
	delete[] kxyc;																	//
	////////// memory release ////////////////////////////////////////////////////////
	
	return ans;
}

//interpolation by using Lagrange formula
double Fmv0::ipol( double rr,double *xx,double *yy,int order )
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

//interpolation for x-derivative of bv0
double Fmv0::bv0_ipol_diff_x(int jc,int kc,int lc,double xc,double yc,double zc,
double facx,double qh,int order,int number)
{
	////////// variables definition //////////////////////////////////////////////////
	double ans=0.;									// final answer
	double dxi=get_dxi();							// 1/dx
	double fac1=facx;								// deviation from nearest smaller grid / dx
	double fac2=1.-facx;							// 1- fac1
	double mindis=1.0e-2;							// minimum acceptted distance to neares grid

	int eoro=(order % 2);							// even or odd
	int lorr=int(fac2+0.5)*eoro;					// left or right
	
	int buf=int(order/2)-1+lorr;					// buffer for smaller grid
	int bufdiff=buf;								// buffer shift for smaller grid(?)
	int midend=buf+1;								// grid number for the end of first half
	int midsta=midend;								// grid number for starting after half
	int fin=order;									// grid number for the end
	int shift=0;									// shift for larger grid(?)
	////////// variables definition //////////////////////////////////////////////////
	

	////////// shift of grid points and corresponding necessary redefs ///////////////
	if(fac1<mindis)																	//
	{
		if(eoro==0)
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}
	if(fac2<mindis)
	{
		if(eoro==0)
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}																				//
	////////// shift of grid points and corresponding necessary redefs ///////////////
	
	
	////////// memory allocation /////////////////////////////////////////////////////
	double *rr,***kkk,**kyc, *kzyc;													//
	
	rr  = new double[order];							// coordinate values of grid points
	kkk  = new double**[order];							// values at points in the cubic region
	kyc = new double*[order];							// interpolated for y
	kzyc= new double[order];							// interpolated for y and z
	////////// memory allocation /////////////////////////////////////////////////////
	

	////////// preparation of values of variables to be used /////////////////////////////////
	//in the interpolation																	//
	for(int j=0;j<order;j++)
	{
		kkk[j] = new double*[order];
		kyc[j]= new double[order];
		kzyc[j]=0.;
		for(int l=0;l<order;l++)
		{
			kkk[j][l] = new double[order];
			kyc[j][l] = 0.;
			rr[l]  =0.;
		}
	}
	
	
	for(int l=0;l<order;l++)
	{
		for(int k=0;k<order;k++)
		{
			for(int j=0;j<midend;j++)
			 kkk[j][l][k]=get_bv0(lc-buf+l,kc-buf+k,jc-bufdiff+j,number);
			
			for(int j=midsta;j<fin;j++)
			 kkk[j-shift][l][k]=get_bv0(lc-buf+l,kc-buf+k,jc-bufdiff+j,number);
			
		}
	}																						//
	////////// preparation of values of variables to be used /////////////////////////////////
	

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++){
		rr[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int j=0;j<order;j++){
		for(int l=0;l<order;l++){
			kyc[j][l] = ipol( yc, rr, kkk[j][l], order );
		}
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++){
		rr[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	for(int j=0;j<order;j++){
		kzyc[j] = ipol( zc, rr, kyc[j], order );
	}																				//
	////////// z-coordinate interpolation ////////////////////////////////////////////

	////////// x-derivative calculation with interpolation ///////////////////////////////////////////////////////////
	if(order==4)																									//
	{
		//4th order finite dif
		double A=(fac1*(fac2 + pow(fac2,2)))
				/((1 + fac1)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double B=-(((1 + fac1)*(fac2 + pow(fac2,2)))
				/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double C=(fac1*(1 + fac1)*(1 + fac2))
				/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double D=-((fac1*(1 + fac1)*fac2)
				/((1 + fac2)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double E=(-fac1 - pow(fac1,2) + fac2 - 2*pow(fac1,2)*fac2 + pow(fac2,2) 
				+ 2*fac1*pow(fac2,2))/(fac1*(1 + fac1)*fac2*(1 + fac2));

		ans = (A*kzyc[0] +B*kzyc[1] +C*kzyc[2] +D*kzyc[3] +E*qh)*dxi;
	}
	else if(order==3)
	{
		if(lorr==1)
		{
			double A=(fac1*fac2)/((1 + fac1)*(1 + fac1 + fac2));
			double B=-(((1 + fac1)*fac2)/(fac1*(fac1 + fac2)));
			double C=(fac1*(1 + fac1))/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
			double D=-((fac1 + pow(fac1,2) - fac2 - 2*fac1*fac2)/(fac1*fac2 + pow(fac1,2)*fac2));
			
			ans=(A*kzyc[0] +B*kzyc[1] +C*kzyc[2] +D*qh)*dxi;
		}
		else
		{
			
			double A=-((fac2*(1 + fac2))/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
			double B=(fac1*(1 + fac2))/(fac2*(fac1 + fac2));
			double C=-((fac1*fac2)/((1 + fac2)*(1 + fac1 + fac2)));
			double D=(-fac1 + fac2 - 2*fac1*fac2 + pow(fac2,2))/(fac1*fac2*(1 + fac2));
			
			ans=(A*kzyc[0] +B*kzyc[1] +C*kzyc[2] +D*qh)*dxi;
		}
	}
	else if(order==2)
	{
		double sumfac=fac1+fac2;
		double A=-fac2/(fac1*sumfac);
		double B=fac1/(fac2*sumfac);
		double C=(-fac1+fac2)/(fac1*fac2);

		ans = (A*kzyc[0] +B*kzyc[1] +C*qh)*dxi;
	}
	else
	{
		cout << "Interpolation Error 2 ..." << endl;
		return 0;
	}																												//
	////////// x-derivative calculation with interpolation ///////////////////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int j=0;j<order;j++){														//
		for(int l=0;l<order;l++){
			delete[] kkk[j][l];
		}
		delete[] kkk[j];
		delete[] kyc[j];
	}
	delete[] kkk;
	delete[] kyc;
	delete[] kzyc;
	delete[] rr;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for y-derivative of bv0
double Fmv0::bv0_ipol_diff_y(int jc,int kc,int lc,double xc,double yc,double zc,
double facy,double qh,int order,int number)
{
	////////// variables definition //////////////////////////////////////////////////
	double ans=0.;																	//
	double dyi=get_dyi();
	double fac1=facy;
	double fac2=1.-facy;
	double mindis=1.0e-2;

	int eoro=(order % 2);
	int lorr=int(fac2+0.5)*eoro;
	
	int buf=int(order/2)-1+lorr;//buffer for smaller grid
	int bufdiff=buf;
	int midend=buf+1;//grid number for the end of first half
	int midsta=midend;//grid number for starting after half
	int fin=order;//grid number for the end
	int shift=0;																	//
	////////// variables definition //////////////////////////////////////////////////
	

	////////// shift of grid points and corresponding necessary redefs ///////////////
	if(fac1<mindis)																	//
	{
		if(eoro==0)
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}
	if(fac2<mindis)
	{
		if(eoro==0)
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}																				//
	////////// shift of grid points and corresponding necessary redefs ///////////////
	

	////////// memory allocation /////////////////////////////////////////////////////
	double *rr,***kkk,**kxc, *kzxc;													//
	
	rr  = new double[order];
	kkk  = new double**[order];
	kxc = new double*[order];
	kzxc= new double[order];														//
	////////// memory allocation /////////////////////////////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////////////
	//in the interpolation																	//
	for(int k=0;k<order;k++)
	{
		kkk[k] = new double*[order];
		kxc[k]= new double[order];
		kzxc[k]=0.;
		for(int l=0;l<order;l++)
		{
			kkk[k][l] = new double[order];
			kxc[k][l] = 0.;
			rr[l]  = 0.;
		}
	}
	
	for(int l=0;l<order;l++)
	{
		for(int j=0;j<order;j++)
		{
			for(int k=0;k<midend;k++)
			 kkk[k][l][j]=get_bv0(lc-buf+l,kc-bufdiff+k,jc-buf+j,number);
			for(int k=midsta;k<fin;k++)
			 kkk[k-shift][l][j]=get_bv0(lc-buf+l,kc-bufdiff+k,jc-buf+j,number);
		}
	}																						//
	////////// preparation of values of variables to be used /////////////////////////////////


	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++){
		rr[j] = get_x(jc-buf+j);
	}

	// x-interpolation
	for(int k=0;k<order;k++){
		for(int l=0;l<order;l++){
			kxc[k][l] = ipol( xc, rr, kkk[k][l], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++){
		rr[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	for(int k=0;k<order;k++){
		kzxc[k] = ipol( zc, rr, kxc[k], order );
	}																				//
	////////// z-coordinate interpolation ////////////////////////////////////////////

	
	////////// y-derivative calculation with interpolation ///////////////////////////////////////////////////////////
	if(order==4)																									//
	{
		//4th order finite dif
		double A=(fac1*(fac2 + pow(fac2,2)))
				/((1 + fac1)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double B=-(((1 + fac1)*(fac2 + pow(fac2,2)))
				/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double C=(fac1*(1 + fac1)*(1 + fac2))
				/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double D=-((fac1*(1 + fac1)*fac2)
				/((1 + fac2)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double E=(-fac1 - pow(fac1,2) + fac2 - 2*pow(fac1,2)*fac2 + pow(fac2,2) 
				+ 2*fac1*pow(fac2,2))/(fac1*(1 + fac1)*fac2*(1 + fac2));

		ans = (A*kzxc[0] +B*kzxc[1] +C*kzxc[2] +D*kzxc[3] +E*qh)*dyi;
	}
	else if(order==3)
	{
		if(lorr==1)
		{
			double A=(fac1*fac2)/((1 + fac1)*(1 + fac1 + fac2));
			double B=-(((1 + fac1)*fac2)/(fac1*(fac1 + fac2)));
			double C=(fac1*(1 + fac1))/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
			double D=-((fac1 + pow(fac1,2) - fac2 - 2*fac1*fac2)/(fac1*fac2 + pow(fac1,2)*fac2));
			
			ans=(A*kzxc[0] +B*kzxc[1] +C*kzxc[2] +D*qh)*dyi;
		}
		else
		{
			
			double A=-((fac2*(1 + fac2))/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
			double B=(fac1*(1 + fac2))/(fac2*(fac1 + fac2));
			double C=-((fac1*fac2)/((1 + fac2)*(1 + fac1 + fac2)));
			double D=(-fac1 + fac2 - 2*fac1*fac2 + pow(fac2,2))/(fac1*fac2*(1 + fac2));
			
			ans=(A*kzxc[0] +B*kzxc[1] +C*kzxc[2] +D*qh)*dyi;
		}
	}
	else if(order==2)
	{
		double sumfac=fac1+fac2;
		double A=-fac2/(fac1*sumfac);
		double B=fac1/(fac2*sumfac);
		double C=(-fac1+fac2)/(fac1*fac2);
		
		ans = (A*kzxc[0] +B*kzxc[1] +C*qh)*dyi;
	}
	else
	{
		cout << "Interpolation Error 2 ..." << endl;
		return 0;
	}																												//
	////////// y-derivative calculation with interpolation ///////////////////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int k=0;k<order;k++){														//
		for(int l=0;l<order;l++){
			delete[] kkk[k][l];
		}
		delete[] kkk[k];
		delete[] kxc[k];
	}
	delete[] kkk;
	delete[] kxc;
	delete[] kzxc;
	delete[] rr;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for z-derivative of bv0
double Fmv0::bv0_ipol_diff_z(int jc,int kc,int lc,double xc,double yc,double zc,
double facz,double qh,int order,int number)
{
	////////// variables definition //////////////////////////////////////////////////
	double ans=0.;																	//
	double dzi=get_dzi();
	double fac1=facz;
	double fac2=1.-facz;
	double mindis=1.0e-2;
	
	int eoro=(order % 2);
	int lorr=int(fac2+0.5)*eoro;

	int buf=int(order/2)-1;//buffer for smaller grid
	int bufdiff=buf;
	int midend=buf+1;//grid number for the end of first half
	int midsta=midend;//grid number for starting after half
	int fin=order;//grid number for the end
	int shift=0;																	//
	////////// variables definition //////////////////////////////////////////////////

	////////// shift of grid points and corresponding necessary redefs ///////////////
	if(fac1<mindis)																	//
	{
		if(eoro==0)
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}
	if(fac2<mindis)
	{
		if(eoro==0)
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}																				//
	////////// shift of grid points and corresponding necessary redefs ///////////////
	

	////////// memory allocation /////////////////////////////////////////////////////
	double *rr,***kkk,**kxc, *kyxc;													//

	rr  = new double[order];
	kkk  = new double**[order];
	kxc = new double*[order];
	kyxc= new double[order];														//
	////////// memory allocation /////////////////////////////////////////////////////


	////////// preparation of values of variables to be used /////////////////////////////////
	//in the interpolation																	//
	for(int l=0;l<order;l++)
	{
		kkk[l] = new double*[order];
		kxc[l]= new double[order];
		kyxc[l]=0.;
		for(int k=0;k<order;k++)
		{
			kkk[l][k] = new double[order];
			kxc[l][k] = 0.;
			rr[k]  =0.;
		}
	}
	
	for(int j=0;j<order;j++)
	{
		for(int k=0;k<order;k++)
		{
			for(int l=0;l<midend;l++)
			 kkk[l][k][j]=get_bv0(lc-bufdiff+l,kc-buf+k,jc-buf+j,number);
			for(int l=midsta;l<fin;l++)
			 kkk[l-shift][k][j]=get_bv0(lc-bufdiff+l,kc-buf+k,jc-buf+j,number);
		}
	}																						//
	////////// preparation of values of variables to be used /////////////////////////////////
	
	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++){
		rr[j] = get_x(jc-buf+j);
	}

	// x-interpolation
	for(int l=0;l<order;l++){
		for(int k=0;k<order;k++){
			kxc[l][k] = ipol( xc, rr, kkk[l][k], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++){
		rr[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int l=0;l<order;l++){
		kyxc[l] = ipol( yc, rr, kxc[l], order );
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-derivative calculation with interpolation ///////////////////////////////////////////////////////////
	if(order==4)																									//
	{
		//4th order finite dif
		double A=(fac1*(fac2 + pow(fac2,2)))
				/((1 + fac1)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double B=-(((1 + fac1)*(fac2 + pow(fac2,2)))
				/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double C=(fac1*(1 + fac1)*(1 + fac2))
				/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double D=-((fac1*(1 + fac1)*fac2)
				/((1 + fac2)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double E=(-fac1 - pow(fac1,2) + fac2 - 2*pow(fac1,2)*fac2 + pow(fac2,2) 
				+ 2*fac1*pow(fac2,2))/(fac1*(1 + fac1)*fac2*(1 + fac2));

		ans = (A*kyxc[0] +B*kyxc[1] +C*kyxc[2] +D*kyxc[3] +E*qh)*dzi;
	}
	else if(order==3)
	{
		if(lorr==1)
		{
			double A=(fac1*fac2)/((1 + fac1)*(1 + fac1 + fac2));
			double B=-(((1 + fac1)*fac2)/(fac1*(fac1 + fac2)));
			double C=(fac1*(1 + fac1))/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
			double D=-((fac1 + pow(fac1,2) - fac2 - 2*fac1*fac2)/(fac1*fac2 + pow(fac1,2)*fac2));
			
			ans=(A*kyxc[0] +B*kyxc[1] +C*kyxc[2] +D*qh)*dzi;
		}
		else
		{
			
			double A=-((fac2*(1 + fac2))/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
			double B=(fac1*(1 + fac2))/(fac2*(fac1 + fac2));
			double C=-((fac1*fac2)/((1 + fac2)*(1 + fac1 + fac2)));
			double D=(-fac1 + fac2 - 2*fac1*fac2 + pow(fac2,2))/(fac1*fac2*(1 + fac2));
			
			ans=(A*kyxc[0] +B*kyxc[1] +C*kyxc[2] +D*qh)*dzi;
		}
	}
	else if(order==2)
	{
		double sumfac=fac1+fac2;
		double A=-fac2/(fac1*sumfac);
		double B=fac1/(fac2*sumfac);
		double C=(-fac1+fac2)/(fac1*fac2);
		
		ans = (A*kyxc[0] +B*kyxc[1] +C*qh)*dzi;
	}
	else
	{
		cout << "Interpolation Error 2 ..." << endl;
		return 0;
	}																												//
	////////// z-derivative calculation with interpolation ///////////////////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int l=0;l<order;l++){														//
		for(int k=0;k<order;k++){
			delete[] kkk[l][k];
		}
		delete[] kkk[l];
		delete[] kxc[l];
	}
	delete[] kkk;
	delete[] kxc;
	delete[] kyxc;
	delete[] rr;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for bv
double Fmv0::bv_ipol(int jc,int kc,int lc,double xc,double yc,double zc,
				int order,int number)
//jc,kc,lc:gridnumber of first smaller point
//xc,yc,zc:coordinate value of the point
//order:order of the interpolation
//number:specifies the variable for interpolation
{
	////////// variables definition and memory allocation ////////////////////////////
	double ans=0.;										// final answer 
	double *xx, ***kkk, **kxc, *kxyc;
	xx  = new double[order];							// coordinate values of grid points
	kkk  = new double**[order];							// values at points in the cubic region
	kxc = new double*[order];							// interpolated for x
	kxyc= new double[order];							// interpolated for x and y
	int buf=order/2-1;									// number of shift to smaller grids
	////////// variables definition and memory allocation ////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////
	//in the pnterpolation															//
	for(int l=0;l<order;l++)
	{
		kkk[l] = new double*[order];
		kxc[l]= new double[order];
		xx[l]  =0.;
		kxyc[l]=0.;
		int lt=lc-buf+l;
		for(int k=0;k<order;k++)
		{
			kkk[l][k] = new double[order];
			kxc[l][k] = 0.;
			int kt=kc-buf+k;
			for(int j=0;j<order;j++)
			{
				int jt=jc-buf+j;
				kkk[l][k][j]=get_bv(lt,kt,jt,number);
			}
		}
	}																				//
	////////// preparation of values of variables to be used /////////////////////////

	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++)
	{
		xx[j] = get_x(jc-buf+j);
	}

	// x-interpolation
	for(int l=0;l<order;l++)
	{
		for(int k=0;k<order;k++)
		{
			kxc[l][k] = ipol( xc, xx, kkk[l][k], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++)
	{
		xx[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int l=0;l<order;l++)
	{
		kxyc[l] = ipol( yc, xx, kxc[l], order );
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++)
	{
		xx[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	ans = ipol( zc, xx, kxyc, order );												//
	////////// z-coordinate interpolation ////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int l=0;l<order;l++)														//
	{
		for(int k=0;k<order;k++)
		{
			delete[] kkk[l][k];
		}
		delete[] kxc[l];
		delete[] kkk[l];
	}
	delete[] xx;
	delete[] kkk;
	delete[] kxc;
	delete[] kxyc;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for x-derivative of bv0
double Fmv0::bv_ipol_diff_x(int jc,int kc,int lc,double xc,double yc,double zc,
double facx,double qh,int order,int number)
{
	////////// variables definition //////////////////////////////////////////////////
	double ans=0.;																	//
	double dxi=get_dxi();
	double fac1=facx;
	double fac2=1.-facx;
	double mindis=1.0e-2;

	int eoro=(order % 2);
	int lorr=int(fac2+0.5)*eoro;
	
	int buf=int(order/2)-1+lorr;				// buffer for smaller grid
	int bufdiff=buf;
	int midend=buf+1;							// grid number for the end of first half
	int midsta=midend;							// grid number for starting after half
	int fin=order;								// grid number for the end
	int shift=0;																	//
	////////// variables definition //////////////////////////////////////////////////

	////////// shift of grid points and corresponding necessary redefs ///////////////
	if(fac1<mindis)																	//
	{
		if(eoro==0)
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}
	if(fac2<mindis)
	{
		if(eoro==0)
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}																				//
	////////// shift of grid points and corresponding necessary redefs ///////////////
	
	////////// memory allocation /////////////////////////////////////////////////////
	double *rr,***kkk,**kyc, *kzyc;													//
	
	rr  = new double[order];
	kkk  = new double**[order];
	kyc = new double*[order];
	kzyc= new double[order];														//
	////////// memory allocation /////////////////////////////////////////////////////


	////////// preparation of values of variables to be used /////////////////////////////////
	//in the interpolation																	//
	for(int j=0;j<order;j++)
	{
		kkk[j] = new double*[order];
		kyc[j]= new double[order];
		kzyc[j]=0.;
		for(int l=0;l<order;l++)
		{
			kkk[j][l] = new double[order];
			kyc[j][l] = 0.;
			rr[l]  =0.;
		}
	}
	
	
	for(int l=0;l<order;l++)
	{
		for(int k=0;k<order;k++)
		{
			for(int j=0;j<midend;j++)
			 kkk[j][l][k]=get_bv(lc-buf+l,kc-buf+k,jc-bufdiff+j,number);
			
			for(int j=midsta;j<fin;j++)
			 kkk[j-shift][l][k]=get_bv(lc-buf+l,kc-buf+k,jc-bufdiff+j,number);
			
		}
	}																						//
	////////// preparation of values of variables to be used /////////////////////////////////
	

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++){
		rr[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int j=0;j<order;j++){
		for(int l=0;l<order;l++){
			kyc[j][l] = ipol( yc, rr, kkk[j][l], order );
		}
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++){
		rr[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	for(int j=0;j<order;j++){
		kzyc[j] = ipol( zc, rr, kyc[j], order );
	}																				//
	////////// z-coordinate interpolation ////////////////////////////////////////////

	////////// x-derivative calculation with interpolation ///////////////////////////////////////////////////////////
	if(order==4)																									//
	{
		//4th order finite dif
		double A=(fac1*(fac2 + pow(fac2,2)))
				/((1 + fac1)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double B=-(((1 + fac1)*(fac2 + pow(fac2,2)))
				/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double C=(fac1*(1 + fac1)*(1 + fac2))
				/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double D=-((fac1*(1 + fac1)*fac2)
				/((1 + fac2)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double E=(-fac1 - pow(fac1,2) + fac2 - 2*pow(fac1,2)*fac2 + pow(fac2,2) 
				+ 2*fac1*pow(fac2,2))/(fac1*(1 + fac1)*fac2*(1 + fac2));

		ans = (A*kzyc[0] +B*kzyc[1] +C*kzyc[2] +D*kzyc[3] +E*qh)*dxi;
	}
	else if(order==3)
	{
		if(lorr==1)
		{
			double A=(fac1*fac2)/((1 + fac1)*(1 + fac1 + fac2));
			double B=-(((1 + fac1)*fac2)/(fac1*(fac1 + fac2)));
			double C=(fac1*(1 + fac1))/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
			double D=-((fac1 + pow(fac1,2) - fac2 - 2*fac1*fac2)/(fac1*fac2 + pow(fac1,2)*fac2));
			
			ans=(A*kzyc[0] +B*kzyc[1] +C*kzyc[2] +D*qh)*dxi;
		}
		else
		{
			
			double A=-((fac2*(1 + fac2))/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
			double B=(fac1*(1 + fac2))/(fac2*(fac1 + fac2));
			double C=-((fac1*fac2)/((1 + fac2)*(1 + fac1 + fac2)));
			double D=(-fac1 + fac2 - 2*fac1*fac2 + pow(fac2,2))/(fac1*fac2*(1 + fac2));
			
			ans=(A*kzyc[0] +B*kzyc[1] +C*kzyc[2] +D*qh)*dxi;
		}
	}
	else if(order==2)
	{
		double sumfac=fac1+fac2;
		double A=-fac2/(fac1*sumfac);
		double B=fac1/(fac2*sumfac);
		double C=(-fac1+fac2)/(fac1*fac2);

		ans = (A*kzyc[0] +B*kzyc[1] +C*qh)*dxi;
	}
	else
	{
		cout << "Interpolation Error 2 ..." << endl;
		return 0;
	}																												//
	////////// x-derivative calculation with interpolation ///////////////////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int j=0;j<order;j++){														//
		for(int l=0;l<order;l++){
			delete[] kkk[j][l];
		}
		delete[] kkk[j];
		delete[] kyc[j];
	}
	delete[] kkk;
	delete[] kyc;
	delete[] kzyc;
	delete[] rr;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for y-derivative of bv
double Fmv0::bv_ipol_diff_y(int jc,int kc,int lc,double xc,double yc,double zc,
double facy,double qh,int order,int number)
{
	////////// variables definition //////////////////////////////////////////////////
	double ans=0.;																	//
	double dyi=get_dyi();
	double fac1=facy;
	double fac2=1.-facy;
	double mindis=1.0e-2;

	int eoro=(order % 2);
	int lorr=int(fac2+0.5)*eoro;
	
	int buf=int(order/2)-1+lorr;					// buffer for smaller grid
	int bufdiff=buf;
	int midend=buf+1;								// grid number for the end of first half
	int midsta=midend;								// grid number for starting after half
	int fin=order;									// grid number for the end
	int shift=0;																	//
	////////// variables definition //////////////////////////////////////////////////
	

	////////// shift of grid points and corresponding necessary redefs ///////////////
	if(fac1<mindis)																	//
	{
		if(eoro==0)
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}
	if(fac2<mindis)
	{
		if(eoro==0)
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}																				//
	////////// shift of grid points and corresponding necessary redefs ///////////////
	
	////////// memory allocation /////////////////////////////////////////////////////
	double *rr,***kkk,**kxc, *kzxc;													//
	
	rr  = new double[order];
	kkk  = new double**[order];
	kxc = new double*[order];
	kzxc= new double[order];														//
	////////// memory allocation /////////////////////////////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////////////
	//in the interpolation																	//
	for(int k=0;k<order;k++)
	{
		kkk[k] = new double*[order];
		kxc[k]= new double[order];
		kzxc[k]=0.;
		for(int l=0;l<order;l++)
		{
			kkk[k][l] = new double[order];
			kxc[k][l] = 0.;
			rr[l]  = 0.;
		}
	}
	
	for(int l=0;l<order;l++)
	{
		for(int j=0;j<order;j++)
		{
			for(int k=0;k<midend;k++)
			 kkk[k][l][j]=get_bv(lc-buf+l,kc-bufdiff+k,jc-buf+j,number);
			for(int k=midsta;k<fin;k++)
			 kkk[k-shift][l][j]=get_bv(lc-buf+l,kc-bufdiff+k,jc-buf+j,number);
		}
	}																						//
	////////// preparation of values of variables to be used /////////////////////////////////

	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++){
		rr[j] = get_x(jc-buf+j);
	}

	// x-interpolation
	for(int k=0;k<order;k++){
		for(int l=0;l<order;l++){
			kxc[k][l] = ipol( xc, rr, kkk[k][l], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// z-coordinate interpolation ////////////////////////////////////////////
	// z-grids																		//
	for(int l=0;l<order;l++){
		rr[l] = get_z(lc-buf+l);
	}

	// z-interpolation
	for(int k=0;k<order;k++){
		kzxc[k] = ipol( zc, rr, kxc[k], order );
	}																				//
	////////// z-coordinate interpolation ////////////////////////////////////////////

	
	////////// y-derivative calculation with interpolation ///////////////////////////////////////////////////////////
	if(order==4)																									//
	{
		//4th order finite dif
		double A=(fac1*(fac2 + pow(fac2,2)))
				/((1 + fac1)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double B=-(((1 + fac1)*(fac2 + pow(fac2,2)))
				/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double C=(fac1*(1 + fac1)*(1 + fac2))
				/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double D=-((fac1*(1 + fac1)*fac2)
				/((1 + fac2)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double E=(-fac1 - pow(fac1,2) + fac2 - 2*pow(fac1,2)*fac2 + pow(fac2,2) 
				+ 2*fac1*pow(fac2,2))/(fac1*(1 + fac1)*fac2*(1 + fac2));

		ans = (A*kzxc[0] +B*kzxc[1] +C*kzxc[2] +D*kzxc[3] +E*qh)*dyi;
	}
	else if(order==3)
	{
		if(lorr==1)
		{
			double A=(fac1*fac2)/((1 + fac1)*(1 + fac1 + fac2));
			double B=-(((1 + fac1)*fac2)/(fac1*(fac1 + fac2)));
			double C=(fac1*(1 + fac1))/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
			double D=-((fac1 + pow(fac1,2) - fac2 - 2*fac1*fac2)/(fac1*fac2 + pow(fac1,2)*fac2));
			
			ans=(A*kzxc[0] +B*kzxc[1] +C*kzxc[2] +D*qh)*dyi;
		}
		else
		{
			
			double A=-((fac2*(1 + fac2))/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
			double B=(fac1*(1 + fac2))/(fac2*(fac1 + fac2));
			double C=-((fac1*fac2)/((1 + fac2)*(1 + fac1 + fac2)));
			double D=(-fac1 + fac2 - 2*fac1*fac2 + pow(fac2,2))/(fac1*fac2*(1 + fac2));
			
			ans=(A*kzxc[0] +B*kzxc[1] +C*kzxc[2] +D*qh)*dyi;
		}
	}
	else if(order==2)
	{
		double sumfac=fac1+fac2;
		double A=-fac2/(fac1*sumfac);
		double B=fac1/(fac2*sumfac);
		double C=(-fac1+fac2)/(fac1*fac2);
		
		ans = (A*kzxc[0] +B*kzxc[1] +C*qh)*dyi;
	}
	else
	{
		cout << "Interpolation Error 2 ..." << endl;
		return 0;
	}																												//
	////////// y-derivative calculation with interpolation ///////////////////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int k=0;k<order;k++){														//
		for(int l=0;l<order;l++){
			delete[] kkk[k][l];
		}
		delete[] kkk[k];
		delete[] kxc[k];
	}
	delete[] kkk;
	delete[] kxc;
	delete[] kzxc;
	delete[] rr;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for z-derivative of bv
double Fmv0::bv_ipol_diff_z(int jc,int kc,int lc,double xc,double yc,double zc,
double facz,double qh,int order,int number)
{
	////////// variables definition //////////////////////////////////////////////////
	double ans=0.;																	//
	double dzi=get_dzi();
	double fac1=facz;
	double fac2=1.-facz;
	double mindis=1.0e-2;

	int eoro=(order % 2);
	int lorr=int(fac2+0.5)*eoro;

	int buf=int(order/2)-1;							// buffer for smaller grid
	int bufdiff=buf;
	int midend=buf+1;								// grid number for the end of first half
	int midsta=midend;								// grid number for starting after half
	int fin=order;									// grid number for the end
	int shift=0;																	//
	////////// variables definition //////////////////////////////////////////////////

	////////// shift of grid points and corresponding necessary redefs ///////////////
	if(fac1<mindis)																	//
	{
		if(eoro==0)
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}
	if(fac2<mindis)
	{
		if(eoro==0)
		{
			fac2+=1.;
			midsta++;
			fin++;
			shift++;
		}
		else
		{
			bufdiff++;
			fac1+=1.;
			midsta++;
			fin++;
			shift++;
		}
	}																				//
	////////// shift of grid points and corresponding necessary redefs ///////////////
	
	////////// memory allocation /////////////////////////////////////////////////////
	double *rr,***kkk,**kxc, *kyxc;													//
	
	rr  = new double[order];
	kkk  = new double**[order];
	kxc = new double*[order];
	kyxc= new double[order];														//
	////////// memory allocation /////////////////////////////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////////////
	//in the interpolation																	//
	for(int l=0;l<order;l++)
	{
		kkk[l] = new double*[order];
		kxc[l]= new double[order];
		kyxc[l]=0.;
		for(int k=0;k<order;k++)
		{
			kkk[l][k] = new double[order];
			kxc[l][k] = 0.;
			rr[k]  =0.;
		}
	}
	
	for(int j=0;j<order;j++)
	{
		for(int k=0;k<order;k++)
		{
			for(int l=0;l<midend;l++)
			 kkk[l][k][j]=get_bv(lc-bufdiff+l,kc-buf+k,jc-buf+j,number);
			for(int l=midsta;l<fin;l++)
			 kkk[l-shift][k][j]=get_bv(lc-bufdiff+l,kc-buf+k,jc-buf+j,number);
		}
	}																						//
	////////// preparation of values of variables to be used /////////////////////////////////
	
	////////// x-coordinate interpolation ////////////////////////////////////////////
	// x-grids																		//
	for(int j=0;j<order;j++){
		rr[j] = get_x(jc-buf+j);
	}

	// x-interpolation
	for(int l=0;l<order;l++){
		for(int k=0;k<order;k++){
			kxc[l][k] = ipol( xc, rr, kkk[l][k], order );
		}
	}																				//
	////////// x-coordinate interpolation ////////////////////////////////////////////

	////////// y-coordinate interpolation ////////////////////////////////////////////
	// y-grids																		//
	for(int k=0;k<order;k++){
		rr[k] = get_y(kc-buf+k);
	}

	// y-interpolation
	for(int l=0;l<order;l++){
		kyxc[l] = ipol( yc, rr, kxc[l], order );
	}																				//
	////////// y-coordinate interpolation ////////////////////////////////////////////

	////////// z-derivative calculation with interpolation ///////////////////////////////////////////////////////////
	if(order==4)																									//
	{
		//4th order finite dif
		double A=(fac1*(fac2 + pow(fac2,2)))
				/((1 + fac1)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double B=-(((1 + fac1)*(fac2 + pow(fac2,2)))
				/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double C=(fac1*(1 + fac1)*(1 + fac2))
				/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
		double D=-((fac1*(1 + fac1)*fac2)
				/((1 + fac2)*(2 + 3*fac1 + pow(fac1,2) + 3*fac2 + 2*fac1*fac2 + pow(fac2,2))));
		double E=(-fac1 - pow(fac1,2) + fac2 - 2*pow(fac1,2)*fac2 + pow(fac2,2) 
				+ 2*fac1*pow(fac2,2))/(fac1*(1 + fac1)*fac2*(1 + fac2));

		ans = (A*kyxc[0] +B*kyxc[1] +C*kyxc[2] +D*kyxc[3] +E*qh)*dzi;
	}
	else if(order==3)
	{
		if(lorr==1)
		{
			double A=(fac1*fac2)/((1 + fac1)*(1 + fac1 + fac2));
			double B=-(((1 + fac1)*fac2)/(fac1*(fac1 + fac2)));
			double C=(fac1*(1 + fac1))/(fac2*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2)));
			double D=-((fac1 + pow(fac1,2) - fac2 - 2*fac1*fac2)/(fac1*fac2 + pow(fac1,2)*fac2));
			
			ans=(A*kyxc[0] +B*kyxc[1] +C*kyxc[2] +D*qh)*dzi;
		}
		else
		{
			
			double A=-((fac2*(1 + fac2))/(fac1*(fac1 + pow(fac1,2) + fac2 + 2*fac1*fac2 + pow(fac2,2))));
			double B=(fac1*(1 + fac2))/(fac2*(fac1 + fac2));
			double C=-((fac1*fac2)/((1 + fac2)*(1 + fac1 + fac2)));
			double D=(-fac1 + fac2 - 2*fac1*fac2 + pow(fac2,2))/(fac1*fac2*(1 + fac2));
			
			ans=(A*kyxc[0] +B*kyxc[1] +C*kyxc[2] +D*qh)*dzi;
		}
	}
	else if(order==2)
	{
		double sumfac=fac1+fac2;
		double A=-fac2/(fac1*sumfac);
		double B=fac1/(fac2*sumfac);
		double C=(-fac1+fac2)/(fac1*fac2);
		
		ans = (A*kyxc[0] +B*kyxc[1] +C*qh)*dzi;
	}
	else
	{
		cout << "Interpolation Error 2 ..." << endl;
		return 0;
	}																												//
	////////// z-derivative calculation with interpolation ///////////////////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	for(int l=0;l<order;l++){														//
		for(int k=0;k<order;k++){
			delete[] kkk[l][k];
		}
		delete[] kkk[l];
		delete[] kxc[l];
	}
	delete[] kkk;
	delete[] kxc;
	delete[] kyxc;
	delete[] rr;																	//
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}


