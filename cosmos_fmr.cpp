/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR FIXED MESH REFINENMENT :: BSSN evolution Class of COSMOS                                */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

//time evolution of this layer
void Fmv1::onestep(int btype)
{
	set01();						//xxx1=xxx0 
	#pragma omp barrier
	setv0();						//xxx0=xxx 
	#pragma omp barrier
	
	set_zero_r();					////xxxr=0 
	#pragma omp barrier
	
	//each step start
	BSSN(1);
	#pragma omp barrier
	set_boundary(btype,intp);
	#pragma omp barrier
	
	//checking constraint 
	check_const();
	#pragma omp barrier

	BSSN(2);
	#pragma omp barrier
	set_boundary(btype,intp);
	#pragma omp barrier
    BSSN(3);
	#pragma omp barrier
    set_boundary(btype+1,intp);
	#pragma omp barrier
	BSSN(4);
	#pragma omp barrier
	set_boundary(btype+1,intp);
	#pragma omp barrier
	//each step end

	set_dtpp(get_dtp());			//one before previous time step
	set_dtp(get_dt0());				//previous time step

	//upper layers 
	if(mrf)
    {
        ulay->evolve();
        #pragma omp barrier
        boundary_quarter();
    }
	
	t+=dt0;							//time forward
	
	return;
}

//refinenment of the lower layer
void Fmv1::refine_llay()
{
	#pragma omp parallel for 
	for(int lj=ljfmrl;lj<=ljfmru;lj++)
	{
		int j=2*lj;
		for(int ll=llfmrl;ll<=llfmru;ll++)
		{
			int l=2*ll;
			for(int lk=lkfmrl;lk<=lkfmru;lk++)
			{
				int k=2*lk;
				for(int i=0;i<nn;i++)
				{
					llay->set_bv(ll,lk,lj,i)=get_bv(l,k,j,i);
				}
			}
		}
	}
	
	return;
}

//time evolution
void Fmv1::evolve()
{
	dt0=0.5*llay->get_dt0();

    set_boundary(1,-2);
	#pragma omp barrier
	onestep(2);
	#pragma omp barrier
	onestep(4);
	#pragma omp barrier
	refine_llay();
}

//inital setting when this layer starts
void Fmv1::set_fmr_initial()
{
	//initial parameter setting start
	cfl=llay->get_cfl();
	etaa=llay->get_etaa();
	etab=llay->get_etab();
	etabb=llay->get_etabb();
	lambda=llay->get_lambda();
	dt0=0.5*llay->get_dt0();
	dtp=dt0;
	dtpp=dt0;
	fluidw=llay->get_fluidw();
	scalarm=llay->get_fluidw();
	t=llay->get_t();
	Hb=llay->get_Hb();
	tini=llay->get_tini();
	KOep=llay->get_KOep();
	exg=llay->get_exg();
	kap_MUSCL=llay->get_Mkap();
	b_minmod=llay->get_b();
	mrf=false;

	llli=llay->get_lli();
	llui=llay->get_lui();
	lkli=llay->get_kli();
	lkui=llay->get_kui();
	ljli=llay->get_jli();
	ljui=llay->get_jui();
	//initial parameter setting end
	
	set_zero_all();
	#pragma omp barrier
	
	//initial setting for inhomogeneous coordinate
	set_flat();

	//setting values with interpolation from the lower layer
	#pragma omp parallel for
	for(int j=jl[tab];j<=ju[tab];j++)
	{
		int lj=int(j/2);
		for(int l=lli;l<=lu[tab];l++)
		{
			int ll=int(l/2);
			for(int k=kli;k<=ku[tab];k++)
			{
				int lk=int(k/2);
				for(int i=0;i<nn;i++)
				{
					set_bv(l,k,j,i)=llay->bv_ipol(lj,lk,ll,get_x(j),get_y(k),get_z(l),ipolo,i);
				}
			}
		}
	}
	#pragma omp barrier
	
	boundary_quarter();			//boundary setting
	#pragma omp barrier
	setv0();					//xxxv=xxx0
	#pragma omp barrier
	dyntoprim();				//fluid primitive variables from dynamical variables
	
	return;
}

//boundary setting from the lower layer
void Fmv1::set_boundary(int itype,int mm)
{
	double delt=llay->get_dt();				//time interval in the lower layer
	double deltp=llay->get_dtpp();			//previous step time interval in the lower layer
	
	double aa,bb,cc;
	
	//setting coefficients for the interpolation
	switch(itype){
	case 1:
		aa=0.;
		bb=1.;
		cc=0.;
		break;
	case 2:
		aa=(-3*pow(delt,2)*pow(deltp,-1)*pow(delt + deltp,-1))/16.;
		bb=(3*(delt + 4*deltp)*pow(deltp,-1))/16.;
		cc=((delt + 4*deltp)*pow(delt + deltp,-1))/16.;
		break;
	case 3:
		aa=-(pow(delt,2)*pow(deltp,-1)*pow(delt + deltp,-1))/4.;
		bb=((delt + 2*deltp)*pow(deltp,-1))/4.;
		cc=((delt + 2*deltp)*pow(delt + deltp,-1))/4.;
		break;
	case 4:
		aa=(-3*pow(delt,2)*pow(deltp,-1)*pow(delt + deltp,-1))/16.;
		bb=0.25 + (3*delt*pow(deltp,-1))/16.;
		cc=(3*(3*delt + 4*deltp)*pow(delt + deltp,-1))/16.;
		break;
	case 5:
		aa=0.;
		bb=0.;
		cc=1.;
		break;
	default:
		aa=0.;
		bb=1.;
		cc=0.;
	}

	/////////////////////////////////////////////////////////////
	int km=kui+mm;
    int lm=lui+mm;
    int jm=jui+mm;
    int jmm=jli-mm;
    
	//run in x direction
	#pragma omp parallel for
	for(int j=jmin;j<=jmax;j++)
	{
		int jj=int(j/2);
		for(int k=km;k<=kmax;k++)
		{
			int kk=int(k/2);
			for(int l=lli;l<=lmax;l++)
			{
				int ll=int(l/2);
				for(int i=0;i<nn;i++)
				tstep_ipol(l,k,j,ll,kk,jj,i,aa,bb,cc);
			}
		}

		for(int k=kli;k<km;k++)
		{
			int kk=int(k/2);
			for(int l=lm;l<=lmax;l++)
			{
				int ll=int(l/2);
				for(int i=0;i<nn;i++)
				tstep_ipol(l,k,j,ll,kk,jj,i,aa,bb,cc);
			}
		}
	}

	//run in z direction
	#pragma omp parallel for 
	for(int l=lli;l<lm;l++)
	{	
		int ll=int(l/2);
		for(int k=kli;k<km;k++)
		{
			int kk=int(k/2);
			for(int j=jmin;j<=jmm;j++)
			{
				int jj=int(j/2);
				for(int i=0;i<nn;i++)
				tstep_ipol(l,k,j,ll,kk,jj,i,aa,bb,cc);
			}
			for(int j=jm;j<=jmax;j++)
			{
				int jj=int(j/2);
				for(int i=0;i<nn;i++)
				tstep_ipol(l,k,j,ll,kk,jj,i,aa,bb,cc);
			}
		}
	}

	#pragma omp barrier
	boundary_quarter();
	
	return;
}

//interpolation function
void Fmv1::tstep_ipol(int l,int k,int j,int ll, int kk, int jj, int i,double aa,double bb,double cc)
{
	if(l%2==0 && j%2==0 && k%2==0)
	set_bv(l,k,j,i)=aa*llay->get_bv1(ll,kk,jj,i)+bb*llay->get_bv0(ll,kk,jj,i)+cc*llay->get_bv(ll,kk,jj,i);
	else
	set_bv(l,k,j,i)=aa*(llay->bv1_ipol(jj,kk,ll,get_x(j),get_y(k),get_z(l),ipolo,i))
			+bb*(llay->bv0_ipol(jj,kk,ll,get_x(j),get_y(k),get_z(l),ipolo,i))
			+cc*(llay->bv_ipol(jj,kk,ll,get_x(j),get_y(k),get_z(l),ipolo,i));

}

