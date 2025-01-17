/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* OUTPUT FUNCTIONS :: BSSN evolution Class of COSMOS                                                    */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"
#include <iomanip>
#include <fstream>

//output as functions of x, k:ygrid num, l:z grid num
void Fmv0::print_x(ofstream& fout, int k, int l)
{
	if(l<lmin)
	{
		fout << "error : l<lmin" << endl;
		return;

	}
	else if(l>lmax)
	{
		fout << "error : l>lmax" << endl;
		return;
	}
	else if(k<kmin)
	{
		fout << "error : k>kmin" << endl;
		return;
	}
	else if(k>kmax)
	{
		fout << "error : k>kmax" << endl;
		return;
	}

	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " k=" << k
		<< " l=" << l
		<< endl;

	for(int j=jmin;j<=jmax;j++)
	{
		double x = get_x(j);
		
		/////////////// output for geometry start ////////////////////////////////////////////////
		double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double M_x,M_y,M_z,nM_x,nM_y,nM_z,dGamx,dGamy,dGamz;
		double alp,hamn,ham;
		double Kinv;
		
		fxx=get_flat_df2x(j);
		fyy=get_flat_df2y(k);
		fzz=get_flat_df2z(l);

		alp= get_bv(l,k,j, 0);
		bx =  get_bv(l,k,j, 1);
		by =  get_bv(l,k,j, 2);
		bz =  get_bv(l,k,j, 3);
		bbx=  get_bv(l,k,j, 4);
		bby=  get_bv(l,k,j, 5);
		bbz=  get_bv(l,k,j, 6);
		gxx=  fxx+get_bv(l,k,j, 7);
		gyy=  fyy+get_bv(l,k,j, 8);
		gzz=  fzz+get_bv(l,k,j, 9);
		gxy=  get_bv(l,k,j,10);
		gxz=  get_bv(l,k,j,11);
		gyz=  get_bv(l,k,j,12);
		wa=   get_bv(l,k,j,13);
		akxx= get_bv(l,k,j,14);
		akyy= get_bv(l,k,j,15);
		akzz= get_bv(l,k,j,16);
		akxy= get_bv(l,k,j,17);
		akxz= get_bv(l,k,j,18);
		akyz= get_bv(l,k,j,19);
		ek=   get_bv(l,k,j,20);
		zgx=  get_bv(l,k,j,21);
		zgy=  get_bv(l,k,j,22);
		zgz=  get_bv(l,k,j,23);
		
		hamn  = get_con(l,k,j,0);
		ham = get_con(l,k,j,1);
		nM_x=get_con(l,k,j,2);
		M_x=get_con(l,k,j,3);
		nM_y=get_con(l,k,j,4);
		M_y=get_con(l,k,j,5);
		nM_z=get_con(l,k,j,6);
		M_z=get_con(l,k,j,7);
		dGamx=get_con(l,k,j,8);
		dGamy=get_con(l,k,j,9);
		dGamz=get_con(l,k,j,10);
		
		Kinv=get_outv(l,k,j,0);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(19) << x 				//1
		<< " "  << setw(20) << alp  			//2
		<< " "  << setw(20) << bx   			//3
		<< " "  << setw(20) << by			//4
		<< " "  << setw(20) << bz			//5
		<< " "  << setw(20) << bbx  			//6
		<< " "  << setw(20) << bby			//7
		<< " "  << setw(20) << bbz			//8
		<< " "  << setw(20) << gxx			//9
		<< " "  << setw(20) << gyy			//10
		<< " "  << setw(20) << gzz			//11
		<< " "  << setw(20) << gxy			//12
		<< " "  << setw(20) << gxz			//13
		<< " "  << setw(20) << gyz			//14
		<< " "  << setw(20) << wa			//15
		<< " "  << setw(20) << akxx			//16
		<< " "  << setw(20) << akyy			//17
		<< " "  << setw(20) << akzz			//18
		<< " "  << setw(20) << akxy			//19
		<< " "  << setw(20) << akxz			//20
		<< " "  << setw(20) << akyz			//21
		<< " "  << setw(20) << ek			//22
		<< " "  << setw(20) << zgx			//23
		<< " "  << setw(20) << zgy			//24
		<< " "  << setw(20) << zgz			//25
		<< " "  << setw(20) << hamn			//26
		<< " "  << setw(20) << ham			//27
		<< " "  << setw(20) << nM_x			//28
		<< " "  << setw(20) << M_x			//29
		<< " "  << setw(20) << nM_y			//30
		<< " "  << setw(20) << M_y			//31
		<< " "  << setw(20) << nM_z			//32
		<< " "  << setw(20) << M_z			//33
		<< " "  << setw(20) << dGamx			//34
		<< " "  << setw(20) << dGamy			//35
		<< " "  << setw(20) << dGamz			//36
		<< " "  << setw(20) << Kinv;			//37
													//
		/////////////// output for geometry end //////////////////////////////////////////////////

		/////////////// output for fluid start ///////////////////////////////////////////////////
		if(fluidevo)
		{
			double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
			
			Ene=  get_bv(l,k,j,24);
			px=   get_bv(l,k,j,25);
			py=   get_bv(l,k,j,26);
			pz=   get_bv(l,k,j,27);
			Den=  get_bv(l,k,j,28);
			
			rho=  get_primv(l,k,j,0);
			Vx=   get_primv(l,k,j,1);
			Vy=   get_primv(l,k,j,2);
			Vz=   get_primv(l,k,j,3);
			eps=  get_primv(l,k,j,4);
			
			fout 
			<< " "  << setw(20) << Ene		//38
			<< " "  << setw(20) << px		//39
			<< " "  << setw(20) << py		//40
			<< " "  << setw(20) << pz		//41
			<< " "  << setw(20) << Den		//42
			<< " "  << setw(20) << rho		//43
			<< " "  << setw(20) << Vx		//44
			<< " "  << setw(20) << Vy		//45
			<< " "  << setw(20) << Vz		//46
			<< " "  << setw(20) << eps;		//47
		}											//
		/////////////// output for fluid end   ///////////////////////////////////////////////////
		
		/////////////// output for scalar start///////////////////////////////////////////////////
		if(scalarevo)										//
		{
			double phii=get_bv(l,k,j,nsc);
			double Pi=get_bv(l,k,j,nscp);

			fout
			<< " "  << setw(20) << phii		//48	//38
			<< " "  << setw(20) << Pi;		//49	//39
		}											//
		/////////////// output for scalar end  ///////////////////////////////////////////////////
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output as functions of y, j:xgrid num, l:z grid num
void Fmv0::print_y(ofstream& fout, int j, int l)
{
	if(j<jmin)
	{
		fout << "error : j<jmin" << endl;
		return;

	}
	else if(j>jmax)
	{
		fout << "error : j>jmax" << endl;
		return;
	}
	else if(l<lmin)
	{
		fout << "error : l<lmin" << endl;
		return;
	}
	else if(l>lmax)
	{
		fout << "error : l>lmax" << endl;
		return;
	}

	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " j="  << j
		<< " l="  << l
		<< endl;

	for(int k=kmin;k<=kmax;k++)
	{
		double y = get_y(k);
		
		/////////////// output for geometry start ////////////////////////////////////////////////
		double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double M_x,M_y,M_z,nM_x,nM_y,nM_z,dGamx,dGamy,dGamz;
		double alp,hamn,ham;
		double Kinv;
		
		fxx=get_flat_df2x(j);
		fyy=get_flat_df2y(k);
		fzz=get_flat_df2z(l);

		alp= get_bv(l,k,j, 0);
		bx =  get_bv(l,k,j, 1);
		by =  get_bv(l,k,j, 2);
		bz =  get_bv(l,k,j, 3);
		bbx=  get_bv(l,k,j, 4);
		bby=  get_bv(l,k,j, 5);
		bbz=  get_bv(l,k,j, 6);
		gxx=  fxx+get_bv(l,k,j, 7);
		gyy=  fyy+get_bv(l,k,j, 8);
		gzz=  fzz+get_bv(l,k,j, 9);
		gxy=  get_bv(l,k,j,10);
		gxz=  get_bv(l,k,j,11);
		gyz=  get_bv(l,k,j,12);
		wa=   get_bv(l,k,j,13);
		akxx= get_bv(l,k,j,14);
		akyy= get_bv(l,k,j,15);
		akzz= get_bv(l,k,j,16);
		akxy= get_bv(l,k,j,17);
		akxz= get_bv(l,k,j,18);
		akyz= get_bv(l,k,j,19);
		ek=   get_bv(l,k,j,20);
		zgx=  get_bv(l,k,j,21);
		zgy=  get_bv(l,k,j,22);
		zgz=  get_bv(l,k,j,23);
		
		hamn  = get_con(l,k,j,0);
		ham = get_con(l,k,j,1);
		nM_x=get_con(l,k,j,2);
		M_x=get_con(l,k,j,3);
		nM_y=get_con(l,k,j,4);
		M_y=get_con(l,k,j,5);
		nM_z=get_con(l,k,j,6);
		M_z=get_con(l,k,j,7);
		dGamx=get_con(l,k,j,8);
		dGamy=get_con(l,k,j,9);
		dGamz=get_con(l,k,j,10);
		
		Kinv=get_outv(l,k,j,0);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(19) << y 				//1
		<< " "  << setw(20) << alp  			//2
		<< " "  << setw(20) << bx   			//3
		<< " "  << setw(20) << by			//4
		<< " "  << setw(20) << bz			//5
		<< " "  << setw(20) << bbx  			//6
		<< " "  << setw(20) << bby			//7
		<< " "  << setw(20) << bbz			//8
		<< " "  << setw(20) << gxx			//9
		<< " "  << setw(20) << gyy			//10
		<< " "  << setw(20) << gzz			//11
		<< " "  << setw(20) << gxy			//12
		<< " "  << setw(20) << gxz			//13
		<< " "  << setw(20) << gyz			//14
		<< " "  << setw(20) << wa			//15
		<< " "  << setw(20) << akxx			//16
		<< " "  << setw(20) << akyy			//17
		<< " "  << setw(20) << akzz			//18
		<< " "  << setw(20) << akxy			//19
		<< " "  << setw(20) << akxz			//20
		<< " "  << setw(20) << akyz			//21
		<< " "  << setw(20) << ek			//22
		<< " "  << setw(20) << zgx			//23
		<< " "  << setw(20) << zgy			//24
		<< " "  << setw(20) << zgz			//25
		<< " "  << setw(20) << hamn			//26
		<< " "  << setw(20) << ham			//27
		<< " "  << setw(20) << nM_x			//28
		<< " "  << setw(20) << M_x			//29
		<< " "  << setw(20) << nM_y			//30
		<< " "  << setw(20) << M_y			//31
		<< " "  << setw(20) << nM_z			//32
		<< " "  << setw(20) << M_z			//33
		<< " "  << setw(20) << dGamx			//34
		<< " "  << setw(20) << dGamy			//35
		<< " "  << setw(20) << dGamz			//36
		<< " "  << setw(20) << Kinv;			//37
													//
		/////////////// output for geometry end //////////////////////////////////////////////////

		/////////////// output for fluid start ///////////////////////////////////////////////////
		if(fluidevo)
		{
			double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
			
			Ene=  get_bv(l,k,j,24);
			px=   get_bv(l,k,j,25);
			py=   get_bv(l,k,j,26);
			pz=   get_bv(l,k,j,27);
			Den=  get_bv(l,k,j,28);
			
			rho=  get_primv(l,k,j,0);
			Vx=   get_primv(l,k,j,1);
			Vy=   get_primv(l,k,j,2);
			Vz=   get_primv(l,k,j,3);
			eps=  get_primv(l,k,j,4);
			
			fout 
			<< " "  << setw(20) << Ene		//38
			<< " "  << setw(20) << px		//39
			<< " "  << setw(20) << py		//40
			<< " "  << setw(20) << pz		//41
			<< " "  << setw(20) << Den		//42
			<< " "  << setw(20) << rho		//43
			<< " "  << setw(20) << Vx		//44
			<< " "  << setw(20) << Vy		//45
			<< " "  << setw(20) << Vz		//46
			<< " "  << setw(20) << eps;		//47
		}											//
		/////////////// output for fluid end   ///////////////////////////////////////////////////
		
		/////////////// output for scalar start///////////////////////////////////////////////////
		if(scalarevo)										//
		{
			double phii=get_bv(l,k,j,nsc);
			double Pi=get_bv(l,k,j,nscp);

			fout
			<< " "  << setw(20) << phii		//48	//38
			<< " "  << setw(20) << Pi;		//49	//39
		}											//
		/////////////// output for scalar end  ///////////////////////////////////////////////////
		
		fout << endl;
	}
	fout << endl << endl;
return;
}

//output as functions of z, j:xgrid num, k:y grid num
void Fmv0::print_z(ofstream& fout, int j, int k)
{

	if(j<jmin)
	{
		fout << "error : j<jmin" << endl;
		return;

	}
	else if(j>jmax)
	{
		fout << "error : j>jmax" << endl;
		return;
	}
	else if(k<kmin)
	{
		fout << "error : k<kmin" << endl;
		return;
	}
	else if(k>kmax)
	{
		fout << "error : k>kmax" << endl;
		return;
	}


	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " j="  << j
		<< " k="  << k
		<< endl;

	for(int l=lmin;l<=lmax;l++)
	{
		double z = get_z(l);
		
		/////////////// output for geometry start ////////////////////////////////////////////////
		double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double M_x,M_y,M_z,nM_x,nM_y,nM_z,dGamx,dGamy,dGamz;
		double alp,hamn,ham;
		double Kinv;

		fxx=get_flat_df2x(j);
		fyy=get_flat_df2y(k);
		fzz=get_flat_df2z(l);
		
		alp= get_bv(l,k,j, 0);
		bx =  get_bv(l,k,j, 1);
		by =  get_bv(l,k,j, 2);
		bz =  get_bv(l,k,j, 3);
		bbx=  get_bv(l,k,j, 4);
		bby=  get_bv(l,k,j, 5);
		bbz=  get_bv(l,k,j, 6);
		gxx=  fxx+get_bv(l,k,j, 7);
		gyy=  fyy+get_bv(l,k,j, 8);
		gzz=  fzz+get_bv(l,k,j, 9);
		gxy=  get_bv(l,k,j,10);
		gxz=  get_bv(l,k,j,11);
		gyz=  get_bv(l,k,j,12);
		wa=   get_bv(l,k,j,13);
		akxx= get_bv(l,k,j,14);
		akyy= get_bv(l,k,j,15);
		akzz= get_bv(l,k,j,16);
		akxy= get_bv(l,k,j,17);
		akxz= get_bv(l,k,j,18);
		akyz= get_bv(l,k,j,19);
		ek=   get_bv(l,k,j,20);
		zgx=  get_bv(l,k,j,21);
		zgy=  get_bv(l,k,j,22);
		zgz=  get_bv(l,k,j,23);
		
		hamn  = get_con(l,k,j,0);
		ham = get_con(l,k,j,1);
		nM_x=get_con(l,k,j,2);
		M_x=get_con(l,k,j,3);
		nM_y=get_con(l,k,j,4);
		M_y=get_con(l,k,j,5);
		nM_z=get_con(l,k,j,6);
		M_z=get_con(l,k,j,7);
		dGamx=get_con(l,k,j,8);
		dGamy=get_con(l,k,j,9);
		dGamz=get_con(l,k,j,10);
		
		Kinv=get_outv(l,k,j,0);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(19) << z 				//1
		<< " "  << setw(20) << alp  			//2
		<< " "  << setw(20) << bx   			//3
		<< " "  << setw(20) << by			//4
		<< " "  << setw(20) << bz			//5
		<< " "  << setw(20) << bbx  			//6
		<< " "  << setw(20) << bby			//7
		<< " "  << setw(20) << bbz			//8
		<< " "  << setw(20) << gxx			//9
		<< " "  << setw(20) << gyy			//10
		<< " "  << setw(20) << gzz			//11
		<< " "  << setw(20) << gxy			//12
		<< " "  << setw(20) << gxz			//13
		<< " "  << setw(20) << gyz			//14
		<< " "  << setw(20) << wa			//15
		<< " "  << setw(20) << akxx			//16
		<< " "  << setw(20) << akyy			//17
		<< " "  << setw(20) << akzz			//18
		<< " "  << setw(20) << akxy			//19
		<< " "  << setw(20) << akxz			//20
		<< " "  << setw(20) << akyz			//21
		<< " "  << setw(20) << ek			//22
		<< " "  << setw(20) << zgx			//23
		<< " "  << setw(20) << zgy			//24
		<< " "  << setw(20) << zgz			//25
		<< " "  << setw(20) << hamn			//26
		<< " "  << setw(20) << ham			//27
		<< " "  << setw(20) << nM_x			//28
		<< " "  << setw(20) << M_x			//29
		<< " "  << setw(20) << nM_y			//30
		<< " "  << setw(20) << M_y			//31
		<< " "  << setw(20) << nM_z			//32
		<< " "  << setw(20) << M_z			//33
		<< " "  << setw(20) << dGamx			//34
		<< " "  << setw(20) << dGamy			//35
		<< " "  << setw(20) << dGamz			//36
		<< " "  << setw(20) << Kinv;			//37
													//
		/////////////// output for geometry end //////////////////////////////////////////////////

		/////////////// output for fluid start ///////////////////////////////////////////////////
		if(fluidevo)
		{
			double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
			
			Ene=  get_bv(l,k,j,24);
			px=   get_bv(l,k,j,25);
			py=   get_bv(l,k,j,26);
			pz=   get_bv(l,k,j,27);
			Den=  get_bv(l,k,j,28);
			
			rho=  get_primv(l,k,j,0);
			Vx=   get_primv(l,k,j,1);
			Vy=   get_primv(l,k,j,2);
			Vz=   get_primv(l,k,j,3);
			eps=  get_primv(l,k,j,4);
			
			fout 
			<< " "  << setw(20) << Ene		//38
			<< " "  << setw(20) << px		//39
			<< " "  << setw(20) << py		//40
			<< " "  << setw(20) << pz		//41
			<< " "  << setw(20) << Den		//42
			<< " "  << setw(20) << rho		//43
			<< " "  << setw(20) << Vx		//44
			<< " "  << setw(20) << Vy		//45
			<< " "  << setw(20) << Vz		//46
			<< " "  << setw(20) << eps;		//47
		}											//
		/////////////// output for fluid end   ///////////////////////////////////////////////////
		
		/////////////// output for scalar start///////////////////////////////////////////////////
		if(scalarevo)										//
		{
			double phii=get_bv(l,k,j,nsc);
			double Pi=get_bv(l,k,j,nscp);
			
			fout
			<< " "  << setw(20) << phii		//48	//38
			<< " "  << setw(20) << Pi;		//49	//39
		}											//
		/////////////// output for scalar end  ///////////////////////////////////////////////////
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output as functions of x and y, l:z grid num
void Fmv0::print_xy(ofstream& fout, int l)
{

	if(l<lmin)
	{
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout << "l<lmin" << endl;
		return;
	}
	else if(l>lmax)
	{
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout << "l>lmax" << endl;
		return;
	}

	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " l=" << l 
		<< endl;

	for(int k=kmin;k<=kmax;k++)
	{
		for(int j=jmin;j<=jmax;j++)
		{
			double x = get_x(j);
			double y = get_y(k);
			
			/////////////// output for geometry start ////////////////////////////////////////////////
			double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
			double akxx,akyy,akzz,akxy,akxz,akyz,ek;
			double zgx,zgy,zgz;
			double bx,by,bz,bbx,bby,bbz;
			double M_x,M_y,M_z,nM_x,nM_y,nM_z,dGamx,dGamy,dGamz;
			double alp,hamn,ham;
			double Kinv;

			fxx=get_flat_df2x(j);
			fyy=get_flat_df2y(k);
			fzz=get_flat_df2z(l);
			
			alp= get_bv(l,k,j, 0);
			bx =  get_bv(l,k,j, 1);
			by =  get_bv(l,k,j, 2);
			bz =  get_bv(l,k,j, 3);
			bbx=  get_bv(l,k,j, 4);
			bby=  get_bv(l,k,j, 5);
			bbz=  get_bv(l,k,j, 6);
			gxx=  fxx+get_bv(l,k,j, 7);
			gyy=  fyy+get_bv(l,k,j, 8);
			gzz=  fzz+get_bv(l,k,j, 9);
			gxy=  get_bv(l,k,j,10);
			gxz=  get_bv(l,k,j,11);
			gyz=  get_bv(l,k,j,12);
			wa=   get_bv(l,k,j,13);
			akxx= get_bv(l,k,j,14);
			akyy= get_bv(l,k,j,15);
			akzz= get_bv(l,k,j,16);
			akxy= get_bv(l,k,j,17);
			akxz= get_bv(l,k,j,18);
			akyz= get_bv(l,k,j,19);
			ek=   get_bv(l,k,j,20);
			zgx=  get_bv(l,k,j,21);
			zgy=  get_bv(l,k,j,22);
			zgz=  get_bv(l,k,j,23);
			
			hamn  = get_con(l,k,j,0);
			ham = get_con(l,k,j,1);
			nM_x=get_con(l,k,j,2);
			M_x=get_con(l,k,j,3);
			nM_y=get_con(l,k,j,4);
			M_y=get_con(l,k,j,5);
			nM_z=get_con(l,k,j,6);
			M_z=get_con(l,k,j,7);
			dGamx=get_con(l,k,j,8);
			dGamy=get_con(l,k,j,9);
			dGamz=get_con(l,k,j,10);
			
			Kinv=get_outv(l,k,j,0);
			
			fout.setf(ios_base::fixed, ios_base::floatfield);
			fout.precision(16);
			fout << setw(19) << x 				//1
			<< " "  << setw(20) << y  			//2
			<< " "  << setw(20) << alp  			//3
			<< " "  << setw(20) << bx   			//4
			<< " "  << setw(20) << by			//5
			<< " "  << setw(20) << bz			//6
			<< " "  << setw(20) << bbx  			//7
			<< " "  << setw(20) << bby			//8
			<< " "  << setw(20) << bbz			//9
			<< " "  << setw(20) << gxx			//10
			<< " "  << setw(20) << gyy			//11
			<< " "  << setw(20) << gzz			//12
			<< " "  << setw(20) << gxy			//13
			<< " "  << setw(20) << gxz			//14
			<< " "  << setw(20) << gyz			//15
			<< " "  << setw(20) << wa			//16
			<< " "  << setw(20) << akxx			//17
			<< " "  << setw(20) << akyy			//18
			<< " "  << setw(20) << akzz			//19
			<< " "  << setw(20) << akxy			//20
			<< " "  << setw(20) << akxz			//21
			<< " "  << setw(20) << akyz			//22
			<< " "  << setw(20) << ek			//23
			<< " "  << setw(20) << zgx			//24
			<< " "  << setw(20) << zgy			//25
			<< " "  << setw(20) << zgz			//26
			<< " "  << setw(20) << hamn			//27
			<< " "  << setw(20) << ham			//28
			<< " "  << setw(20) << nM_x			//29
			<< " "  << setw(20) << M_x			//30
			<< " "  << setw(20) << nM_y			//31
			<< " "  << setw(20) << M_y			//32
			<< " "  << setw(20) << nM_z			//33
			<< " "  << setw(20) << M_z			//34
			<< " "  << setw(20) << dGamx			//35
			<< " "  << setw(20) << dGamy			//36
			<< " "  << setw(20) << dGamz			//37
			<< " "  << setw(20) << Kinv;			//38
														//
			/////////////// output for geometry end //////////////////////////////////////////////////

			/////////////// output for fluid start ///////////////////////////////////////////////////
			if(fluidevo)
			{
				double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
				
				Ene=  get_bv(l,k,j,24);
				px=   get_bv(l,k,j,25);
				py=   get_bv(l,k,j,26);
				pz=   get_bv(l,k,j,27);
				Den=  get_bv(l,k,j,28);
				
				rho=  get_primv(l,k,j,0);
				Vx=   get_primv(l,k,j,1);
				Vy=   get_primv(l,k,j,2);
				Vz=   get_primv(l,k,j,3);
				eps=  get_primv(l,k,j,4);
				
				fout 
				<< " "  << setw(20) << Ene		//39
				<< " "  << setw(20) << px		//40
				<< " "  << setw(20) << py		//41
				<< " "  << setw(20) << pz		//42
				<< " "  << setw(20) << Den		//43
				<< " "  << setw(20) << rho		//44
				<< " "  << setw(20) << Vx		//45
				<< " "  << setw(20) << Vy		//46
				<< " "  << setw(20) << Vz		//47
				<< " "  << setw(20) << eps;		//48
			}											//
			/////////////// output for fluid end   ///////////////////////////////////////////////////
			
			/////////////// output for scalar start///////////////////////////////////////////////////
			if(scalarevo)										//
			{
				double phii=get_bv(l,k,j,nsc);
				double Pi=get_bv(l,k,j,nscp);

				fout
				<< " "  << setw(20) << phii		//49	//39
				<< " "  << setw(20) << Pi;		//50	//40
			}											//
			/////////////// output for scalar end  ///////////////////////////////////////////////////
			
			fout << endl;
		}
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output as functions of x and z, k:y grid num
void Fmv0::print_xz(ofstream& fout, int k)
{
	if(k<kmin)
	{
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout << "k<kmin" << endl;
		return;
	}
	else if(k>kmax)
	{
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout << "k>kmax" << endl;
		return;
	}
	
	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " k=" << k 
		<< endl;

	for(int l=lmin;l<=lmax;l++)
	{
		for(int j=jmin;j<=jmax;j++)
		{
			double x = get_x(j);
			double z = get_z(l);
			
			/////////////// output for geometry start ////////////////////////////////////////////////
			double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
			double akxx,akyy,akzz,akxy,akxz,akyz,ek;
			double zgx,zgy,zgz;
			double bx,by,bz,bbx,bby,bbz;
			double M_x,M_y,M_z,nM_x,nM_y,nM_z,dGamx,dGamy,dGamz;
			double alp,hamn,ham;
			double Kinv;

			fxx=get_flat_df2x(j);
			fyy=get_flat_df2y(k);
			fzz=get_flat_df2z(l);
			
			alp= get_bv(l,k,j, 0);
			bx =  get_bv(l,k,j, 1);
			by =  get_bv(l,k,j, 2);
			bz =  get_bv(l,k,j, 3);
			bbx=  get_bv(l,k,j, 4);
			bby=  get_bv(l,k,j, 5);
			bbz=  get_bv(l,k,j, 6);
			gxx=  fxx+get_bv(l,k,j, 7);
			gyy=  fyy+get_bv(l,k,j, 8);
			gzz=  fzz+get_bv(l,k,j, 9);
			gxy=  get_bv(l,k,j,10);
			gxz=  get_bv(l,k,j,11);
			gyz=  get_bv(l,k,j,12);
			wa=   get_bv(l,k,j,13);
			akxx= get_bv(l,k,j,14);
			akyy= get_bv(l,k,j,15);
			akzz= get_bv(l,k,j,16);
			akxy= get_bv(l,k,j,17);
			akxz= get_bv(l,k,j,18);
			akyz= get_bv(l,k,j,19);
			ek=   get_bv(l,k,j,20);
			zgx=  get_bv(l,k,j,21);
			zgy=  get_bv(l,k,j,22);
			zgz=  get_bv(l,k,j,23);
			
			hamn  = get_con(l,k,j,0);
			ham = get_con(l,k,j,1);
			nM_x=get_con(l,k,j,2);
			M_x=get_con(l,k,j,3);
			nM_y=get_con(l,k,j,4);
			M_y=get_con(l,k,j,5);
			nM_z=get_con(l,k,j,6);
			M_z=get_con(l,k,j,7);
			dGamx=get_con(l,k,j,8);
			dGamy=get_con(l,k,j,9);
			dGamz=get_con(l,k,j,10);
			
			Kinv=get_outv(l,k,j,0);
			
			fout.setf(ios_base::fixed, ios_base::floatfield);
			fout.precision(16);
			fout << setw(19) << x 				//1
			<< " "  << setw(20) << z  			//2
			<< " "  << setw(20) << alp  			//3
			<< " "  << setw(20) << bx   			//4
			<< " "  << setw(20) << by			//5
			<< " "  << setw(20) << bz			//6
			<< " "  << setw(20) << bbx  			//7
			<< " "  << setw(20) << bby			//8
			<< " "  << setw(20) << bbz			//9
			<< " "  << setw(20) << gxx			//10
			<< " "  << setw(20) << gyy			//11
			<< " "  << setw(20) << gzz			//12
			<< " "  << setw(20) << gxy			//13
			<< " "  << setw(20) << gxz			//14
			<< " "  << setw(20) << gyz			//15
			<< " "  << setw(20) << wa			//16
			<< " "  << setw(20) << akxx			//17
			<< " "  << setw(20) << akyy			//18
			<< " "  << setw(20) << akzz			//19
			<< " "  << setw(20) << akxy			//20
			<< " "  << setw(20) << akxz			//21
			<< " "  << setw(20) << akyz			//22
			<< " "  << setw(20) << ek			//23
			<< " "  << setw(20) << zgx			//24
			<< " "  << setw(20) << zgy			//25
			<< " "  << setw(20) << zgz			//26
			<< " "  << setw(20) << hamn			//27
			<< " "  << setw(20) << ham			//28
			<< " "  << setw(20) << nM_x			//29
			<< " "  << setw(20) << M_x			//30
			<< " "  << setw(20) << nM_y			//31
			<< " "  << setw(20) << M_y			//32
			<< " "  << setw(20) << nM_z			//33
			<< " "  << setw(20) << M_z			//34
			<< " "  << setw(20) << dGamx			//35
			<< " "  << setw(20) << dGamy			//36
			<< " "  << setw(20) << dGamz			//37
			<< " "  << setw(20) << Kinv;			//38
														//
			/////////////// output for geometry end //////////////////////////////////////////////////

			/////////////// output for fluid start ///////////////////////////////////////////////////
			if(fluidevo)
			{
				double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
				
				Ene=  get_bv(l,k,j,24);
				px=   get_bv(l,k,j,25);
				py=   get_bv(l,k,j,26);
				pz=   get_bv(l,k,j,27);
				Den=  get_bv(l,k,j,28);
				
				rho=  get_primv(l,k,j,0);
				Vx=   get_primv(l,k,j,1);
				Vy=   get_primv(l,k,j,2);
				Vz=   get_primv(l,k,j,3);
				eps=  get_primv(l,k,j,4);
				
				fout 
				<< " "  << setw(20) << Ene		//39
				<< " "  << setw(20) << px		//40
				<< " "  << setw(20) << py		//41
				<< " "  << setw(20) << pz		//42
				<< " "  << setw(20) << Den		//43
				<< " "  << setw(20) << rho		//44
				<< " "  << setw(20) << Vx		//45
				<< " "  << setw(20) << Vy		//46
				<< " "  << setw(20) << Vz		//47
				<< " "  << setw(20) << eps;		//48
			}											//
			/////////////// output for fluid end   ///////////////////////////////////////////////////
			
			/////////////// output for scalar start///////////////////////////////////////////////////
			if(scalarevo)										//
			{
				double phii=get_bv(l,k,j,nsc);
				double Pi=get_bv(l,k,j,nscp);

				fout
				<< " "  << setw(20) << phii		//49	//39
				<< " "  << setw(20) << Pi;		//50	//40
			}											//
			/////////////// output for scalar end  ///////////////////////////////////////////////////
			
			fout << endl;
		}
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output as functions of y and z, j:x grid num
void Fmv0::print_yz(ofstream& fout, int j)
{
	
	if(j<jmin)
	{
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout << "j<jmin" << endl;
		return;
	}
	else if(j>jmax)
	{
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout << "j>jmax" << endl;
		return;
	}
	
	double tt = get_t();
	
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " j="  << j
		<< endl;

	for(int l=lmin;l<=lmax;l++)
	{
		for(int k=kmin;k<=kmax;k++)
		{
			double y = get_y(k);
			double z = get_z(l);
			
			/////////////// output for geometry start ////////////////////////////////////////////////
			double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
			double akxx,akyy,akzz,akxy,akxz,akyz,ek;
			double zgx,zgy,zgz;
			double bx,by,bz,bbx,bby,bbz;
			double M_x,M_y,M_z,nM_x,nM_y,nM_z,dGamx,dGamy,dGamz;
			double alp,hamn,ham;
			double Kinv;

			fxx=get_flat_df2x(j);
			fyy=get_flat_df2y(k);
			fzz=get_flat_df2z(l);
			
			alp= get_bv(l,k,j, 0);
			bx =  get_bv(l,k,j, 1);
			by =  get_bv(l,k,j, 2);
			bz =  get_bv(l,k,j, 3);
			bbx=  get_bv(l,k,j, 4);
			bby=  get_bv(l,k,j, 5);
			bbz=  get_bv(l,k,j, 6);
			gxx=  fxx+get_bv(l,k,j, 7);
			gyy=  fyy+get_bv(l,k,j, 8);
			gzz=  fzz+get_bv(l,k,j, 9);
			gxy=  get_bv(l,k,j,10);
			gxz=  get_bv(l,k,j,11);
			gyz=  get_bv(l,k,j,12);
			wa=   get_bv(l,k,j,13);
			akxx= get_bv(l,k,j,14);
			akyy= get_bv(l,k,j,15);
			akzz= get_bv(l,k,j,16);
			akxy= get_bv(l,k,j,17);
			akxz= get_bv(l,k,j,18);
			akyz= get_bv(l,k,j,19);
			ek=   get_bv(l,k,j,20);
			zgx=  get_bv(l,k,j,21);
			zgy=  get_bv(l,k,j,22);
			zgz=  get_bv(l,k,j,23);
			
			hamn  = get_con(l,k,j,0);
			ham = get_con(l,k,j,1);
			nM_x=get_con(l,k,j,2);
			M_x=get_con(l,k,j,3);
			nM_y=get_con(l,k,j,4);
			M_y=get_con(l,k,j,5);
			nM_z=get_con(l,k,j,6);
			M_z=get_con(l,k,j,7);
			dGamx=get_con(l,k,j,8);
			dGamy=get_con(l,k,j,9);
			dGamz=get_con(l,k,j,10);
			
			Kinv=get_outv(l,k,j,0);
			
			fout.setf(ios_base::fixed, ios_base::floatfield);
			fout.precision(16);
			fout << setw(19) << y 				//1
			<< " "  << setw(20) << z  			//2
			<< " "  << setw(20) << alp  			//3
			<< " "  << setw(20) << bx   			//4
			<< " "  << setw(20) << by			//5
			<< " "  << setw(20) << bz			//6
			<< " "  << setw(20) << bbx  			//7
			<< " "  << setw(20) << bby			//8
			<< " "  << setw(20) << bbz			//9
			<< " "  << setw(20) << gxx			//10
			<< " "  << setw(20) << gyy			//11
			<< " "  << setw(20) << gzz			//12
			<< " "  << setw(20) << gxy			//13
			<< " "  << setw(20) << gxz			//14
			<< " "  << setw(20) << gyz			//15
			<< " "  << setw(20) << wa			//16
			<< " "  << setw(20) << akxx			//17
			<< " "  << setw(20) << akyy			//18
			<< " "  << setw(20) << akzz			//19
			<< " "  << setw(20) << akxy			//20
			<< " "  << setw(20) << akxz			//21
			<< " "  << setw(20) << akyz			//22
			<< " "  << setw(20) << ek			//23
			<< " "  << setw(20) << zgx			//24
			<< " "  << setw(20) << zgy			//25
			<< " "  << setw(20) << zgz			//26
			<< " "  << setw(20) << hamn			//27
			<< " "  << setw(20) << ham			//28
			<< " "  << setw(20) << nM_x			//29
			<< " "  << setw(20) << M_x			//30
			<< " "  << setw(20) << nM_y			//31
			<< " "  << setw(20) << M_y			//32
			<< " "  << setw(20) << nM_z			//33
			<< " "  << setw(20) << M_z			//34
			<< " "  << setw(20) << dGamx			//35
			<< " "  << setw(20) << dGamy			//36
			<< " "  << setw(20) << dGamz			//37
			<< " "  << setw(20) << Kinv;			//38
														//
			/////////////// output for geometry end //////////////////////////////////////////////////

			/////////////// output for fluid start ///////////////////////////////////////////////////
			if(fluidevo)
			{
				double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
				
				Ene=  get_bv(l,k,j,24);
				px=   get_bv(l,k,j,25);
				py=   get_bv(l,k,j,26);
				pz=   get_bv(l,k,j,27);
				Den=  get_bv(l,k,j,28);
				
				rho=  get_primv(l,k,j,0);
				Vx=   get_primv(l,k,j,1);
				Vy=   get_primv(l,k,j,2);
				Vz=   get_primv(l,k,j,3);
				eps=  get_primv(l,k,j,4);
				
				fout 
				<< " "  << setw(20) << Ene		//39
				<< " "  << setw(20) << px		//40
				<< " "  << setw(20) << py		//41
				<< " "  << setw(20) << pz		//42
				<< " "  << setw(20) << Den		//43
				<< " "  << setw(20) << rho		//44
				<< " "  << setw(20) << Vx		//45
				<< " "  << setw(20) << Vy		//46
				<< " "  << setw(20) << Vz		//47
				<< " "  << setw(20) << eps;		//48
			}											//
			/////////////// output for fluid end   ///////////////////////////////////////////////////
			
			/////////////// output for scalar start///////////////////////////////////////////////////
			if(scalarevo)										//
			{
				double phii=get_bv(l,k,j,nsc);
				double Pi=get_bv(l,k,j,nscp);

				fout
				<< " "  << setw(20) << phii		//49	//39
				<< " "  << setw(20) << Pi;		//50	//40
			}											//
			/////////////// output for scalar end  ///////////////////////////////////////////////////
			
			fout << endl;
		}
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

void Fmv0::print_Kremax(ofstream& fout)
//printing curvature invariants
{
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << setw(20) << get_t()					//1
	<< " "  << setw(20) << get_Kremax()				//2
	<< " "  << setw(20) << get_x(jkm)				//3
	<< " "  << setw(20) << get_y(kkm)				//4
	<< " "  << setw(20) << get_z(lkm)				//5
	<< " "  << setw(20) << get_Weylmax()				//6
	<< " "  << setw(20) << get_x(jwm)				//7
	<< " "  << setw(20) << get_y(kwm)				//8
	<< " "  << setw(20) << get_z(lwm)				//9
	<< endl;

}

void Fmv0::print_const(ofstream& fout)
//printing constraints
{
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << setw(20) << get_t()					//1
	<< " "  << setw(20) << get_ham()				//2
	<< " "  << setw(20) << get_hammax()				//3
	<< " "  << setw(20) << get_mom()				//4
	<< " "  << setw(20) << get_mommax()				//5
	<< " "  << setw(20) << get_x(jhm)				//6
	<< " "  << setw(20) << get_y(khm)				//7
	<< " "  << setw(20) << get_z(lhm)				//8
	<< endl;

}

//output all variables
void Fmv0::print_all(ofstream& fout)
{
	double tt = get_t();
	cout << "print_all:time=" << tt 
		<< endl;
	
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << "##" 
		<< "time=" << tt 
		<< " dtp=" << dtp 
		<< " dtpp=" << dtpp 
		<< endl;

	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jli;j<=jui;j++)
			{
				double hxx,hyy,hzz,gxy,gxz,gyz,wa;
				double akxx,akyy,akzz,akxy,akxz,akyz,ek;
				double zgx,zgy,zgz;
				double bx,by,bz,bbx,bby,bbz;
				double alp;
				
				alp= get_bv(l,k,j, 0);
				bx =  get_bv(l,k,j, 1);
				by =  get_bv(l,k,j, 2);
				bz =  get_bv(l,k,j, 3);
				bbx=  get_bv(l,k,j, 4);
				bby=  get_bv(l,k,j, 5);
				bbz=  get_bv(l,k,j, 6);
				hxx=  get_bv(l,k,j, 7);
				hyy=  get_bv(l,k,j, 8);
				hzz=  get_bv(l,k,j, 9);
				gxy=  get_bv(l,k,j,10);
				gxz=  get_bv(l,k,j,11);
				gyz=  get_bv(l,k,j,12);
				wa=   get_bv(l,k,j,13);
				akxx= get_bv(l,k,j,14);
				akyy= get_bv(l,k,j,15);
				akzz= get_bv(l,k,j,16);
				akxy= get_bv(l,k,j,17);
				akxz= get_bv(l,k,j,18);
				akyz= get_bv(l,k,j,19);
				ek=   get_bv(l,k,j,20);
				zgx=  get_bv(l,k,j,21);
				zgy=  get_bv(l,k,j,22);
				zgz=  get_bv(l,k,j,23);
				
				fout.setf(ios_base::fixed, ios_base::floatfield);
				fout.precision(16);
				fout << setw(20) << alp 	 		//1
				<< " "  << setw(20) << bx   			//2
				<< " "  << setw(20) << by			//3
				<< " "  << setw(20) << bz			//4
				<< " "  << setw(20) << bbx  			//5
				<< " "  << setw(20) << bby			//6
				<< " "  << setw(20) << bbz			//7
				<< " "  << setw(20) << hxx			//8
				<< " "  << setw(20) << hyy			//9
				<< " "  << setw(20) << hzz			//10
				<< " "  << setw(20) << gxy			//11
				<< " "  << setw(20) << gxz			//12
				<< " "  << setw(20) << gyz			//13
				<< " "  << setw(20) << wa			//14
				<< " "  << setw(20) << akxx			//15
				<< " "  << setw(20) << akyy			//16
				<< " "  << setw(20) << akzz			//17
				<< " "  << setw(20) << akxy			//18
				<< " "  << setw(20) << akxz			//19
				<< " "  << setw(20) << akyz			//20
				<< " "  << setw(20) << ek			//21
				<< " "  << setw(20) << zgx			//22
				<< " "  << setw(20) << zgy			//23
				<< " "  << setw(20) << zgz;			//24
				
				if(fluidevo)
				{
					double Ene,px,py,pz,Den;

					Ene=  get_bv(l,k,j,24);
					px=   get_bv(l,k,j,25);
					py=   get_bv(l,k,j,26);
					pz=   get_bv(l,k,j,27);
					Den=  get_bv(l,k,j,28);
					
					fout 
					<< " "  << setw(20) << Ene		//25
					<< " "  << setw(20) << px		//26
					<< " "  << setw(20) << py		//27
					<< " "  << setw(20) << pz		//28
					<< " "  << setw(20) << Den;		//29
				}	
				
				if(scalarevo)
				{
					
					double phii,Pi;

					phii=get_bv(l,k,j,nsc);
					Pi=get_bv(l,k,j,nscp);
					
					fout
					<< " "  << setw(20) << phii		//30	//25
					<< " "  << setw(20) << Pi;		//31	//26
				}
				
				fout << endl;
			}
			fout << endl;
		}
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output for 3d animation
void Fmv0::print_3d(ofstream& fout)
{
	double tt = get_t();
	
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##" 
		<< "time=" << tt 
		<< endl;

	for(int l=lli;l<=lui;l++)
	{
		double z = get_z(l);

		for(int k=kli;k<=kui;k++)
		{
			double y = get_y(k);

			for(int j=jli;j<=jui;j++)
			{
				double x = get_x(j);
				
				double fxx,fyy,fzz,gxx,gyy,gzz,gxy,gxz,gyz,wa;
				double akxx,akyy,akzz,akxy,akxz,akyz,ek;
				double zgx,zgy,zgz;
				double bx,by,bz,bbx,bby,bbz;
				double alp;

				fxx=get_flat_df2x(j);
				fyy=get_flat_df2y(k);
				fzz=get_flat_df2z(l);
			
				alp= get_bv(l,k,j, 0);
				bx =  get_bv(l,k,j, 1);
				by =  get_bv(l,k,j, 2);
				bz =  get_bv(l,k,j, 3);
				bbx=  get_bv(l,k,j, 4);
				bby=  get_bv(l,k,j, 5);
				bbz=  get_bv(l,k,j, 6);
				gxx=  fxx+get_bv(l,k,j, 7);
				gyy=  fyy+get_bv(l,k,j, 8);
				gzz=  fzz+get_bv(l,k,j, 9);
				gxy=  get_bv(l,k,j,10);
				gxz=  get_bv(l,k,j,11);
				gyz=  get_bv(l,k,j,12);
				wa=   get_bv(l,k,j,13);
				akxx= get_bv(l,k,j,14);
				akyy= get_bv(l,k,j,15);
				akzz= get_bv(l,k,j,16);
				akxy= get_bv(l,k,j,17);
				akxz= get_bv(l,k,j,18);
				akyz= get_bv(l,k,j,19);
				ek=   get_bv(l,k,j,20);
				zgx=  get_bv(l,k,j,21);
				zgy=  get_bv(l,k,j,22);
				zgz=  get_bv(l,k,j,23);

				fout.setf(ios_base::fixed, ios_base::floatfield);
				fout.precision(16);
				fout << setw(20) << x			//1
				<< " "  << setw(20) << y		//2
				<< " "  << setw(20) << z		//3
				<< " "  << setw(20) << alp 	 	//4
				<< " "  << setw(20) << bx   		//5
				<< " "  << setw(20) << by		//6
				<< " "  << setw(20) << bz		//7
				<< " "  << setw(20) << bbx  		//8
				<< " "  << setw(20) << bby		//9
				<< " "  << setw(20) << bbz		//10
				<< " "  << setw(20) << gxx		//11
				<< " "  << setw(20) << gyy		//12
				<< " "  << setw(20) << gzz		//13
				<< " "  << setw(20) << gxy		//14
				<< " "  << setw(20) << gxz		//15
				<< " "  << setw(20) << gyz		//16
				<< " "  << setw(20) << wa		//17
				<< " "  << setw(20) << akxx		//18
				<< " "  << setw(20) << akyy		//19
				<< " "  << setw(20) << akzz		//20
				<< " "  << setw(20) << akxy		//21
				<< " "  << setw(20) << akxz		//22
				<< " "  << setw(20) << akyz		//23
				<< " "  << setw(20) << ek		//24
				<< " "  << setw(20) << zgx		//25
				<< " "  << setw(20) << zgy		//26
				<< " "  << setw(20) << zgz;		//27
			
				if(fluidevo)
				{
					double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;

					Ene=  get_bv(l,k,j,24);
					px=   get_bv(l,k,j,25);
					py=   get_bv(l,k,j,26);
					pz=   get_bv(l,k,j,27);
					Den=  get_bv(l,k,j,28);
					
					rho=  get_primv(l,k,j,0);
					Vx=   get_primv(l,k,j,1);
					Vy=   get_primv(l,k,j,2);
					Vz=   get_primv(l,k,j,3);
					eps=  get_primv(l,k,j,4);
					
					fout
					<< " "  << setw(20) << Ene		//28
					<< " "  << setw(20) << px		//29
					<< " "  << setw(20) << py		//30
					<< " "  << setw(20) << pz		//31
					<< " "  << setw(20) << Den		//32
					<< " "  << setw(20) << rho		//33
					<< " "  << setw(20) << Vx		//34
					<< " "  << setw(20) << Vy		//35
					<< " "  << setw(20) << Vz		//36
					<< " "  << setw(20) << eps;		//37
					
				}
				
				if(scalarevo)
				{
					double phii,Pi;

					phii=get_bv(l,k,j,nsc);
					Pi=get_bv(l,k,j,nscp);
					
					fout
					<< " "  << setw(20) << phii		//38	//28
					<< " "  << setw(20) << Pi;		//39	//29
				}	
					
				fout << endl;
			}
			fout << endl;
		}
		fout << endl;
	}
	fout << endl << endl << endl;
	return;
}
