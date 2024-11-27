/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* BOUNDARY CONDITIONS :: BSSN evolution Class of COSMOS                                                 */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos.h"

//reflection boundary for scalar field
void Fmv::boundary_reflection_scalar()
{	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(l,k,jl[j],i)=get_bv(l,k,jli+j,i);
					set_bv(l,k,ju[j],i)=get_bv(l,k,jui-j,i);
				}
			}
		}
	}
	
	// y-fixed 
	//
	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(l,kl[k],j,i)=get_bv(l,kli+k,-j,i);
					set_bv(l,ku[k],j,i)=get_bv(l,kui-k,j,i);
				}
			}
		}
	}
	
	// z-fixed
	//
	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(ll[l],k,j,i)=get_bv(lli+l,k,j,i);
					set_bv(lu[l],k,j,i)=get_bv(lui-l,k,j,i);
				}
			}
		}
	}
	//x,y-fixed
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(l,kl[k],jl[j],i)=get_bv(l,kli+k,-jli-j,i);
					set_bv(l,ku[k],jl[j],i)=get_bv(l,kui-k,jli+j,i);
					set_bv(l,kl[k],ju[j],i)=get_bv(l,kli+k,-jui+j,i);
					set_bv(l,ku[k],ju[j],i)=get_bv(l,kui-k,jui-j,i);
				}
			}
		}
	}

	// x,z-fixed

	#pragma omp parallel for
	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(ll[l],k,jl[j],i)=get_bv(lli+l,k,jli+j,i);
					set_bv(ll[l],k,ju[j],i)=get_bv(lli+l,k,jui-j,i);
					set_bv(lu[l],k,jl[j],i)=get_bv(lui-l,k,jli+j,i);
					set_bv(lu[l],k,ju[j],i)=get_bv(lui-l,k,jui-j,i);
				}
			}
		}
	}
		
	// y,z-fixed

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(ll[l],kl[k],j,i)=get_bv(lli+l,kli+k,-j,i);
					set_bv(ll[l],ku[k],j,i)=get_bv(lli+l,kui-k,j,i);
					set_bv(lu[l],kl[k],j,i)=get_bv(lui-l,kli+k,-j,i);
					set_bv(lu[l],ku[k],j,i)=get_bv(lui-l,kui-k,j,i);
				}
			}
		}
	}

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(ll[l],kl[k],jl[j],i)=get_bv(lli+l,kli+k,-jli-j,i);
					set_bv(ll[l],kl[k],ju[j],i)=get_bv(lli+l,kli+k,-jui+j,i);
					set_bv(ll[l],ku[k],jl[j],i)=get_bv(lli+l,kui-k,jli+j,i);
					set_bv(ll[l],ku[k],ju[j],i)=get_bv(lli+l,kui-k,jui-j,i);
					set_bv(lu[l],kl[k],jl[j],i)=get_bv(lui-l,kli+k,-jli-j,i);
					set_bv(lu[l],kl[k],ju[j],i)=get_bv(lui-l,kli+k,-jui+j,i);
					set_bv(lu[l],ku[k],jl[j],i)=get_bv(lui-l,kui-k,jli+j,i);
					set_bv(lu[l],ku[k],ju[j],i)=get_bv(lui-l,kui-k,jui-j,i);
				}
			}
		}
	}
	
	return;
}

//reflection boundary for trK
void Fmv::boundary_reflection_even(int i)
{
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				set_bv(l,k,jl[j],i)=get_bv(l,k,jli+j,i);
				set_bv(l,k,ju[j],i)=get_bv(l,k,jui-j,i);
			}
		}
	}
	
	// y-fixed 
	//

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=1;k<=tab;k++)
			{
				set_bv(l,kl[k],j,i)=get_bv(l,kli+k,-j,i);
				set_bv(l,ku[k],j,i)=get_bv(l,kui-k,j,i);
			}
		}
	}
	
	// z-fixed
	//

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_bv(ll[l],k,j,i)=get_bv(lli+l,k,j,i);
				set_bv(lu[l],k,j,i)=get_bv(lui-l,k,j,i);
			}
		}
	}
	//x,y-fixed
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				set_bv(l,kl[k],jl[j],i)=get_bv(l,kli+k,-jli-j,i);
				set_bv(l,ku[k],jl[j],i)=get_bv(l,kui-k,jli+j,i);
				set_bv(l,kl[k],ju[j],i)=get_bv(l,kli+k,-jui+j,i);
				set_bv(l,ku[k],ju[j],i)=get_bv(l,kui-k,jui-j,i);
			}
		}
	}

	// x,z-fixed

	#pragma omp parallel for
	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_bv(ll[l],k,jl[j],i)=get_bv(lli+l,k,jli+j,i);
				set_bv(ll[l],k,ju[j],i)=get_bv(lli+l,k,jui-j,i);
				set_bv(lu[l],k,jl[j],i)=get_bv(lui-l,k,jli+j,i);
				set_bv(lu[l],k,ju[j],i)=get_bv(lui-l,k,jui-j,i);
			}
		}
	}
		
	// y,z-fixed

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_bv(ll[l],kl[k],j,i)=get_bv(lli+l,kli+k,-j,i);
				set_bv(ll[l],ku[k],j,i)=get_bv(lli+l,kui-k,j,i);
				set_bv(lu[l],kl[k],j,i)=get_bv(lui-l,kli+k,-j,i);
				set_bv(lu[l],ku[k],j,i)=get_bv(lui-l,kui-k,j,i);
			}
		}
	}

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_bv(ll[l],kl[k],jl[j],i)=get_bv(lli+l,kli+k,-jli-j,i);
				set_bv(ll[l],kl[k],ju[j],i)=get_bv(lli+l,kli+k,-jui+j,i);
				set_bv(ll[l],ku[k],jl[j],i)=get_bv(lli+l,kui-k,jli+j,i);
				set_bv(ll[l],ku[k],ju[j],i)=get_bv(lli+l,kui-k,jui-j,i);
				set_bv(lu[l],kl[k],jl[j],i)=get_bv(lui-l,kli+k,-jli-j,i);
				set_bv(lu[l],kl[k],ju[j],i)=get_bv(lui-l,kli+k,-jui+j,i);
				set_bv(lu[l],ku[k],jl[j],i)=get_bv(lui-l,kui-k,jli+j,i);
				set_bv(lu[l],ku[k],ju[j],i)=get_bv(lui-l,kui-k,jui-j,i);
			}
		}
	}
	
	return;
}



//reflection boundary for fluid
void Fmv::boundary_reflection_fluid()
{
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

	// x-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28
	//
	//		vector scalar:	  2,3,    5,6,     22,23,   26,27
	//		vector vector:	1,      4,      21,      25,
	//
	//		tensor scalar:	      12,        19
	//		tensor vector:	10,11,     17,18,
	
	for(int i=24;i<29;i++)
	 refsign[i]=1; 

	refsign[25]=-1;
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(l,k,jl[j],i)=refsign[i]*get_bv(l,k,jli+j,i);
					set_bv(l,k,ju[j],i)=refsign[i]*get_bv(l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
	// y-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28
	//
	//		vector scalar:	1,  3, 4,  6,  21,   23, 25,  27
	//		vector vector:	  2,     5,       22,       26,
	//
	//		tensor scalar:	   11,           18,
	//		tensor vector:	10,   12,     17,   19
	
	for(int i=24;i<29;i++)
    refsign[i]=1; 

	refsign[26]=-1;
    
	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-j,i);
					set_bv(l,ku[k],j,i)=refsign[i]*get_bv(l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// z-fixed
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28
	//
	//		vector scalar:	1,2,     4,5,      21,22,       25,26,  
	//		vector vector:	    3,       6,          23,          27
	//
	//		tensor scalar:	10,           17,
	//		tensor vector:	   11,12,        18,19

	for(int i=24;i<29;i++)
	 refsign[i]=1; 

	refsign[27]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(ll[l],k,j,i)=refsign[i]*get_bv(lli+l,k,j,i);
					set_bv(lu[l],k,j,i)=refsign[i]*get_bv(lui-l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//x,y-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28
	//
	//		vector positive:	    3,      6,        23,       27
	//		vector negative:	1,2,    4,5,    21,22,    25,26,
	//
	//		tensor positive:	10,         17,  
	//		tensor negative:	   11,12,      18,19

	for(int i=24;i<29;i++)
	refsign[i]=1;
    
	refsign[25]=-1;
	refsign[26]=-1;
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(l,kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-jli-j,i);
					set_bv(l,ku[k],jl[j],i)=refsign[i]*get_bv(l,kui-k,jli+j,i);
					set_bv(l,kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-jui+j,i);
					set_bv(l,ku[k],ju[j],i)=refsign[i]*get_bv(l,kui-k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// x,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28 
	//
	//		vector positive:	  2,       5,         22,       26,
	//		vector negative:	1,  3,   4,  6,    21,   23, 25,   27
	//
	//		tensor positive:	   11,         18,  
	//		tensor negative:	10,   12,   17,   19
	for(int i=24;i<29;i++)
	 refsign[i]=1; 

	refsign[25]=-1;
	refsign[27]=-1;

	#pragma omp parallel for
	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(ll[l],k,jl[j],i)=refsign[i]*get_bv(lli+l,k,jli+j,i);
					set_bv(ll[l],k,ju[j],i)=refsign[i]*get_bv(lli+l,k,jui-j,i);
					set_bv(lu[l],k,jl[j],i)=refsign[i]*get_bv(lui-l,k,jli+j,i);
					set_bv(lu[l],k,ju[j],i)=refsign[i]*get_bv(lui-l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
		
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28 
	//
	//		vector positive:	1,      4,      21,       25,
	//		vector negative:	  2,3,    5,6,     22,23,    26,27
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=24;i<29;i++)
    refsign[i]=1; 
    
	refsign[26]=-1;
	refsign[27]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-j,i);
					set_bv(ll[l],ku[k],j,i)=refsign[i]*get_bv(lli+l,kui-k,j,i);
					set_bv(lu[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(lui-l,kli+k,-j,i);
					set_bv(lu[l],ku[k],j,i)=refsign[i]*get_bv(lui-l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//		scalar(positive):	0,7,8,9,13,14,15,16,20,24,28 
	//
	//		vector(negative):	1,2,3,	4,5,6,	21,22,23, 25,26,27
	//
	//		tensor(positive):	10,11,12,	17,18,19
	for(int i=24;i<29;i++)
    refsign[i]=1; 

	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge alpha, \beta^i
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma_{ij}
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK 
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma^i

	refsign[25]=-1;
	refsign[26]=-1;
	refsign[27]=-1;

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=24;i<29;i++)
				{
					set_bv(ll[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-jli-j,i);
					set_bv(ll[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-jui+j,i);
					set_bv(ll[l],ku[k],jl[j],i)=refsign[i]*get_bv(lli+l,kui-k,jli+j,i);
					set_bv(ll[l],ku[k],ju[j],i)=refsign[i]*get_bv(lli+l,kui-k,jui-j,i);
					set_bv(lu[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_bv(lui-l,kli+k,-jli-j,i);
					set_bv(lu[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_bv(lui-l,kli+k,-jui+j,i);
					set_bv(lu[l],ku[k],jl[j],i)=refsign[i]*get_bv(lui-l,kui-k,jli+j,i);
					set_bv(lu[l],ku[k],ju[j],i)=refsign[i]*get_bv(lui-l,kui-k,jui-j,i);
				}
			}
		}
	}
	
	return;
}

//reflection boundary for geometry + fluid (+ scalar)
void Fmv::boundary_reflection()
{
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
	
	// x-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	  2,3,    5,6,     22,23,   26,27
	//		vector vector:	1,      4,      21,      25,
	//
	//		tensor scalar:	      12,        19
	//		tensor vector:	10,11,     17,18,
	
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[4]=-1;
	refsign[10]=-1;
	refsign[11]=-1;
	refsign[17]=-1;
	refsign[18]=-1;
	refsign[21]=-1;

	if(fluidevo)
	refsign[25]=-1;
	
	#pragma omp parallel for
	for(int j=1;j<=tab;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=kli;k<=kui;k++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(l,k,jl[j],i)=refsign[i]*get_bv(l,k,jli+j,i);
					set_bv(l,k,ju[j],i)=refsign[i]*get_bv(l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier


	// y-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	1,  3, 4,  6,  21,   23, 25,  27
	//		vector vector:	  2,     5,       22,       26,
	//
	//		tensor scalar:	   11,           18,
	//		tensor vector:	10,   12,     17,   19
	
	for(int i=0;i<nn;i++)
	refsign[i]=1; 
    
	refsign[2]=-1;
	refsign[5]=-1;
	refsign[22]=-1;
	refsign[10]=-1;
	refsign[12]=-1;
	refsign[17]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[26]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-j,i);
					set_bv(l,ku[k],j,i)=refsign[i]*get_bv(l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// z-fixed
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	1,2,     4,5,      21,22,       25,26,  
	//		vector vector:	    3,       6,          23,          27
	//
	//		tensor scalar:	10,           17,
	//		tensor vector:	   11,12,        18,19

	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[3]=-1;
	refsign[6]=-1;
	refsign[23]=-1;
	refsign[11]=-1;
	refsign[12]=-1;
	refsign[18]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[27]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(ll[l],k,j,i)=refsign[i]*get_bv(lli+l,k,j,i);
					set_bv(lu[l],k,j,i)=refsign[i]*get_bv(lui-l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier

	//x,y-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	    3,      6,        23,       27
	//		vector negative:	1,2,    4,5,    21,22,    25,26,
	//
	//		tensor positive:	10,         17,  
	//		tensor negative:	   11,12,      18,19

	for(int i=0;i<nn;i++)
    refsign[i]=1; 
    
	refsign[1]=-1;
	refsign[2]=-1;
	refsign[4]=-1;
	refsign[5]=-1;
	refsign[21]=-1;
	refsign[22]=-1;
	refsign[11]=-1;
	refsign[12]=-1;
	refsign[18]=-1;
	refsign[19]=-1;

	if(fluidevo)
	{
		refsign[25]=-1;
		refsign[26]=-1;
	}

	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(l,kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-jli-j,i);
					set_bv(l,ku[k],jl[j],i)=refsign[i]*get_bv(l,kui-k,jli+j,i);
					set_bv(l,kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-jui+j,i);
					set_bv(l,ku[k],ju[j],i)=refsign[i]*get_bv(l,kui-k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// x,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	  2,       5,         22,       26,
	//		vector negative:	1,  3,   4,  6,    21,   23, 25,   27
	//
	//		tensor positive:	   11,         18,  
	//		tensor negative:	10,   12,   17,   19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[3]=-1;
	refsign[4]=-1;
	refsign[6]=-1;
	refsign[21]=-1;
	refsign[23]=-1;
	refsign[10]=-1;
	refsign[12]=-1;
	refsign[17]=-1;
	refsign[19]=-1;

	if(fluidevo)
	{
		refsign[25]=-1;
		refsign[27]=-1;
	}

	#pragma omp parallel for
	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(ll[l],k,jl[j],i)=refsign[i]*get_bv(lli+l,k,jli+j,i);
					set_bv(ll[l],k,ju[j],i)=refsign[i]*get_bv(lli+l,k,jui-j,i);
					set_bv(lu[l],k,jl[j],i)=refsign[i]*get_bv(lui-l,k,jli+j,i);
					set_bv(lu[l],k,ju[j],i)=refsign[i]*get_bv(lui-l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
		
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	1,      4,      21,       25,
	//		vector negative:	  2,3,    5,6,     22,23,    26,27
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nn;i++)
	{
        refsign[i]=1;
    }

	refsign[2]=-1;
	refsign[3]=-1;
	refsign[5]=-1;
	refsign[6]=-1;
	refsign[22]=-1;
	refsign[23]=-1;
	refsign[10]=-1;
	refsign[11]=-1;
	refsign[17]=-1;
	refsign[18]=-1;

	if(fluidevo)
	{
		refsign[26]=-1;
		refsign[27]=-1;
	}

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-j,i);
					set_bv(ll[l],ku[k],j,i)=refsign[i]*get_bv(lli+l,kui-k,j,i);
					set_bv(lu[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(lui-l,kli+k,-j,i);
					set_bv(lu[l],ku[k],j,i)=refsign[i]*get_bv(lui-l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//		scalar(positive):	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector(negative):	1,2,3,	4,5,6,	21,22,23, 25,26,27
	//
	//		tensor(positive):	10,11,12,	17,18,19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge alpha, \beta^i
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma_{ij}
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK 
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma^i

	refsign[1]=-1;
	refsign[2]=-1;
	refsign[3]=-1;
	refsign[4]=-1;
	refsign[5]=-1;
	refsign[6]=-1;
	refsign[21]=-1;
	refsign[22]=-1;
	refsign[23]=-1;

	if(fluidevo)
	{
		refsign[25]=-1;
		refsign[26]=-1;
		refsign[27]=-1;
	}

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_bv(ll[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-jli-j,i);
					set_bv(ll[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-jui+j,i);
					set_bv(ll[l],ku[k],jl[j],i)=refsign[i]*get_bv(lli+l,kui-k,jli+j,i);
					set_bv(ll[l],ku[k],ju[j],i)=refsign[i]*get_bv(lli+l,kui-k,jui-j,i);
					set_bv(lu[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_bv(lui-l,kli+k,-jli-j,i);
					set_bv(lu[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_bv(lui-l,kli+k,-jui+j,i);
					set_bv(lu[l],ku[k],jl[j],i)=refsign[i]*get_bv(lui-l,kui-k,jli+j,i);
					set_bv(lu[l],ku[k],ju[j],i)=refsign[i]*get_bv(lui-l,kui-k,jui-j,i);
				}
			}
		}
	}
	
	return;
}


//reflection boundary for fluid primitive variable
void Fmv::boundary_prim_reflection()
{
	// @@@@@@@@@@@@   primitive variables for fluid   @@@@@@@@@@@
	//  0:rho , 1:V^x  ,  2:V^y  ,  3:V^z   ,  4:epsilon

	// x-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(l,k,jl[j],i)=refsign[i]*get_primv(l,k,jli+j,i);
					set_primv(l,k,ju[j],i)=refsign[i]*get_primv(l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
	
	// y-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[2]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_primv(l,kli+k,-j,i);
					set_primv(l,ku[k],j,i)=refsign[i]*get_primv(l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	
	// z-fixed
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		

	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[3]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(ll[l],k,j,i)=refsign[i]*get_primv(lli+l,k,j,i);
					set_primv(lu[l],k,j,i)=refsign[i]*get_primv(lui-l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//x,y-fixed
	//

	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[2]=-1;
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(l,kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_primv(l,kli+k,-jli-j,i);
					set_primv(l,ku[k],jl[j],i)=refsign[i]*get_primv(l,kui-k,jli+j,i);
					set_primv(l,kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_primv(l,kli+k,-jui+j,i);
					set_primv(l,ku[k],ju[j],i)=refsign[i]*get_primv(l,kui-k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// x,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	  2,       5,         22,       28,
	//		vector negative:	1,  3,   4,  6,    21,   23, 27,   29
	//
	//		tensor positive:	   11,         18,  
	//		tensor negative:	10,   12,   17,   19
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[3]=-1;

	#pragma omp parallel for
	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(ll[l],k,jl[j],i)=refsign[i]*get_primv(lli+l,k,jli+j,i);
					set_primv(ll[l],k,ju[j],i)=refsign[i]*get_primv(lli+l,k,jui-j,i);
					set_primv(lu[l],k,jl[j],i)=refsign[i]*get_primv(lui-l,k,jli+j,i);
					set_primv(lu[l],k,ju[j],i)=refsign[i]*get_primv(lui-l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
		
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	1,      4,      21,       27,
	//		vector negative:	  2,3,    5,6,     22,23,    28,29
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[3]=-1;
	
	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_primv(lli+l,kli+k,-j,i);
					set_primv(ll[l],ku[k],j,i)=refsign[i]*get_primv(lli+l,kui-k,j,i);
					set_primv(lu[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_primv(lui-l,kli+k,-j,i);
					set_primv(lu[l],ku[k],j,i)=refsign[i]*get_primv(lui-l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//		scalar(positive):	0,7,8,9,13,14,15,16,20 
	//
	//		vector(negative):	1,2,3,	4,5,6,	21,22,23, 27,28,29
	//
	//		tensor(positive):	10,11,12,	17,18,19
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge alpha, \beta^i
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma_{ij}
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK 
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma^i

	refsign[1]=-1;
	refsign[2]=-1;
	refsign[3]=-1;

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<npr;i++)
				{
					set_primv(ll[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_primv(lli+l,kli+k,-jli-j,i);
					set_primv(ll[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_primv(lli+l,kli+k,-jui+j,i);
					set_primv(ll[l],ku[k],jl[j],i)=refsign[i]*get_primv(lli+l,kui-k,jli+j,i);
					set_primv(ll[l],ku[k],ju[j],i)=refsign[i]*get_primv(lli+l,kui-k,jui-j,i);
					set_primv(lu[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_primv(lui-l,kli+k,-jli-j,i);
					set_primv(lu[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_primv(lui-l,kli+k,-jui+j,i);
					set_primv(lu[l],ku[k],jl[j],i)=refsign[i]*get_primv(lui-l,kui-k,jli+j,i);
					set_primv(lu[l],ku[k],ju[j],i)=refsign[i]*get_primv(lui-l,kui-k,jui-j,i);
				}
			}
		}
	}

	return;
}

//reflection boundary for dbv of geometry + fluid (+ scalar)
void Fmv::boundary_d_reflection()
{
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

	// x-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	  2,3,    5,6,     22,23,    26,27
	//		vector vector:	1,      4,      21,       25,
	//
	//		tensor scalar:	      12,        19
	//		tensor vector:	10,11,     17,18,
	
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[4]=-1;
	refsign[10]=-1;
	refsign[11]=-1;
	refsign[17]=-1;
	refsign[18]=-1;
	refsign[21]=-1;

	if(fluidevo)
	refsign[25]=-1;
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(l,k,jl[j],i)=refsign[i]*get_dbv(l,k,jli+j,i);
					set_dbv(l,k,ju[j],i)=refsign[i]*get_dbv(l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
	
	// y-fixed 
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector scalar:	1,  3,   4,  6,    21,   23, 25,   27
	//		vector vector:	  2,       5,         22,       26,
	//
	//		tensor scalar:	   11,         18,  
	//		tensor vector:	10,   12,   17,   19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[5]=-1;
	refsign[22]=-1;
	refsign[10]=-1;
	refsign[12]=-1;
	refsign[17]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[26]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int l=lli;l<=lui;l++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_dbv(l,kli+k,-j,i);
					set_dbv(l,ku[k],j,i)=refsign[i]*get_dbv(l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	
	// z-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	1,2,    4,5,    21,22,    25,26,
	//		vector vector:	    3,      6,        23,       27
	//
	//		tensor scalar:	10,         17,  
	//		tensor vector:	   11,12,      18,19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[3]=-1;
	refsign[6]=-1;
	refsign[23]=-1;
	refsign[11]=-1;
	refsign[12]=-1;
	refsign[18]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[27]=-1;

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(ll[l],k,j,i)=refsign[i]*get_dbv(lli+l,k,j,i);
					set_dbv(lu[l],k,j,i)=refsign[i]*get_dbv(lui-l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//x,y-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	    3,      6,        23,       27
	//		vector negative:	1,2,    4,5,    21,22,    25,26,
	//
	//		tensor positive:	10,         17,  
	//		tensor negative:	   11,12,      18,19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[2]=-1;
	refsign[4]=-1;
	refsign[5]=-1;
	refsign[21]=-1;
	refsign[22]=-1;
	refsign[11]=-1;
	refsign[12]=-1;
	refsign[18]=-1;
	refsign[19]=-1;

	if(fluidevo)
	{
		refsign[25]=-1;
		refsign[26]=-1;
	}
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(l,kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_dbv(l,kli+k,-jli-j,i);
					set_dbv(l,ku[k],jl[j],i)=refsign[i]*get_dbv(l,kui-k,jli+j,i);
					set_dbv(l,kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_dbv(l,kli+k,-jui+j,i);
					set_dbv(l,ku[k],ju[j],i)=refsign[i]*get_dbv(l,kui-k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// x,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	  2,       5,         22,       26,
	//		vector negative:	1,  3,   4,  6,    21,   23, 25,   27
	//
	//		tensor positive:	   11,         18,  
	//		tensor negative:	10,   12,   17,   19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[3]=-1;
	refsign[4]=-1;
	refsign[6]=-1;
	refsign[21]=-1;
	refsign[23]=-1;
	refsign[10]=-1;
	refsign[12]=-1;
	refsign[17]=-1;
	refsign[19]=-1;

	if(fluidevo)
	{
		refsign[25]=-1;
		refsign[27]=-1;
	}

	#pragma omp parallel for
	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(ll[l],k,jl[j],i)=refsign[i]*get_dbv(lli+l,k,jli+j,i);
					set_dbv(ll[l],k,ju[j],i)=refsign[i]*get_dbv(lli+l,k,jui-j,i);
					set_dbv(lu[l],k,jl[j],i)=refsign[i]*get_dbv(lui-l,k,jli+j,i);
					set_dbv(lu[l],k,ju[j],i)=refsign[i]*get_dbv(lui-l,k,jui-j,i);
				}
			}
		}
	}
	#pragma omp barrier
		
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	1,      4,      21,       25,
	//		vector negative:	  2,3,    5,6,     22,23,    26,27
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[3]=-1;
	refsign[5]=-1;
	refsign[6]=-1;
	refsign[22]=-1;
	refsign[23]=-1;
	refsign[10]=-1;
	refsign[11]=-1;
	refsign[17]=-1;
	refsign[18]=-1;

	if(fluidevo)
	{
		refsign[26]=-1;
		refsign[27]=-1;
	}

	#pragma omp parallel for
	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_dbv(lli+l,kli+k,-j,i);
					set_dbv(ll[l],ku[k],j,i)=refsign[i]*get_dbv(lli+l,kui-k,j,i);
					set_dbv(lu[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_dbv(lui-l,kli+k,-j,i);
					set_dbv(lu[l],ku[k],j,i)=refsign[i]*get_dbv(lui-l,kui-k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	//		scalar(positive):	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector(negative):	1,2,3,	4,5,6,	21,22,23, 25,26,27
	//
	//		tensor(positive):	10,11,12,	17,18,19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[2]=-1;
	refsign[3]=-1;
	refsign[4]=-1;
	refsign[5]=-1;
	refsign[6]=-1;
	refsign[21]=-1;
	refsign[22]=-1;
	refsign[23]=-1;

	if(fluidevo)
	{
		refsign[25]=-1;
		refsign[26]=-1;
		refsign[27]=-1;
	}

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nn;i++)
				{
					set_dbv(ll[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_dbv(lli+l,kli+k,-jli-j,i);
					set_dbv(ll[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_dbv(lli+l,kli+k,-jui+j,i);
					set_dbv(ll[l],ku[k],jl[j],i)=refsign[i]*get_dbv(lli+l,kui-k,jli+j,i);
					set_dbv(ll[l],ku[k],ju[j],i)=refsign[i]*get_dbv(lli+l,kui-k,jui-j,i);
					set_dbv(lu[l],kl[k],jl[j],i)=rotsign[i]*refsign[i]*get_dbv(lui-l,kli+k,-jli-j,i);
					set_dbv(lu[l],kl[k],ju[j],i)=rotsign[i]*refsign[i]*get_dbv(lui-l,kli+k,-jui+j,i);
					set_dbv(lu[l],ku[k],jl[j],i)=refsign[i]*get_dbv(lui-l,kui-k,jli+j,i);
					set_dbv(lu[l],ku[k],ju[j],i)=refsign[i]*get_dbv(lui-l,kui-k,jui-j,i);
				}
			}
		}
	}
	//cout << "boundary reflection" << endl;

	return;
}


//periodic boundary condition
void Fmv::boundary_periodic()
{

	//BSSN
	for(int i=0;i<nn;i++){ 

		// x-fixed
		/////////////////////////////////////
		//     0 1 2 3 4 5 6 7 8 9 
		//     B B B I I I I B B B
		//     3 <=> 6
		//     2 <= 5, 1 <= 4, 0 <= 3  || jl[j] <= jui-j
		//     7 <= 4, 8 <= 5, 9 <= 6  || ju[j] <= jli+j

		for(int l=lli;l<=lui;l++){
			for(int k=kli;k<=kui;k++){
				for(int j=1;j<=tab;j++){
					set_bv(l,k,jl[j],i)=get_bv(l,k,jui-j,i);
					set_bv(l,k,ju[j],i)=get_bv(l,k,jli+j,i);
				}
			}
		}
		// y-fixed
		for(int l=lli;l<=lui;l++){
			for(int j=jli;j<=jui;j++){
				for(int k=1;k<=tab;k++){
					set_bv(l,kl[k],j,i)=get_bv(l,kui-k,j,i);
					set_bv(l,ku[k],j,i)=get_bv(l,kli+k,j,i);
				}
			}
		}
		// z-fixed
		for(int k=kli;k<=kui;k++){
			for(int j=jli;j<=jui;j++){
				for(int l=1;l<=tab;l++){
					set_bv(ll[l],k,j,i)=get_bv(lui-l,k,j,i);
					set_bv(lu[l],k,j,i)=get_bv(lli+l,k,j,i);
				}
			}
		}

		for(int l=lli;l<=lui;l++){
			for(int k=1;k<=tab;k++){
				for(int j=1;j<=tab;j++){
					set_bv(l,kl[k],jl[j],i)=get_bv(l,kui-k,jui+j,i);
					set_bv(l,ku[k],jl[j],i)=get_bv(l,kli+k,jui+j,i);
					set_bv(l,kl[k],ju[j],i)=get_bv(l,kui-k,jli-j,i);
					set_bv(l,ku[k],ju[j],i)=get_bv(l,kli+k,jli-j,i);
				}
			}
		}

		// z-fixed
		for(int k=kli;k<=kui;k++){
			for(int j=1;j<=tab;j++){
				for(int l=1;l<=tab;l++){
					set_bv(ll[l],k,jl[j],i)=get_bv(lui+l,k,jui-j,i);
					set_bv(ll[l],k,ju[j],i)=get_bv(lui+l,k,jli+j,i);
					set_bv(lu[l],k,jl[j],i)=get_bv(lli-l,k,jui-j,i);
					set_bv(lu[l],k,ju[j],i)=get_bv(lli-l,k,jli+j,i);
				}
			}
		}
		for(int j=jli;j<=jui;j++){
			for(int k=1;k<=tab;k++){
				for(int l=1;l<=tab;l++){
					set_bv(ll[l],kl[k],j,i)=get_bv(lui+l,kui-k,j,i);
					set_bv(ll[l],ku[k],j,i)=get_bv(lui+l,kli+k,j,i);
					set_bv(lu[l],kl[k],j,i)=get_bv(lli-l,kui-k,j,i);
					set_bv(lu[l],ku[k],j,i)=get_bv(lli-l,kli+k,j,i);
				}
			}
		}

		// z-fixed
		for(int j=1;j<=tab;j++){
			for(int k=1;k<=tab;k++){
				for(int l=1;l<=tab;l++){
					set_bv(ll[l],kl[k],jl[j],i)=get_bv(lui-l,kui-k,jui-j,i);
					set_bv(ll[l],kl[k],ju[j],i)=get_bv(lui-l,kui-k,jli+j,i);
					set_bv(ll[l],ku[k],jl[j],i)=get_bv(lui-l,kli+k,jui-j,i);
					set_bv(ll[l],ku[k],ju[j],i)=get_bv(lui-l,kli+k,jli+j,i);
					set_bv(lu[l],kl[k],jl[j],i)=get_bv(lli+l,kui-k,jui-j,i);
					set_bv(lu[l],kl[k],ju[j],i)=get_bv(lli+l,kui-k,jli+j,i);
					set_bv(lu[l],ku[k],jl[j],i)=get_bv(lli+l,kli+k,jui-j,i);
					set_bv(lu[l],ku[k],ju[j],i)=get_bv(lli+l,kli+k,jli+j,i);
				}
			}
		}
	}

	return;
}

//reflection boundary condition for w in momentum constraint solver
void Fmv::boundary_w()
{
	short int *wrefsign;
	
	// @@@@@@@@@@@@   variables (BSSN +GAUGE)   @@@@@@@@@@@
	//  0:wx  ,  1:wy  ,  2:wz 
	
	int nw=3;
	
	wrefsign = new short int[nw];
	
	//cout << "boundary reflection" << endl;
	
	// x-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//          |      |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j-1
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	  2,3,    5,6,     22,23
	//		vector vector:	1,      4,      21,
	//
	//		tensor scalar:	      12,        19
	//		tensor vector:	10,11,     17,18,
	
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[0]=-1;
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(l,k,jl[j],i)=wrefsign[i]*get_wvec(l,k,jli+j,i);
					set_wvec(l,k,ju[j],i)=wrefsign[i]*get_wvec(l,k,jui-j,i);
				}
			}
		}
	}
	
	// y-fixed 
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	1,  3,   4,  6,    21,   23
	//		vector vector:	  2,       5,         22,
	//
	//		tensor scalar:	   11,         18,  
	//		tensor vector:	10,   12,   17,   19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[1]=-1;
	wrefsign[0]=-1;

	for(int l=lli;l<=lui;l++)
	{
		for(int j=jli;j<=jui;j++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(l,kl[k],j,i)=wrefsign[i]*get_wvec(l,kli+k,-j,i);
					set_wvec(l,ku[k],j,i)=wrefsign[i]*get_wvec(l,kui-k,j,i);
				}
			}
		}
	}
	
	// z-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	1,2,    4,5,    21,22,
	//		vector vector:	    3,      6,        23
	//
	//		tensor scalar:	10,         17,  
	//		tensor vector:	   11,12,      18,19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[2]=-1;

	for(int k=kli;k<=kui;k++)
	{
		for(int j=jli;j<=jui;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(ll[l],k,j,i)=wrefsign[i]*get_wvec(lli+l,k,j,i);
					set_wvec(lu[l],k,j,i)=wrefsign[i]*get_wvec(lui-l,k,j,i);
				}
			}
		}
	}
	//x,y-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	    3,      6,        23
	//		vector negative:	1,2,    4,5,    21,22,
	//
	//		tensor positive:	10,         17,  
	//		tensor negative:	   11,12,      18,19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[1]=-1;
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(l,kl[k],jl[j],i)=wrefsign[i]*get_wvec(l,kli+k,-jli-j,i);
					set_wvec(l,ku[k],jl[j],i)=wrefsign[i]*get_wvec(l,kui-k,jli+j,i);
					set_wvec(l,kl[k],ju[j],i)=wrefsign[i]*get_wvec(l,kli+k,-jui+j,i);
					set_wvec(l,ku[k],ju[j],i)=wrefsign[i]*get_wvec(l,kui-k,jui-j,i);
				}
			}
		}
	}

	// x,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	  2,       5,         22,
	//		vector negative:	1,  3,   4,  6,    21,   23
	//
	//		tensor positive:	   11,         18,  
	//		tensor negative:	10,   12,   17,   19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[0]=-1;
	wrefsign[2]=-1;

	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(ll[l],k,jl[j],i)=wrefsign[i]*get_wvec(lli+l,k,jli+j,i);
					set_wvec(ll[l],k,ju[j],i)=wrefsign[i]*get_wvec(lli+l,k,jui-j,i);
					set_wvec(lu[l],k,jl[j],i)=wrefsign[i]*get_wvec(lui-l,k,jli+j,i);
					set_wvec(lu[l],k,ju[j],i)=wrefsign[i]*get_wvec(lui-l,k,jui-j,i);
				}
			}
		}
	}
		
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	1,      4,      21,
	//		vector negative:	  2,3,    5,6,     22,23
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[0]=-1;
	wrefsign[1]=-1;
	wrefsign[2]=-1;

	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(ll[l],kl[k],j,i)=wrefsign[i]*get_wvec(lli+l,kli+k,-j,i);
					set_wvec(ll[l],ku[k],j,i)=wrefsign[i]*get_wvec(lli+l,kui-k,j,i);
					set_wvec(lu[l],kl[k],j,i)=wrefsign[i]*get_wvec(lui-l,kli+k,-j,i);
					set_wvec(lu[l],ku[k],j,i)=wrefsign[i]*get_wvec(lui-l,kui-k,j,i);
				}
			}
		}
	}
	//		scalar(positive):	0,7,8,9,13,14,15,16,20 
	//
	//		vector(negative):	1,2,3,	4,5,6,	21,22,23
	//
	//		tensor(positive):	10,11,12,	17,18,19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge alpha, \beta^i
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma_{ij}
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK 
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma^i

	wrefsign[1]=-1;
	wrefsign[2]=-1;

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_wvec(ll[l],kl[k],jl[j],i)=wrefsign[i]*get_wvec(lli+l,kli+k,-jli-j,i);
					set_wvec(ll[l],kl[k],ju[j],i)=wrefsign[i]*get_wvec(lli+l,kli+k,-jui+j,i);
					set_wvec(ll[l],ku[k],jl[j],i)=wrefsign[i]*get_wvec(lli+l,kui-k,jli+j,i);
					set_wvec(ll[l],ku[k],ju[j],i)=wrefsign[i]*get_wvec(lli+l,kui-k,jui-j,i);
					set_wvec(lu[l],kl[k],jl[j],i)=wrefsign[i]*get_wvec(lui-l,kli+k,-jli-j,i);
					set_wvec(lu[l],kl[k],ju[j],i)=wrefsign[i]*get_wvec(lui-l,kli+k,-jui+j,i);
					set_wvec(lu[l],ku[k],jl[j],i)=wrefsign[i]*get_wvec(lui-l,kui-k,jli+j,i);
					set_wvec(lu[l],ku[k],ju[j],i)=wrefsign[i]*get_wvec(lui-l,kui-k,jui-j,i);
				}
			}
		}
	}
	//cout << "boundary reflection" << endl;

	delete[] wrefsign;	
	return;
}

//reflection boundary condition for shift
void Fmv::boundary_beta()
{
	short int *wrefsign;
	
	// @@@@@@@@@@@@   variables (BSSN +GAUGE)   @@@@@@@@@@@
	//  0:wx  ,  1:wy  ,  2:wz 
	
	int nw=3;
	
	wrefsign = new short int[nw];
	
	//cout << "boundary reflection" << endl;
	
	// x-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//          |      |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j-1
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	  2,3,    5,6,     22,23
	//		vector vector:	1,      4,      21,
	//
	//		tensor scalar:	      12,        19
	//		tensor vector:	10,11,     17,18,
	
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[0]=-1;
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(l,k,jl[j],i+1)=wrefsign[i]*get_bv(l,k,jli+j,i+1);
					set_bv(l,k,ju[j],i+1)=wrefsign[i]*get_bv(l,k,jui-j,i+1);
				}
			}
		}
	}
	
	// y-fixed 
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	1,  3,   4,  6,    21,   23
	//		vector vector:	  2,       5,         22,
	//
	//		tensor scalar:	   11,         18,  
	//		tensor vector:	10,   12,   17,   19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[1]=-1;

	for(int l=lli;l<=lui;l++)
	{
		for(int j=jli;j<=jui;j++)
		{
			for(int k=1;k<=tab;k++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(l,kl[k],j,i+1)=wrefsign[i]*get_bv(l,kli+k,-j,i+1);
					set_bv(l,ku[k],j,i+1)=wrefsign[i]*get_bv(l,kui-k,j,i+1);
				}
			}
		}
	}
	
	// z-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	1,2,    4,5,    21,22,
	//		vector vector:	    3,      6,        23
	//
	//		tensor scalar:	10,         17,  
	//		tensor vector:	   11,12,      18,19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[2]=-1;

	for(int k=kli;k<=kui;k++)
	{
		for(int j=jli;j<=jui;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(ll[l],k,j,i+1)=wrefsign[i]*get_bv(lli+l,k,j,i+1);
					set_bv(lu[l],k,j,i+1)=wrefsign[i]*get_bv(lui-l,k,j,i+1);
				}
			}
		}
	}
	//x,y-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	    3,      6,        23
	//		vector negative:	1,2,    4,5,    21,22,
	//
	//		tensor positive:	10,         17,  
	//		tensor negative:	   11,12,      18,19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[0]=-1;
	wrefsign[1]=-1;
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(l,kl[k],jl[j],i+1)=wrefsign[i]*get_bv(l,kli+k,-jli-j,i+1);
					set_bv(l,ku[k],jl[j],i+1)=wrefsign[i]*get_bv(l,kui-k,jli+j,i+1);
					set_bv(l,kl[k],ju[j],i+1)=wrefsign[i]*get_bv(l,kli+k,-jui+j,i+1);
					set_bv(l,ku[k],ju[j],i+1)=wrefsign[i]*get_bv(l,kui-k,jui-j,i+1);
				}
			}
		}
	}

	// x,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	  2,       5,         22,
	//		vector negative:	1,  3,   4,  6,    21,   23
	//
	//		tensor positive:	   11,         18,  
	//		tensor negative:	10,   12,   17,   19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[0]=-1;
	wrefsign[2]=-1;

	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(ll[l],k,jl[j],i+1)=wrefsign[i]*get_bv(lli+l,k,jli+j,i+1);
					set_bv(ll[l],k,ju[j],i+1)=wrefsign[i]*get_bv(lli+l,k,jui-j,i+1);
					set_bv(lu[l],k,jl[j],i+1)=wrefsign[i]*get_bv(lui-l,k,jli+j,i+1);
					set_bv(lu[l],k,ju[j],i+1)=wrefsign[i]*get_bv(lui-l,k,jui-j,i+1);
				}
			}
		}
	}
		
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	1,      4,      21,
	//		vector negative:	  2,3,    5,6,     22,23
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	wrefsign[1]=-1;
	wrefsign[2]=-1;

	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(ll[l],kl[k],j,i+1)=wrefsign[i]*get_bv(lli+l,kli+k,-j,i+1);
					set_bv(ll[l],ku[k],j,i+1)=wrefsign[i]*get_bv(lli+l,kui-k,j,i+1);
					set_bv(lu[l],kl[k],j,i+1)=wrefsign[i]*get_bv(lui-l,kli+k,-j,i+1);
					set_bv(lu[l],ku[k],j,i+1)=wrefsign[i]*get_bv(lui-l,kui-k,j,i+1);
				}
			}
		}
	}
	//		scalar(positive):	0,7,8,9,13,14,15,16,20 
	//
	//		vector(negative):	1,2,3,	4,5,6,	21,22,23
	//
	//		tensor(positive):	10,11,12,	17,18,19
	for(int i=0;i<nw;i++)
	 wrefsign[i]=1; 

	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge alpha, \beta^i
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma_{ij}
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK 
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma^i

	wrefsign[0]=-1;
	wrefsign[1]=-1;
	wrefsign[2]=-1;

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				for(int i=0;i<nw;i++)
				{
					set_bv(ll[l],kl[k],jl[j],i+1)=wrefsign[i]*get_bv(lli+l,kli+k,-jli-j,i+1);
					set_bv(ll[l],kl[k],ju[j],i+1)=wrefsign[i]*get_bv(lli+l,kli+k,-jui+j,i+1);
					set_bv(ll[l],ku[k],jl[j],i+1)=wrefsign[i]*get_bv(lli+l,kui-k,jli+j,i+1);
					set_bv(ll[l],ku[k],ju[j],i+1)=wrefsign[i]*get_bv(lli+l,kui-k,jui-j,i+1);
					set_bv(lu[l],kl[k],jl[j],i+1)=wrefsign[i]*get_bv(lui-l,kli+k,-jli-j,i+1);
					set_bv(lu[l],kl[k],ju[j],i+1)=wrefsign[i]*get_bv(lui-l,kli+k,-jui+j,i+1);
					set_bv(lu[l],ku[k],jl[j],i+1)=wrefsign[i]*get_bv(lui-l,kui-k,jli+j,i+1);
					set_bv(lu[l],ku[k],ju[j],i+1)=wrefsign[i]*get_bv(lui-l,kui-k,jui-j,i+1);
				}
			}
		}
	}
	//cout << "boundary reflection" << endl;

	delete[] wrefsign;	
	return;
}

//reflection boundary condition for psi
void Fmv::boundary_psi_initial()
{
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=kli;k<=kui;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				set_psi(l,k,jl[j])=get_psi(l,k,jli+j);
				set_psi(l,k,ju[j])=get_psi(l,k,jui-j);
			}
		}
	}
	
	// y-fixed 
	//

	for(int l=lli;l<=lui;l++)
	{
		for(int j=jli;j<=jui;j++)
		{
			for(int k=1;k<=tab;k++)
			{
				set_psi(l,kl[k],j)=get_psi(l,kli+k,-j);
				set_psi(l,ku[k],j)=get_psi(l,kui-k,j);
			}
		}
	}
	
	// z-fixed
	//

	for(int k=kli;k<=kui;k++)
	{
		for(int j=jli;j<=jui;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_psi(ll[l],k,j)=get_psi(lli+l,k,j);
				set_psi(lu[l],k,j)=get_psi(lui-l,k,j);
			}
		}
	}
	//x,y-fixed
	
	for(int l=lli;l<=lui;l++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int j=1;j<=tab;j++)
			{
				set_psi(l,kl[k],jl[j])=get_psi(l,kli+k,-jli-j);
				set_psi(l,ku[k],jl[j])=get_psi(l,kui-k,jli+j);
				set_psi(l,kl[k],ju[j])=get_psi(l,kli+k,-jui+j);
				set_psi(l,ku[k],ju[j])=get_psi(l,kui-k,jui-j);
			}
		}
	}

	// x,z-fixed

	for(int k=kli;k<=kui;k++)
	{
		for(int j=1;j<=tab;j++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_psi(ll[l],k,jl[j])=get_psi(lli+l,k,jli+j);
				set_psi(ll[l],k,ju[j])=get_psi(lli+l,k,jui-j);
				set_psi(lu[l],k,jl[j])=get_psi(lui-l,k,jli+j);
				set_psi(lu[l],k,ju[j])=get_psi(lui-l,k,jui-j);
			}
		}
	}
		
	// y,z-fixed

	for(int j=jli;j<=jui;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_psi(ll[l],kl[k],j)=get_psi(lli+l,kli+k,-j);
				set_psi(ll[l],ku[k],j)=get_psi(lli+l,kui-k,j);
				set_psi(lu[l],kl[k],j)=get_psi(lui-l,kli+k,-j);
				set_psi(lu[l],ku[k],j)=get_psi(lui-l,kui-k,j);
			}
		}
	}

	// vertex
	for(int j=1;j<=tab;j++)
	{
		for(int k=1;k<=tab;k++)
		{
			for(int l=1;l<=tab;l++)
			{
				set_psi(ll[l],kl[k],jl[j])=get_psi(lli+l,kli+k,-jli-j);
				set_psi(ll[l],kl[k],ju[j])=get_psi(lli+l,kli+k,-jui+j);
				set_psi(ll[l],ku[k],jl[j])=get_psi(lli+l,kui-k,jli+j);
				set_psi(ll[l],ku[k],ju[j])=get_psi(lli+l,kui-k,jui-j);
				set_psi(lu[l],kl[k],jl[j])=get_psi(lui-l,kli+k,-jli-j);
				set_psi(lu[l],kl[k],ju[j])=get_psi(lui-l,kli+k,-jui+j);
				set_psi(lu[l],ku[k],jl[j])=get_psi(lui-l,kui-k,jli+j);
				set_psi(lu[l],ku[k],ju[j])=get_psi(lui-l,kui-k,jui-j);
			}
		}
	}
	//cout << "boundary reflection" << endl;
	
	return;
}

//extrinsic curvature reflection boundary
void Fmv0::boundary_quarter_scalar()
{
    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(l,kl[k],j,i)=get_bv(l,kli+k,-j,i);
				}
			}
		}
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(ll[l],k,j,i)=get_bv(lli+l,k,j,i);
				}
			}
		}
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=nsc;i<=nscp;i++)
				{
					set_bv(ll[l],kl[k],j,i)=get_bv(lli+l,kli+k,-j,i);
				}
			}
		}
	}
	
	return;
}

//even variable reflection boundary
void Fmv0::boundary_quarter_even(int i)
{
    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				set_bv(l,kl[k],j,i)=get_bv(l,kli+k,-j,i);
			}
		}
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_bv(ll[l],k,j,i)=get_bv(lli+l,k,j,i);
			}
		}
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_bv(ll[l],kl[k],j,i)=get_bv(lli+l,kli+k,-j,i);
			}
		}
	}
	
	return;
}


//reflection boundary for fluid
void Fmv0::boundary_quarter_fluid()
{
	// y-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28
	//
	//		vector scalar:	1,  3, 4,  6,  21,   23, 25,  27
	//		vector vector:	  2,     5,       22,       26,
	//
	//		tensor scalar:	   11,           18,
	//		tensor vector:	10,   12,     17,   19
	
	for(int i=24;i<29;i++)
	 refsign[i]=1; 

	refsign[26]=-1;

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				for(int i=24;i<29;i++)
				set_bv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-j,i);
			}
		}
	}
	#pragma omp barrier
	
	// z-fixed
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28
	//
	//		vector scalar:	1,2,     4,5,      21,22,       25,26,  
	//		vector vector:	    3,       6,          23,          27
	//
	//		tensor scalar:	10,           17,
	//		tensor vector:	   11,12,        18,19

	for(int i=24;i<29;i++)
	 refsign[i]=1; 

	refsign[27]=-1;

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=24;i<29;i++)
				{
					set_bv(ll[l],k,j,i)=refsign[i]*get_bv(lli+l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
    
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28 
	//
	//		vector positive:	1,      4,      21,       25,
	//		vector negative:	  2,3,    5,6,     22,23,    26,27
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=24;i<29;i++)
	 refsign[i]=1; 

	refsign[26]=-1;
	refsign[27]=-1;

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=24;i<29;i++)
				{
					set_bv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-j,i);
				}
			}
		}
	}
	#pragma omp barrier

	return;
}

//reflection boundary for geometry + fluid (+ scalar)
void Fmv0::boundary_quarter()
{
	
	// y-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	1,  3, 4,  6,  21,   23, 25,  27
	//		vector vector:	  2,     5,       22,       26,
	//
	//		tensor scalar:	   11,           18,
	//		tensor vector:	10,   12,     17,   19
	
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[5]=-1;
	refsign[22]=-1;
	refsign[10]=-1;
	refsign[12]=-1;
	refsign[17]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[26]=-1;

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				for(int i=0;i<nn;i++)
				set_bv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(l,kli+k,-j,i);
			}
		}
	}
	#pragma omp barrier
	
	// z-fixed
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	1,2,     4,5,      21,22,       25,26,  
	//		vector vector:	    3,       6,          23,          27
	//
	//		tensor scalar:	10,           17,
	//		tensor vector:	   11,12,        18,19

	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[3]=-1;
	refsign[6]=-1;
	refsign[23]=-1;
	refsign[11]=-1;
	refsign[12]=-1;
	refsign[18]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[27]=-1;

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<nn;i++)
				{
					set_bv(ll[l],k,j,i)=refsign[i]*get_bv(lli+l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	1,      4,      21,       25,
	//		vector negative:	  2,3,    5,6,     22,23,    26,27
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[3]=-1;
	refsign[5]=-1;
	refsign[6]=-1;
	refsign[22]=-1;
	refsign[23]=-1;
	refsign[10]=-1;
	refsign[11]=-1;
	refsign[17]=-1;
	refsign[18]=-1;

	if(fluidevo)
	{
		refsign[26]=-1;
		refsign[27]=-1;
	}

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<nn;i++)
				{
					set_bv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_bv(lli+l,kli+k,-j,i);
				}
			}
		}
	}
	
	return;
}

//boundary for fluid primitive variable
void Fmv0::boundary_prim_quarter()
{
	// y-fixed 
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		
	
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[2]=-1;

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				for(int i=0;i<npr;i++)
				{
					set_primv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_primv(l,kli+k,-j,i);
				}
			}
		}
	}
	#pragma omp barrier
	
	// z-fixed
	/////////////////////////////////////
	//     0 1 2 3 4 5 6 7 8 9 
	//     B B B I I I I B B B --> inner and buffer
	//           |     |       --> boundary position
	//     2 <= 3, 1 <= 4, 0 <= 5  || jl[j] <= jli+j
	//     7 <= 5, 8 <= 4, 9 <= 3  || ju[j] <= jui-j  for scalar
	//		

	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[3]=-1;

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<npr;i++)
				{
					set_primv(ll[l],k,j,i)=refsign[i]*get_primv(lli+l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier

	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	1,      4,      21,       27,
	//		vector negative:	  2,3,    5,6,     22,23,    28,29
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<npr;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[3]=-1;
	
	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<npr;i++)
				{
					set_primv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_primv(lli+l,kli+k,-j,i);
				}
			}
		}
	}

	return;
}

//reflection boundary for dbv of geometry + fluid (+ scalar)
void Fmv0::boundary_d_quarter()
{
	
	// y-fixed 
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector scalar:	1,  3,   4,  6,    21,   23, 25,   27
	//		vector vector:	  2,       5,         22,       26,
	//
	//		tensor scalar:	   11,         18,  
	//		tensor vector:	10,   12,   17,   19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[5]=-1;
	refsign[22]=-1;
	refsign[10]=-1;
	refsign[12]=-1;
	refsign[17]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[26]=-1;

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				for(int i=0;i<nn;i++)
				{
					set_dbv(l,kl[k],j,i)=rotsign[i]*refsign[i]*get_dbv(l,kli+k,-j,i);
				}
			}
		}
	}
	#pragma omp barrier
	// z-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30
	//
	//		vector scalar:	1,2,    4,5,    21,22,    25,26,
	//		vector vector:	    3,      6,        23,       27
	//
	//		tensor scalar:	10,         17,  
	//		tensor vector:	   11,12,      18,19
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[3]=-1;
	refsign[6]=-1;
	refsign[23]=-1;
	refsign[11]=-1;
	refsign[12]=-1;
	refsign[18]=-1;
	refsign[19]=-1;

	if(fluidevo)
	refsign[27]=-1;

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<nn;i++)
				{
					set_dbv(ll[l],k,j,i)=refsign[i]*get_dbv(lli+l,k,j,i);
				}
			}
		}
	}
	#pragma omp barrier
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20,24,28,29,30 
	//
	//		vector positive:	1,      4,      21,       25,
	//		vector negative:	  2,3,    5,6,     22,23,    26,27
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nn;i++)
	 refsign[i]=1; 

	refsign[2]=-1;
	refsign[3]=-1;
	refsign[5]=-1;
	refsign[6]=-1;
	refsign[22]=-1;
	refsign[23]=-1;
	refsign[10]=-1;
	refsign[11]=-1;
	refsign[17]=-1;
	refsign[18]=-1;

	if(fluidevo)
	{
		refsign[26]=-1;
		refsign[27]=-1;
	}

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<nn;i++)
				{
					set_dbv(ll[l],kl[k],j,i)=rotsign[i]*refsign[i]*get_dbv(lli+l,kli+k,-j,i);
				}
			}
		}
	}

	return;
}


//reflection boundary cond. for excision flags
void Fmv0::boundary_quarter_excflags()
{

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				set_bflag(l,kl[k],j)=get_bflag(l,kli+k,-j);
			}
		}
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_bflag(ll[l],k,j)=get_bflag(lli+l,k,j);
			}
		}
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_bflag(ll[l],kl[k],j)=get_bflag(lli+l,kli+k,-j);
			}
		}
	}


	return;
}

//reflection boundary condition for horizon flags
void Fmv0::boundary_quarter_hflags()
{

    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				set_hflag(l,kl[k],j)=get_hflag(l,kli+k,-j);
			}
		}
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_hflag(ll[l],k,j)=get_hflag(lli+l,k,j);
			}
		}
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_hflag(ll[l],kl[k],j)=get_hflag(lli+l,kli+k,-j);
			}
		}
	}

	return;
}

//reflection boundary condition for shift
void Fmv0::boundary_beta_quarter()
{
	// @@@@@@@@@@@@   variables (BSSN +GAUGE)   @@@@@@@@@@@
	//  0:wx  ,  1:wy  ,  2:wz 
	
	int nw=3;
	// y-fixed 
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	1,  3,   4,  6,    21,   23
	//		vector vector:	  2,       5,         22,
	//
	//		tensor scalar:	   11,         18,  
	//		tensor vector:	10,   12,   17,   19
	for(int i=0;i<nw;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
    #pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				for(int i=0;i<nw;i++)
				{
					set_bv(l,kl[k],j,i+1)=refsign[i]*get_bv(l,kli+k,j,i+1);
				}
			}
		}
	}
	#pragma omp barrier
	// z-fixed
	//
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector scalar:	1,2,    4,5,    21,22,
	//		vector vector:	    3,      6,        23
	//
	//		tensor scalar:	10,         17,  
	//		tensor vector:	   11,12,      18,19
	for(int i=0;i<nw;i++)
	 refsign[i]=1; 

	refsign[2]=-1;

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<nw;i++)
				{
					set_bv(ll[l],k,j,i+1)=refsign[i]*get_bv(lli+l,k,j,i+1);
				}
			}
		}
	}
	#pragma omp barrier
	// y,z-fixed
	//		pure scalar:	0,7,8,9,13,14,15,16,20 
	//
	//		vector positive:	1,      4,      21,
	//		vector negative:	  2,3,    5,6,     22,23
	//
	//		tensor positive:	      12,        19
	//		tensor negative:	10,11,     17,18,
	for(int i=0;i<nw;i++)
	 refsign[i]=1; 

	refsign[1]=-1;
	refsign[2]=-1;

	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				for(int i=0;i<nw;i++)
				{
					set_bv(ll[l],kl[k],j,i+1)=refsign[i]*get_bv(lli+l,kli+k,j,i+1);
				}
			}
		}
	}
	return;
}

void Fmv0::boundary_psi_initial_quarter()
{
	#pragma omp parallel for
    for(int j=jmin;j<=jmax;j++)
    {
        for(int l=lli;l<=lu[tab];l++)
	    {
            for(int k=1;k<=tab;k++)
            {
				set_psi(l,kl[k],j)=get_psi(l,kli+k,-j);
			}
		}
        for(int k=kli;k<=ku[tab];k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_psi(ll[l],k,j)=get_psi(lli+l,k,j);
			}
		}
        for(int k=1;k<=tab;k++)
	    {
            for(int l=1;l<=tab;l++)
            {
				set_psi(ll[l],kl[k],j)=get_psi(lli+l,kli+k,-j);
			}
		}
	}

	return;
}

//asymptotically flat boundary condition 
//Assumptions:
//geometry + scalar system with Cartesian coord. 
//scalar field is initially localized -> 0 in the background

void Fmv::asymcond(int l,int k,int j,int i,double bgv)
{
	//Coordinate should be Cartesian 
	double R=sqrt(pow(get_x(j),2)+pow(get_y(k),2)+pow(get_z(l),2));
	
	double delR=sqrt(dx*dx+dy*dy+dz*dz);
	double contrfac=(R-delR)/R;

	double zc=contrfac*get_z(l);
	double yc=contrfac*get_y(k);
	double xc=contrfac*get_x(j);
	
	int lc=int(zc*dzi);
	int kc=int(yc*dyi);
	int jc=int(xc*dxi);

	if(lc+2>lu[0])
	lc=lu[0]-2;
	if(kc+2>ku[0])
	kc=ku[0]-2;
	if(jc+2>ju[0])
	jc=ju[0]-2;
	if(jc-2<jl[0])
	jc=jl[0]+2;
	if(kc-2<kl[0])
	kc=kl[0]+2;
	if(lc-2<ll[0])
	lc=ll[0]+2;

	double Qv=bv_ipol(jc,kc,lc,xc,yc,zc,4,i)-bgv;

	set_bv(l,k,j,i)=bgv+(1.-delR/R)*Qv;
	
	return;
	
}

void Fmv::asymcond(int l,int k,int j,int i,double bgv0,double bgv,double dtime,int itype)
{
	//Coordinate should be Cartesian 
	double R=sqrt(pow(get_x(j),2)+pow(get_y(k),2)+pow(get_z(l),2));
	
	// double alpha=get_bv0(l,k,j,0);
	// double wa=get_bv0(l,k,j,13);
	// double delR=dtime*alpha*exp(-2.*wa);
	double delR=dtime;
	double contrfac=(R-delR)/R;

	double zc=contrfac*get_z(l);
	double yc=contrfac*get_y(k);
	double xc=contrfac*get_x(j);
	
	int lc=int(zc*dzi);
	int kc=int(yc*dyi);
	int jc=int(xc*dxi);

	if(lc+2>lu[tab])
	lc=lu[0];
	if(kc+2>ku[tab])
	kc=ku[0];
	if(jc+2>ju[tab])
	jc=ju[0];
	if(jc-2<jl[tab])
	jc=jl[0];
	if(kc-2<kl[tab])
	kc=kl[0];
	if(lc-2<ll[tab])
	lc=ll[0];

	double Qv;
	Qv=bv0_ipol(jc,kc,lc,xc,yc,zc,4,i)-bgv0;
	
	set_bv(l,k,j,i)=bgv+(1.-delR/R)*Qv;
	
	return;
	
}

void Fmv::boundary_asym(int itype)
{
	double dtime;
	double bgv0[nn],bgv[nn];

	if(itype==1 || itype==2)
	dtime=0.5*dt0;
	else
	dtime=dt0;

	bgv0[0]=1.;
	bgv[0]=1.;

	for(int i=1;i<nn;i++)
	{
		bgv0[i]=0.;
		bgv[i]=0.;
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

	// if(fluidevo)
	// {
	// 	double psibg1,trkbg1;
	// 	double psibg,trkbg;
	// 	double tt1=t-dtp,tt;

	// 	psibg1=1./(3.*(1.+fluidw))*log(tt1/tini);
	// 	trkbg1=-2./(1.+fluidw)/tt1;

	// 	if(itype==1 || itype==2)
	// 	tt=t+0.5*dt0;
	// 	else
	// 	tt=t+dt0;

	// 	psibg=1./(3.*(1.+fluidw))*log(tt/tini);
	// 	trkbg=-2./(1.+fluidw)/tt;

	// 	bgv1[13]=psibg1;
	// 	bgv[13]=psibg;
	// 	bgv1[20]=trkbg1;
	// 	bgv[20]=trkbg;
	// }
		
	//run in x direction
	#pragma omp parallel for
	for(int j=jmin;j<=jmax;j++)
	{
		for(int k=kui+1;k<=kmax;k++)
		{
			for(int l=lli;l<=lmax;l++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv0[i],bgv[i],dtime,itype);
			}
		}

		for(int k=kli;k<=kui;k++)
		{
			for(int l=lui+1;l<=lmax;l++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv0[i],bgv[i],dtime,itype);
			}
		}
	}

	//run in z direction
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{	
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jmin;j<jli;j++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv0[i],bgv[i],dtime,itype);
			}
			for(int j=jui+1;j<=jmax;j++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv0[i],bgv[i],dtime,itype);
			}
		}
	}

	#pragma omp barrier
	boundary_quarter();

	return;
}

void Fmv::boundary_asym0()
{
	double bgv[nn];

	bgv[0]=1.;

	for(int i=1;i<nn;i++)
	{
		bgv[i]=0.;
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

	// if(fluidevo)
	// {
	// 	double psibg,trkbg;

	// 	psibg=1./(3.*(1.+fluidw))*log(t/tini);
	// 	trkbg=-2./(1.+fluidw)/t;

	// 	bgv[13]=psibg;
	// 	bgv[20]=trkbg;
	// }
	
	//run in x direction
	#pragma omp parallel for
	for(int j=jmin;j<=jmax;j++)
	{
		for(int k=kui+1;k<=kmax;k++)
		{
			for(int l=lli;l<=lmax;l++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv[i]);
			}
		}

		for(int k=kli;k<=kui;k++)
		{
			for(int l=lui+1;l<=lmax;l++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv[i]);
			}
		}
	}

	//run in z direction
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{	
		for(int k=kli;k<=kui;k++)
		{
			for(int j=jmin;j<jli;j++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv[i]);
			}
			for(int j=jui+1;j<=jmax;j++)
			{
				for(int i=0;i<nn;i++)
				asymcond(l,k,j,i,bgv[i]);
			}
		}
	}

	#pragma omp barrier
	boundary_quarter();

	return;
}

