##################################################################################
##   parameters for COSMOS code            #######################################
##   ver1.00 by Chulmoon Yoo               #######################################
##################################################################################
9999999         # maximum step of the main loop 
400.	        # maximum time to evolve
3	            # tab number of the bufer grids
10.	            # amp
-40	            # minimum grid number of x =-nmax-1
40	            # maximum grid number of x =imax/2-1
0	            # minimum grid number of y
40	            # maximum grid number of y
0	            # minimum grid number of z
40	            # maximum grid number of z
-1.	            # minimum coordinate of x
1.	            # maximum coordinate of x
0.	            # minimum coordinate of y
1.	            # maximum coordinate of y
0.	            # minimum coordinate of z
1.	            # maximum coordinate of z

##################################################################################
###  parameters for Evolution
##################################################################################
0.05	        # CFL number
0.1	            # cdt number
2.0	            # etaa
5.0	            # etab(eta)
0.75	        # etabb(k)
0.	            # KO dissipation
5	            # excision grid

##################################################################################
###  initial data parameter
##################################################################################
0	            # 0:no continue 1:continue
ini_all.dat	    # continue file
0.83	        # amplitude 
10.	            # wave number 
0.	            # xi2 nonsphericity parameter 1
0.	            # xi3 nonsphericity parameter 2
0.	            # w3  alignment angle 
0.	            # amplitude for the scalar field
10.	            # wave number for the scalar field
15.	            # xi2s
0.	            # xi3s
50.0	        # Hubble

##################################################################################
###  fluid parameters
##################################################################################
0.3333333333333	# fluidw
-1.0 	        # Mkap
2.0	            # bminmod

##################################################################################
###  parameters for output
##################################################################################
10.0	            #1st part print interval boundary time
10.0	            #2nd part
100.	        #changing time for print interval
