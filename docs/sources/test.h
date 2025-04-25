/**
 * @file cosmos.h
 * @brief Header file for the BSSN evolution classes of the COSMOS code.
 *
 * Defines the base class Fmv0 for grid management, variable storage,
 * finite differencing, and basic evolution utilities. Also defines derived
 * classes Fmv (for specific initial data and boundary conditions) and
 * Fmv1 (for adaptive mesh refinement).
 *
 * @version 1.00
 * @author Chulmoon Yoo
 */

 #ifndef _COSMOS_H_
 #define _COSMOS_H_
 // ... (includes and definitions) ...
 
 using namespace std;
 
 /**
  * @class Fmv0
  * @brief Base class for Finite Mesh Variables version 0.
  *
  * This class manages the computational grid, stores dynamical variables
  * (BSSN, gauge, fluid, scalar fields), provides methods for finite
  * differencing, basic time stepping support, constraint checking,
  * and memory management. It forms the foundation for specific simulation
  * setups and mesh refinement.
  */
 class Fmv0{
 
 protected:
     // ... (protected member variables - not typically part of public API docs) ...
 
 public:
     /**
      * @brief Constructor for the Fmv0 class.
      * @param tabs Number of buffer zones (ghost zones).
      * @param jupper Upper index boundary in x for the physical domain.
      * @param jlower Lower index boundary in x for the physical domain.
      * @param kupper Upper index boundary in y for the physical domain.
      * @param klower Lower index boundary in y for the physical domain.
      * @param lupper Upper index boundary in z for the physical domain.
      * @param llower Lower index boundary in z for the physical domain.
      * @param xupper Upper coordinate boundary in x.
      * @param xlower Lower coordinate boundary in x.
      * @param yupper Upper coordinate boundary in y.
      * @param ylower Lower coordinate boundary in y.
      * @param zupper Upper coordinate boundary in z.
      * @param zlower Lower coordinate boundary in z.
      * @param am Amplitude for inhomogeneous grid mapping (if used).
      * @param fld Boolean flag to enable fluid evolution.
      * @param scl Boolean flag to enable scalar field evolution.
      * @param cuev Boolean flag to enable curvature evaluation.
      *
      * Initializes grid parameters, allocates memory for all variables,
      * coordinate arrays, flags, and temporary storage. Sets up basic
      * parameters and constants.
      */
     Fmv0(int tabs,int jupper,int jlower,int kupper,int klower,int lupper,int llower,
     double xupper,double xlower,double yupper,double ylower,double zupper,double zlower,double am,bool fld, bool scl, bool cuev);
 
     /**
      * @brief Virtual destructor for the Fmv0 class.
      *
      * Deallocates all dynamically allocated memory associated with the grid,
      * variables, flags, and coordinate arrays.
      */
     virtual ~Fmv0();
 
     //--------------------------------------
     // GET functions (Accessors)
     //--------------------------------------
 
     /** @brief Check if fluid evolution is enabled. */
     bool get_fluidevo() const;
     /** @brief Check if horizon formation has been detected (flag). */
     bool get_hform() const;
     /** @brief Check if excision is enabled. */
     bool get_exc() const;
     /** @brief Check if mesh refinement is active (flag, relevant for derived classes). */
     bool get_mrf() const;
 
     /** @brief Get the maximum grid index in x (including buffer). */
     int get_jmax() const;
     /** @brief Get the minimum grid index in x (including buffer). */
     int get_jmin() const;
     /** @brief Get the maximum grid index in y (including buffer). */
     int get_kmax() const;
     /** @brief Get the minimum grid index in y (including buffer). */
     int get_kmin() const;
     /** @brief Get the maximum grid index in z (including buffer). */
     int get_lmax() const;
     /** @brief Get the minimum grid index in z (including buffer). */
     int get_lmin() const;
 
     /** @brief Get the upper index boundary in x for the evolved physical domain. */
     int get_jui() const;
     /** @brief Get the lower index boundary in x for the evolved physical domain. */
     int get_jli() const;
     /** @brief Get the upper index boundary in y for the evolved physical domain. */
     int get_kui() const;
     /** @brief Get the lower index boundary in y for the evolved physical domain. */
     int get_kli() const;
     /** @brief Get the upper index boundary in z for the evolved physical domain. */
     int get_lui() const;
     /** @brief Get the lower index boundary in z for the evolved physical domain. */
     int get_lli() const;
 
     /** @brief Get the total number of grid points in x (including buffer). */
     int get_nx() const;
     /** @brief Get the total number of grid points in y (including buffer). */
     int get_ny() const;
     /** @brief Get the total number of grid points in z (including buffer). */
     int get_nz() const;
 
     /** @brief Get the x-index of the maximum Hamiltonian constraint violation. */
     int get_jhm() const;
     /** @brief Get the y-index of the maximum Hamiltonian constraint violation. */
     int get_khm() const;
     /** @brief Get the z-index of the maximum Hamiltonian constraint violation. */
     int get_lhm() const;
     /** @brief Get the x-index of the maximum Momentum constraint violation component. */
     int get_jmm() const;
     /** @brief Get the y-index of the maximum Momentum constraint violation component. */
     int get_kmm() const;
     /** @brief Get the z-index of the maximum Momentum constraint violation component. */
     int get_lmm() const;
     /** @brief Get the x-index of the maximum Kretschmann invariant. */
     int get_jkm() const;
     /** @brief Get the y-index of the maximum Kretschmann invariant. */
     int get_kkm() const;
     /** @brief Get the z-index of the maximum Kretschmann invariant. */
     int get_lkm() const;
     /** @brief Get the x-index of the maximum Weyl invariant. */
     int get_jwm() const;
     /** @brief Get the y-index of the maximum Weyl invariant. */
     int get_kwm() const;
     /** @brief Get the z-index of the maximum Weyl invariant. */
     int get_lwm() const;
 
     /** @brief Get the size (radius in grid points) of the excision region. */
     int get_exg() const;
     /** @brief Get the layer number (for AMR, 0 for base layer). */
     int get_layn() const;
     /** @brief Get the minimum x-index for constraint checking. */
     int get_jjmin() const;
     /** @brief Get the minimum y-index for constraint checking. */
     int get_kkmin() const;
     /** @brief Get the minimum z-index for constraint checking. */
     int get_llmin() const;
 
     /** @brief Get the current simulation time. */
     double get_t() const;
     /** @brief Get the current time step size (dt). */
     double get_dt() const;
     /** @brief Get the initial time step size (dt0). */
     double get_dt0() const;
     /** @brief Get the previous time step size (dtp). */
     double get_dtp() const;
     /** @brief Get the time step size from two steps ago (dtpp). */
     double get_dtpp() const;
     /** @brief Get the maximum simulation time. */
     double get_tmax() const;
     /** @brief Get the CFL factor. */
     double get_cfl() const;
 
     /** @brief Get the grid spacing in x (dx). */
     double get_dx() const;
     /** @brief Get the grid spacing in y (dy). */
     double get_dy() const;
     /** @brief Get the grid spacing in z (dz). */
     double get_dz() const;
     /** @brief Get the grid cell volume (dx*dy*dz). */
     double get_dvol() const;
     /** @brief Get the inverse grid spacing in x (1/dx). */
     double get_dxi() const;
     /** @brief Get the inverse grid spacing in y (1/dy). */
     double get_dyi() const;
     /** @brief Get the inverse grid spacing in z (1/dz). */
     double get_dzi() const;
     /** @brief Get 1/(2*dx). */
     double get_dxi2() const;
     /** @brief Get 1/(2*dy). */
     double get_dyi2() const;
     /** @brief Get 1/(2*dz). */
     double get_dzi2() const;
     /** @brief Get 1/(4*dx). */
     double get_dxi4() const;
     /** @brief Get 1/(4*dy). */
     double get_dyi4() const;
     /** @brief Get 1/(4*dz). */
     double get_dzi4() const;
 
     /** @brief Get the upper coordinate boundary in x. */
     double get_xu() const;
     /** @brief Get the lower coordinate boundary in x. */
     double get_xl() const;
     /** @brief Get the upper coordinate boundary in y. */
     double get_yu() const;
     /** @brief Get the lower coordinate boundary in y. */
     double get_yl() const;
     /** @brief Get the upper coordinate boundary in z. */
     double get_zu() const;
     /** @brief Get the lower coordinate boundary in z. */
     double get_zl() const;
 
     /** @brief Get the average Hamiltonian constraint violation. */
     double get_ham() const;
     /** @brief Get the maximum Hamiltonian constraint violation. */
     double get_hammax() const;
     /** @brief Get the maximum Kretschmann invariant value. */
     double get_Kremax() const;
     /** @brief Get the maximum Weyl invariant value. */
     double get_Weylmax() const;
     /** @brief Get the average Momentum constraint violation. */
     double get_mom() const;
     /** @brief Get the maximum Momentum constraint violation component. */
     double get_mommax() const;
 
     /** @brief Get the gauge parameter eta for lapse evolution (1+log slicing). */
     double get_etaa() const;
     /** @brief Get the gauge parameter eta for shift evolution (Gamma driver). */
     double get_etab() const;
     /** @brief Get the gauge parameter eta_b for shift evolution (Gamma driver damping term). */
     double get_etabb() const;
     /** @brief Get the cosmological constant lambda. */
     double get_lambda() const;
     /** @brief Get the initial time tini (often for analytic solutions). */
     double get_tini() const;
     /** @brief Get the Kreiss-Oliger dissipation coefficient epsilon. */
     double get_KOep() const;
     /** @brief Get the MUSCL interpolation parameter kappa. */
     double get_Mkap() const;
     /** @brief Get the minmod limiter parameter b. */
     double get_b() const;
     /** @brief Get the initial Hubble parameter Hb (for initial data). */
     double get_Hb() const;
     /** @brief Get the fluid equation of state parameter w (e.g., P = w * rho * epsilon for ideal fluid). */
     double get_fluidw() const;
     /** @brief Get the scalar field mass parameter. */
     double get_scalarm() const;
 
     /**
      * @brief Get the boundary flag at a specific grid point.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Boundary flag value (e.g., 0=physical, -1=excision, >0 buffer/boundary type).
      */
     int get_bflag(int l,int k,int j) const;
     /**
      * @brief Get the horizon flag at a specific grid point.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Horizon flag value (e.g., 0=outside, 1=inside apparent horizon).
      */
     int get_hflag(int l,int k,int j) const;
     /**
      * @brief Get the fine mesh refinement flag at a specific grid point.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Flag indicating if this point is covered by a finer grid (1 if yes, 0 if no).
      */
     int get_fmrflag(int l,int k,int j) const;
 
     /** @brief Get the x-coordinate at index j. */
     double get_x(int j)const;
     /** @brief Get the y-coordinate at index k. */
     double get_y(int k)const;
     /** @brief Get the z-coordinate at index l. */
     double get_z(int l)const;
     /** @brief Calculate the x-coordinate corresponding to index j (even outside stored range). */
     double get_ext_x(int j)const;
     /** @brief Calculate the y-coordinate corresponding to index k (even outside stored range). */
     double get_ext_y(int k)const;
     /** @brief Calculate the z-coordinate corresponding to index l (even outside stored range). */
     double get_ext_z(int l)const;
 
     /**
      * @brief Get the value of a dynamical variable.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Value of bv[i][l][k][j].
      */
     double get_bv(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a time derivative variable (RHS).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Value of dbv[i][l][k][j].
      */
     double get_dbv(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a dynamical variable at the previous time step.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Value of bv0[i][l][k][j].
      */
     double get_bv0(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a dynamical variable at two time steps ago.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Value of bv1[i][l][k][j].
      */
     double get_bv1(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of the Runge-Kutta sum variable.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Value of bvr[i][l][k][j].
      */
     double get_bvr(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a constraint variable.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Constraint index (0-nc).
      * @return Value of con[i][l][k][j].
      */
     double get_con(int l,int k,int j,int i) const;
 
     /** @brief Get the flat metric coordinate mapping function derivative d^2f/dx^2 at index j. */
     double get_flat_df2x(int j)const;
     /** @brief Get the flat metric coordinate mapping function derivative d^2f/dy^2 at index k. */
     double get_flat_df2y(int k)const;
     /** @brief Get the flat metric coordinate mapping function derivative d^2f/dz^2 at index l. */
     double get_flat_df2z(int l)const;
     /** @brief Get the flat metric Christoffel symbol Gamma^x_xx at index j. */
     double get_flat_Gamx(int j)const;
     /** @brief Get the flat metric Christoffel symbol Gamma^y_yy at index k. */
     double get_flat_Gamy(int k)const;
     /** @brief Get the flat metric Christoffel symbol Gamma^z_zz at index l. */
     double get_flat_Gamz(int l)const;
     /** @brief Get the derivative of the flat metric Christoffel symbol d(Gamma^x_xx)/dx at index j. */
     double get_flat_dGamx(int j)const;
     /** @brief Get the derivative of the flat metric Christoffel symbol d(Gamma^y_yy)/dy at index k. */
     double get_flat_dGamy(int k)const;
     /** @brief Get the derivative of the flat metric Christoffel symbol d(Gamma^z_zz)/dz at index l. */
     double get_flat_dGamz(int l)const;
 
     /**
      * @brief Get the value of an output variable (e.g., curvature invariants).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Output variable index (0-nv).
      * @return Value of outv[i][l][k][j].
      */
     double get_outv(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of the conformal factor psi (temporary storage, often from bv[13]).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Value of psi[l][k][j].
      */
     double get_psi(int l,int k,int j) const;
     /**
      * @brief Get the value of a primitive fluid variable.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Primitive variable index (0-npr).
      * @return Value of primv[i][l][k][j].
      */
     double get_primv(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a fluid flux component in the x-direction.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index (often at i-1/2 interface).
      * @param i Conserved variable index corresponding to the flux (0-npr).
      * @return Value of flux_x[i][l][k][j].
      */
     double get_flux_x(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a fluid flux component in the y-direction.
      * @param l z-index.
      * @param k y-index (often at k-1/2 interface).
      * @param j x-index.
      * @param i Conserved variable index corresponding to the flux (0-npr).
      * @return Value of flux_y[i][l][k][j].
      */
     double get_flux_y(int l,int k,int j,int i) const;
     /**
      * @brief Get the value of a fluid flux component in the z-direction.
      * @param l z-index (often at l-1/2 interface).
      * @param k y-index.
      * @param j x-index.
      * @param i Conserved variable index corresponding to the flux (0-npr).
      * @return Value of flux_z[i][l][k][j].
      */
     double get_flux_z(int l,int k,int j,int i) const;
 
     /** @brief Get a pointer to the 3D array for dynamical variable i. */
     double*** get_bv(int i) const;
     /** @brief Get a pointer to the 3D array for time derivative variable i. */
     double*** get_dbv(int i) const;
     /** @brief Get a pointer to the 3D array for previous time step variable i. */
     double*** get_bv0(int i) const;
     /** @brief Get a pointer to the 3D array for variable i from two steps ago. */
     double*** get_bv1(int i) const;
     /** @brief Get a pointer to the 3D array for Runge-Kutta sum variable i. */
     double*** get_bvr(int i) const;
 
     //--------------------------------------
     // GET derivative functions
     //--------------------------------------
 
     /** @brief Calculate 4th order finite difference derivative d(bv[i])/dx at (l,k,j). */
     double get_f_x(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d(bv[i])/dy at (l,k,j). */
     double get_f_y(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d(bv[i])/dz at (l,k,j). */
     double get_f_z(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d^2(bv[i])/dx^2 at (l,k,j). */
     double get_f_xx(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d^2(bv[i])/dy^2 at (l,k,j). */
     double get_f_yy(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d^2(bv[i])/dz^2 at (l,k,j). */
     double get_f_zz(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d^2(bv[i])/dxdy at (l,k,j). */
     double get_f_xy(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d^2(bv[i])/dxdz at (l,k,j). */
     double get_f_xz(int l,int k,int j,int i) const;
     /** @brief Calculate 4th order finite difference derivative d^2(bv[i])/dydz at (l,k,j). */
     double get_f_yz(int l,int k,int j,int i) const;
 
     /** @brief Calculate 3rd order interpolation of bv[i] at the midpoint (j-1/2, k, l). */
     double get_ipol_x_lower_mid(int l,int k,int j,int i) const;
     /** @brief Calculate 3rd order interpolation of bv[i] at the midpoint (j, k-1/2, l). */
     double get_ipol_y_lower_mid(int l,int k,int j,int i) const;
     /** @brief Calculate 3rd order interpolation of bv[i] at the midpoint (j, k, l-1/2). */
     double get_ipol_z_lower_mid(int l,int k,int j,int i) const;
 
     /** @brief Calculate 2nd order finite difference derivative d(bv[i])/dx at (l,k,j). */
     double get_2ndf_x(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d(bv[i])/dy at (l,k,j). */
     double get_2ndf_y(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d(bv[i])/dz at (l,k,j). */
     double get_2ndf_z(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d^2(bv[i])/dx^2 at (l,k,j). */
     double get_2ndf_xx(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d^2(bv[i])/dy^2 at (l,k,j). */
     double get_2ndf_yy(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d^2(bv[i])/dz^2 at (l,k,j). */
     double get_2ndf_zz(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d^2(bv[i])/dxdy at (l,k,j). */
     double get_2ndf_xy(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d^2(bv[i])/dxdz at (l,k,j). */
     double get_2ndf_xz(int l,int k,int j,int i) const;
     /** @brief Calculate 2nd order finite difference derivative d^2(bv[i])/dydz at (l,k,j). */
     double get_2ndf_yz(int l,int k,int j,int i) const;
 
     //--------------------------------------
     // SET functions (Mutators)
     //--------------------------------------
 
     /** @brief Enable/disable fluid evolution. */
     void set_fluidevo(bool f);
     /** @brief Enable/disable scalar field evolution. */
     void set_scalarevo(bool s);
     /** @brief Set the current simulation time. */
     void set_t(double time);
     /** @brief Set the current time step size (dt). */
     void set_dt(double time);
     /** @brief Set the initial time step size (dt0). */
     void set_dt0(double time);
     /** @brief Set the previous time step size (dtp). */
     void set_dtp(double time);
     /** @brief Set the time step size from two steps ago (dtpp). */
     void set_dtpp(double time);
     /** @brief Set the maximum simulation time. */
     void set_tmax(double time);
     /** @brief Set the CFL factor. */
     void set_cfl(double c);
     /** @brief Set the gauge parameter etaa. */
     void set_etaa(double e);
     /** @brief Set the gauge parameter etab. */
     void set_etab(double e);
     /** @brief Set the gauge parameter etabb. */
     void set_etabb(double e);
     /** @brief Set the cosmological constant lambda. */
     void set_lambda(double l);
     /** @brief Set the initial Hubble parameter Hb. */
     void set_Hb(double hb);
     /** @brief Set the initial time tini. */
     void set_tini(double t);
     /** @brief Set the Kreiss-Oliger dissipation coefficient epsilon. */
     void set_KOep(double l);
     /** @brief Set the size (radius in grid points) of the excision region. */
     void set_exg(int eg);
     /** @brief Set the amplitude for inhomogeneous grid mapping. */
     void set_amp(double a);
     /** @brief Set the fluid equation of state parameter w. */
     void set_fluidw(double fw);
     /** @brief Set the scalar field mass parameter. */
     void set_scalarm(double sm);
     /** @brief Set the MUSCL interpolation parameter kappa. */
     void set_Mkap(double k);
     /** @brief Set the minmod limiter parameter b. */
     void set_b(double b);
     /** @brief Enable/disable excision. */
     void set_exc(bool e);
     /** @brief Enable/disable mesh refinement flag (used by derived classes). */
     void set_mrf(bool s);
     /** @brief Set the minimum x-index for constraint checking. */
     void set_jjmin(int jjj);
     /** @brief Set the minimum y-index for constraint checking. */
     void set_kkmin(int kkk);
     /** @brief Set the minimum z-index for constraint checking. */
     void set_llmin(int lll);
 
     /**
      * @brief Set the boundary flag at a specific grid point. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to bflag[l-lmin][k-kmin][j-jmin].
      */
     int& set_bflag(int l,int k,int j);
     /**
      * @brief Set the horizon flag at a specific grid point. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to hflag[l-lmin][k-kmin][j-jmin].
      */
     int& set_hflag(int l,int k,int j);
     /**
      * @brief Set the fine mesh refinement flag at a specific grid point. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to fmrflag[l-lmin][k-kmin][j-jmin].
      */
     int& set_fmrflag(int l,int k,int j);
     /**
      * @brief Set the fine mesh refinement flag to 1 within a specified region.
      * @param lmi Minimum z-index of the region.
      * @param lma Maximum z-index of the region.
      * @param kmi Minimum y-index of the region.
      * @param kma Maximum y-index of the region.
      * @param jmi Minimum x-index of the region.
      * @param jma Maximum x-index of the region.
      */
     void set_fmrregion(int lmi,int lma,int kmi,int kma,int jmi,int jma);
 
     /**
      * @brief Set the value of a dynamical variable. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Reference to bv[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_bv(int l,int k,int j,int i);
     /**
      * @brief Set the value of a time derivative variable (RHS). Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Reference to dbv[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_dbv(int l,int k,int j,int i);
     /**
      * @brief Set the value of a dynamical variable at the previous time step. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Reference to bv0[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_bv0(int l,int k,int j,int i);
     /**
      * @brief Set the value of a dynamical variable at two time steps ago. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Reference to bv1[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_bv1(int l,int k,int j,int i);
     /**
      * @brief Set the value of the Runge-Kutta sum variable. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (0-nn).
      * @return Reference to bvr[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_bvr(int l,int k,int j,int i);
     /**
      * @brief Set the value of a constraint variable. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Constraint index (0-nc).
      * @return Reference to con[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_con(int l,int k,int j,int i);
 
     /** @brief Set the flat metric coordinate mapping function derivative d^2f/dx^2 at index j. Returns a reference. */
     double& set_flat_df2x(int j);
     /** @brief Set the flat metric coordinate mapping function derivative d^2f/dy^2 at index k. Returns a reference. */
     double& set_flat_df2y(int k);
     /** @brief Set the flat metric coordinate mapping function derivative d^2f/dz^2 at index l. Returns a reference. */
     double& set_flat_df2z(int l);
     /** @brief Set the flat metric Christoffel symbol Gamma^x_xx at index j. Returns a reference. */
     double& set_flat_Gamx(int j);
     /** @brief Set the flat metric Christoffel symbol Gamma^y_yy at index k. Returns a reference. */
     double& set_flat_Gamy(int k);
     /** @brief Set the flat metric Christoffel symbol Gamma^z_zz at index l. Returns a reference. */
     double& set_flat_Gamz(int l);
     /** @brief Set the derivative of the flat metric Christoffel symbol d(Gamma^x_xx)/dx at index j. Returns a reference. */
     double& set_flat_dGamx(int j);
     /** @brief Set the derivative of the flat metric Christoffel symbol d(Gamma^y_yy)/dy at index k. Returns a reference. */
     double& set_flat_dGamy(int k);
     /** @brief Set the derivative of the flat metric Christoffel symbol d(Gamma^z_zz)/dz at index l. Returns a reference. */
     double& set_flat_dGamz(int l);
 
     /**
      * @brief Set the value of an output variable. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Output variable index (0-nv).
      * @return Reference to outv[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_outv(int l,int k,int j,int i);
     /**
      * @brief Set the value of a primitive fluid variable. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Primitive variable index (0-npr).
      * @return Reference to primv[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_primv(int l,int k,int j,int i);
     /**
      * @brief Set the value of a fluid flux component in the x-direction. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Conserved variable index corresponding to the flux (0-npr).
      * @return Reference to flux_x[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_flux_x(int l,int k,int j,int i);
     /**
      * @brief Set the value of a fluid flux component in the y-direction. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Conserved variable index corresponding to the flux (0-npr).
      * @return Reference to flux_y[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_flux_y(int l,int k,int j,int i);
     /**
      * @brief Set the value of a fluid flux component in the z-direction. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Conserved variable index corresponding to the flux (0-npr).
      * @return Reference to flux_z[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_flux_z(int l,int k,int j,int i);
 
     /** @brief Set pointer to the 3D array for dynamical variable i (use with caution). */
     double*** set_bv(int i) const;
     /** @brief Set pointer to the 3D array for time derivative variable i (use with caution). */
     double*** set_dbv(int i) const;
     /** @brief Set pointer to the 3D array for previous time step variable i (use with caution). */
     double*** set_bv0(int i) const;
     /** @brief Set pointer to the 3D array for variable i from two steps ago (use with caution). */
     double*** set_bv1(int i) const;
     /** @brief Set pointer to the 3D array for Runge-Kutta sum variable i (use with caution). */
     double*** set_bvr(int i) const;
 
     /**
      * @brief Set the value of the conformal factor psi (temporary storage). Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to psi[l-lmin][k-kmin][j-jmin].
      */
     double& set_psi(int l,int k,int j);
 
     //--------------------------------------
     // SET zero functions
     //--------------------------------------
 
     /** @brief Set all dynamical variables (bv, dbv, bv0, bv1), constraints (con), output (outv), flat metric helpers, psi, primitive variables (primv), and fluxes to zero. */
     void set_zero_all();
     /** @brief Set time derivatives (dbv), constraints (con), and output (outv) to zero only within the excised region (bflag == -1). */
     void set_zero_all_exc();
     /** @brief Set all primitive fluid variables (primv) to zero. */
     void set_zero_primv();
     /** @brief Set current dynamical variables (bv) to zero. */
     void set_zero();
     /** @brief Set previous time step dynamical variables (bv0) to zero. */
     void set_zero_0();
     /** @brief Set dynamical variables from two steps ago (bv1) to zero. */
     void set_zero_1();
     /** @brief Set time derivative variables (dbv) to zero. */
     void set_zero_d();
     /** @brief Set Runge-Kutta sum variables (bvr) to zero. */
     void set_zero_r();
     /** @brief Set all boundary flags (bflag) to zero. */
     void set_bflag_zero();
     /** @brief Set all horizon flags (hflag) to zero. */
     void set_hflag_zero();
     /** @brief Set all fine mesh refinement flags (fmrflag) to zero. */
     void set_fmrflag_zero();
 
     //--------------------------------------
     // UPDATE functions (Time stepping support)
     //--------------------------------------
 
     /** @brief Copy current dynamical variables (bv) to previous time step storage (bv0). */
     void setv0();
     /** @brief Copy previous time step variables (bv0) to two steps ago storage (bv1). */
     void set01();
     /** @brief Copy previous time step variables (bv0) back to current variables (bv). (Used for RK substeps). */
     void set0v();
     /**
      * @brief Perform one Runge-Kutta accumulation step: bvr = bvr + dbv * dt.
      * @param dt Time step increment for this RK stage.
      */
     void runge_kutta(double dt);
     /**
      * @brief Update variables using Euler step or final RK step: bv = bv0 + dbv * dt.
      * @param dt Full time step size.
      */
     void new_bv(double dt);
     /** @brief Update variables using accumulated RK sum: bv = bv0 + bvr. (Assumes bvr contains sum of k_i/b_i). */
     void new_bv4();
     /** @brief Copy the conformal factor (bv[13]) into the temporary psi array. */
     void set_psi();
 
     //--------------------------------------
     // ERROR evaluation functions
     //--------------------------------------
 
     /**
      * @brief Calculate Hamiltonian and Momentum constraint violations.
      * Computes average and maximum values over the physical domain
      * (excluding excised/refined regions) and stores them in member variables
      * (ham, hammax, mom, mommax, dGam, dGammax) along with the location of maxima.
      * Uses values stored in the `con` array.
      */
     void check_const();
     /**
      * @brief Calculate the maximum Kretschmann scalar.
      * Finds the maximum value over the physical domain (excluding excised/refined regions)
      * and stores it in Kremax and the location in (lkm, kkm, jkm).
      * Uses values stored in the `outv` array (index 0).
      */
     void check_Kremax();
     /**
      * @brief Calculate the maximum Weyl scalar invariant I.
      * Finds the maximum value over the physical domain (excluding excised/refined regions)
      * and stores it in Weylmax and the location in (lwm, kwm, jwm).
      * Uses values stored in the `outv` array (index 3).
      */
     void check_Weylmax();
 
     //--------------------------------------
     // Excision setup
     //--------------------------------------
 
     /**
      * @brief Sets up boundary flags (bflag) for a square excision region.
      * Marks points inside `exg` distance from origin as excised (-1) or
      * buffer layers (1, 2, 3). Calls boundary routines to fill flags in other quadrants.
      * Zeros out variables in the newly excised region.
      */
     void set_excflags_square();
 
     //--------------------------------------
     // Fluid helper functions (Declarations)
     //--------------------------------------
 
     /**
      * @brief Calculate fluid pressure from density (Equation of State).
      * @param rho Rest mass density.
      * @return Pressure P. (Implementation specific)
      */
     double pres(double rho);
     /**
      * @brief Calculate derivative of pressure with respect to density dP/d(rho*epsilon).
      * @param rho Rest mass density.
      * @return Derivative dP/d(rho*epsilon). (Implementation specific)
      */
     double dpres(double rho);
     /**
      * @brief Convert conserved fluid variables (E, S_i, D) to primitive (rho, Gamma*v^i).
      * @param Ene Conserved energy density.
      * @param S Conserved momentum density magnitude |S_i|.
      * @param[out] rho Resulting rest mass density.
      * @param[out] Gam Resulting Lorentz factor Gamma. (Implementation specific)
      */
     void get_rhoGam(double Ene, double S,double& rho,double& Gam);
     /**
      * @brief Minmod limiter function.
      * @param a First input value.
      * @param b Second input value.
      * @return minmod(a, b). (Implementation specific)
      */
     double minmod(double a,double b);
     /**
      * @brief Sign function.
      * @param A Input value.
      * @return Sign of A (+1, -1, or 0).
      */
     double sign(double A);
     /**
      * @brief Convert conserved dynamical fluid variables (bv[24]..bv[28]) to primitive variables (primv) across the grid.
      * (Implementation specific)
      */
     void dyntoprim();
 
     //--------------------------------------
     // Interpolation functions (Declarations)
     //--------------------------------------
 
     /**
      * @brief General 1D polynomial interpolation.
      * @param rr Point to interpolate at.
      * @param xx Array of x-coordinates of known points.
      * @param yy Array of y-coordinates (function values) of known points.
      * @param order Order of interpolation.
      * @return Interpolated value at rr. (Implementation specific)
      */
     double ipol( double rr,double *xx,double *yy,int order );
     /**
      * @brief 3D interpolation of previous time step variable bv0.
      * @param jc, kc, lc Indices of the 'center' grid point near the interpolation target.
      * @param xc, yc, zc Coordinates of the target point.
      * @param order Order of interpolation.
      * @param number Index of the variable to interpolate (0-nn).
      * @return Interpolated value of bv0[number] at (xc, yc, zc). (Implementation specific)
      */
     double bv0_ipol(int jc,int kc,int lc,double xc,double yc,double zc,int order,int number);
     /**
      * @brief 3D interpolation of variable bv1 (two steps ago).
      * @param jc, kc, lc Indices of the 'center' grid point near the interpolation target.
      * @param xc, yc, zc Coordinates of the target point.
      * @param order Order of interpolation.
      * @param number Index of the variable to interpolate (0-nn).
      * @return Interpolated value of bv1[number] at (xc, yc, zc). (Implementation specific)
      */
     double bv1_ipol(int jc,int kc,int lc,double xc,double yc,double zc,int order,int number);
     // ... (similar declarations for bv0_ipol_diff_x, _y, _z) ...
     /**
      * @brief 3D interpolation of current time step variable bv.
      * @param jc, kc, lc Indices of the 'center' grid point near the interpolation target.
      * @param xc, yc, zc Coordinates of the target point.
      * @param order Order of interpolation.
      * @param number Index of the variable to interpolate (0-nn).
      * @return Interpolated value of bv[number] at (xc, yc, zc). (Implementation specific)
      */
     double bv_ipol(int jc,int kc,int lc,double xc,double yc,double zc,int order,int number);
     // ... (similar declarations for bv_ipol_diff_x, _y, _z) ...
 
     //--------------------------------------
     // BSSN Evolution functions (Declarations)
     //--------------------------------------
 
     /** @brief Calculate advection terms for BSSN evolution. (Implementation specific) */
     void BSSN_adv();
     /**
      * @brief Calculate the right-hand sides (RHS) of the BSSN evolution equations.
      * @param itype Type parameter (e.g., for different RK stages or boundary handling). (Implementation specific)
      */
     void BSSN(int itype);
     /** @brief Enforce algebraic constraints (e.g., determinant of tilde_gamma = 1). (Implementation specific) */
     void enforce_const();
     /**
      * @brief Enforce algebraic constraints at a single grid point.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index. (Implementation specific)
      */
     void enforce_const_gp(int l,int k,int j);
     /** @brief Add Kreiss-Oliger dissipation to the RHS (dbv). (Implementation specific) */
     void KOdiss();
     /** @brief Apply excision boundary conditions (zeroing inside). (Implementation specific) */
     void excision();
     /** @brief Fill flux values, potentially at boundaries or refinement interfaces. (Implementation specific) */
     void flux_fill();
 
     //--------------------------------------
     // BOUNDARY condition functions (Declarations - assuming quadrant symmetry)
     //--------------------------------------
 
     /** @brief Apply boundary conditions to dynamical variables (bv) assuming quadrant symmetry. (Implementation specific) */
     void boundary_quarter();
     /** @brief Apply boundary conditions to RHS variables (dbv) assuming quadrant symmetry. (Implementation specific) */
     void boundary_d_quarter();
     /** @brief Apply even reflection boundary conditions to a specific variable (bv[i]) assuming quadrant symmetry. (Implementation specific) */
     void boundary_quarter_even(int i);
     /** @brief Apply boundary conditions to primitive fluid variables (primv) assuming quadrant symmetry. (Implementation specific) */
     void boundary_prim_quarter();
     /** @brief Apply boundary conditions to boundary flags (bflag) assuming quadrant symmetry. (Implementation specific) */
     void boundary_quarter_excflags();
     /** @brief Apply boundary conditions to horizon flags (hflag) assuming quadrant symmetry. (Implementation specific) */
     void boundary_quarter_hflags();
     /** @brief Apply boundary conditions to initial psi array assuming quadrant symmetry. (Implementation specific) */
     void boundary_psi_initial_quarter();
     /** @brief Apply boundary conditions to shift vector components (beta^i) assuming quadrant symmetry. (Implementation specific) */
     void boundary_beta_quarter();
     /** @brief Apply boundary conditions to fluid variables assuming quadrant symmetry. (Implementation specific) */
     void boundary_quarter_fluid();
     /** @brief Apply boundary conditions to scalar field variables assuming quadrant symmetry. (Implementation specific) */
     void boundary_quarter_scalar();
 
     //--------------------------------------
     // Initial Data functions (Declarations)
     //--------------------------------------
 
     /**
      * @brief Read simulation state from a continuation file.
      * @param fcontinue Input file stream opened to the continuation data file. (Implementation specific)
      */
     void initial_continue(ifstream& fcontinue);
     /** @brief Coordinate mapping function f(X) for inhomogeneous grids. (Implementation specific) */
     double funcf(double X);
     /** @brief Inverse coordinate mapping function f^{-1}(x). (Implementation specific) */
     double ifuncf(double X);
     /** @brief First derivative of coordinate mapping function df/dX. (Implementation specific) */
     double df(double X);
     /** @brief Second derivative of coordinate mapping function d^2f/dX^2. (Implementation specific) */
     double ddf(double X);
     /** @brief Third derivative of coordinate mapping function d^3f/dX^3. (Implementation specific) */
     double dddf(double X);
     /** @brief Set up flat metric components and related quantities for inhomogeneous coordinates. (Implementation specific) */
     void set_flat();
     /** @brief Set up Christoffel symbols for the flat metric in inhomogeneous coordinates. (Implementation specific) */
     void set_Gam();
     /** @brief Set initial energy-momentum tensor components based on initial data setup. (Implementation specific) */
     void set_enemomini();
     /**
      * @brief Set various simulation parameters read from input.
      * @param cfli CFL factor.
      * @param etaai Gauge parameter etaa.
      * @param etabi Gauge parameter etab.
      * @param etabbi Gauge parameter etabb.
      * @param lambdai Cosmological constant.
      * @param dt0i Initial time step.
      * @param dtpi Previous time step (for restart).
      * @param dtppi Time step two steps ago (for restart).
      * @param ti Current time (for restart).
      * @param tinii Initial time tini.
      * @param Hbi Initial Hubble parameter.
      * @param KOepi Kreiss-Oliger dissipation coefficient.
      * @param exgi Excision radius.
      * @param fluidwi Fluid EOS parameter w.
      * @param scalarmi Scalar field mass.
      * @param kap_MUSCLi MUSCL parameter kappa.
      * @param b_minmodi Minmod parameter b.
      */
     void initial_params(double cfli,double etaai,double etabi,double etabbi,double lambdai,double dt0i,double dtpi,double dtppi,double ti,double tinii,double Hbi,double KOepi,int exgi,double fluidwi,double scalarmi,double kap_MUSCLi,double b_minmodi);
 
     //--------------------------------------
     // OUTPUT functions (Declarations)
     //--------------------------------------
 
     /** @brief Print data along an x-line at fixed y=k, z=l to a file. (Implementation specific) */
     void print_x(ofstream& fn, int k,int l);
     /** @brief Print data along a y-line at fixed x=j, z=l to a file. (Implementation specific) */
     void print_y(ofstream& fn, int j,int l);
     /** @brief Print data along a z-line at fixed x=j, y=k to a file. (Implementation specific) */
     void print_z(ofstream& fn, int j,int k);
     /** @brief Print data on the xy-plane at fixed z=l to a file. (Implementation specific) */
     void print_xy(ofstream& fn, int l);
     /** @brief Print data on the xz-plane at fixed y=k to a file. (Implementation specific) */
     void print_xz(ofstream& fn, int k);
     /** @brief Print data on the yz-plane at fixed x=j to a file. (Implementation specific) */
     void print_yz(ofstream& fn, int j);
     /** @brief Print the maximum Kretschmann value to a file. (Implementation specific) */
     void print_Kremax(ofstream& fout);
     /** @brief Print constraint violation information to a file. (Implementation specific) */
     void print_const(ofstream& fout);
     /** @brief Print all relevant simulation data (e.g., for restart) to a file. (Implementation specific) */
     void print_all(ofstream& fout);
     /** @brief Print selected 3D data to a file. (Implementation specific) */
     void print_3d(ofstream& fout);
 
     //--------------------------------------
     // Potential functions (Scalar Field)
     //--------------------------------------
 
     /**
      * @brief Calculate the scalar field potential V(phi).
      * @param p Scalar field value phi.
      * @return Potential value V(p) = 0.5 * (scalarm * p)^2.
      */
     double funcV(double p);
     /**
      * @brief Calculate the derivative of the scalar field potential dV/dphi.
      * @param p Scalar field value phi.
      * @return Derivative value dV/dp = scalarm^2 * p.
      */
     double funcdV(double p);
 
 }; // End class Fmv0
 
 /**
  * @class Fmv
  * @brief Derived class from Fmv0, implementing specific boundary conditions
  *        and initial data setups for a single, non-refined grid.
  */
 class Fmv : public Fmv0{
 private:
     // No private members specific to Fmv shown
 public:
     /**
      * @brief Constructor for the Fmv class.
      * @param tabs Number of buffer zones (ghost zones).
      * @param jupper Upper index boundary in x for the physical domain.
      * @param jlower Lower index boundary in x for the physical domain.
      * @param kupper Upper index boundary in y for the physical domain.
      * @param klower Lower index boundary in y for the physical domain.
      * @param lupper Upper index boundary in z for the physical domain.
      * @param llower Lower index boundary in z for the physical domain.
      * @param xupper Upper coordinate boundary in x.
      * @param xlower Lower coordinate boundary in x.
      * @param yupper Upper coordinate boundary in y.
      * @param ylower Lower coordinate boundary in y.
      * @param zupper Upper coordinate boundary in z.
      * @param zlower Lower coordinate boundary in z.
      * @param am Amplitude for inhomogeneous grid mapping (if used).
      * @param fld Boolean flag to enable fluid evolution.
      * @param scl Boolean flag to enable scalar field evolution.
      * @param cuev Boolean flag to enable curvature evaluation.
      *
      * Calls the Fmv0 constructor and sets the layer number `layn` to 0.
      */
     Fmv(int tabs,int jupper,int jlower,int kupper,int klower,int lupper,int llower,
     double xupper,double xlower,double yupper,double ylower,double zupper,double zlower,double am,bool fld, bool scl, bool cuev);
 
     //--------------------------------------
     // BOUNDARY condition functions (Declarations - Specific Implementations)
     //--------------------------------------
 
     /** @brief Apply periodic boundary conditions. (Implementation specific) */
     void boundary_periodic();
     /** @brief Apply reflection boundary conditions (typically across coordinate planes). (Implementation specific) */
     void boundary_reflection();
     /** @brief Apply reflection boundary conditions to RHS variables (dbv). (Implementation specific) */
     void boundary_d_reflection();
     /** @brief Apply reflection boundary conditions to RHS variables (dbv) using ghost zone values. (Implementation specific) */
     void boundary_d_gs_reflection();
     /** @brief Apply even reflection boundary conditions to a specific variable (bv[i]). (Implementation specific) */
     void boundary_reflection_even(int i);
     /** @brief Apply reflection boundary conditions to primitive fluid variables (primv). (Implementation specific) */
     void boundary_prim_reflection();
     /** @brief Apply reflection boundary conditions to horizon flags (hflag). (Implementation specific) */
     void boundary_reflection_hflags();
     /** @brief Apply boundary conditions (type depends on implementation, potentially outflow/Sommerfeld). (Implementation specific) */
     void boundary_w();
     /** @brief Apply boundary conditions to initial psi array. (Implementation specific) */
     void boundary_psi_initial();
     /** @brief Apply boundary conditions to shift vector components (beta^i). (Implementation specific) */
     void boundary_beta();
     /** @brief Apply reflection boundary conditions to fluid variables. (Implementation specific) */
     void boundary_reflection_fluid();
     /** @brief Apply reflection boundary conditions to scalar field variables. (Implementation specific) */
     void boundary_reflection_scalar();
 
     /**
      * @brief Apply asymptotic boundary condition to a single variable at a point (for time level n).
      * @param l, k, j Grid point indices.
      * @param i Variable index.
      * @param bgv Background (asymptotic) value. (Implementation specific)
      */
     void asymcond(int l,int k,int j,int i,double bgv);
     /** @brief Apply asymptotic boundary conditions to all variables at time level n. (Implementation specific) */
     void boundary_asym0();
     /**
      * @brief Apply asymptotic boundary condition to a single variable at a point (for time level n+1, using previous levels).
      * @param l, k, j Grid point indices.
      * @param i Variable index.
      * @param bgv1 Background value at time n-1.
      * @param bgv Background value at time n.
      * @param dt Time step.
      * @param itype Type parameter (e.g., RK stage). (Implementation specific)
      */
     void asymcond(int l,int k,int j,int i,double bgv1,double bgv,double dt,int itype);
     /**
      * @brief Apply asymptotic boundary conditions to all variables at time level n+1.
      * @param itype Type parameter (e.g., RK stage). (Implementation specific)
      */
     void boundary_asym(int itype);
 
     //--------------------------------------
     // Initial Condition functions (Declarations - Specific Setups)
     //--------------------------------------
 
     /**
      * @brief Set conformal factor Psi for non-spherical initial data (static part).
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter.
      * @param xi2, xi3 Non-spherical deformation parameters. (Implementation specific)
      */
     void set_Psi_nonsph(double mu,double k,double xi2,double xi3);
     /**
      * @brief Set conformal factor Psi for non-spherical initial data (including time-dependent part).
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter.
      * @param xi2, xi3 Non-spherical deformation parameters.
      * @param xit2, xit3 Time derivatives of deformation parameters.
      * @param w Frequency parameter. (Implementation specific)
      */
     void set_Psi_nonsph(double mu,double k,double xi2,double xi3,double xit2,double xit3,double w);
     /**
      * @brief Set up full initial data for non-spherical case (including time-dependent part).
      * Calls set_Psi_nonsph and set_ini_from_Psi.
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter.
      * @param xi2, xi3 Non-spherical deformation parameters.
      * @param xit2, xit3 Time derivatives of deformation parameters.
      * @param w Frequency parameter. (Implementation specific)
      */
     void initial_nonsph(double mu,double k,double xi2,double xi3,double xit2,double xit3,double w);
     /**
      * @brief Set up full initial data for non-spherical case (static part).
      * Calls set_Psi_nonsph and set_ini_from_Psi.
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter.
      * @param xi2, xi3 Non-spherical deformation parameters. (Implementation specific)
      */
     void initial_nonsph(double mu,double k,double xi2,double xi3);
     /**
      * @brief Set initial data specifically for the scalar field part in a non-spherical setup.
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter.
      * @param xi2, xi3 Non-spherical deformation parameters. (Implementation specific)
      */
     void set_initial_scalar(double mu,double k,double xi2,double xi3);
     /**
      * @brief Set initial data specifically for the fluid part in a non-spherical setup.
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter.
      * @param xi2, xi3 Non-spherical deformation parameters. (Implementation specific)
      */
     void set_initial_fluid(double mu,double k,double xi2,double xi3);
     /**
      * @brief Set up initial data based on a mass parameter (likely spherical).
      * @param mu Mass parameter. (Implementation specific)
      */
     void initial(double mu);
     /**
      * @brief Set up initial data based on mass and decay/wave parameters (likely spherical).
      * @param mu Mass parameter.
      * @param k Wave number/decay parameter. (Implementation specific)
      */
     void initial(double mu,double k);
     /**
      * @brief Set BSSN variables based on a pre-calculated conformal factor Psi (stored in the psi array).
      * Solves Hamiltonian/Momentum constraints assuming conformal flatness and K=const. (Implementation specific)
      */
     void set_ini_from_Psi();
 
 }; // End class Fmv
 
 /**
  * @class Fmv1
  * @brief Derived class from Fmv0, implementing Adaptive Mesh Refinement (AMR).
  *
  * Manages a refinement level grid, interacting with a coarser grid (`llay`)
  * and potentially a finer grid (`ulay`). Handles interpolation for boundary
  * conditions at the coarse-fine interface and restriction of data back to
  * the coarse grid.
  */
 class Fmv1 : public Fmv0{
 private:
     int ljli,ljui,lkli,lkui,llui,llli,intp,llfmrl,lkfmrl,ljfmrl,llfmru,lkfmru,ljfmru;
     Fmv0* llay; // Pointer to the lower (coarser) layer
     Fmv1* ulay; // Pointer to the upper (finer) layer (if exists)
 
 public:
     /**
      * @brief Constructor for the Fmv1 (AMR layer) class.
      * @param tabs Number of buffer zones (ghost zones).
      * @param jupper Upper index boundary in x for the physical domain of this layer.
      * @param jlower Lower index boundary in x for the physical domain of this layer.
      * @param kupper Upper index boundary in y for the physical domain of this layer.
      * @param klower Lower index boundary in y for the physical domain of this layer.
      * @param lupper Upper index boundary in z for the physical domain of this layer.
      * @param llower Lower index boundary in z for the physical domain of this layer.
      * @param xupper Upper coordinate boundary in x for this layer.
      * @param xlower Lower coordinate boundary in x for this layer.
      * @param yupper Upper coordinate boundary in y for this layer.
      * @param ylower Lower coordinate boundary in y for this layer.
      * @param zupper Upper coordinate boundary in z for this layer.
      * @param zlower Lower coordinate boundary in z for this layer.
      * @param am Amplitude for inhomogeneous grid mapping (if used).
      * @param fld Boolean flag to enable fluid evolution.
      * @param scl Boolean flag to enable scalar field evolution.
      * @param cuev Boolean flag to enable curvature evaluation.
      * @param lolay Pointer to the Fmv0 object representing the coarser layer below this one.
      *
      * Calls the Fmv0 constructor, sets the layer number `layn`, stores the pointer
      * to the lower layer, sets AMR parameters (like interpolation order), and
      * marks the region covered by this fine grid on the lower layer using `llay->set_fmrregion`.
      */
     Fmv1(int tabs,int jupper,int jlower,int kupper,int klower,int lupper,int llower,
     double xupper,double xlower,double yupper,double ylower,double zupper,double zlower,double am, bool fld, bool scl, bool cuev, Fmv0* lolay);
 
     /**
      * @brief Destructor for the Fmv1 class.
      */
      ~Fmv1();
 
     //--------------------------------------
     // AMR specific Getters/Setters
     //--------------------------------------
     /** @brief Get lower x-index boundary of the physical domain on this layer. */
     int get_ljli() const;
     /** @brief Set lower x-index boundary of the physical domain on this layer. */
     void set_ljli(int p);
     /** @brief Get lower y-index boundary of the physical domain on this layer. */
     int get_lkli() const;
     /** @brief Set lower y-index boundary of the physical domain on this layer. */
     void set_lkli(int p);
     /** @brief Get lower z-index boundary of the physical domain on this layer. */
     int get_llli() const;
     /** @brief Set lower z-index boundary of the physical domain on this layer. */
     void set_llli(int p);
     /** @brief Set the pointer to the upper (finer) refinement layer. */
     void set_ulay(Fmv1* uplay);
 
     //--------------------------------------
     // AMR specific functions (Declarations)
     //--------------------------------------
 
     /**
      * @brief Set boundary conditions for this refinement level.
      * Handles physical boundaries (if this level extends to them) and
      * coarse-fine interface boundaries (interpolating from `llay`).
      * @param btype Type of physical boundary condition.
      * @param mm Parameter related to boundary condition type or stage. (Implementation specific)
      */
     void set_boundary(int btype,int mm);
     /**
      * @brief Initialize data on this fine grid by interpolating from the coarse grid (`llay`).
      * Called when the refinement level is first created. (Implementation specific)
      */
     void set_fmr_initial();
     /**
      * @brief Evolve this refinement level, potentially involving sub-cycling in time.
      * Coordinates boundary data exchange with `llay` and `ulay`. (Implementation specific)
      */
     void evolve();
     /**
      * @brief Restrict data from this fine grid back to the underlying coarse grid (`llay`).
      * Averages or injects data into the `fmrflag` region of `llay`. (Implementation specific)
      */
     void refine_llay();
     /**
      * @brief Perform one time step evolution on this refinement level.
      * Includes calculating RHS, applying boundary conditions, and updating variables.
      * @param btype Type of physical boundary condition. (Implementation specific)
      */
     void onestep(int btype);
     /**
      * @brief Perform time interpolation for boundary conditions at coarse-fine interfaces.
      * Used during sub-cycling when the fine grid needs data from the coarse grid
      * at intermediate times. Interpolates between bv0 and bv1 of the coarse grid.
      * @param l, k, j Fine grid point indices needing boundary data.
      * @param ll, kk, jj Corresponding coarse grid indices.
      * @param i Variable index.
      * @param aa, bb, cc Time interpolation coefficients. (Implementation specific)
      */
     void tstep_ipol(int l,int k,int j,int ll, int kk, int jj, int i,double aa,double bb,double cc);
 
 }; // End class Fmv1
 
 
 #endif // _COSMOS_H_
 