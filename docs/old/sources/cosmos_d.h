/**
 * @file cosmos.h
 * @brief Header file defining the core BSSN evolution classes (Fmv0, Fmv, Fmv1) for the COSMOS code.
 * @details Contains base class Fmv0 for grid variables and evolution,
 *          derived class Fmv for specific boundary conditions/initial data,
 *          and derived class Fmv1 for Fixed Mesh Refinement (FMR).
 * @version 1.00 // Or your current version
 * @date 2023-10-27 // Or the relevant date
 * @author Chulmoon Yoo // Or the primary author(s)
 */

 #ifndef COSMOS_H // Ensure include guards are present
 #define COSMOS_H
 
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <string.h>
 // Include other necessary headers like complex.h if used
 
 // --- Constants and Macros (Example) ---
 /** @brief Gravitational constant G (often set to 1 in simulations). */
 #define G_CONST 1.0
 /** @brief Speed of light c (often set to 1 in simulations). */
 #define C_LIGHT 1.0
 // Add other constants/macros if present
 
 // --- Enums (Example if you have them) ---
 /** @brief Enum defining possible boundary condition types. */
 typedef enum {
     BOUNDARY_PERIODIC,
     BOUNDARY_REFLECTING,
     BOUNDARY_OUTFLOW
     // Add other boundary types
 } BoundaryType;
 
 
 // --- Fmv0 Class ---
 /**
  * @class Fmv0
  * @brief Base class for Finite Mesh Variables version 0.
  * @details Handles the fundamental grid structure, memory allocation,
  *          coordinate mapping, basic BSSN evolution steps (RHS calculation),
  *          and Kreiss-Oliger dissipation. Does not implement specific
  *          boundary conditions or mesh refinement.
  */
 class Fmv0 {
 public: // Or protected, depending on your design
     // --- Member Variables (Document important public/protected ones) ---
     /** @brief Number of buffer/ghost zones. */
     int tab;
     /** @brief Minimum physical grid index in x. */
     int jmin;
     /** @brief Maximum physical grid index in x. */
     int jmax;
     // ... (similarly for kmin, kmax, lmin, lmax)
     /** @brief Minimum coordinate value in x. */
     double xmin;
     /** @brief Maximum coordinate value in x. */
     double xmax;
     // ... (similarly for ymin, ymax, zmin, zmax)
     /** @brief Grid spacing in x (uniform logical grid). */
     double dx;
     // ... (similarly for dy, dz)
     /** @brief CFL factor for time step calculation. */
     double cfl;
     /** @brief Gauge parameter eta for 1+log lapse. */
     double etaa;
     /** @brief Gauge parameter eta for Gamma-driver shift. */
     double etab;
     /** @brief Gauge parameter eta_b (or similar) for Gamma-driver shift. */
     double etabb;
     /** @brief Kreiss-Oliger dissipation coefficient epsilon. */
     double KOep;
     /** @brief Initial Hubble parameter. */
     double Hb;
     /** @brief Fluid equation of state parameter w (if fluid enabled). */
     double fluidw;
     /** @brief MUSCL reconstruction parameter kappa (if fluid enabled). */
     double kap_MUSCL;
     /** @brief Minmod limiter parameter b (if fluid enabled). */
     double b_minmod;
     /** @brief Grid size for excision region (if enabled). */
     int exg;
     /** @brief Maximum simulation time. */
     double tmax;
     /** @brief Grid stretching parameter for coordinate mapping. */
     double amp;
 
     // Pointers to grid function arrays (Example for lapse 'alp')
     /** @brief Pointer to the lapse function grid data. */
     double ***alp;
     /** @brief Pointer to the shift vector component beta^x grid data. */
     double ***betax;
     // ... (Document ALL BSSN variables: betay, betaz, Gamh1, Gamh2, Gamh3, trK, A11, A12, ...)
     // ... (Document scalar field variables: sf, Pi)
     // ... (Document fluid variables: rho, Sx, Sy, Sz, tau)
     // ... (Document auxiliary variables: Bx, By, Bz)
     // ... (Document boundary flags: bflag)
 
     // --- Constructor & Destructor ---
     /**
      * @brief Constructor for the Fmv0 class.
      * @param tabs Number of buffer zones (ghost zones).
      * @param jupper Upper index boundary in x for the physical domain.
      * @param jlower Lower index boundary in x for the physical domain.
      * @param kupper Upper index boundary in y.
      * @param klower Lower index boundary in y.
      * @param lupper Upper index boundary in z.
      * @param llower Lower index boundary in z.
      * @param xupper Upper coordinate boundary in x.
      * @param xlower Lower coordinate boundary in x.
      * @param yupper Upper coordinate boundary in y.
      * @param ylower Lower coordinate boundary in y.
      * @param zupper Upper coordinate boundary in z.
      * @param zlower Lower coordinate boundary in z.
      * @param cfl_in CFL factor.
      * @param etaa_in Gauge parameter eta for lapse.
      * @param etab_in Gauge parameter eta for shift damping (B^i).
      * @param etabb_in Gauge parameter eta_b for shift evolution.
      * @param KOep_in Kreiss-Oliger dissipation coefficient.
      * @param exg_in Excision region grid size.
      * @param Hb_in Initial Hubble parameter.
      * @param fluidw_in Fluid equation of state parameter w.
      * @param kap_MUSCL_in MUSCL kappa parameter.
      * @param b_minmod_in Minmod b parameter.
      * @param amp_in Grid stretching parameter.
      */
     Fmv0(int tabs, int jupper, int jlower, int kupper, int klower, int lupper, int llower,
          double xupper, double xlower, double yupper, double ylower, double zupper, double zlower,
          double cfl_in, double etaa_in, double etab_in, double etabb_in, double KOep_in, int exg_in,
          double Hb_in, double fluidw_in, double kap_MUSCL_in, double b_minmod_in, double amp_in);
 
     /**
      * @brief Destructor for the Fmv0 class. Frees allocated memory.
      */
     virtual ~Fmv0(); // Make destructor virtual if inheritance is used
 
     // --- Core Methods ---
     /**
      * @brief Calculates the right-hand sides (RHS) of the BSSN evolution equations.
      * @param dt Current time step size (may be needed for some terms).
      * @param t Current simulation time.
      * @note This is typically the most complex function, calculating spatial derivatives
      *       and combining terms for each evolved variable.
      * @warning Ensure boundary conditions are applied *before* calling this for points
      *          needing boundary data.
      */
     virtual void rhs(double dt, double t); // Make virtual if overridden
 
     /**
      * @brief Applies Kreiss-Oliger numerical dissipation to evolved variables.
      * @param q Pointer to the grid function data array to apply dissipation to.
      * @param f Pointer to the RHS array corresponding to q (dissipation added here).
      */
     void KOdiss(double ***q, double ***f);
 
     /**
      * @brief Coordinate mapping function from logical coordinate (e.g., j*dx) to physical coordinate.
      * @param x Logical coordinate.
      * @return Physical coordinate.
      */
     double funcf(double x);
     /** @brief First derivative of the coordinate mapping function funcf. */
     double df(double x);
     /** @brief Second derivative of the coordinate mapping function funcf. */
     double ddf(double x);
     /** @brief Third derivative of the coordinate mapping function funcf. */
     double dddf(double x);
 
     // --- Getters/Setters (Examples) ---
     /** @brief Get the number of buffer zones. */
     int get_tab() const;
     /** @brief Get the CFL factor. */
     double get_cfl() const;
     /** @brief Check if fluid evolution is enabled (based on fluidw perhaps). */
     bool get_fluidevo() const; // Implement logic based on your setup
     /**
      * @brief Set the boundary flag at a specific grid point. Returns a reference.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to bflag[l-lmin][k-kmin][j-jmin].
      */
     int& set_bflag(int l, int k, int j);
     /**
      * @brief Get the boundary flag at a specific grid point.
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Value of bflag[l-lmin][k-kmin][j-jmin].
      */
     int get_bflag(int l, int k, int j) const;
 
     // --- Utility Functions (Examples) ---
     /** @brief Allocates memory for a 3D grid function array. */
     double*** matrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
     /** @brief Frees memory allocated by matrix(). */
     void free_matrix(double ***m, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
     // Add similar functions for int*** if needed (like for bflag)
 
 protected: // Or private
     // --- Helper functions for RHS (Examples) ---
     /** @brief Calculates spatial derivatives (e.g., first derivative in x). */
     double D1x(double ***f, int l, int k, int j); // Add similar for y, z, higher orders
     /** @brief Calculates Christoffel symbols or connection functions. */
     void calculate_christoffel(int l, int k, int j /*, output parameters */);
     /** @brief Calculates Ricci tensor components. */
     void calculate_ricci(int l, int k, int j /*, output parameters */);
     /** @brief Calculates matter source terms (Stress-Energy Tensor T_munu). */
     void calculate_matter_sources(int l, int k, int j /*, output parameters */);
 
     // Add other protected members/methods if necessary
 };
 
 
 // --- Fmv Class ---
 /**
  * @class Fmv
  * @brief Derived class from Fmv0, implementing specific boundary conditions
  *        and initial data routines for the COSMOS simulation.
  * @details Overrides virtual functions like `set_boundary` and provides
  *          methods like `initial_data` or `initial_nonsph`.
  */
 class Fmv : public Fmv0 {
 public:
     /**
      * @brief Constructor for the Fmv class.
      * @details Passes most parameters up to the Fmv0 constructor. May add
      *          specific parameters for this class if needed.
      * @param [in] params Parameters identical to Fmv0 constructor.
      */
     Fmv(/* Parameters matching Fmv0 constructor */); // List all params again
 
     /**
      * @brief Destructor for the Fmv class.
      */
     virtual ~Fmv();
 
     // --- Overridden/Specific Methods ---
     /**
      * @brief Sets the boundary conditions for all evolved grid functions.
      * @details Implements the specific boundary logic (e.g., periodic, reflection, outflow)
      *          for the COSMOS simulation setup. Called typically after each RK step.
      */
     virtual void set_boundary(); // Override if Fmv0::set_boundary is virtual
 
     /**
      * @brief Sets up the initial data for the simulation.
      * @param t Initial time (usually 0).
      * @param mu Initial geometric perturbation amplitude.
      * @param kk Initial geometric perturbation scale/wavenumber.
      * @param xi2 Initial geometric perturbation non-spherical parameter.
      * @param xi3 Initial geometric perturbation non-spherical parameter.
      * @param w3 Initial geometric perturbation additional parameter.
      * @note Calls helper functions like initial_nonsph, set_initial_scalar.
      */
     void initial_data(double t, double mu, double kk, double xi2, double xi3, double w3);
 
     /**
      * @brief Sets the initial data for the geometric (BSSN) variables with non-spherical perturbations.
      * @param t Initial time.
      * @param mu Amplitude parameter.
      * @param kk Scale/wavenumber parameter.
      * @param xi2 Non-spherical parameter.
      * @param xi3 Non-spherical parameter.
      * @param w3 Additional parameter.
      */
     void initial_nonsph(double t, double mu, double kk, double xi2, double xi3, double w3);
 
     /**
      * @brief Sets the initial data for the scalar field variables.
      * @param mus Amplitude parameter for scalar field perturbation.
      * @param kks Scale/wavenumber parameter for scalar field perturbation.
      * @param xi2s Non-spherical parameter for scalar field perturbation.
      * @param xi3s Non-spherical parameter for scalar field perturbation.
      */
     void set_initial_scalar(double mus, double kks, double xi2s, double xi3s);
 
     /**
      * @brief Sets the initial data for the fluid variables (if enabled).
      * @param t Initial time.
      * // Add parameters relevant to fluid initial data if needed
      */
     void set_initial_fluid(double t /*, ... */);
 
     // Add other specific methods for this class
 };
 
 
 // --- Fmv1 Class ---
 /**
  * @class Fmv1
  * @brief Derived class from Fmv0, implementing Fixed Mesh Refinement (FMR).
  * @details Manages multiple refinement levels, interpolation between levels,
  *          and potentially modified evolution steps considering refinement boundaries.
  * @note This class structure might be significantly more complex depending on the
  *       FMR implementation details (e.g., using pointers to Fmv0/Fmv objects for each level).
  *       The documentation below is a basic placeholder.
  */
 class Fmv1 : public Fmv0 {
 public:
     /** @brief Maximum number of refinement levels. */
     int laymax;
     /** @brief Array of pointers to Fmv0/Fmv objects representing each refinement level. */
     Fmv0 **fmv; // Or Fmv **fmv;
     /** @brief Array storing x-boundary indices for refinement regions on coarser levels. */
     int *jbs;
     /** @brief Array storing y-boundary indices for refinement regions on coarser levels. */
     int *kbs;
     /** @brief Array storing z-boundary indices for refinement regions on coarser levels. */
     int *lbs;
     /** @brief Array storing lapse threshold values for triggering refinement. */
     double *alp_fmr;
 
     /**
      * @brief Constructor for the Fmv1 class (FMR manager).
      * @param laymax_in Maximum number of refinement levels.
      * @param jbs_in Array of x-boundary indices.
      * @param kbs_in Array of y-boundary indices.
      * @param lbs_in Array of z-boundary indices.
      * @param alp_fmr_in Array of lapse trigger thresholds.
      * @param [in] base_params Parameters needed to construct the base level (level 0),
      *                       likely similar to Fmv0/Fmv constructor.
      */
     Fmv1(int laymax_in, int *jbs_in, int *kbs_in, int *lbs_in, double *alp_fmr_in /*, base level params */);
 
     /**
      * @brief Destructor for the Fmv1 class. Deallocates levels.
      */
     virtual ~Fmv1();
 
     /**
      * @brief Performs one full evolution step across all refinement levels.
      * @details Handles Runge-Kutta steps, boundary conditions, inter-level
      *          interpolation (prolongation/restriction), and potentially
      *          checks for creating new refinement levels.
      * @param dt Current time step size for the finest level.
      * @param t Current simulation time.
      */
     void evolution(double dt, double t);
 
     /**
      * @brief Checks if new refinement levels need to be created based on triggers.
      * @param t Current simulation time.
      */
     void check_refinement_trigger(double t);
 
     /**
      * @brief Creates a new refinement level.
      * @param level Index of the new level to create.
      */
     void create_refinement_level(int level);
 
     /**
      * @brief Interpolates data from a coarser level to a finer level (Prolongation).
      * @param coarse_level Index of the coarser level.
      * @param fine_level Index of the finer level.
      */
     void prolongate(int coarse_level, int fine_level);
 
     /**
      * @brief Averages data from a finer level to a coarser level (Restriction).
      * @param fine_level Index of the finer level.
      * @param coarse_level Index of the coarser level.
      */
     void restrict(int fine_level, int coarse_level);
 
     // Add other FMR-specific methods
 };
 
 
 #endif // COSMOS_H
 