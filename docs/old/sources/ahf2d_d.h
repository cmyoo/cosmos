/**
 * @file ahf2d.h
 * @brief Header file for the 2D Apparent Horizon Finder (AHF) class.
 * @details Implements methods to locate apparent horizons on spatial slices
 *          using a 2D spherical coordinate grid (theta, phi) and an
 *          iterative elliptic solver.
 * @version 1.0 // Or your current version
 * @date 2023-10-27 // Or the relevant date
 * @author Chulmoon Yoo // Or the primary author(s)
 */

 #ifndef AHF2D_H
 #define AHF2D_H
 
 #include "cosmos.h" // Include if it needs access to Fmv0/Fmv grid data/methods
 #include <stdio.h>
 #include <math.h>
 // Include other necessary headers
 
 /**
  * @class Ahf2d
  * @brief Implements a 2D Apparent Horizon Finder.
  * @details Solves the apparent horizon equation Theta = 0, where Theta is the
  *          expansion of outgoing null geodesics, on a spherical coordinate grid.
  *          Uses an iterative method involving a Poisson-like solver.
  */
 class Ahf2d {
 public:
     // --- Member Variables ---
     /** @brief Number of grid points in theta direction (excluding buffers). */
     int ntheta;
     /** @brief Number of grid points in phi direction (excluding buffers). */
     int nphi;
     /** @brief Grid spacing in theta. */
     double dth;
     /** @brief Grid spacing in phi. */
     double dphi;
     /** @brief Initial guess for the horizon radius (coordinate radius). */
     double hini;
     /** @brief Acceleration parameter eta for the elliptic solver. Also distinguishes AH/CH. */
     double var; // Corresponds to etaah in par_ahf.d
     /** @brief Mixing factor for updating horizon radius guess during iteration. */
     double fac; // Corresponds to facn in par_ahf.d
     /** @brief Tolerance for the Poisson solver part. */
     double err_p;
     /** @brief Tolerance for the main AH equation solver part. */
     double err_e;
     /** @brief Maximum number of iterations allowed for the solver. */
     int ahloop;
 
     // Pointers to AHF internal grid arrays (Example for horizon shape 'h')
     /** @brief Horizon shape function h(theta, phi). */
     double **h;
     /** @brief Expansion Theta(theta, phi). */
     double **ths;
     // Add other internal arrays if needed (e.g., for solver residuals, coefficients)
 
     /** @brief Pointer to the Fmv0/Fmv object containing the spacetime grid data. */
     Fmv0 *fm; // Use Fmv0* or Fmv* depending on what data is needed
 
     // --- Constructor & Destructor ---
     /**
      * @brief Constructor for the Ahf2d class.
      * @param ntheta_in Number of theta grid points.
      * @param nphi_in Number of phi grid points.
      * @param hini_in Initial guess for horizon radius.
      * @param var_in Solver acceleration parameter eta (etaah).
      * @param fac_in Radius update mixing factor (facn).
      * @param err_p_in Poisson solver tolerance.
      * @param err_e_in AH equation solver tolerance.
      * @param ahloop_in Maximum solver iterations.
      * @param fm_in Pointer to the Fmv0/Fmv object with spacetime data.
      */
     Ahf2d(int ntheta_in, int nphi_in, double hini_in, double var_in, double fac_in,
           double err_p_in, double err_e_in, int ahloop_in, Fmv0 *fm_in);
 
     /**
      * @brief Destructor for the Ahf2d class. Frees allocated memory.
      */
     ~Ahf2d();
 
     // --- Core Methods ---
     /**
      * @brief Main function to find the apparent horizon.
      * @details Iteratively solves for the horizon shape function h(theta, phi)
      *          until the expansion Theta is sufficiently close to zero or
      *          the maximum number of iterations is reached.
      * @param[out] ah_mass Pointer to store the calculated irreducible mass (optional).
      * @param[out] ah_spin Pointer to store the calculated dimensionless spin (optional).
      * @return 0 on success, non-zero on failure (e.g., no convergence).
      */
     int find_ah(double *ah_mass = nullptr, double *ah_spin = nullptr); // Example output params
 
     /**
      * @brief Calculates the expansion Theta at all points on the AHF grid.
      * @details Requires interpolating spacetime variables (metric, extrinsic curvature)
      *          from the Cartesian grid (fm) to the current guess of the horizon
      *          surface defined by h(theta, phi).
      */
     void calc_expansion();
 
     /**
      * @brief Solves the elliptic (Poisson-like) equation arising in the iteration.
      * @details This is likely the core iterative solver step (e.g., SOR, CG).
      * @param[in] source The source term for the elliptic equation (related to Theta).
      * @param[out] solution The solution array (update to h or correction).
      * @return Number of iterations taken, or -1 on failure.
      */
     int solve_elliptic(double **source, double **solution);
 
     /**
      * @brief Interpolates a spacetime grid function from the Cartesian grid to a point on the AH surface.
      * @param f Pointer to the 3D Cartesian grid function data (e.g., fm->alp).
      * @param r Coordinate radius of the point.
      * @param th Theta angle of the point.
      * @param ph Phi angle of the point.
      * @return Interpolated value of f at (r, th, ph).
      * @note Needs careful implementation of coordinate transformation and interpolation (e.g., trilinear).
      */
     double interpolate_to_ah(double ***f, double r, double th, double ph);
 
     /**
      * @brief Calculates the irreducible mass of the found horizon.
      * @return Irreducible mass M_irr = sqrt(Area / (16 * PI)). Returns -1 if horizon not found.
      */
     double calculate_ah_mass();
 
     /**
      * @brief Calculates the dimensionless spin of the found horizon.
      * @return Dimensionless spin parameter J/M^2. Returns -1 if horizon not found or spin calculation fails.
      * @note Spin calculation often requires integrating the extrinsic curvature over the horizon surface.
      */
     double calculate_ah_spin();
 
     /**
      * @brief Prints the current horizon shape h(theta, phi) to a file.
      * @param filename Name of the output file.
      * @param t Current simulation time to include in the output.
      */
     void print_horizon_shape(const char* filename, double t);
 
     // --- Utility Functions ---
     /** @brief Allocates memory for a 2D AHF grid array. */
     double** matrix2d(long nrl, long nrh, long ncl, long nch);
     /** @brief Frees memory allocated by matrix2d(). */
     void free_matrix2d(double **m, long nrl, long nrh, long ncl, long nch);
 
 protected: // Or private
     // Add helper functions if needed, e.g., for calculating derivatives on the spherical grid
     /** @brief Calculates the theta derivative of a 2D AHF grid function. */
     double D1th(double **f, int i, int j);
     /** @brief Calculates the phi derivative of a 2D AHF grid function. */
     double D1ph(double **f, int i, int j);
     // Add second derivatives if needed (D2th, D2ph, Dthph)
 };
 
 
 #endif // AHF2D_H
 