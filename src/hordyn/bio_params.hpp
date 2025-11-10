#pragma once

// 2012
struct BioParams {
		// Aging function parameters
		const double gamma1 = 4.0;
		const double gamma2 = 1.2;

		// Maturation function parameters
		const double tau_h = 1.2;   // maturation velocity rate 
		const double c1 = 0.68;
		const double c2 = 0.08;
		// basal FSH scale
		const double u_bar = 0.02; 

		// Apoptosis function parameter
		const double Lambda_bar = 6.0;
		const double gamma_bar = 0.02; 

		// plasma FSH function parameter
		const double U_min = 0.075;
		const double c = 2.0;
		const double M_ref = 4.5;

		// local FSH function parameter
		const double b1 = 0.08;
		const double b2 = 2.25;
		const double b3 = 1450.0;

		// Uninitialized
		double y_s; // Threshold in maturity axis (y_s)
		BioParams(double y_s_);
};

/*
struct BioParams {
		// Aging function parameters
		const double gamma1 = 2.0;
		const double gamma2 = 2.0;

		// Maturation function parameters
		const double tau_h = 0.7;   // maturation velocity rate 
		const double c1 = 0.68;
		const double c2 = 0.08;
		// basal FSH scale
		const double u_bar = 0.02; 

		// Apoptosis function parameter
		const double Lambda_bar = 0.1;
		const double gamma_bar = 0.01; 

		// plasma FSH function parameter
		const double U_min = 0.5;
		const double c = 2.0;
		const double M_ref = 4.5;

		// local FSH function parameter
		const double b1 = 0.08;
		const double b2 = 2.25;
		const double b3 = 1450.0;

		// Uninitialized
		double y_s; // Threshold in maturity axis (y_s)
		BioParams(double y_s_);
};
*/
