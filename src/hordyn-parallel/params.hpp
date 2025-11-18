// 2012
#pragma once
#include <cmath>

struct GlobalBioParams {
		// Apoptosis function parameter
		double Lambda_bar = 6.0;
		double gamma_bar = 0.02; 

		// plasma FSH function parameter
		double U_min = 0.075;
		double c = 2.0;
		double M_ref = 4.5;

		// Threshold in maturity axis (y_s)
		double y_s = 0.3;
};

struct FollicleParams {
		// Location
		double cx = 0.25;
		double cy = 0.15;
		double sigma = std::sqrt(0.002); 

		// Aging function parameters
		double gamma1 = 4.0;
		double gamma2 = 1.2;

		// Maturation function parameters
	 	// maturation velocity rate 
		double tau_h = 1.0;
		double c1 = 0.68;
		double c2 = 0.08;
		// basal FSH scale
		double u_bar = 0.02; 

		// local FSH function parameter
		double b1 = 0.08;
		double b2 = 2.25;
		double b3 = 1450.0;

};
