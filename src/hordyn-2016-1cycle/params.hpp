// 2012
#pragma once
#include <cmath>

struct GlobalBioParams {
		// apoptosis function parameter
		double Lambda_bar = 0.2;//6.0;
		double gamma_bar = 0.1;//0.02; 

		// plasma fsh function parameter
        double U_max = 0.15;//1.0;
		double U_min = 0.075;
		double c = 2.0;
		double M_ref = 4.5;

		// threshold in maturity axis (y_s)
		double ys = 0.3;
};

struct FollicleParams {
		// location
		double cx = 0.25;
		double cy = 0.15;
		double sigma = std::sqrt(0.002); 
        // 
        double mu1 = 0.1;
        double mu2 = 0.4;

		// aging function parameters
		double gamma1 = 4.0;
		double gamma2 = 1.2;

		// maturation function parameters
	 	// maturation velocity rate 
		double tau_h = 0.4;//1.0;
		double c1 = 0.01;//0.68;
		double c2 = 1.2;//0.08;
		// basal fsh scale
		double u_bar = 0.1;//0.02; 

		// local fsh function parameter
		double b1 = 0.08;
		double b2 = 2.25;
		double b3 = 1450.0;

};
