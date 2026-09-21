#pragma once
#include "domain.hpp"
#include <vector>
#include "params.hpp"

void compute_phi(const Domain& domain,
			   	std::vector<double>& phi,
			   	double cx, double cy, double sigma);
double compute_mass(const Domain& domain,
			   	const std::vector<double>& phi);
double compute_local_maturity(const Domain& domain,
			   	const std::vector<double>& phi);
double compute_global_fsh(const GlobalBioParams& global_bio_params,
			   	const double global_maturity);
double compute_local_fsh(const FollicleParams& follicle_params,
			   	const double local_maturity,
				const double global_fsh);
void compute_g_h_source(const Domain& domain, 
				std::vector<double>& g,
			   	std::vector<double>& h,
				std::vector<double>& source,
			   	const GlobalBioParams& global_bio_params, 
				const FollicleParams& follicle_params,
				const double local_fsh,
			   	const double global_fsh);
