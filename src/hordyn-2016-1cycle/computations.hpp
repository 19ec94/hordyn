#pragma once
#include "domain.hpp"
#include "params.hpp"

#include <vector>
#include <cmath>
#include <numbers>

void compute_phi(const Domain& domain, std::vector<double>& phi,
        const FollicleParams& follicle_params);
double compute_mass(const Domain& domain, const std::vector<double>& phi);
double compute_mf(const Domain& domain, const std::vector<double>& phi);
double compute_U(const GlobalBioParams& global_bio_params, const double M);
double compute_uf(const FollicleParams& follicle_params, const double mf, const double U);
void compute_g_h(const Domain& domain,
        std::vector<double>& g, std::vector<double>& h,
        const FollicleParams& follicle_params,
        const double uf) ;
void compute_source(const Domain& domain,
        std::vector<double>& source,
        const GlobalBioParams& global_bio_params, const double U);
