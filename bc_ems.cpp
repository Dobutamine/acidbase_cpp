#include <iostream>
#include <limits>    // for std::numeric_limits
#include <cmath>     // for std::abs, std::log10, std::exp, std::tanh
#include <functional>

// Normal neonatal p50_0: https://pmc.ncbi.nlm.nih.gov/articles/PMC6626707/
/*). 
Hgb F, the predominant hemoglobin in newborns, has a high oxygen affinity and, as a result, a lower P50 
(partial pressure of oxygen at which 50% of Hgb is saturated by oxygen) of ≈ 18–19 mm Hg, 
whereas adult Hgb A has a P50 of ≈ 26–27 mm Hg (Emond et al. 1993; Maurer et al. 1970). 
*/
// include emscripten dependencies
#include <emscripten.h>
#include <emscripten/bind.h>

// -----------------------------------------------------------------------------
// Compile with command : 
// emcc bc_ems.cpp -o bc_ems.js  -O3  -s ENVIRONMENT="web,worker" --std=c++17  --bind 
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Struct for input parameters
// -----------------------------------------------------------------------------
struct Input {
    double tco2;
    double to2;
    double temp;   
    double hemoglobin;
    double na;
    double k;
    double ca;
    double mg;
    double cl;
    double lact;
    double albumin;
    double phosphates;
    double dpg;
    double uma;
    double prev_ph;
    double prev_po2;
    double p50_0; // P50_0 = 18.8 for fetal hemoglobin, 26.7 for adult hemoglobin
};

// -----------------------------------------------------------------------------
// Struct for computed outputs
// -----------------------------------------------------------------------------
struct Output {
    double ph;
    double pco2;
    double hco3;
    double be;
    double po2;
    double so2;
    int iterations;
    bool error;
};

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
constexpr double gas_constant    = 62.36367; // L·mmHg/(K·mol)
constexpr double kw              = 2.5119e-11; // water dissociation constant
constexpr double kc              = 7.94328235e-4; // carbonic acid dissociation constant
constexpr double kd              = 6.0255959e-8; // bicarbonate dissociation constant
constexpr double alpha_co2p      = 0.03067; // CO2 solubility coefficient
constexpr double left_hp_wide    = 5.848931925e-6; // lower bound for H⁺ concentration
constexpr double right_hp_wide   = 3.16227766017e-4; // upper bound for H⁺ concentration
constexpr double delta_ph_limits = 0.1; // delta for pH limits
constexpr double n               = 2.7; // Hill coefficient
constexpr double alpha_o2        = 1.38e-5; // O2 solubility coefficient
constexpr double left_o2_wide    = 0; // lower bound for pO2
constexpr double right_o2_wide   = 800.0; // upper bound for pO2
constexpr double delta_o2_limits = 10.0; // delta for pO2 limits
constexpr double brent_accuracy  = 1e-8;
constexpr int    max_iterations  = 100;

// -----------------------------------------------------------------------------
// Brent’s method root‐finder
// -----------------------------------------------------------------------------
template<typename F>
double brent_root_fast(F f, double a, double b, double tolerance, int max_iterations);

// -----------------------------------------------------------------------------
// Core computation functions
// -----------------------------------------------------------------------------
double net_charge_plasma(double hp_estimate);
double calc_so2(double po2_estimate);
double dO2_content(double po2_estimate);

Output calc_blood_composition(Input _input);

// initialize state variables
double tco2{};
double to2{};
double sid{};
double albumin{};
double phosphates{};
double uma{};
double hemoglobin{};
double temp{};
double dpg{};
double ph{};
double pco2{};
double hco3{};
double be{};
double po2{};
double so2{};
double prev_po2{};
double prev_ph{};
bool error{};
int iterations{};
bool error_ab{};
bool error_oxy{};
int iterations_ab{};
int iterations_oxy{};

double P50_0 = 20.0; // PO2 at which 50% of Hgb is saturated by O2 (fetal = 18.8 (high Hb O2 affinity), neonatal = 20.0, adult = 26.7)
double P50{};
double log10_p50{};
double P50_n{};
double left_o2 = 0; // lower bound for pO2
double right_o2 = 800.0; // upper bound for pO2
double left_hp = 5.848931925e-6; // lower bound for H⁺ concentration
double right_hp = 3.16227766017e-4; // upper bound for H⁺ concentration

// main function to calculate blood gas parameters
Output calc_blood_composition(Input _input) {

    // initialize state from inputs
    tco2       = _input.tco2;
    to2        = _input.to2;
    sid        = _input.na + _input.k + 2*_input.ca + 2*_input.mg - _input.cl - _input.lact;
    albumin    = _input.albumin;
    phosphates = _input.phosphates;
    uma        = _input.uma;
    hemoglobin = _input.hemoglobin;
    temp       = _input.temp;
    dpg        = _input.dpg;
    prev_po2   = _input.prev_po2;
    prev_ph    = _input.prev_ph;
    P50_0      = _input.p50_0;

    // initialize the output struct
    Output out; 

    // reset the iteration and error counters
    iterations = 0;
    error = false;

    // default to -1 if root-finding fails
    out.ph = out.pco2 = out.hco3 = out.be = out.po2  = out.so2 = -1;

    // ACIDBASE
    bool limits_used_ab = false;
    // set the limits based
    left_hp = left_hp_wide; // lower bound for H⁺ concentration
    right_hp = right_hp_wide; // upper bound for H⁺ concentration
    // set the limits based on the previous calculations if available
    if (prev_ph > 0) {
        left_hp = std::pow(10.0, -(prev_ph + delta_ph_limits)) * 1000.0;
        right_hp = std::pow(10.0, -(prev_ph - delta_ph_limits)) * 1000.0;
        limits_used_ab = true;
    }

    // 1) Solve for [H⁺]
    double hp_est = brent_root_fast([](double x){ return net_charge_plasma(x); }, left_hp, right_hp, brent_accuracy, max_iterations);

    // check if the root-finding was successful
    iterations_ab = iterations;
    error_ab = error;

    // if not successful and limits were used, try again with wider limits
    if (error_ab && limits_used_ab) {
        left_hp = left_hp_wide; // wide lower bound for H⁺ concentration
        if (left_hp < 0) left_hp = 0; // ensure lower bound is not negative
        right_hp = right_hp_wide; // wide upper bound for H⁺ concentration
        hp_est = brent_root_fast([](double x){ return net_charge_plasma(x); }, left_hp, right_hp, brent_accuracy, max_iterations);
        iterations_ab = iterations;
        error_ab = error;
    }

    // std::cout << "iterations acid base: " << iterations_ab << "\n";
    // std::cout << "error acid base: " << error_ab << "\n";

    if (hp_est > 0) {
        // compute base excess
        be = (hco3 - 25.1 + (2.3*hemoglobin + 7.7)*(ph - 7.4)) * (1.0 - 0.023*hemoglobin);
        out.ph   = ph;
        out.pco2 = pco2;
        out.hco3 = hco3;
        out.be   = be;
    }

    // reset the iteration and error counters
    iterations = 0;
    error = false;

    // OXYGENATION
     // Compute shifts from standard conditions
    double dpH = ph - 7.40;             //Bohr effect: ↓pH → right shift → ↑P₅₀)
    double dpCO2 = pco2 - 40.0;         // Haldane effect: ↑pCO2 → right shift → ↑P₅₀
    double dT = temp - 37.0;            // ↑T → right shift → ↑P₅₀
    double dDPG = dpg - 5.0;            // ↑DPG → right shift → ↑P₅₀

    // Adjust P50 on log10 scale
    log10_p50 = std::log10(P50_0) - 0.48 * dpH + 0.014 * dpCO2 + 0.024 * dT + 0.051 * dDPG;
    P50 = std::pow(10.0, log10_p50);
    P50_n = std::pow(P50, n);

    // set the limits based on the previous calculations
    bool limits_used_oxy = false;
    // reset the iteration and error counters
    iterations = 0;
    error = false;
    // set the wide limits as the previous limits were not successful
    left_o2 = left_o2_wide; // lower bound po2
    right_o2 = right_o2_wide; // upper bound po2
    // set the limits based on the previous calculations if available
    if (prev_po2 > 0) {
        left_o2 = prev_po2 - delta_o2_limits;
        if (left_o2 < 0) left_o2 = 0; // ensure lower bound is not negative
        right_o2 = prev_po2 + delta_o2_limits;
        limits_used_oxy = true;
    }

    // 2) Solve for pO2
    double po2_est = brent_root_fast([](double x){ return dO2_content(x); }, left_o2, right_o2, brent_accuracy, max_iterations);

    // check if the root-finding was successful
    iterations_oxy = iterations;
    error_oxy = error;

    // if not successful and limits were used, try again with wider limits
    if (error_oxy && limits_used_oxy) {
        // reset the iteration and error counters
        iterations = 0;
        error = false;
        // set the wide limits as the previous limits were not successful
        left_o2 = left_o2_wide; // wide lower bound po2
        right_o2 = right_o2_wide; // wide upper bound po2
        po2_est = brent_root_fast([](double x){ return dO2_content(x); }, left_o2, right_o2, brent_accuracy, max_iterations);
        iterations_oxy = iterations;
        error_oxy = error;
    }

    // std::cout << "iterations oxygenation: " << iterations_oxy << "\n";
    // std::cout << "error oxygenation: " << error_oxy << "\n";

    if (po2_est > -1) {
        out.po2 = po2_est;
        out.so2 = so2 * 100.0;
    }
    
    // check for errors
    out.iterations = iterations_ab + iterations_oxy;
    out.error = error_ab || error_oxy;

    return out;
}

// declare function to calculate net charge of plasma
double net_charge_plasma(double hp_estimate) {
    ph = -std::log10(hp_estimate / 1000.0);
    double cco2p = tco2 / (1.0 + kc/hp_estimate + (kc*kd)/(hp_estimate*hp_estimate));
    hco3       = (kc * cco2p) / hp_estimate;
    double co3p = (kd * hco3) / hp_estimate;
    double ohp  = kw / hp_estimate;

    pco2 = cco2p / alpha_co2p;

    double a_base = albumin*(0.123*ph - 0.631) + phosphates*(0.309*ph - 0.469);

    return hp_estimate + sid - hco3 - 2.0*co3p - ohp - a_base - uma;
}

// calculate the O2 saturation from the po2 estimate
double calc_so2(double po2_estimate){
    double po2_n = std::pow(po2_estimate, n);
    double denom = po2_n + P50_n;
    return po2_n / denom;
}

// calculate the difference between the TO2 and the calculated TO2 from the po2 estimate
double dO2_content(double po2_estimate){
    // calculate the O2 saturation
    so2 = calc_so2(po2_estimate);
    // calculate the difference between the TO2 and the calculated TO2
    return hemoglobin * so2 + alpha_o2 * po2_estimate - to2;
}

// Brent’s method root‐finder
template<typename F>
double brent_root_fast(F f, double a, double b, double tolerance, int max_iterations) {
    double fa = f(a);
    double fb = f(b);
    
    if (std::abs(fb) < tolerance) return b;
    if (std::abs(fa) < tolerance) return a;
    if (fa * fb > 0.0) {
        error = true;
        return -1;
    }
    
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    
    double c = a, fc = fa, d = b - a, e = d;
    const double eps2 = 2.0 * std::numeric_limits<double>::epsilon();
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        iterations++;
        if (std::abs(fb) < tolerance) return b;
        
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b); std::swap(fa, fb);
            c = a; fc = fa; e = d = b - a;
        }
        
        const double tol = eps2 * std::abs(b) + tolerance;
        const double m = 0.5 * (c - b);
        if (std::abs(m) <= tol) return b;
        
        if (std::abs(e) < tol || std::abs(fa) <= std::abs(fb)) {
            d = e = m;
        } else {
            const double s = fb / fa;
            double p, q;
            
            if (a == c) {
                p = m * s;
                q = 0.5 - 0.5 * s;
            } else {
                const double qa = fa / fc, r = fb / fc;
                p = s * (m * qa * (qa - r) - (b - a) * (r - 1.0));
                q = (qa - 1.0) * (r - 1.0) * (s - 1.0);
            }
            
            if (p > 0.0) q = -q; else p = -p;
            
            const double old_e = e;
            e = d;
            
            if ((p + p < 3.0 * m * q - std::abs(tol * q)) && 
                (p + p < std::abs(old_e * q))) {
                d = p / q;
            } else {
                d = e = m;
            }
        }
        
        a = b; fa = fb;
        b += (std::abs(d) > tol) ? d : ((m > 0.0) ? tol : -tol);
        fb = f(b);
        
        if ((fb > 0.0) == (fc > 0.0)) {
            c = a; fc = fa; e = d = b - a;
        }
    }
    error = true;
    return b;
}

// set the emscripten bindings
EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::function("calc_blood_composition", &calc_blood_composition);
    emscripten::value_object<Input>("Input")
        .field("tco2", &Input::tco2)
        .field("to2", &Input::to2)
        .field("temp", &Input::temp)
        .field("hemoglobin", &Input::hemoglobin)
        .field("na", &Input::na)
        .field("k", &Input::k)
        .field("ca", &Input::ca)
        .field("mg", &Input::mg)
        .field("cl", &Input::cl)
        .field("lact", &Input::lact)
        .field("albumin", &Input::albumin)
        .field("phosphates", &Input::phosphates)
        .field("dpg", &Input::dpg)
        .field("uma", &Input::uma)
        .field("prev_ph", &Input::prev_ph)
        .field("prev_po2", &Input::prev_po2)
        .field("p50_0", &Input::p50_0);
    emscripten::value_object<Output>("Output")
        .field("ph", &Output::ph)
        .field("pco2", &Output::pco2)
        .field("hco3", &Output::hco3)
        .field("be", &Output::be)
        .field("po2", &Output::po2)
        .field("so2", &Output::so2)
        .field("error", &Output::error)
        .field("iterations", &Output::iterations);
};


// int main() {
//     // Example usage
//     Input input;
//     input.tco2 = 23.5;
//     input.to2 = 4.02;
//     input.temp = 37.0;
//     input.hemoglobin = 8.0;
//     input.na = 138.0;
//     input.k = 3.5;
//     input.ca = 1.0;
//     input.mg = 0.75;
//     input.cl = 108.0;
//     input.lact = 1.0;
//     input.albumin = 25.0;
//     input.phosphates = 1.64;
//     input.dpg = 5.0;
//     input.uma = 3.8;
//     input.prev_po2 = 18.7;
//     input.prev_ph = 7.37;
//     input.p50_0 = 18.8; // P50_0 = 18.8 for fetal hemoglobin, 26.7 for adult hemoglobin

//     Output result = calc_blood_composition(input);
//     std::cout << "pH: " << result.ph << "\n";
//     std::cout << "pco2: " << result.pco2 << "\n";
//     std::cout << "hco3: " << result.hco3 << "\n";
//     std::cout << "be: " << result.be << "\n";
//     std::cout << "po2: " << result.po2 << "\n";
//     std::cout << "so2: " << result.so2 << "\n";

// }
