#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>

double monte_carlo_digital_option(bool is_call, double S0, double K, double r, double sigma, double T, int num_simulations){
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::normal_distribution<> norm(0.0, 1.0);
    double sum_payoffs = 0.0;
    for (int i = 0; i < num_simulations; i++){
        double Z = norm(rng);
        double ST = S0 * std::exp((r-0.5*sigma*sigma)*T + sigma * std::sqrt(T)*Z);
        double payoff = (is_call ? (ST > K) : (ST < K)) ? 1.0:0.0;
        sum_payoffs += payoff;
    }
    double discounted_mean = std::exp(-r*T)*(sum_payoffs/num_simulations);
    return discounted_mean;
}

int main(){
    double S0 = 100.0;
    double K = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;
    int num_simulations = 1000000;
    auto start = std::chrono::high_resolution_clock::now();
    double put_price = monte_carlo_digital_option(false, S0, K, r, sigma, T, num_simulations);
    double call_price = monte_carlo_digital_option(true, S0, K, r, sigma, T, num_simulations);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Digital Put Price: " << put_price << std::endl;
    std::cout << "Digital Call Price: " << call_price << std::endl;
    std::cout << "Simulation Time: " << duration.count() << " seconds" << std::endl;
    return 0;
}