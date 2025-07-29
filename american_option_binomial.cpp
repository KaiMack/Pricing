#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>

double monte_carlo_put_price(double S0, double K, double r, double sigma, double T, int num_simulations){
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::normal_distribution<> norm(0.0, 1.0);

    double sum_payoffs = 0.0;
    for (int i = 0; i < num_simulations; i++){
        double Z = norm(rng);
        double ST = S0*std::exp((r-0.5*sigma*sigma) * T + sigma * std::sqrt(T)*Z);
        double payoff = std::max(K-ST,0.0);
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
    double price = monte_carlo_put_price(S0,K,r,sigma, T, num_simulations);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Monte Carlo Estimated Put Price: " << price << std::endl;
    std::cout << "Simulation Time: " << duration.count() << "seconds" << std::endl;
    return 0; 
}