#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>

double monte_carlo_asian_put(double S0, double K, double r, double sigma, double T, int num_paths, int num_steps){
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::normal_distribution<> norm(0.0, 1.0);

    double dt = T / num_steps;
    double sum_payoffs = 0.0;
    for (int i = 0; i < num_paths; ++i){
        double S = S0;
        double path_sum = 0.0;
        for (int j = 0; j < num_steps; ++j){
            double Z = norm(rng);
            S *= std::exp((r-0.5*sigma*sigma) * dt + sigma * std::sqrt(dt) * Z);
            path_sum +=S;
        }
        double average_price = path_sum/num_steps;
        double payoff = std::max(K - average_price, 0.0);
        sum_payoffs += payoff;  
    }
    double discounted_mean = std::exp(-r * T) * (sum_payoffs/num_paths);
    return discounted_mean;
}

int main(){
    double S0 = 100.0;
    double K = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;
    int num_paths = 100000;
    int num_steps = 50;
    auto start = std::chrono::high_resolution_clock::now();
    double price = monte_carlo_asian_put(S0,K,r,sigma, T, num_paths, num_steps);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Monte Carlo Asian Put Price: " << price << std::endl;
    std::cout << "Simulation Time: " << duration.count() << "seconds" << std::endl;
    return 0; 
}