#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>
#include <numeric>

void monte_carlo_worker(bool is_call, double S0, double K, double r, double sigma, double T, int num_sims, double& result){
    std::mt19937_64 rng(std::random_device{}());
    std::normal_distribution<> norm(0.0, 1.0);
    double local_sum = 0.0;
    for (int i = 0; i < num_sims; ++i){
        double Z = norm(rng);
        double ST = S0 * std::exp((r-0.5*sigma*sigma)*T + sigma * std::sqrt(T)*Z);
        double payoff = is_call ? std::max(ST-K,0.0) : std::max(K-ST,0.0);
        local_sum += payoff;
    }
    result = local_sum;
}

double monte_carlo_option_std_thread(bool is_call, double S0, double K, double r, double sigma, double T, int num_sims, int num_threads){
    std::vector<std::thread> threads;
    std::vector<double> results(num_threads, 0.0);
    int sims_per_thread = num_sims/num_threads;
    for (int i = 0; i < num_threads; ++i){
        threads.emplace_back(monte_carlo_worker, is_call, S0, K, r, sigma, T, sims_per_thread, std::ref(results[i]));
    }
    for (auto& t: threads){
        t.join();
    }
    double total_payoff = std::accumulate(results.begin(), results.end(), 0.0);
    return std::exp(-r*T)*(total_payoff/num_sims);
}

int main(){
    double S0 = 100.0;
    double K = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;
    int num_sims = 1'000'000;
    int num_threads = std::thread::hardware_concurrency(); 

    double call_price = monte_carlo_option_std_thread(true, S0, K, r, sigma, T, num_sims, num_threads);
    double put_price  = monte_carlo_option_std_thread(false, S0, K, r, sigma, T, num_sims, num_threads);

    std::cout << "European Call Price (std::thread): " << call_price << "\n";
    std::cout << "European Put Price  (std::thread): " << put_price  << "\n";
    return 0;
}