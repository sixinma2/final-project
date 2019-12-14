//
//  continuous_bar(pt1).cpp
//  Final
//
//  Created by SIXIN MA.
//

#include <iostream>
#include <cmath>

using namespace std;

double risk_free_rate;
double strike_price;
double initial_stock_price;
double expiration_time;
double volatility;
double bar_price;
double no_of_divs;

int no_of_trials;

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

double max(double a, double b) {
    return (b < a )? a:b;
}

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z)*(z / 2.0));
    double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
};

double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1 - (sigma*time_sqrt);
    return K * exp(-r * time)*N(-d2) - S * N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time)      // time to maturity
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1 - (sigma*time_sqrt);
    return S * N(d1) - K * exp(-r * time)*N(d2);
};

double closed_form_down_and_out_european_call_option()
{
    // I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
    double K = (2 * risk_free_rate) / (volatility*volatility);
    double A = option_price_call_black_scholes(initial_stock_price, strike_price,
                                               risk_free_rate, volatility, expiration_time);
    double B = (bar_price*bar_price) / initial_stock_price;
    double C = pow(initial_stock_price / bar_price, -(K - 1));
    double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
    return (A - D * C);
}

double closed_form_down_and_in_european_put_option()
{
    // just making it easier by renaming the global variables locally
    double S = initial_stock_price;
    double r = risk_free_rate;
    double T = expiration_time;
    double sigma = volatility;
    double H = bar_price;
    double X = strike_price;
    
    // Took these formulae from some online reference
    double lambda = (r + ((sigma*sigma) / 2)) / (sigma*sigma);
    double temp = 2 * lambda - 2.0;
    double x1 = (log(S / H) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    double y = (log(H*H / (S*X)) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    double y1 = (log(H / S) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    return (-S * N(-x1) + X * exp(-r * T)*N(-x1 + sigma * sqrt(T)) +
            S * pow(H / S, 2 * lambda)*(N(y) - N(y1)) -
            X * exp(-r * T)*pow(H / S, temp)*(N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

double closed_form_down_and_out_european_put_option()
{
    double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
                                                        risk_free_rate, volatility, expiration_time);
    double put_down_in = closed_form_down_and_in_european_put_option();
    return (vanilla_put - put_down_in);
}

pair<double,double> continuous_sim_MCS(){
    double dT;
    double dR;
    double dSD;
    
    dT = expiration_time/no_of_divs;
    dR = (risk_free_rate - pow(volatility,2)/2)*dT;
    dSD = volatility*sqrt(dT);
    
    double call_price = 0;
    double put_price = 0;
    
    for (int counter = 0; counter<no_of_trials; counter++){
        double S1[2], S2[2], S3[2], S4[2];
        
        S1[1] = initial_stock_price;
        S2[1] = initial_stock_price;
        S3[1] = initial_stock_price;
        S4[1] = initial_stock_price;
        
        S1[0] = 1;          //toggle triggering of barrier, 1 as not triggered
        S2[0] = 1;
        S3[0] = 1;
        S4[0] = 1;
        
        for (int counter_2 = 0; counter_2<no_of_divs; counter_2++){
            double uniform1 = get_uniform();
            double uniform2 = get_uniform();
            double a = sqrt(-2.0*log(uniform1)) * cos(6.283185307999998*uniform2);
            double b = sqrt(-2.0*log(uniform1)) * sin(6.283185307999998*uniform2);
            
            S1[1] = S1[1] * exp(dR + dSD*a);
            S2[1] = S2[1] * exp(dR - dSD*a);
            S3[1] = S3[1] * exp(dR + dSD*b);
            S4[1] = S4[1] * exp(dR - dSD*b);
            
            if (S1[1]<=bar_price)
                S1[0] = 0;
            if (S2[1]<=bar_price)
                S2[0] = 0;
            if (S3[1]<=bar_price)
                S3[0] = 0;
            if (S4[1]<=bar_price)
                S4[0] = 0;
        }
        
        double call_payoff_1 = max(0, S1[1] - strike_price);
        call_price += call_payoff_1 * S1[0];
        double call_payoff_2 = max(0, S2[1] - strike_price);
        call_price += call_payoff_2 * S2[0];
        double call_payoff_3 = max(0, S3[1] - strike_price);
        call_price += call_payoff_3 * S3[0];
        double call_payoff_4 = max(0, S4[1] - strike_price);
        call_price += call_payoff_4 * S4[0];
        
        
        double put_payoff_1 = max(0, strike_price - S1[1]);
        put_price += put_payoff_1*S1[0];
        double put_payoff_2 = max(0, strike_price - S2[1]);
        put_price += put_payoff_2*S2[0];
        double put_payoff_3 = max(0, strike_price - S3[1]);
        put_price += put_payoff_3*S3[0];
        double put_payoff_4 = max(0, strike_price - S4[1]);
        put_price += put_payoff_4*S4[0];
    }
    double final_call_price = call_price/4/(double)no_of_trials;
    double final_put_price = put_price/4/(double)no_of_trials;
    
    final_call_price *= exp(-risk_free_rate*expiration_time);
    final_put_price *= exp(-risk_free_rate*expiration_time);
    
    return make_pair(final_call_price, final_put_price);
};

double pc(double cur_price, double ini_price, double bar)       //pc formula
{
    if ((ini_price <= bar) || (cur_price <= bar))
        return 1;
    else
        return exp(-(2*log(ini_price/bar)*log(cur_price/bar))/(pow(volatility,2)*expiration_time));
};

pair<double, double> continuous_sim_pc(){
    double dT;
    double dR;
    double dSD;
    
    dT = expiration_time/no_of_divs;
    dR = (risk_free_rate - pow(volatility,2)/2)*dT;
    dSD = volatility*sqrt(dT);
    
    double call_price = 0;
    double put_price = 0;
    
    for (int counter = 0; counter<no_of_trials; counter++){
        double S1, S2, S3, S4;
        
        S1 = initial_stock_price;
        S2 = S1;
        S3 = S2;
        S4 = S3;
        
        for (int counter_2 = 0; counter_2<no_of_divs; counter_2++){
            double uniform1 = get_uniform();
            double uniform2 = get_uniform();
            double a = sqrt(-2.0*log(uniform1)) * cos(6.283185307999998*uniform2);
            double b = sqrt(-2.0*log(uniform1)) * sin(6.283185307999998*uniform2);
            
            S1 = S1 * exp(dR + dSD*a);
            S2 = S2 * exp(dR - dSD*a);
            S3 = S3 * exp(dR + dSD*b);
            S4 = S4 * exp(dR - dSD*b);
        }
        
        double call_payoff_1 = max(0, S1 - strike_price);
        call_price += call_payoff_1*(1-pc(S1, initial_stock_price, bar_price));
        double call_payoff_2 = max(0, S2 - strike_price);
        call_price += call_payoff_2*(1-pc(S2, initial_stock_price, bar_price));
        double call_payoff_3 = max(0, S3 - strike_price);
        call_price += call_payoff_3*(1-pc(S3, initial_stock_price, bar_price));
        double call_payoff_4 = max(0, S4 - strike_price);
        call_price += call_payoff_4*(1-pc(S4, initial_stock_price, bar_price));
        
        
        double put_payoff_1 = max(0, strike_price - S1);
        put_price += put_payoff_1*(1-pc(S1, initial_stock_price, bar_price));
        double put_payoff_2 = max(0, strike_price - S2);
        put_price += put_payoff_2*(1-pc(S2, initial_stock_price, bar_price));
        double put_payoff_3 = max(0, strike_price - S3);
        put_price += put_payoff_3*(1-pc(S3, initial_stock_price, bar_price));
        double put_payoff_4 = max(0, strike_price - S4);
        put_price += put_payoff_4*(1-pc(S4, initial_stock_price, bar_price));
    }
    double final_call_price = call_price/4/(double)no_of_trials;
    double final_put_price = put_price/4/(double)no_of_trials;
    
    final_call_price *= exp(-risk_free_rate*expiration_time);
    final_put_price *= exp(-risk_free_rate*expiration_time);
    
    return make_pair(final_call_price, final_put_price);
};


int main(int argc, char* argv[])
{
    
    sscanf(argv[1], "%lf", &expiration_time);
    sscanf(argv[2], "%lf", &risk_free_rate);
    sscanf(argv[3], "%lf", &volatility);
    sscanf(argv[4], "%lf", &initial_stock_price);
    sscanf(argv[5], "%lf", &strike_price);
    sscanf(argv[6], "%d", &no_of_trials);
    sscanf(argv[7], "%lf", &no_of_divs);
    sscanf(argv[8], "%lf", &bar_price);
    
    cout << "European Down-and-Out Continuous Barrier Option Pricing via Monte Carlo Simulaiton" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interst Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Barrier Price = " << bar_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Points = " << no_of_divs << endl;
    
    pair<double,double> output_MCS = continuous_sim_MCS();
    pair<double,double> output_pc = continuous_sim_pc();
    double theo_call = closed_form_down_and_out_european_call_option();
    double theo_put = closed_form_down_and_out_european_put_option();
    
    cout<<"---------------------------------------------------------"<<endl;
    cout<<"Theoretical Call Price:                "<<theo_call<<endl;
    cout<<"Average Call Price (Explicit MCS):     "<<get<0>(output_MCS)<<endl;
    cout<<"Average Call Price ((1-p)-adjustment): "<<get<0>(output_pc)<<endl;
    cout<<"---------------------------------------------------------"<<endl;
    cout<<"Theoretical Put Price:                 "<<theo_put<<endl;
    cout<<"Average Put Price (Explicit MCS):      "<<get<1>(output_MCS)<<endl;
    cout<<"Average Put Price ((1-p)-adjustment):  "<<get<1>(output_pc)<<endl;
    
}
