#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>

using namespace std;

double expiration_time, risk_free_rate, volatility, initial_stock_price, strike_price, bar_price, no_of_bars;

int no_of_trials;

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

double max(double a, double b)
{
    return (a > b) ? a : b;
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

double mu(double t, double a, double b, double t1, double t2)
{
    return a + (t-t1) / (t2-t1)*(b - a);
}

double sigma(double t, double a, double b, double t1, double t2)
{
    return sqrt((t-t1)*(t2-t)/(t2-t1));
}

double pd(double a, double b, double t1, double t2)
{
    double result = 1;
    
    for (int i = 0; i < no_of_bars; i++)
    {
        result *= (1 - N((bar_price - mu(i, a, b, t1, t2)) / sigma(i, a, b, t1, t2)));
    }
    
    return result;
}


pair<double,double> discret_sim(){
    double dT;
    double dR;
    double dSD;
    
    dT = expiration_time/(double)no_of_bars;
    dR = (risk_free_rate - pow(volatility,2)/2)*dT;
    dSD = volatility*sqrt(dT);
    
    double call_price = 0;
    double put_price = 0;
    
    for (int counter = 0; counter<no_of_trials; counter++){
        double S1, S2, S3, S4;
        
        S1 = initial_stock_price;
        S2 = initial_stock_price;
        S3 = initial_stock_price;
        S4 = initial_stock_price;
        
        for (int counter_2 = 0; counter_2<no_of_bars; counter_2++){
            double uniform1 = get_uniform();
            double uniform2 = get_uniform();
            double a = sqrt(-2.0*log(uniform1)) * cos(6.283185307999998*uniform2);
            double b = sqrt(-2.0*log(uniform1)) * sin(6.283185307999998*uniform2);
            
            S1 = S1 * exp(dR + dSD*a);
            S2 = S2 * exp(dR - dSD*a);
            S3 = S3 * exp(dR + dSD*b);
            S4 = S4 * exp(dR - dSD*b);
        }
        
        
        double call_1 = max(S1 - strike_price, 0) * (pd(initial_stock_price, S1, 0, no_of_bars));
        call_price += call_1;
        
        double call_2 = max(S2 - strike_price, 0) * (pd(initial_stock_price, S2, 0, no_of_bars));
        call_price += call_2;
        
        double call_3 = max(S3 - strike_price, 0) * (pd(initial_stock_price, S3, 0, no_of_bars));
        call_price += call_3;
        
        double call_4 = max(S4 - strike_price, 0) * (pd(initial_stock_price, S4, 0, no_of_bars));
        call_price += call_4;
        
        double put_1 = max(strike_price - S1, 0) * (pd(initial_stock_price, S1, 0, no_of_bars));
        double put_2 = max(strike_price - S2, 0) * (pd(initial_stock_price, S2, 0, no_of_bars));
        double put_3 = max(strike_price - S3, 0) * (pd(initial_stock_price, S3, 0, no_of_bars));
        double put_4 = max(strike_price - S4, 0) * (pd(initial_stock_price, S4, 0, no_of_bars));
        put_price = put_price + put_1 + put_2 + put_3 + put_4;
        
    }
    
    double final_call = call_price/(double)no_of_trials/4.0*exp(-risk_free_rate*expiration_time);
    double final_put = put_price/(double)no_of_trials/4.0*exp(-risk_free_rate*expiration_time);
    
    return make_pair(final_call, final_put);
}

pair<double,double> continuous_sim_MCS(){
    double dT;
    double dR;
    double dSD;
    
    dT = expiration_time/no_of_bars;
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
        
        for (int counter_2 = 0; counter_2<no_of_bars; counter_2++){
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

int main(int argc, char* argv[])
{
    
    sscanf(argv[1], "%lf", &expiration_time);
    sscanf(argv[2], "%lf", &risk_free_rate);
    sscanf(argv[3], "%lf", &volatility);
    sscanf(argv[4], "%lf", &initial_stock_price);
    sscanf(argv[5], "%lf", &strike_price);
    sscanf(argv[6], "%d", &no_of_trials);
    sscanf(argv[7], "%lf", &no_of_bars);
    sscanf(argv[8], "%lf", &bar_price);
    
    cout << "European Down-and-Out Discret Barrier Option Pricing via Monte Carlo Simulaiton" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interst Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Barrier Price = " << bar_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Barriers = " << no_of_bars << endl;
    
    pair<double,double> output_MCS = continuous_sim_MCS();
    pair<double,double> output_Dis = discret_sim();
    
    cout<<"---------------------------------------------------------"<<endl;
    cout<<"Average Call Price (Explicit MCS):     "<<get<0>(output_MCS)<<endl;
    cout<<"Average Call Price (Discret): "<<get<0>(output_Dis)<<endl;
    cout<<"---------------------------------------------------------"<<endl;
    cout<<"Average Put Price (Explicit MCS):      "<<get<1>(output_MCS)<<endl;
    cout<<"Average Put Price (Discret):  "<<get<1>(output_Dis)<<endl;
    
}
