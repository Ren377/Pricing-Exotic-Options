//
//  Deng_Renren.cpp
//  IE523_Final_Assignment
//
//  Created by Renren on 12/7/17.
//  Copyright © 2017 Renren Deng. All rights reserved.
//  Reference by IE523 Class

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <vector>
using namespace std;

double  risk_free_rate, strike_price, barrier_price, R, SD;
double initial_stock_price, expiration_time, volatility;
int no_of_trials,no_of_sample_points;

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
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
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};
void discrete_simulation()
{
    R = (risk_free_rate - 0.5*pow(volatility,2))*(expiration_time/(double)no_of_sample_points);
    SD = volatility*sqrt(expiration_time/(double)no_of_sample_points);
    
    double mean_i_1, mean_i_2, mean_i_3, mean_i_4, variance_i_1;
    double touch1, touch2, touch3, touch4;
    double st1, st2, st3, st4;
    vector< double > prob;
    
    
    double put_option_price_1 = 0.0;
    double call_option_price_1 = 0.0;
    double put_option_price_2 = 0.0;
    double call_option_price_2 = 0.0;
    for (int i = 0; i < no_of_trials; i++) {
        touch1 = 1,
        touch2 = 1,
        touch3 = 1,
        touch4 = 1;
        st1 = initial_stock_price;
        st2 = initial_stock_price;
        st3 = initial_stock_price;
        st4 = initial_stock_price;
        
        for (int j = 0; j < 4; j++)
        {
            prob.push_back(1.0);
        }
        
        for (int k = 0; k < no_of_sample_points; k++)
        {
            // generate unit-normals using Box-Muller Transform
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            st1 = st1 * exp(R + SD*a);
            st2 = st2 * exp(R - SD*a);
            st3 = st3 * exp(R + SD*b);
            st4 = st4 * exp(R - SD*b);
            
            if(st1 <= barrier_price || st2 <= barrier_price || st3 <= barrier_price || st4 <= barrier_price)
            {
                if (st1 <= barrier_price)
                    touch1 = 0;
                if (st2 <= barrier_price)
                    touch2 = 0;
                if (st3 <= barrier_price)
                    touch3 = 0;
                if (st4 <= barrier_price)
                    touch4 = 0;
            }
        }
        
        for (int j = 1; j <= no_of_sample_points; j++)
        {
            mean_i_1 = initial_stock_price +( ((double) j)/((double) no_of_sample_points)*(st1 - initial_stock_price));
            mean_i_2 = initial_stock_price +( ((double) j)/((double) no_of_sample_points)*(st2 - initial_stock_price));
            mean_i_3 = initial_stock_price +( ((double) j)/((double) no_of_sample_points)*(st3 - initial_stock_price));
            mean_i_4 = initial_stock_price +( ((double) j)/((double) no_of_sample_points)*(st4 - initial_stock_price));
            
            variance_i_1 = (j*(expiration_time/(double)no_of_sample_points)* (expiration_time - j*(expiration_time/(double)no_of_sample_points)))/(expiration_time);
            if (j == 1)
            {
                prob[0] = (1.0 - N((barrier_price - mean_i_1)/sqrt(variance_i_1)));
                prob[1] = (1.0 - N((barrier_price - mean_i_2)/sqrt(variance_i_1)));
                prob[2] = (1.0 - N((barrier_price - mean_i_3)/sqrt(variance_i_1)));
                prob[3] = (1.0 - N((barrier_price - mean_i_4)/sqrt(variance_i_1)));
            }
            else
            {
                prob[0] = prob[0] * (1 - N((barrier_price - mean_i_1)/sqrt(variance_i_1)));
                prob[1] = prob[1] * (1 - N((barrier_price - mean_i_2)/sqrt(variance_i_1)));
                prob[2] = prob[2] * (1 - N((barrier_price - mean_i_3)/sqrt(variance_i_1)));
                prob[3] = prob[3] * (1 - N((barrier_price - mean_i_4)/sqrt(variance_i_1)));
                
            }
        }
        
        
        if (st1 <= barrier_price)
            prob[0] = 1.0;
        if (st2 <= barrier_price)
            prob[1] = 1.0;
        if (st3 <= barrier_price)
            prob[2] = 1.0;
        if (st4 <= barrier_price)
            prob[3] = 1.0;
        
        
        
        call_option_price_1 += (max(0.0, st1 - strike_price)*touch1 +
                                max(0.0, st2 - strike_price)*touch2 +
                                max(0.0, st3 - strike_price)*touch3 +
                                max(0.0, st4 - strike_price)*touch4)/4.0;
        put_option_price_1 += (max(0.0, strike_price - st1)*touch1 +
                               max(0.0, strike_price - st2)*touch2 +
                               max(0.0, strike_price - st3)*touch3 +
                               max(0.0, strike_price - st4)*touch4)/4.0;
        
        call_option_price_2 += (max(0.0, st1 - strike_price)* prob[0] +
                                max(0.0, st2 - strike_price)* prob[1] +
                                max(0.0, st3 - strike_price)* prob[2] +
                                max(0.0, st4 - strike_price)* prob[3])/4.0;
        put_option_price_2 += (max(0.0, strike_price - st1)* prob[0] +
                               max(0.0, strike_price - st2)* prob[1] +
                               max(0.0, strike_price - st3)* prob[2] +
                               max(0.0, strike_price - st4)* prob[3])/4.0;
    }
    call_option_price_1 = exp(-risk_free_rate*expiration_time)*(call_option_price_1/((double) no_of_trials));
    put_option_price_1 = exp(-risk_free_rate*expiration_time)*(put_option_price_1/((double) no_of_trials));
    call_option_price_2 = exp(-risk_free_rate*expiration_time)*(call_option_price_2/((double) no_of_trials));
    put_option_price_2 = exp(-risk_free_rate*expiration_time)*(put_option_price_2/((double) no_of_trials));
    
    cout << "The average Call Price via explicit simulation of price paths = " << call_option_price_1<< endl;
    cout << "The average Call Price with Brownian-Bridge correction on the final price = " << call_option_price_2 << endl;
    cout << "The average Put Price via explicit simulation of price paths = " << put_option_price_1 << endl;
    cout << "The average Put Price with Brownian-Bridge correction on the final price = " << put_option_price_2 << endl;
}


int main (int argc, char* argv[])
{
    
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &no_of_trials);
    sscanf (argv[7], "%d", &no_of_sample_points);
    sscanf (argv[8], "%lf", &barrier_price);
    
    cout << "--------------------------------" << endl;
    cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "barrier price = " << barrier_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of discrete barriers = " << no_of_sample_points << endl;
    cout << "--------------------------------" << endl;
    
    
    discrete_simulation();
    
}


