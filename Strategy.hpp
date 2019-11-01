//
//  MMStrategy.hpp
//  HFT_FIX_C
//
//  Created by Chris Chen on 9/25/19.
//  Copyright © 2019 Chris Chen. All rights reserved.
//

#ifndef MMStrategy_hpp
#define MMStrategy_hpp

#include <stdio.h>
#include <iostream>
#include <boost/circular_buffer.hpp>
#include <algorithm>
#include <future>


#include "Tick_Listener.h"
#include "Instrument.h"
#include "ExecutionReport_Listener.h"
#include "PNL_struct.h"


class Application;


#define MAX_LEN 300
#define WIDTH 21
#define MULTIPLIER 10000

class Strategy:
public Tick_Listener,
public ExecutionReport_Listener{
    
public:
    
    std::shared_ptr<Application> app_ptr; //the app;
    Instrument ins; // Instrument to be traded on;
    PNL pnl; // object to keep track of PNL of the strategy;
    double unit;
    
    //files and variables for admin uses.
    std::string strategy_name;
    std::string ini_file;
    std::mutex *p_mutex;
    bool permit_to_trade = true;
    std::vector<std::string> openordertable;
    boost::circular_buffer<std::shared_ptr<std::array<double, WIDTH>>> buffer;
   
    std::future<void> config_handle;
    std::future<void> partial_handle;
    
private:
    //trading variables
    double* coef;
    double intercept;
    double* mean;
    double* scale;

    
    long utc_start;
    long utc_end;
    double totalPNL = 0;

    double gamma;
    double eta;
    double kappa;
    int period = 5;
    double max_q = 20;
    
    //admin variables
    long tickcount = 1;
    int tagcounter = 1;
    int wday = 0;
    std::string prev_bid_label = "STARTBID";
    std::string prev_ask_label = "STARTASK";
    double prev_bid = 0;
    double prev_ask = 0;
    
    
public:
    Strategy(Instrument& inst, double _unit, std::mutex &m);
    
    void register_app(std::shared_ptr<Application>);
    
    void onTick(std::shared_ptr<std::array<double, WIDTH>>, int size);
   
    void onExecutionReport(const FIX44::ExecutionReport& message);
    
    std::string create_label();
    
    void load_config();
    
    int timecheck(double timestamp);
    
    void configcheck(const std::string &path);
    
    void partialordercheck();
    
private:
    void clear_pos();
};





#endif /* Strategy_hpp */
