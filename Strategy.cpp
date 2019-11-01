//
//  MMStrategy.cpp
//  HFT_FIX_C
//
//  Created by Chris Chen on 9/25/19.
//  Copyright © 2019 Chris Chen. All rights reserved.
//

#include "Strategy.hpp"
#include "Application.h"
#include "Ctools.hpp"

#include <time.h>
#include <cstdlib>
#include <future>
#include <boost/filesystem/operations.hpp>


Strategy::Strategy(Instrument& inst, double _unit, std::mutex &m){
   
    this->ins = inst;
    this->unit = _unit;
    this->strategy_name = ins.Ticker + "MM";
    this->ini_file = "Resources//" + ins.Ticker + "_MM_L2.ini";
    this->buffer = boost::circular_buffer<std::shared_ptr<std::array<double, WIDTH>>> (MAX_LEN);
    this->p_mutex = &m;
    
//    config_handle = std::async(std::launch::deferred, &Strategy::configcheck, this, ini_file);
    
    Strategy::load_config();
    
    time_t now = time(0);
    auto day = gmtime(&now);
    wday = (*day).tm_wday;
    
    std::cout<< "Seccessfually initialised - " << ins.Ticker << " " << utc_end << "\n";
}

void Strategy::register_app(std::shared_ptr<Application> app){
    this->app_ptr = app;
}

void Strategy::onTick(std::shared_ptr<std::array<double, WIDTH>> ptr, int size){
    
    buffer.push_back(ptr);
    buffer.linearize();
    
    if ((timecheck((*ptr)[0]) == 1) &&  permit_to_trade){
        if ((buffer.array_one().second >= MAX_LEN) && (tickcount%period == 1)){
            std::shared_ptr<std::array<double, 21>>* x = buffer.array_one().first;
            std::array<double, 50> feature = Ctools::MM2_alpha(x, static_cast<int>(buffer.array_one().second));
            
            std::array<double, 50> f = Ctools::transform_(feature, mean, scale);
            double pred = Ctools::predict_(f, coef, intercept);
            
            double t = ((*ptr)[0] - utc_start)/(utc_end - utc_start);
            
            std::vector<double> mid (period);
            for (int i = (buffer.array_one().second - period); i <(buffer.array_one().second); i++){
                mid[i - (buffer.array_one().second - period)] = ((*(x[i]))[1] +(*(x[i]))[11]) / 2;
            }
            double sigma = Ctools::vector_std_(mid);
            
            std::array<double, 3> hjb_arry = Ctools::HJB_ARM(((*ptr)[1] + (*ptr)[11])/2, gamma, kappa, eta, sigma, pred, t, pnl.hjb_q);
            
            int decimal = 5;
            
            if (ins.Ticker == "USDJPY"){
                decimal = 3;
            }
            
            double bid = Ctools::rounding(hjb_arry[0] - hjb_arry[1], decimal);
            double ask = Ctools::rounding(hjb_arry[0] + hjb_arry[2], decimal);
            
            std::cout<<"pred: " << pred << " bid: " << bid << " ask: "<< ask << "\n";
            
            if (bid > (*ptr)[1]){
                bid = (*ptr)[1];
            }
            
            if (ask < (*ptr)[11]){
                ask = (*ptr)[11];
            }
              
            if (pnl.hjb_q < max_q){
                if (bid == prev_bid){
                    if (!Ctools::checkexistence(openordertable, prev_bid_label)){
                        std::string label = create_label();
                        app_ptr->NEW_ORDER(ins, label, ordertype::LIMIT, side::BID, bid, unit);
                        prev_bid_label = label;
                    }
                }else{
                    if (!Ctools::checkexistence(openordertable, prev_bid_label)){
                        std::string label = create_label();
                        app_ptr->NEW_ORDER(ins, label, ordertype::LIMIT, side::BID, bid, unit);
                        prev_bid_label = label;
                    }else{
                        std::string label = create_label();
                        app_ptr->REPLACE_ORDER_REQUEST(ins, prev_bid_label, label, side::BID, ordertype::LIMIT, bid, unit);
                        prev_bid_label = label;
                    }
                }
            }else{
                if (!Ctools::checkexistence(openordertable, prev_bid_label)){
                    std::string label = create_label();
                    app_ptr->CANCEL_ORDER(ins, prev_bid_label, label, side::BID, unit);
                }
            }
            
            if (pnl.hjb_q > -max_q){
                if (ask == prev_ask){
                    if (!Ctools::checkexistence(openordertable, prev_ask_label)){
                        std::string label = create_label();
                        app_ptr->NEW_ORDER(ins, label, ordertype::LIMIT, side::ASK, ask, unit);
                        prev_ask_label = label;
                    }
                }else{
                    if (!Ctools::checkexistence(openordertable, prev_ask_label)){
                        std::string label = create_label();
                        app_ptr->NEW_ORDER(ins, label, ordertype::LIMIT, side::ASK, ask, unit);
                        prev_ask_label = label;
                    }else{
                        std::string label = create_label();
                        app_ptr->REPLACE_ORDER_REQUEST(ins, prev_ask_label, label, side::ASK, ordertype::LIMIT, ask, unit);
                        prev_ask_label = label;
                    }
                }
            }else{
                if (!Ctools::checkexistence(openordertable, prev_ask_label)){
                   std::string label = create_label();
                   app_ptr->CANCEL_ORDER(ins, prev_ask_label, label, side::BID, unit);
                }
            }
            prev_bid = bid;
            prev_ask = ask;
        }
            tickcount++;
        
    }else if (timecheck((*ptr)[0]) == 2){
        clear_pos();
    }else{
        permit_to_trade = true;
    }
    
    //PNL update//
    if (pnl.cur_px == 0){
        pnl.contract_num > 0 ? pnl.cur_px = (*ptr)[1] : pnl.cur_px = (*ptr)[11];
        pnl.offset = -(pnl.cur_px * pnl.contract_num * MULTIPLIER + pnl.cash);
    }else{
        pnl.contract_num > 0 ? pnl.cur_px = (*ptr)[1] : pnl.cur_px = (*ptr)[11];
    }
    
    if (ins.needinverse == '1'){
       pnl.pnl_ = (pnl.cur_px * pnl.contract_num * MULTIPLIER + pnl.cash + pnl.offset)/pnl.cur_px;
    }else{
      pnl.pnl_ = pnl.cur_px * pnl.contract_num * MULTIPLIER + pnl.cash + pnl.offset;
    }
    
    p_mutex->lock();
    std::cout << (ins.Ticker + " PNL: ") << pnl.pnl_ << " current quantity: " << pnl.contract_num <<"\n";
    p_mutex->unlock();
}

void Strategy::onExecutionReport(const FIX44::ExecutionReport &message){
    
    FIX::ExecType execType;
    FIX::OrdStatus ordstatus;
    
    switch (message.get(execType).getValue()) {
        case '0':{
            FIX::ClOrdID clOrdID;
            std::string orderid = message.get(clOrdID).getValue();
            openordertable.push_back(orderid);
            break;}
        case '4':{
            FIX::OrigClOrdID origClOrdID;
            std::string origorderid = message.get(origClOrdID).getValue();
            if (Ctools::checkexistence(openordertable, origorderid)){
                 Ctools::Removefromvect<std::string>(openordertable, origorderid);
            }
            break;}
        case '5':{
            FIX::ClOrdID clOrdID;
            FIX::OrigClOrdID origClOrdID;
            std::string orderid = message.get(clOrdID).getValue();
            std::string origorderid = message.get(origClOrdID).getValue();
            openordertable.push_back(orderid);
            if (Ctools::checkexistence(openordertable, origorderid)){
                 Ctools::Removefromvect<std::string>(openordertable, origorderid);
            }
            break;}
        case '8':{
            break;}
        case 'F':{
            FIX::ClOrdID clOrdID;
            FIX::Side side;
            FIX::LastQty lastQty;
            FIX::LeavesQty leavesQty;
            FIX::LastPx lastPx;
            std::string orderid = message.get(clOrdID).getValue();
            char side_ = message.get(side).getValue();
            double last_Q = message.get(lastQty).getValue();
            double leaves_Q = message.get(leavesQty).getValue();
            double px = message.get(lastPx).getValue();
            
            //update PNL and inventory Q//
            side_ == '1'? pnl.contract_num += last_Q : pnl.contract_num -=last_Q;
            side_ == '1'? pnl.cash -= (MULTIPLIER * last_Q*px) : pnl.cash += (MULTIPLIER *last_Q * px);
            pnl.hjb_q = pnl.contract_num/unit;
        
            //update oustanding ordertable - add check to make sure to take care of partial fills//
            if ((leaves_Q == 0) && Ctools::checkexistence(openordertable, orderid)){
                Ctools::Removefromvect<std::string>(openordertable, orderid);
            }else if (leaves_Q != 0){
                if (side_ == '1'){
                    app_ptr->CANCEL_ORDER(ins, orderid, create_label(), side::BID, leaves_Q);
                }else{
                    app_ptr->CANCEL_ORDER(ins, orderid, create_label(), side::ASK, leaves_Q);
                }
            }
        }
        default:
            break;
    }
}

std::string Strategy::create_label(){
    return ins.SecurityID + "-" + std::to_string(tagcounter++);
    }
    

void Strategy::clear_pos(){
    for (std::string openorder : openordertable){
       app_ptr->CANCEL_ORDER(ins, openorder, create_label(), side::BID, 1);
    };
        
    permit_to_trade = false;

    if ((wday == 5) && (pnl.contract_num !=0)){
        pnl.contract_num > 0 ? app_ptr->NEW_ORDER(ins, "CLOSE"+create_label(), ordertype::MARKET, side::ASK, 0, pnl.contract_num) :
        app_ptr->NEW_ORDER(ins, "CLOSE"+create_label(), ordertype::MARKET, side::BID, 0, -pnl.contract_num);
    }
}

void Strategy::load_config(){
    
     this->coef = Ctools::ini_parsing(ini_file, "regressors", "coef");
     this->intercept = *Ctools::ini_parsing(ini_file, "regressors", "intercept");
     this->mean = Ctools::ini_parsing(ini_file, "scale", "mean");
     this->scale = Ctools::ini_parsing(ini_file, "scale", "scale");
     this->gamma = *Ctools::ini_parsing(ini_file, "HJB", "gamma");
     this->kappa = *Ctools::ini_parsing(ini_file, "HJB", "kappa");
     this->eta = *Ctools::ini_parsing(ini_file, "HJB", "eta");
     this->utc_start = *Ctools::ini_parsing(ini_file, "Hour", "utc_start");
     this->utc_end = *Ctools::ini_parsing(ini_file, "Hour", "utc_end");
    
}

int Strategy::timecheck(double timestamp){
    
    if (timestamp < utc_start){
        return 0;
    }else if ((timestamp > utc_start)&&(timestamp<utc_end)){
        return 1;
    }else{
        return 2;
    }
}

void Strategy::configcheck(const std::string &path){
    boost::filesystem::path p (path);
    std::cout<<"path - " << p << "\n";
    std::time_t t = boost::filesystem::last_write_time(path);
    std::cout<< "initial read: " << t << "\n";
    
    while (true) {
        sleep(60);
        std::time_t t_1 = boost::filesystem::last_write_time(path);
        
        if (t_1 != t){
           std::cout<< "Config file changed!\n";
            Strategy::load_config();
            permit_to_trade = true;
            t = t_1;
        }
    }
}

void Strategy::partialordercheck(){
    while (true) {
        std::cout<<"checking partial\n";
        if (openordertable.size()>2){
            for (int i =0 ; i < openordertable.size(); i++){
                std::cout << openordertable[i] << std::endl;
                if ((openordertable[i] != prev_ask_label) && (openordertable[i] != prev_ask_label)){
                    app_ptr->CANCEL_ORDER(ins, openordertable[i], create_label(), side::BID, 0);
                }
            }
        }
        sleep(20);
    }
}
