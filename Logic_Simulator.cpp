/*
 * Logic_Simulator.cpp
 *
 *  Created on: Dec 1, 2019
 *      Author: chris
 *
 * Introduction:
 * This file is part of the implementation of the backtest system.
 * The C++ system consists mainly 4 parts:
 *
 * 1. The server object which reads recorded market data from database and feeds into the trading logic.
 *    The server side is also responsible for order processing for the orders from the trading logic.
 *
 * 2. Logic Simulator: receives incoming market data feed from the server object, generates signals and trade actions, and conduct other administrative order management functions.
 *
 * 3. Ctools library: a collection of alpha signals
 *
 * 4. Machine learning object: an object that supports real time learning and update model parameters.
 *
 * This file implements part 2 of the overal system.
 */

#include "Logic_Simulator.hpp"
#include "App_Sim.hpp"
#include "ThostFtdcUserApiStruct.h"
#include "Ctools.hpp"


/*
 Constuctor function
 */
Signal_Generator::Signal_Generator(std::string mode_ , std::string model_name_, std::string tradingday_, double threshold_ ,int corr_lookback_, std::string trade_, std::string pair_1_, std::set<std::string> ins_vect_){
	this->corr_lookback = corr_lookback_;
	this->trade_ins = trade_;
	this->pair_1 = pair_1_;
	if (mode_ == "B") {
		this->ini_file = "Resources/model_files/" + trade_ + "/" + trade_ + "_CC.ini";
        std::cout << ini_file << "\n";
	}
    
	for (std::string ins : ins_vect_) {
		mkt_set.insert(ins);
		size_map[ins] = false;
		if ((ins != pair_1) && (ins != trade_ins)) {
			corr_map[ins] = -10;
		}
	}
        
    load_config();
	register_func_map();
    
	///update header for X file out put///
	if (mode_ == "M") {
		header.push_back(header_[0]);
		if (trade_ins != pair_1) {
			for (int i = 0; i < pre_fix.size(); i++) {
				for (int k = 1; k < header_.size(); k++) {
					std::string temp = pre_fix[i] + "_" + header_[k];
					header.push_back(temp);
				}
			}
			if (ins_vect_.size() > 2) {
				for (int i = 0; i < header_p.size(); i++) {
					std::string temp = "3_" + header_p[i];
					header.push_back(temp);
				}
			}
			if (ins_vect_.size() > 2) {
				for (int i = 0; i < pre_fix_r.size(); i++) {
					for (int k = 0; k < header_r.size(); k++) {
						std::string temp = pre_fix_r[i] + "_" + header_r[k];
						header.push_back(temp);
					}
				}
			}
			else {
				for (int k = 0; k < header_r.size(); k++) {
					std::string temp = "4_" + header_r[k];
					header.push_back(temp);
				}
			}
		}
		else {
			for (int i = 0; i < pre_fix.size(); i++) {
				if (pre_fix[i] != "2") {
					for (int k = 1; k < header_.size(); k++) {
						std::string temp = pre_fix[i] + "_" + header_[k];
						header.push_back(temp);
					}
				}
			}
			if (ins_vect_.size() > 2) {
				for (int i = 0; i < header_p.size(); i++) {
					std::string temp = "3_" + header_p[i];
					header.push_back(temp);
				}
			}
			if (ins_vect_.size() > 2) {
				for (int i = 0; i < pre_fix_r.size(); i++) {
					if (pre_fix_r[i] != "4") {
						for (int k = 1; k < header_r.size(); k++) {
							std::string temp = pre_fix_r[i] + "_" + header_r[k];
							header.push_back(temp);
						}
					}
				}
			}
		}
	}
	else {
		header.push_back(header_[0]);
		for (int i = 0; i < main_feature.size(); i++) {
				std::string temp = "1_" + main_feature[i];
				header.push_back(temp);
		}
		if (p1_feature[0].length()!=0) {
			for (int i = 0; i < p1_feature.size(); i++) {
				std::string temp = "2_" + p1_feature[i];
				header.push_back(temp);
			}
		}
		if (p2_feature[0].length()!=0) {
			for (int i = 0; i < p2_feature.size(); i++) {
				std::string temp = "3_" + p2_feature[i];
				header.push_back(temp);
			}
		}

		if (r1_feature[0].length()!=0) {
			for (int i = 0; i < r1_feature.size(); i++) {
				std::string temp = "4_" + r1_feature[i];
				header.push_back(temp);
			}
		}

		if (r2_feature[0].length()!=0) {
			for (int i = 0; i < r2_feature.size(); i++) {
				std::string temp = "5_" + r2_feature[i];
				header.push_back(temp);
			}
		}
	}
    
    std::cout << "ini_file threshold: " << threshold << "\n";
    if (mode_ != "B"){
        this->model = ML(mode_, model_name_, threshold_, corr_lookback_, static_cast<int>(header.size()-1));
    }else{
        this->model = ML(mode_, model_name_, threshold, corr_lookback_, static_cast<int>(header.size()-1));
        std::string feature_path = "Resources/model_files/" + trade_ + "/" + trade_ + "_features.csv";
        std::string model_path = "Resources/model_files/" + trade_ + "/" + trade_ + "_" + model_name_ + ".xml";
        model.load_main_feature(feature_path);
        mlpack::data::Load(model_path, model_name_, model.regressor);
        std::cout << "header size: " << header.size() << "\n";
        //std::cout<< " intercept: " << model.regressor.Intercept() << "\n";
    }
}


/*
 Function to read the recorded exchange timestamp and convert into C tm object. Function is used to determine the appropriate time to flat all the outstanding positions.
 */
tm Signal_Generator::make_time(std::string t_str){
	tm tm_;
	int year, month, day, hour, minute, second, micro;
	sscanf(t_str.c_str(),"%d-%d-%d %d:%d:%d.%d", &year, &month, &day, &hour, &minute, &second, &micro);
	tm_.tm_year  = year-1900;
	tm_.tm_mon   = month-1;
	tm_.tm_mday  = day;
	tm_.tm_hour  = hour;
	tm_.tm_min   = minute;
	tm_.tm_sec   = second;
	tm_.tm_isdst = 0;

	return tm_;
}

/*
 Function that processes incoming market data feed when the system is in the "Modeling" mode.
 It generates all the alpha signals that are available in the Ctool library.
 The full signals would be outputed as a csv file.
 This is to facilitate feature selection process which is done in Python.
 */
void Signal_Generator::onModeling(std::string time_, DATA& x, int size) {
	tm t_ = make_time(time_);
	if (mkt_set.find(x.buffer[size - 1].InstrumentID) != mkt_set.end()) {
		if ((size >= MAX_LEN) && (bt->data_map[pair_1].mid.size() >= MAX_LEN)) {
			if ((x.buffer[size - 1].InstrumentID != trade_ins) && (x.buffer[size - 1].InstrumentID != pair_1)) {
				//std::cout << "pass " << x.buffer[size - 1].InstrumentID << "size: " << x.ret.size() << "\n";
				double cor = Ctools::compute_cos(bt->data_map[pair_1], x, corr_lookback);
				if (cor != 0) {
					corr_map[x.buffer[size - 1].InstrumentID] = cor;
				}
			}
			if (!permit_to_trade) {
				size_map[x.buffer[size - 1].InstrumentID] = true;
				permit_check();
			}
		}
	}
	if (time_check(t_) == 1) {
		if ((permit_to_trade) && (bt->data_map[trade_ins].mid.size() >= MAX_LEN) && (x.buffer[size - 1].InstrumentID == trade_ins)){
			std::vector<double> features;
			if (trade_ins == pair_1) {
				std::string p_max_key = Ctools::max_(corr_map);
				std::vector<double> f_1 = alpha(x, size);				
				if (mkt_set.size() > 1) {
					std::vector<double> f_2 = pair_alpha(bt->data_map[p_max_key], size);
					std::vector<double> f_3 = r_alpha(x, bt->data_map[p_max_key]);
					features.insert(features.end(), f_1.begin(), f_1.end());
					features.insert(features.end(), f_2.begin(), f_2.end());
					features.insert(features.end(), f_3.begin(), f_3.end());
				}
				else {
					features.insert(features.end(), f_1.begin(), f_1.end());
				}
			}
			else {
				std::string p_max_key = Ctools::max_(corr_map);
				std::vector<double> f_1 = alpha(x, size);
				std::vector<double> f_2 = alpha(bt->data_map[pair_1], size);
				if (mkt_set.size() > 2) {
					std::vector<double> f_3 = pair_alpha(bt->data_map[p_max_key], size);
					std::vector<double> f_4 = r_alpha(x, bt->data_map[pair_1]);
					std::vector<double> f_5 = r_alpha(x, bt->data_map[p_max_key]);

					features.insert(features.end(), f_1.begin(), f_1.end());
					features.insert(features.end(), f_2.begin(), f_2.end());
					features.insert(features.end(), f_3.begin(), f_3.end());
					features.insert(features.end(), f_4.begin(), f_4.end());
					features.insert(features.end(), f_5.begin(), f_5.end());
				}
				else {
					std::vector<double> f_4 = r_alpha(x, bt->data_map[pair_1]);
					features.insert(features.end(), f_1.begin(), f_1.end());
					features.insert(features.end(), f_2.begin(), f_2.end());
					features.insert(features.end(), f_4.begin(), f_4.end());
				}
			
			}
		
            std::pair<std::string, std::vector<double>> new_pp = std::make_pair(time_, features);
            model.processing_features(new_pp, x.mid[size-1]);
        }
	}
	else if (time_check(t_) == 0) {
		for (std::string ins : mkt_set) {
			size_map[ins] = false;
			if (Ctools::checkkeyexistence(ins, corr_map)) {
				corr_map[ins] = -10;
			}
			
		}
		permit_to_trade = false;
	}
	else {
		for (std::string ins : mkt_set) {
			size_map[ins] = false;
			if (Ctools::checkkeyexistence(ins, corr_map)) {
				corr_map[ins] = -10;
			}
		}
		permit_to_trade = false;
	}
}

/*
 Function to process income market data when the system is in the "Backtest" model.
 Only the selected signals for a specific trading day will be calculated.
 Will trigger long/short action if the prediction is over certain user-defined threshold.
 */
void Signal_Generator::onTick(std::string time_, DATA & x, int size) {
	
	//auto started = std::chrono::high_resolution_clock::now();
	tm t_ = make_time(time_);

	if (mkt_set.find(x.buffer[size - 1].InstrumentID) != mkt_set.end()) {
		if ((size >= MAX_LEN) && (bt->data_map[pair_1].mid.size() >= MAX_LEN)) {
			if ((x.buffer[size - 1].InstrumentID != trade_ins) && (x.buffer[size - 1].InstrumentID != pair_1)) {
				double cor = Ctools::compute_cos(bt->data_map[pair_1], x, corr_lookback);
				if (cor != 0) {
					corr_map[x.buffer[size - 1].InstrumentID] = cor;
				}
			}
			if (!permit_to_trade) {
				size_map[x.buffer[size - 1].InstrumentID] = true;
				permit_check();
			}
		}
	}

	if (time_check(t_) == 1){

		if ((permit_to_trade) && (bt->data_map[trade_ins].mid.size() >= MAX_LEN) && (x.buffer[size - 1].InstrumentID == trade_ins)) {
            
			std::vector<double> features;
			if (trade_ins == pair_1) {
				std::string p_max_key = Ctools::max_(corr_map);
				std::vector<double> f_1 = live_main_alpha(x, size);	
				if (mkt_set.size() > 1) {
					std::vector<double> f_2 = live_p2_alpha(bt->data_map[p_max_key], size);
					std::vector<double> f_3 = live_r2_alpha(x, bt->data_map[p_max_key]);
					features.insert(features.end(), f_1.begin(), f_1.end());
					features.insert(features.end(), f_2.begin(), f_2.end());
					features.insert(features.end(), f_3.begin(), f_3.end());
				}
				else {
					features.insert(features.end(), f_1.begin(), f_1.end());
				}
				
				std::pair<std::string, std::vector<double>> new_pp = std::make_pair(time_, features);
				features_m.push_back(new_pp);
			}
			else {
				std::string p_max_key = Ctools::max_(corr_map);
				std::vector<double> f_1 = live_main_alpha(x, size);
				std::vector<double> f_2 = live_p1_alpha(bt->data_map[pair_1], size);
				if (mkt_set.size() > 2) {
					std::vector<double> f_3 = live_p2_alpha(bt->data_map[p_max_key], size);
					std::vector<double> f_4 = live_r1_alpha(x, bt->data_map[pair_1]);
					std::vector<double> f_5 = live_r2_alpha(x, bt->data_map[p_max_key]);
					features.insert(features.end(), f_1.begin(), f_1.end());
					features.insert(features.end(), f_2.begin(), f_2.end());
					features.insert(features.end(), f_3.begin(), f_3.end());
					features.insert(features.end(), f_4.begin(), f_4.end());
					features.insert(features.end(), f_5.begin(), f_5.end());
				}
				else {
					std::vector<double> f_4 = live_r1_alpha(x, bt->data_map[pair_1]);
					features.insert(features.end(), f_1.begin(), f_1.end());
					features.insert(features.end(), f_2.begin(), f_2.end());
					features.insert(features.end(), f_4.begin(), f_4.end());
				}
				
				std::pair<std::string, std::vector<double>> new_pp = std::make_pair(time_, features);
				features_m.push_back(new_pp);
			}
		
			double pred = predict_(features);
            //std::cout << time_ << " pred: " << pred << "\n";
            
			std::pair<std::string, double> new_p = std::make_pair(time_, pred);
			matrix.push_back(new_p);
            
//			std::cout << time_ << " status check: " << status <<" q: " << ins.q << " pred: " << pred << " outstanding orders: ";
//			for (std::string o : openordertable) {
//				std::cout << o << " ";
//			}
//			std::cout << "\n";

			double bar = Ctools::std_(x.mid, lag, size);
			if (bar < (ask - bid) / 4) {
				bar = (ask - bid) / 2;
			}
            
			if (abs(pred) > (2 * (threshold + bar))) {
				pred = 0;
			}

			if (abs(pred) <= (cancel_threshold + bar)) {
				if (!openordertable.empty()) {
					status = 4;
					CancelOrder(time_);
				}
			}
           
			if (pred > (threshold + bar)) {
				if ((status == 0) || (status == 1)||(status == -1)) {
					if (ins.q < ins.contract_num) {
						status = 2;
						OpenOrder(time_, ins.contract_num - ins.q, 1, x);
						std::cout << time_ << " long: " << pred << " current bid: " << x.buffer[size - 1].BidPrice1 << " current ask: " << x.buffer[size - 1].AskPrice1 << " status: " << status << "\n";
					}
				}
				else if (status == -2) {
					status = 4;
					CancelOrder(time_);
				}
			}
			else if (pred < -(threshold + bar)) {
				if ((status == 0) || (status == 1) || (status == -1)) {
					if (ins.q > -ins.contract_num) {
						status = -2;
						OpenOrder(time_, abs(-ins.contract_num - ins.q), -1, x);
						std::cout << time_ << " short: " << pred << " current bid: " << x.buffer[size - 1].BidPrice1 << " current ask: " << x.buffer[size - 1].AskPrice1 << " status: " << status << "\n";
					}
				}
				else if (status == 2) {
					status = 4;
					CancelOrder(time_);
				}
			}
            
            /// updateing model samples ///
            std::pair<std::string, std::vector<double>> new_pp = std::make_pair(time_, features);
            model.processing_features(new_pp, x.mid[size-1]);
       
            prev_f = features;
            
			///update PNL of the instrument///
			ins.q > 0 ? ins.cur_px = bid : ins.cur_px = ask;
			ins.s_pnl = (ins.cur_px * ins.q * multiplier + ins.cash + ins.offset - ins.comm);
			//std::cout  << "PnL"<< ins.s_pnl <<" cur_px: " << ins.cur_px << " current inventory: "<< ins.q << " cash: " << ins.cash << " commission: "<< ins.comm <<"\n";
			//std::cout << ins.trade_instrument << " PnL: " << ins.s_pnl << " current inventory: "<< ins.q << " trade status: " << this->status << " cancel count: "<< ins.cancel_cout <<"\n";
		}
	}else if(time_check(t_) == 0){
		permit_to_trade = false;
		status = 0;
        model.feature_map[model.model_name].features.clear();
        model.feature_map[model.model_name].y.clear();
        
	}else {
		//std::cout << " time to close: " << t_ << " " << "night end: "<< night_end << " " << ins.tlong << " " << ins.tshort << "\n";
		CancelOrder(time_);
		if (ins.q != 0) {
			if (ins.q < 0) {
				if (abs(ins.q) < (ins.tlong + ins.tshort)) {
					bt->open_order(time_, "close", 1, ask, -ins.q, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, get_reference(), "OPENLONG");
				}
				else {
					if (ins.tshort > 0) {
						bt->open_order(time_, "close", 1, bid, ins.tshort, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, get_reference(), "CLOSESHORT");
					}
					if (ins.tlong > 0) {
						bt->open_order(time_, "close", -1, ask, ins.tlong, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, get_reference(), "CLOSELONG");
					}
				}
			}
			else if (ins.q > 0) {
				if (abs(ins.q) < (ins.tlong + ins.tshort)) {
					bt->open_order(time_, "close", -1, bid, ins.q, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, get_reference(), "OPENSHORT");
				}
				else {
					if (ins.tshort > 0) {
						bt->open_order(time_, "close", 1, bid, ins.tshort, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, get_reference(), "CLOSESHORT");
					}
					if (ins.tlong > 0) {
						bt->open_order(time_, "close", -1, ask, ins.tlong, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, get_reference(), "CLOSELONG");
					}
				}
			}
		}
		for (std::string ins : mkt_set) {
			size_map[ins] = false;
			if (Ctools::checkkeyexistence(ins, corr_map)) {
				corr_map[ins] = -10;
			}
		}
        
        model.feature_map[model.model_name].features.clear();
        model.feature_map[model.model_name].y.clear();
        permit_to_trade = false;
		status = 0;
	}

	/*auto recv_tick = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
	auto epoch = recv_tick.time_since_epoch();
	auto value = std::chrono::duration_cast<std::chrono::microseconds>(epoch);

	std::cout << "finish tick time: " << value.count() << "\n";*/

	/*auto done = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(done - started).count()<<"\n";*/
}


/*
 Function to load model parameters
 */
void Signal_Generator::load_config(){

     this->threshold = Ctools::ini_parsing(ini_file, "Threshold", "threshold")[0];
     this-> multiplier = Ctools::ini_parsing(ini_file, "Multiplier", "multiplier")[0];
     this-> ins.contract_num = Ctools::ini_parsing(ini_file, "Quantity", "quantity")[0];
     this->commission = Ctools::ini_parsing(ini_file, "Commission", "commission")[0];
	 this->cancel_threshold = Ctools::ini_parsing(ini_file, "CancelThreshold", "cancelthreshold")[0];
	 this->lag = Ctools::ini_parsing(ini_file, "Lag", "lag")[0];
	 this->main_feature = Ctools::ini_parsing_string(ini_file, "Feature", "mainfeature");
	 this->p1_feature = Ctools::ini_parsing_string(ini_file, "Feature", "p1feature");
	 this->p2_feature = Ctools::ini_parsing_string(ini_file, "Feature", "p2feature");
	 this->r1_feature = Ctools::ini_parsing_string(ini_file, "Feature", "r1feature");
	 this->r2_feature = Ctools::ini_parsing_string(ini_file, "Feature", "r2feature");
};

/*
 Function to register the server pointer in order to call order requests.
 */
void Signal_Generator::register_backtester(Backtester *bt_){
	this->bt = bt_;
}

/*
 Function to open an position
 */
void Signal_Generator::OpenOrder(std::string time,double position, int direction, DATA &x) {
	int size = static_cast<int> (x.buffer.size());
	std::string reference = get_reference();
	openordertable.insert(reference);
	
	if (direction == 1) {
		if (ins.tshort >= position) {
			std::vector<std::string> trade_log_ = {time, reference, std::to_string(position), std::to_string(ask), "Long", "CLOSESHORT", std::to_string(ins.tlong), std::to_string(ins.tshort)};
			trade_log.push_back(trade_log_);
			bt->open_order(time, "open", direction, ask, position, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, reference, "CLOSESHORT");
		}
		else {
			std::vector<std::string> trade_log_ = { time, reference, std::to_string(position), std::to_string(ask), "Long", "OPENLONG", std::to_string(ins.tlong), std::to_string(ins.tshort) };
			trade_log.push_back(trade_log_);
			bt->open_order(time, "open", direction, ask, position, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, reference, "OPENLONG");
		}
	}
	else {
		if (ins.tlong >= position) {
			std::vector<std::string> trade_log_ = { time, reference, std::to_string(position), std::to_string(bid), "SHORT", "CLOSELONG", std::to_string(ins.tlong), std::to_string(ins.tshort) };
			trade_log.push_back(trade_log_);
			bt->open_order(time, "open", direction, bid, position, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, reference, "CLOSELONG");
		}
		else {
			std::vector<std::string> trade_log_ = { time, reference, std::to_string(position), std::to_string(bid), "SHORT", "OPENSHORT", std::to_string(ins.tlong), std::to_string(ins.tshort) };
			trade_log.push_back(trade_log_);
			bt->open_order(time, "open", direction, bid, position, x.buffer[size - 1].BidVolume1, x.buffer[size - 1].AskVolume1, reference, "OPENSHORT");
		}
	}
};

/*
 Function to initiate an order cancel request.
 */
void Signal_Generator::CancelOrder(std::string time) {
    std::set<std::string>::iterator it = openordertable.begin();
    while (it != openordertable.end()) {
        std::cout << "canceling " << *it << "\n";
		bt->process_cancel(time, *it);
        if (openordertable.size() == 0){
            break;
        }
	}
};


/*
 Function to receive trade confirmation from the server side if the order is succesfully trade.
 Will update position inventory of the logic
 */
void Signal_Generator::onTrade(std::string time, int direction, double price, double quantity, double orig_quant, std::string ref, std::string action_flag){

	if (action_flag == "OPENLONG") {
		ins.tlong += quantity;
	}
	else if (action_flag == "CLOSELONG") {
		ins.tlong -= quantity;
	}
	else if (action_flag == "OPENSHORT") {
		ins.tshort += quantity;
	}
	else {
		ins.tshort -= quantity;
	}

	ins.q = ins.tlong - ins.tshort;
	direction == 1? ins.cash -= (quantity * price * multiplier) : ins.cash += (quantity * price * multiplier);
	ins.comm += commission * quantity;


	if (quantity == orig_quant) {
		openordertable.erase(ref);
		if (ins.q > 0.) {
			status = 1;
		}
		else if (ins.q < 0.) {
			status = -1;
		}
		else {
			status = 0;
		}
	}else{
		if (direction == 1) {
			status = 2;
		}
		else if (direction == -1) {
			status = -2;
		}
	
	}
	
	
	total_trade ++;

	std::string temp;
	direction == 1 ? temp = "Long" : temp = "SHORT";
	std::vector<std::string> trade_log_ = { time, ref, std::to_string(quantity), std::to_string(price), temp, "TRADED", std::to_string(ins.tlong), std::to_string(ins.tshort)};
	trade_log.push_back(trade_log_);

	ins.q > 0 ? ins.cur_px = bid : ins.cur_px = ask;
	ins.s_pnl = (ins.cur_px * ins.q * multiplier + ins.cash + ins.offset - ins.comm);

	if (direction == 1) {
		std::cout << "TRADED status: "<< status <<" buy at: " << price << " pnl: " << ins.s_pnl <<"  fill q: "<< quantity
			<< " current price: " << ins.cur_px <<" current q: " << ins.q << " commission: "<< ins.comm << " cash balance: "<< ins.cash << " current bid: " << bid << " current ask: " << ask
			<< " tlong: "<< ins.tlong << " tshort: "<< ins.tshort <<"\n";
	}
	else {
		std::cout << "TRADED status: " << status << " sell at: " << price << " pnl: " << ins.s_pnl << " fill q: " << quantity
			<< " current price: " << ins.cur_px << " current q: " << ins.q  << " commission: " << ins.comm << " cash balance: " << ins.cash << " current bid: " << bid << " current ask: "<< ask
			<< " tlong: " << ins.tlong << " tshort: " << ins.tshort  << "\n";
	}
};


/*
 Function to process order cancel confirmation from the server side
 */
void Signal_Generator::onCancel(std::string time, std::string ref){
	std::cout << "on cancel. \n";
	if (openordertable.find(ref) != openordertable.end()) {
		openordertable.erase(ref);
		std::vector<std::string> trade_log_ = { time, ref, " ", " ", " ", "CANCELED", std::to_string(ins.tlong), std::to_string(ins.tshort) };
		trade_log.push_back(trade_log_);
	}

	cancel_count++;
	
	if (ins.q > 0) {
		status = 1;
	}
	else if (ins.q < 0) {
		status = -1;
	}
	else {
		status = 0;
	}
    std::cout << "cancel done \n";
}

/*
 Function to check whether current time is under the user-defined allowable trading time period.
 */

int Signal_Generator::time_check(tm &timestamp){
	
	if (((timestamp.tm_hour >= 21) && (timestamp.tm_hour < 23)) || ((timestamp.tm_hour >= 9) && (timestamp.tm_hour < 11)) || ((timestamp.tm_hour >= 13) && (timestamp.tm_hour < 15))) {
		if ((timestamp.tm_hour == 22) || (timestamp.tm_hour == 14)) {
			if (timestamp.tm_min > 58) {
				return 2;
			}
			else {
				return 1;
			}
		}
        else if (timestamp.tm_hour == 11){
           if (timestamp.tm_min > 28) {
                return 2;
            }
            else {
                return 1;
            }
        }
		else if (timestamp.tm_hour == 10) {
			if ((timestamp.tm_min > 13) && (timestamp.tm_min < 35)) {
				return 2;
			}
			else {
				return 1;
			}
		}
		else if (((timestamp.tm_hour == 21) || (timestamp.tm_hour == 9))) {
			if (timestamp.tm_min < 5) {
				return 2;
			}
			else {
				return 1;
			}
		}
		else if (timestamp.tm_hour == 13) {
			if (timestamp.tm_min < 35) {
				return 2;
			}
			else {
				return 1;
			}
		}
	}
	else {
		return 0;
	}
    return 2;
}


/*
 Function to update the main traded product's current bid ask price
 */
void Signal_Generator::updateBidAsk(CThostFtdcDepthMarketDataField &p){
	this->bid = p.BidPrice1;
	this->ask = p.AskPrice1;
}


/*
 Function to generate the full set of alpha signals for its main traded asset
 */
std::vector<double> Signal_Generator::alpha(DATA &data, int size) {
	std::vector<double> f_;
	size_t size_ = header_.size();
	for (size_t i = 1; i < size_; i ++) {
		//std::cout << " " << header[i];
		f_.push_back(func_map[header_[i]](data, size));
	}
	//std::cout << "\n";
	size_t f_size = f_.size();
	for (size_t i = 0; i < f_size; i++) {
		if (isnan(f_[i])) {
			f_[i] = 0;
		}
		else if (!isfinite(f_[i])) {
			f_[i] = 0;
		}
	}
	return f_;
};

/*
Function to generate the full set of alpha signals with respect to the correlation pair asset
*/
std::vector<double> Signal_Generator::pair_alpha(DATA& data, int size) {
	std::vector<double> f_;
	for (std::string ins : header_p) {
		//std::cout << " " << ins;
		f_.push_back(func_map[ins](data, size));
	}
	//std::cout << "\n";
	size_t f_size = f_.size();
	for (size_t i = 0; i < f_size; i++) {
		if (isnan(f_[i])) {
			f_[i] = 0;
		}
		else if (!isfinite(f_[i])) {
			f_[i] = 0;
		}
	}
	return f_;
}

/*
Function to generate the full set of cross-pair alpha signals
*/
std::vector<double> Signal_Generator::r_alpha(DATA& data_1, DATA& data_2) {
	std::vector<double> f_;
	for (std::string ins : header_r) {
		//std::cout << " " << ins;
		f_.push_back(Rfunc_map[ins](data_1, data_2));
	}
	//std::cout << "\n";
	size_t f_size = f_.size();
	for (size_t i = 0; i < f_size; i++) {
		if (isnan(f_[i])) {
			f_[i] = 0;
		}
		else if (!isfinite(f_[i])) {
			f_[i] = 0;
		}
	}
	return f_;
}

/*
Function to generate alpha signals for its main traded asset
*/
std::vector<double> Signal_Generator::live_main_alpha(DATA& data, int size) {
	std::vector<double> f_;
	if (!main_feature.empty()) {
		for (std::string f_name : main_feature) {
			//std::cout << "got in. "<< f_name <<" "<<func_map[f_name](data, size) <<"\n";
			f_.push_back(func_map[f_name](data, size));
		}

		size_t f_size = f_.size();
		for (size_t i = 0; i < f_size; i++) {
			if (isnan(f_[i])) {
				f_[i] = 0;
			}
			else if (!isfinite(f_[i])) {
				f_[i] = 0;
			}
		}
	}
	return f_;
};

/*
Function to generate alpha signals for its pair reference asset
*/
std::vector<double> Signal_Generator::live_p1_alpha(DATA& data, int size) {
	std::vector<double> f_;
	if (!p1_feature.empty()) {
		for (std::string f_name : p1_feature) {
			//std::cout << "got in. " << f_name << " " << func_map[f_name](data, size) << "\n";
			f_.push_back(func_map[f_name](data, size));
		}

		size_t f_size = f_.size();
		for (size_t i = 0; i < f_size; i++) {
			if (isnan(f_[i])) {
				f_[i] = 0;
			}
			else if (!isfinite(f_[i])) {
				f_[i] = 0;
			}
		}
	}
	return f_;
};

/*
Function to generate pair asset alpha signals
*/
std::vector<double> Signal_Generator::live_p2_alpha(DATA& data, int size) {
	std::vector<double> f_;
	if (!p2_feature.empty()) {
		for (std::string f_name : p2_feature) {
			//std::cout << " " << f_name;
			f_.push_back(func_map[f_name](data, size));
		}
		//std::cout << "\n";
		size_t f_size = f_.size();
		for (size_t i = 0; i < f_size; i++) {
			if (isnan(f_[i])) {
				f_[i] = 0;
			}
			else if (!isfinite(f_[i])) {
				f_[i] = 0;
			}
		}
	}
	
	return f_;

};

/*
Function to generate cross-pair alpha signals
*/
std::vector<double> Signal_Generator::live_r1_alpha(DATA& data_1, DATA& data_2) {
	std::vector<double> f_;
	if (!r1_feature.empty()) {
		for (std::string f_name : r1_feature) {
			//std::cout << " " << f_name;
			f_.push_back(Rfunc_map[f_name](data_1, data_2));
		}
		//std::cout << "\n";
		size_t f_size = f_.size();
		for (size_t i = 0; i < f_size; i++) {
			if (isnan(f_[i])) {
				f_[i] = 0;
			}
			else if (!isfinite(f_[i])) {
				f_[i] = 0;
			}
		}
	}
	return f_;
}

/*
Function to generate cross-pair alpha signals
*/
std::vector<double> Signal_Generator::live_r2_alpha(DATA& data_1, DATA& data_2) {
	std::vector<double> f_;
	if (!r2_feature.empty()) {
		for (std::string f_name : r2_feature) {
		//	std::cout << " " <<f_name;
			f_.push_back(Rfunc_map[f_name](data_1, data_2));
		}
		//std::cout << "\n";
		size_t f_size = f_.size();
		for (size_t i = 0; i < f_size; i++) {
			if (isnan(f_[i])) {
				f_[i] = 0;
			}
			else if (!isfinite(f_[i])) {
				f_[i] = 0;
			}
		}
	}
	return f_;
}

/*
Function to perform prediction.
 Parameter: vector of features
*/
double Signal_Generator::predict_(std::vector<double> f_) {
    double pred = 0;
    
//    if (prev_f.size()>2){
//        colvec pre(prev_f);
        colvec cur_predictors(f_);
//       colvec predictors = arma::join_cols(pre, cur_predictors);
        rowvec y;
        model.regressor.Predict(cur_predictors, y);
        pred = y[0];
   // }
    
    return pred;
};

/*
 Admin function that checks whether all referend assets have collected long enough tick data feed.
 Ensure signal processing does not go out of bound
*/
void Signal_Generator::permit_check() {
	int i = 0;
	for (std::pair<std::string, bool> element : size_map) {
		if (!element.second) {
			i++;
		}
	}
	i == 0 ? permit_to_trade = true : permit_to_trade = false;
};

std::string Signal_Generator::get_reference() {
	tagCounter++;
	std::string ref = std::to_string(tagCounter);
	return ref;
};


/*
 Function to register all the alpha processing functions from the Ctools library.
 */
void Signal_Generator::register_func_map() {

	
	Cfunc o010 = &Ctools::o010_;
	func_map["o010"] = o010;
	Cfunc o011 = &Ctools::o011_;
	func_map["o011"] = o011;
	Cfunc o012 = &Ctools::o012_;
	func_map["o012"] = o012;
	Cfunc o013 = &Ctools::o013_;
	func_map["o013"] = o013;
	Cfunc o014 = &Ctools::o014_;
	func_map["o014"] = o014;

	Cfunc f085 = &Ctools::f085_;
	func_map["f085"] = f085;

	Cfunc a1000 = &Ctools::a1000_;
	func_map["a1000"] = a1000;
	Cfunc a1001 = &Ctools::a1001_;
	func_map["a1001"] = a1001;
	Cfunc a1002 = &Ctools::a1002_;
	func_map["a1002"] = a1002;
	Cfunc a1003 = &Ctools::a1003_;
	func_map["a1003"] = a1003;
	Cfunc a1004 = &Ctools::a1004_;
	func_map["a1004"] = a1004;

	Cfunc c100 = &Ctools::c100_;
	func_map["c100"] = c100;
	Cfunc c101 = &Ctools::c101_;
	func_map["c101"] = c101;
	Cfunc c102 = &Ctools::c102_;
	func_map["c102"] = c102;
	Cfunc c103 = &Ctools::c103_;
	func_map["c103"] = c103;
	Cfunc c104 = &Ctools::c104_;
	func_map["c104"] = c104;

	Cfunc c106 = &Ctools::c106_;
	func_map["c106"] = c106;
	Cfunc c107 = &Ctools::c107_;
	func_map["c107"] = c107;
	Cfunc c108 = &Ctools::c108_;
	func_map["c108"] = c108;
	Cfunc c109 = &Ctools::c109_;
	func_map["c109"] = c109;

	Cfunc c2 = &Ctools::c2_;
	func_map["c2"] = c2;
	Cfunc c3 = &Ctools::c3_;
	func_map["c3"] = c3;
	Cfunc c4 = &Ctools::c4_;
	func_map["c4"] = c4;
	Cfunc c5 = &Ctools::c5_;
	func_map["c5"] = c5;
	
	Cfunc c6 = &Ctools::c6_;
	func_map["c6"] = c6;
	Cfunc c7 = &Ctools::c7_;
	func_map["c7"] = c7;
	Cfunc c8 = &Ctools::c8_;
	func_map["c8"] = c8;

	Cfunc c9 = &Ctools::c9_;
	func_map["c9"] = c9;
	Cfunc c10 = &Ctools::c10_;
	func_map["c10"] = c10;
	Cfunc c11 = &Ctools::c11_;
	func_map["c11"] = c11;

	Cfunc f301 = &Ctools::f301_;
	func_map["f301"] = f301;
	Cfunc f304 = &Ctools::f304_;
	func_map["f304"] = f304;
	Cfunc f305 = &Ctools::f305_;
	func_map["f305"] = f305;
	Cfunc f306 = &Ctools::f306_;
	func_map["f306"] = f306;
	Cfunc f307 = &Ctools::f307_;
	func_map["f307"] = f307;
	Cfunc f318 = &Ctools::f318_;
	func_map["f318"] = f318;
	Cfunc f320 = &Ctools::f320_;
	func_map["f320"] = f320;
	Cfunc f321 = &Ctools::f321_;
	func_map["f321"] = f321;
	Cfunc f324 = &Ctools::f324_;
	func_map["f324"] = f324;
	Cfunc f326 = &Ctools::f326_;
	func_map["f326"] = f326;

	Cfunc f319 = &Ctools::f319_;
	func_map["f319"] = f319;
	Cfunc f333 = &Ctools::f333_;
	func_map["f333"] = f333;
	Cfunc f334 = &Ctools::f334_;
	func_map["f334"] = f334;
	Cfunc f335 = &Ctools::f335_;
	func_map["f335"] = f335;

	Cfunc f091 = &Ctools::f091_;
	func_map["f091"] = f091;
	Cfunc f096 = &Ctools::f096_;
	func_map["f096"] = f096;
	Cfunc f097 = &Ctools::f097_;
	func_map["f097"] = f097;
	Cfunc f098 = &Ctools::f098_;
	func_map["f098"] = f098;

	Cfunc f043 = &Ctools::f043_;
	func_map["f043"] = f043;

	Cfunc z_10 = &Ctools::z10_;
	func_map["z10"] = z_10;
	Cfunc z_30 = &Ctools::z30_;
	func_map["z30"] = z_30;
	Cfunc z_60 = &Ctools::z60_;
	func_map["z60"] = z_60;

	Cfunc q = &Ctools::qt;
	func_map["q"] = q;

	Cfunc a500 = &Ctools::a500_;
	func_map["a500"] = a500;
	Cfunc a501 = &Ctools::a501_;
	func_map["a501"] = a501;

	Cfunc a1400 = &Ctools::a1400_;
	func_map["a1400"] = a1400;
	Cfunc a1401 = &Ctools::a1401_;
	func_map["a1401"] = a1401;
	Cfunc a1402 = &Ctools::a1402_;
	func_map["a1402"] = a1402;
	Cfunc a1403 = &Ctools::a1403_;
	func_map["a1403"] = a1403;

	Cfunc b1400 = &Ctools::b1400_;
	func_map["b1400"] = b1400;
	Cfunc b1401 = &Ctools::b1401_;
	func_map["b1401"] = b1401;
	Cfunc b1402 = &Ctools::b1402_;
	func_map["b1402"] = b1402;
	Cfunc b1403 = &Ctools::b1403_;
	func_map["b1403"] = b1403;

	Cfunc a1800 = &Ctools::a1800_;
	func_map["a1800"] = a1800;
	Cfunc a1801 = &Ctools::a1801_;
	func_map["a1801"] = a1801;
	Cfunc a1802 = &Ctools::a1802_;
	func_map["a1802"] = a1802;
	Cfunc a1803 = &Ctools::a1803_;
	func_map["a1803"] = a1803;

	Cfunc a2500 = &Ctools::a2500_;
	func_map["a2500"] = a2500;
	Cfunc a2501 = &Ctools::a2501_;
	func_map["a2501"] = a2501; 
	Cfunc a2502 = &Ctools::a2502_;
	func_map["a2502"] = a2502;

	Cfunc a4200 = &Ctools::a4200_;
	func_map["a4200"] = a4200;
	Cfunc a4201 = &Ctools::a4201_;
	func_map["a4201"] = a4201;
	Cfunc a4202 = &Ctools::a4202_;
	func_map["a4202"] = a4202;
	Cfunc a4203 = &Ctools::a4203_;
	func_map["a4203"] = a4203;

	Cfunc b4201 = &Ctools::b4201_;
	func_map["b4201"] = b4201;
	Cfunc b4202 = &Ctools::b4202_;
	func_map["b4202"] = b4202;
	Cfunc b4203 = &Ctools::b4203_;
	func_map["b4203"] = b4203;

	Cfunc c4201 = &Ctools::c4201_;
	func_map["c4201"] = c4201;
	Cfunc c4202 = &Ctools::c4202_;
	func_map["c4202"] = c4202;

	Cfunc a7800 = &Ctools::a7800_;
	func_map["a7800"] = a7800;
	Cfunc a7801 = &Ctools::a7801_;
	func_map["a7801"] = a7801;
	Cfunc a7802 = &Ctools::a7802_;
	func_map["a7802"] = a7802;

	Cfunc a7803 = &Ctools::a7803_;
	func_map["a7803"] = a7803;
	Cfunc a7804 = &Ctools::a7804_;
	func_map["a7804"] = a7804;
	Cfunc a7805 = &Ctools::a7805_;
	func_map["a7805"] = a7805;

	Cfunc f341 = &Ctools::f341_;
	func_map["f341"] = f341;

	Cfunc f083 = &Ctools::f083_;
	func_map["f083"] = f083;

	Cfunc o000 = &Ctools::o000_;
	func_map["o000"] = o000;
	Cfunc o001 = &Ctools::o001_;
	func_map["o001"] = o001;
	Cfunc o002 = &Ctools::o002_;
	func_map["o002"] = o002;
	Cfunc o003 = &Ctools::o003_;
	func_map["o003"] = o003;
	Cfunc o004 = &Ctools::o004_;
	func_map["o004"] = o004;

	Cfunc o100 = &Ctools::o100_;
	func_map["o100"] = o100;
	Cfunc o101 = &Ctools::o101_;
	func_map["o101"] = o101;
	Cfunc o102 = &Ctools::o102_;
	func_map["o102"] = o102;
	Cfunc o103 = &Ctools::o103_;
	func_map["o103"] = o103;
	Cfunc o104 = &Ctools::o104_;
	func_map["o104"] = o104;

	Cfunc o200 = &Ctools::o200_;
	func_map["o200"] = o200;
	Cfunc o201 = &Ctools::o201_;
	func_map["o201"] = o201;
	Cfunc o202 = &Ctools::o202_;
	func_map["o202"] = o202;
	Cfunc o203 = &Ctools::o203_;
	func_map["o203"] = o203;
	Cfunc o204 = &Ctools::o204_;
	func_map["o204"] = o204;

	Cfunc o300 = &Ctools::o300_;
	func_map["o300"] = o300;
	Cfunc o301 = &Ctools::o301_;
	func_map["o301"] = o301;
	Cfunc o302 = &Ctools::o302_;
	func_map["o302"] = o302;
	Cfunc o303 = &Ctools::o303_;
	func_map["o303"] = o303;
	Cfunc o304 = &Ctools::o304_;
	func_map["o304"] = o304;

	Cfunc o400 = &Ctools::o400_;
	func_map["o400"] = o400;
	Cfunc o401 = &Ctools::o401_;
	func_map["o401"] = o401;
	Cfunc o402 = &Ctools::o402_;
	func_map["o402"] = o402;
	Cfunc o403 = &Ctools::o403_;
	func_map["o403"] = o403;
	Cfunc o404 = &Ctools::o404_;
	func_map["o404"] = o404;

	Cfunc c12 = &Ctools::c12_;
	func_map["c12"] = c12;
	Cfunc c13 = &Ctools::c13_;
	func_map["c13"] = c13;
	Cfunc c14 = &Ctools::c14_;
	func_map["c14"] = c14;
	Cfunc c15 = &Ctools::c15_;
	func_map["c15"] = c15;
	Cfunc c16 = &Ctools::c16_;
	func_map["c16"] = c16;
	Cfunc c17 = &Ctools::c17_;
	func_map["c17"] = c17;
	Cfunc c18 = &Ctools::c18_;
	func_map["c18"] = c18;
	Cfunc c19 = &Ctools::c19_;
	func_map["c19"] = c19;
	Cfunc c20 = &Ctools::c20_;
	func_map["c20"] = c20;
	Cfunc c21 = &Ctools::c21_;
	func_map["c21"] = c21;
	Cfunc c22 = &Ctools::c22_;
	func_map["c22"] = c22;
	Cfunc c23 = &Ctools::c23_;
	func_map["c23"] = c23;
	Cfunc c24 = &Ctools::c24_;
	func_map["c24"] = c24;
	Cfunc c25 = &Ctools::c25_;
	func_map["c25"] = c25;
	Cfunc c26 = &Ctools::c26_;
	func_map["c26"] = c26;
    Cfunc c27 = &Ctools::c27_;
    func_map["c27"] = c27;
    
	Cfunc t100 = &Ctools::t100_;
	func_map["t100"] = t100;
	Cfunc t101 = &Ctools::t101_;
	func_map["t101"] = t101;
	Cfunc t102 = &Ctools::t102_;
	func_map["t102"] = t102;
	Cfunc t103 = &Ctools::t103_;
	func_map["t103"] = t103;
	Cfunc t104 = &Ctools::t104_;
	func_map["t104"] = t104;
	Cfunc t105 = &Ctools::t105_;
	func_map["t105"] = t105;

	Cfunc t200 = &Ctools::t200_;
	func_map["t200"] = t200;
	Cfunc t201 = &Ctools::t201_;
	func_map["t201"] = t201;
	Cfunc t203 = &Ctools::t203_;
	func_map["t203"] = t203;

	Cfunc t303 = &Ctools::t303_;
	func_map["t303"] = t303;
	Cfunc t304 = &Ctools::t304_;
	func_map["t304"] = t304;

	Cfunc e100 = &Ctools::e100_;
	func_map["e100"] = e100;
	Cfunc e101 = &Ctools::e101_;
	func_map["e101"] = e101;
	Cfunc e102 = &Ctools::e102_;
	func_map["e102"] = e102;
	Cfunc e103 = &Ctools::e103_;
	func_map["e103"] = e103;

	Cfunc e200 = &Ctools::e200_;
	func_map["e200"] = e200;
	Cfunc e201 = &Ctools::e201_;
	func_map["e201"] = e201;
	Cfunc e202 = &Ctools::e202_;
	func_map["e202"] = e202;
	Cfunc e203 = &Ctools::e203_;
	func_map["e203"] = e203;


	Cfunc v100 = &Ctools::v100_;
	func_map["v100"] = v100;
	Cfunc v101 = &Ctools::v101_;
	func_map["v101"] = v101;
	Cfunc v102 = &Ctools::v102_;
	func_map["v102"] = v102;
	Cfunc v103 = &Ctools::v103_;
	func_map["v103"] = v103;
	Cfunc v104 = &Ctools::v104_;
	func_map["v104"] = v104;

	Cfunc v200 = &Ctools::v200_;
	func_map["v200"] = v200;
	Cfunc v201 = &Ctools::v201_;
	func_map["v201"] = v201;
	Cfunc v202 = &Ctools::v202_;
	func_map["v202"] = v202;
	Cfunc v203 = &Ctools::v203_;
	func_map["v203"] = v203;
	Cfunc v204 = &Ctools::v204_;
	func_map["v204"] = v204;

	Cfunc v300 = &Ctools::v300_;
	func_map["v300"] = v300;
	Cfunc v301 = &Ctools::v301_;
	func_map["v301"] = v301;
	Cfunc v302 = &Ctools::v302_;
	func_map["v302"] = v302;
	Cfunc v303 = &Ctools::v303_;
	func_map["v303"] = v303;
	Cfunc v304 = &Ctools::v104_;
	func_map["v304"] = v304;





	Rfunc r100 = &Ctools::r100_;
	Rfunc_map["r100"] = r100;
	Rfunc r101 = &Ctools::r101_;
	Rfunc_map["r101"] = r101;
	Rfunc r102 = &Ctools::r102_;
	Rfunc_map["r102"] = r102;
	Rfunc r103 = &Ctools::r103_;
	Rfunc_map["r103"] = r103;

	Rfunc cr100 = &Ctools::cr100_;
	Rfunc_map["cr100"] = cr100;
	Rfunc cr101 = &Ctools::cr101_;
	Rfunc_map["cr101"] = cr101;
	Rfunc cr102 = &Ctools::cr102_;
	Rfunc_map["cr102"] = cr102;
	Rfunc cr103 = &Ctools::cr103_;
	Rfunc_map["cr103"] = cr103;
	Rfunc cr104 = &Ctools::cr104_;
	Rfunc_map["cr104"] = cr104;
	Rfunc cr105 = &Ctools::cr105_;
	Rfunc_map["cr105"] = cr105;
	Rfunc cr106 = &Ctools::cr106_;
	Rfunc_map["cr106"] = cr106;
	Rfunc cr107 = &Ctools::cr107_;
	Rfunc_map["cr107"] = cr107;

	Rfunc cr108 = &Ctools::cr108_;
	Rfunc_map["cr108"] = cr108;
	Rfunc cr109 = &Ctools::cr109_;
	Rfunc_map["cr109"] = cr109;
	Rfunc cr110 = &Ctools::cr110_;
	Rfunc_map["cr110"] = cr110;
	Rfunc cr111 = &Ctools::cr111_;
	Rfunc_map["cr111"] = cr111;
	Rfunc cr112 = &Ctools::cr112_;
	Rfunc_map["cr112"] = cr112;
	Rfunc cr113 = &Ctools::cr113_;
	Rfunc_map["cr113"] = cr113;
	Rfunc cr114 = &Ctools::cr114_;
	Rfunc_map["cr114"] = cr114;
	Rfunc cr115 = &Ctools::cr115_;
	Rfunc_map["cr115"] = cr115;
	Rfunc cr116 = &Ctools::cr116_;
	Rfunc_map["cr116"] = cr116;
	Rfunc cr117 = &Ctools::cr117_;
	Rfunc_map["cr117"] = cr117;

	Rfunc vwr100 = &Ctools::vwr100_;
	Rfunc_map["vwr100"] = vwr100;
	Rfunc vwr101 = &Ctools::vwr101_;
	Rfunc_map["vwr101"] = vwr101;
	Rfunc vwr102 = &Ctools::vwr102_;
	Rfunc_map["vwr102"] = vwr102;
	
}





