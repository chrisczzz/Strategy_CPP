# Strategy_CPP

FURTHER PERFORMANCE OPTIMIZATION NEEDED.

A sample script that is taken from C++ backtest system. 

The full C++ backtest system simulates China commodity futures market trading. 

The sample script processes incoming market data from the market data processor side and takes appropriate actions according to the logic. Class also has several functions that supports administrative tasks such as PNL calculation, order management.

Strategy description:
1. Strategy focuses on signal mining from the main traded product's orderbook dynamics.
2. Signals are further combined with signals from highly-correlated assets. High correlated asset set is analyzed and selected    from statistical studies conducted in Python.
3. Strategy parameters would be real time updated as the new market data feeds into the logic through Machine learning            class object.
4. Strategy is for liquidity taking. 
