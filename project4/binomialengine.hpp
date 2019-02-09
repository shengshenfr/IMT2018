/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003, 2004 Ferdinando Ametrano
 Copyright (C) 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2007 StatPro Italia srl
 Copyright (C) 2007 Affine Group Limited

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file binomialengine.hpp
    \brief Binomial option engine
*/

#ifndef binomial_engine_hpp
#define binomial_engine_hpp

#include <ql/methods/lattices/binomialtree.hpp>
#include <ql/methods/lattices/bsmlattice.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/vanilla/discretizedvanillaoption.hpp>
#include <ql/pricingengines/greeks.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <iostream>
#include <iomanip>
#include <ql/pricingengines/blackscholescalculator.hpp>

namespace QuantLib {

    //! Pricing engine for vanilla options using binomial trees
    /*! \ingroup vanillaengines

        \test the correctness of the returned values is tested by
              checking it against analytic results.

        \todo Greeks are not overly accurate. They could be improved
              by building a tree so that it has three points at the
              current time. The value would be fetched from the middle
              one, while the two side points would be used for
              estimating partial derivatives.
    */
    template <class T>
    class BinomialVanillaEngine_2 : public VanillaOption::engine {
      public:
        BinomialVanillaEngine_2(
             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
             Size timeSteps)
        : process_(process), timeSteps_(timeSteps) {
            QL_REQUIRE(timeSteps >= 2,
                       "at least 2 time steps required, "
                       << timeSteps << " provided");
            registerWith(process_);
        }
        void calculate() const;
      private:
        boost::shared_ptr<GeneralizedBlackScholesProcess> process_;
        Size timeSteps_;
    };
    
    
    
    
	template <class T>
	class BlackScholesLattice_2 : public TreeLattice1D<BlackScholesLattice_2<T>> {
	public:
	BlackScholesLattice_2(const boost::shared_ptr<T>& tree,
										Rate riskFreeRate,
										Rate dividendYield , 
										Time end,
										Size steps , 
										Real volatility, 
										Real strike)
	: TreeLattice1D< BlackScholesLattice_2 <T> >(TimeGrid(end, steps), 2),
	tree_(tree), riskFreeRate_(riskFreeRate), dt_(end/steps),
		discount_(std::exp(-riskFreeRate*(dt_))) {
		pd_ = tree->probability(0, 0, 0);
		pu_ = tree->probability(0, 0, 1);
		t_ = TimeGrid(end, steps) ; 
		volatility_ = volatility ; 
		strike_ =  strike   ; 
		dividendYield_ =  dividendYield ; 
		}
		Size size(Size i) const { return tree_->size(i); }
		DiscountFactor discount(Size,Size) const { return discount_; }
		Real underlying(Size i, Size index) const {
						 return tree_->underlying(i, index);}
		Size descendant(Size i, Size index, Size branch) const {
						 return tree_->descendant(i, index, branch);}
		Real probability(Size i, Size index, Size branch) const {
						 return tree_->probability(i, index, branch);}
		void stepback(Size i, const Array& values, Array& newValues) const {
			for (Size j=0; j<size(i); j++)
			  newValues[j] = (pd_*values[j] + pu_*values[j+1])*discount_;
				
			}
		void partialRollback(DiscretizedAsset & asset, Time to) const {
			 //rollback one step  from the very End of the Lattice 
			// ajustAsset (asset) ;  
			const  TimeGrid& grid = t_;
			Size lastTimeRef =  grid.size() ; 
			Integer iFrom = Integer(t_.index(asset.time())) ;
			Integer iTo = Integer(t_.index(to));
			for (Integer i=iFrom-1; i>=iTo; --i) {
			   Array newValues(this->impl().size(i));
			   if(i == Integer(lastTimeRef -2)){
				  ajustAsset (asset, lastTimeRef); 
			   }else{
			   this->impl().stepback(i, asset.values(), newValues);
				// std::cout << "value of the asset " << asset.values() << std::endl;
			   asset.time() =t_[i];
			   asset.values() = newValues;
				if (i != iTo) // skip the very last adjustment
				asset.adjustValues();
			   }
			  }
		   }
		void  ajustAsset (DiscretizedAsset & asset, Size lastTimeRef)  const {
				//1 - change the value of the of the current data in the Lattice with the blackScholes formula  
				Array coxRoxValues =  asset.values();
				boost::shared_ptr <PlainVanillaPayoff> payOffCall(new  PlainVanillaPayoff (Option::Type::Call,strike_)); 
				Real vol =  volatility_ *std::sqrt(dt_); 
				DiscountFactor growth = std::exp(-(dividendYield_)*dt_);
				for (int i =0 ; i<Integer(coxRoxValues.size()); i++) {	  
				   BlackScholesCalculator bsCalculator(payOffCall, underlying(lastTimeRef -2 , i), growth, vol, discount_) ; // use the Blackschole calculator
				   coxRoxValues[i] = bsCalculator.value() ; 
				   //2 - compute the new values of the with the blackScholes formula  (coxRoxValues(i) = bsformula (coxRoxValues(i) ))
				}
				asset.values() = coxRoxValues;     
				//3 - asset.values() = coxRoxValues
			
		}
	protected:
        boost::shared_ptr<T> tree_;
        Rate riskFreeRate_;
        Time dt_;
        DiscountFactor discount_;
        Real pd_, pu_;
        TimeGrid t_; 
        Real volatility_; 
        Real strike_;  
        Rate dividendYield_;    
	};



    // template definitions

    template <class T>
    void BinomialVanillaEngine_2<T>::calculate() const {

        DayCounter rfdc  = process_->riskFreeRate()->dayCounter();
        DayCounter divdc = process_->dividendYield()->dayCounter();
        DayCounter voldc = process_->blackVolatility()->dayCounter();
        Calendar volcal = process_->blackVolatility()->calendar();

        Real s0 = process_->stateVariable()->value();
        QL_REQUIRE(s0 > 0.0, "negative or null underlying given");
        Volatility v = process_->blackVolatility()->blackVol(
            arguments_.exercise->lastDate(), s0);
        Date maturityDate = arguments_.exercise->lastDate();
        Rate r = process_->riskFreeRate()->zeroRate(maturityDate,
            rfdc, Continuous, NoFrequency);
        Rate q = process_->dividendYield()->zeroRate(maturityDate,
            divdc, Continuous, NoFrequency);
        Date referenceDate = process_->riskFreeRate()->referenceDate();

        // binomial trees with constant coefficient
        Handle<YieldTermStructure> flatRiskFree(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(referenceDate, r, rfdc)));
        Handle<YieldTermStructure> flatDividends(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(referenceDate, q, divdc)));
        Handle<BlackVolTermStructure> flatVol(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(referenceDate, volcal, v, voldc)));

        boost::shared_ptr<PlainVanillaPayoff> payoff =
            boost::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");

        Time maturity = rfdc.yearFraction(referenceDate, maturityDate);

        boost::shared_ptr<StochasticProcess1D> bs(
                         new GeneralizedBlackScholesProcess(
                                      process_->stateVariable(),
                                      flatDividends, flatRiskFree, flatVol));

        TimeGrid grid(maturity, timeSteps_);

        boost::shared_ptr<T> tree(new T(bs, maturity, timeSteps_,
                                        payoff->strike()));
		boost::shared_ptr<BlackScholesLattice_2<T> > lattice(new BlackScholesLattice_2 <T> (tree, r, q, maturity, timeSteps_, v, payoff->strike()));

        DiscretizedVanillaOption option(arguments_, *process_, grid);

        option.initialize(lattice, maturity);

        // Partial derivatives calculated from various points in the
        // binomial tree 
        // (see J.C.Hull, "Options, Futures and other derivatives", 6th edition, pp 397/398)

        // Rollback to third-last step, and get underlying prices (s2) &
        // option values (p2) at this point
        option.rollback(grid[2]);
        Array va2(option.values());
        QL_ENSURE(va2.size() == 3, "Expect 3 nodes in grid at second step");
        Real p2u = va2[2]; // up
        Real p2m = va2[1]; // mid
        Real p2d = va2[0]; // down (low)
        Real s2u = lattice->underlying(2, 2); // up price
        Real s2m = lattice->underlying(2, 1); // middle price
        Real s2d = lattice->underlying(2, 0); // down (low) price

        // calculate gamma by taking the first derivate of the two deltas
        Real delta2u = (p2u - p2m)/(s2u-s2m);
        Real delta2d = (p2m-p2d)/(s2m-s2d);
        Real gamma = (delta2u - delta2d) / ((s2u-s2d)/2);

        // Rollback to second-last step, and get option values (p1) at
        // this point
        option.rollback(grid[1]);
        Array va(option.values());
        QL_ENSURE(va.size() == 2, "Expect 2 nodes in grid at first step");
        Real p1u = va[1];
        Real p1d = va[0];
        Real s1u = lattice->underlying(1, 1); // up (high) price
        Real s1d = lattice->underlying(1, 0); // down (low) price

        Real delta = (p1u - p1d) / (s1u - s1d);

        // Finally, rollback to t=0
        option.rollback(0.0);
        Real p0 = option.presentValue();

        // Store results
        results_.value = p0;
        results_.delta = delta;
        results_.gamma = gamma;
        results_.theta = blackScholesTheta(process_,
                                           results_.value,
                                           results_.delta,
                                           results_.gamma);
    }

}


#endif
