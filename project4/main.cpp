#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include "binomialengine.hpp"
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/methods/lattices/binomialtree.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>

#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

using namespace QuantLib;


#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    Integer sessionId() { return 0; }

}
#endif


int main(int, char* []) {

    try {

        boost::timer timer;
        std::cout << std::endl;

        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(11, January, 2019);
        Date settlementDate(17, January, 2019);
        Settings::instance().evaluationDate() = todaysDate;

        // our options
        Option::Type type(Option::Call);
        Real underlying = 100;
        Real strike = 110;
        Spread dividendYield = 0.00;
        Rate riskFreeRate = 0.03;
        Volatility volatility = 0.20;
        Date maturity(17, February, 2019);
        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Option type = "  << type << std::endl;
        std::cout << "Maturity = "        << maturity << std::endl;
        std::cout << "Underlying price = "        << underlying << std::endl;
        std::cout << "Strike = "                  << strike << std::endl;
        std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
                  << std::endl;
        std::cout << "Dividend yield = " << io::rate(dividendYield)
                  << std::endl;
        std::cout << "Volatility = " << io::volatility(volatility)
                  << std::endl;
        std::cout << std::endl;
        std::string method;
        std::cout << std::endl ;

        // write column headings
        boost::shared_ptr<Exercise> americanExercise(
                                         new AmericanExercise(maturity));


        Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(underlying)));


        // bootstrap the yield/dividend/vol curves
        Handle<YieldTermStructure> flatTermStructure(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
        Handle<YieldTermStructure> flatDividendTS(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter)));
        Handle<BlackVolTermStructure> flatVolTS(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(settlementDate, calendar, volatility,
                                     dayCounter)));
        boost::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));

        boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                new BlackScholesMertonProcess(underlyingH, flatDividendTS,
                                              flatTermStructure, flatVolTS));

        // options
        VanillaOption americanOption(payoff, americanExercise);
        std::ofstream data_result; 
        data_result.open ("data.txt");
        int i;                 
		int i_min = 801 ; 
		int i_max = 1000 ; 	
		//	method = "Binomial Cox-Ross-Rubinstein";
		for ( i = i_min ; i<i_max +1  ; i++) {
			  Size timeSteps (i) ;
			  americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
								new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess, timeSteps)));
			  Real NPV1 = americanOption.NPV();
			  americanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
								new BinomialVanillaEngine_2<CoxRossRubinstein>(bsmProcess, timeSteps)));
		      Real NPV2 = americanOption.NPV();
			  data_result << i << " " << NPV1 << " " << NPV2 << "\n ";

		}
		data_result.close() ; 

		double seconds = timer.elapsed();
        Integer hours = int(seconds/3600);
        seconds -= hours * 3600;
        Integer minutes = int(seconds/60);
        seconds -= minutes * 60;
        std::cout << " \nRun completed in ";
        if (hours > 0)
            std::cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            std::cout << minutes << " m ";
        std::cout << std::fixed << std::setprecision(0)
                  << seconds << " s\n" << std::endl;
        return 0;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}
