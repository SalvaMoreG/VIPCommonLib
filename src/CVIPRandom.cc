
#include "CVIPRandom.h"

#ifdef USING_CLHEP			// USING CLHEP (which is terrible, should be avoided)
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#else
#ifdef USING_BOOST			// USING BOOST

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

boost::minstd_rand intgen;
boost::uniform_01<boost::minstd_rand> boost_dist_flat(intgen);

boost::mt19937 rng; // Forget about seeds, it's impossible to do in a global way...
boost::normal_distribution<> nd(0.0, 1.0);
boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > boost_dist_gauss(rng, nd);

#else						// USING C++11 (new C++), the easiest way...
#include <random>
#endif
#endif

#include <iostream>
using namespace std;

CVIPRandomUniform::CVIPRandomUniform(const double& in_min, const double& in_max, int seed)
    : m_min(in_min)
    , m_max(in_max)
{

//	#ifdef USING_CLHEP
//	    cout << "Using CLHEP Gauss uniform generator" << endl;
//	#elif defined USING_BOOST
//	    cout << "Using BOOST Gauss uniform generator" << endl;
//	#elif not defined USING_BOOST
//	    cout << "Using C++11 uniform generator" << endl;
//	#endif

}

CVIPRandomUniform::~CVIPRandomUniform()
{
}

double
CVIPRandomUniform::GetNewValue()
{
    double X(0.0);
#ifdef USING_CLHEP
    X = CLHEP::RandFlat::shoot(1.0);
    X = m_min + (m_max - m_min) * X;
#elif defined USING_BOOST
    X = boost_dist_flat();
    X = m_min + (m_max - m_min) * X;
#elif not defined USING_BOOST
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(m_min, m_max);
    X = dist(mt);
#endif
    return X;
}

// ==========================================================================

CVIPRandomGauss::CVIPRandomGauss( const double& in_mean, const double& in_sigma, int seed )
    : m_mean(in_mean)
    , m_sigma(in_sigma)
{
/*
#ifdef USING_CLHEP
    cout << "Using CLHEP Gauss random generator" << endl;
#elif defined USING_BOOST
    cout << "Using BOOST Gauss random generator" << endl;
#elif not defined USING_BOOST
    cout << "Using C++11 random generator" << endl;
#endif
*/
#ifdef USING_CLHEP
#elif defined USING_BOOST
#elif not defined USING_BOOST
    m_seed = seed;
    if (m_seed > 0)
        m_engine.seed(m_seed);
#endif
}

CVIPRandomGauss::~CVIPRandomGauss()
{
}

double
CVIPRandomGauss::GetNewValue()
{
    double X(0.0);
#ifdef USING_CLHEP
    X = CLHEP::RandGauss::shoot(m_mean, m_sigma);
#elif defined USING_BOOST
    X = m_mean + (m_sigma * boost_dist_gauss());
#elif not defined USING_BOOST
    if (m_seed > 0)
    {
        std::normal_distribution<double> distribution(m_mean, m_sigma);
        X = distribution( m_engine );
    }
    else
    {
        std::random_device rd;
        std::mt19937 mt(rd());
        std::normal_distribution<double> distribution(m_mean, m_sigma);
        X = distribution( mt );
    }
#endif

    return X;
}


