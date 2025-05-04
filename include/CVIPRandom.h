
#pragma once

#ifndef __CVIP_Random_H___
#define __CVIP_Random_H___

#ifdef USING_BOOST			// USING BOOST

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#else
#include <random>
#endif

class CVIPRandom
{
public:
                        CVIPRandom(){};
    virtual             ~CVIPRandom(){};

    virtual double      GetNewValue() = 0;

private:
};

class CVIPRandomUniform : public CVIPRandom
{
public:
                        CVIPRandomUniform(const double& in_min, const double& in_max, int seed = 0);
    virtual             ~CVIPRandomUniform();

    virtual double      GetNewValue() ;

private:
    double              m_min;
    double              m_max;
};

class CVIPRandomGauss : public CVIPRandom
{
public:
                        CVIPRandomGauss(const double& in_mean, const double& in_sigma, int seed = 0);
    virtual             ~CVIPRandomGauss();

    virtual double      GetNewValue() ;     // function cannot be const because the engine is changing everytime the generator is called
    double              GetMean() const { return m_mean; }
    double              GetSigma() const { return m_sigma; }

private:
    double              m_mean;
    double              m_sigma;

#ifdef USING_CLHEP
    int dummy;
#elif defined USING_BOOST
    int dummy;
#elif not defined USING_BOOST
    std::default_random_engine   m_engine;
    int                 m_seed;
#endif



/*
#if defined USING_BOOST
  	boost::mt19937* 	m_rng;
#endif
*/
};

#endif
