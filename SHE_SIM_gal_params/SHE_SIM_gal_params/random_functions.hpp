/**********************************************************************\
  @file random_functions.hpp

 **********************************************************************

 Copyright (C) 2014, 2015  Bryan R. Gillis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

\**********************************************************************/

// body file: random_functions.cpp

#ifndef _BRG_RANDOM_FUNCTIONS_HPP_INCLUDED_
#define _BRG_RANDOM_FUNCTIONS_HPP_INCLUDED_

#include <SHE_SIM_gal_params/common.hpp>
#include <cassert>
#include <cmath>
#include <random>

#include "SHE_SIM_gal_params/math.hpp"

namespace SHE_SIM {

extern gen_t rng; // Initialised in random_functions.cpp

/** Global function declarations **/
#if (1)

// Generates a random int_t between min and max, inclusive of min, exclusive of max
template< typename T=int_t, typename T_in=int_t >
T irand( T_in && min, T_in && max, gen_t & gen=rng )
{
	return std::uniform_int_distribution<T>(std::forward<T_in>(min),std::forward<T_in>(max)-1)(gen);
} // flt_t drand(flt_t min, flt_t max)

// Generates a random flt_t between min and max
template< typename T=flt_t >
T drand( gen_t & gen=rng )
{
	return std::uniform_real_distribution<T>()(gen);
}
template< typename T=flt_t, typename T_in=flt_t >
T drand( T_in && min, T_in && max, gen_t & gen=rng )
{
	return std::uniform_real_distribution<T>(std::forward<T_in>(min),std::forward<T_in>(max))(gen);
} // flt_t drand(flt_t min, flt_t max)

// Returns a random variable from a Gaussian distribution
template< typename T=flt_t >
T Gaus_rand( gen_t & gen=rng )
{
	return std::normal_distribution<T>()(gen);
} // flt_t Gaus_rand()
template< typename T=flt_t, typename T_in=flt_t >
T Gaus_rand( T_in && mean, T_in && stddev = 1.0, gen_t & gen=rng )
{
	return std::normal_distribution<T>(std::forward<T_in>(mean),std::forward<T_in>(stddev))(gen);

} // flt_t Gaus_rand(flt_t mean, flt_t stddev)

// Returns a random variable from a Gaussian distribution, truncated between min and max
template< typename T=flt_t, typename T_in=flt_t >
T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )
{
	assert(max>min);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		flt_t test_result = Gaus_rand(mean,stddev,gen);
		if((test_result >= min)and(test_result <= max))
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return (min+max)/2.;

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Returns a random variable from a Gaussian distribution in log space
// Note that "mean" here is the desired mean, NOT the peak of the function (which differ in log space). If you want to use
// the peak, simply use the standard Gaus_rand version instead.
template< typename T=flt_t >
T log10Gaus_rand( gen_t & gen=rng )
{
	const flt_t & fact = std::exp( -square( std::log( 10. ) ) / 2 );

	return ( fact * std::pow(10., Gaus_rand<T>(gen) ) );
} // flt_t log10Gaus_rand()

// Returns a random variable from a Gaussian distribution, truncated between min and max
template< typename T=flt_t, typename T_in=flt_t >
T trunc_log10Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )
{
	assert(max>min);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		flt_t test_result = Gaus_rand(mean,stddev,gen);
		if((test_result >= min)and(test_result <= max))
		{
			return std::pow(10.,test_result);
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return std::pow(10.,(min+max)/2.);

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Returns a random variable from a Rayleigh distribution
template< typename T=flt_t >
T Rayleigh_rand( gen_t & gen=rng )
{
	return std::sqrt(-2.*std::log(drand<T>(gen)));
}
template< typename T=flt_t, typename T_in=flt_t >
T Rayleigh_rand( T_in && sigma, gen_t & gen=rng )
{
	return std::forward<T_in>(sigma)*Rayleigh_rand(gen);
}

// Returns a random variable from a Gaussian distribution, truncated between min and max
template< typename T=flt_t, typename T_in=flt_t >
T trunc_Rayleigh_rand( T_in && sigma, T_in && max, gen_t & gen=rng )
{
	assert(max>0.);

	// Try values until one works
	bool good_value = false;
	int attempt_counter = 0;

	while( (!good_value) and (attempt_counter < 1000) )
	{
		flt_t test_result = Rayleigh_rand(sigma,gen);
		if(test_result <= max)
		{
			return test_result;
		}
		else
		{
			++attempt_counter;
		}
	}

	// Failsafe
	return max/2.;

} // T trunc_Gaus_rand( T_in && mean, T_in && stddev, T_in && min, T_in && max, gen_t & gen=rng )

// Get a Rayleigh random variable, smoothly contracted to always be less than max
template< typename T=flt_t, typename T_in=flt_t >
T contracted_Rayleigh_rand( T_in && sigma, T_in && max, T_in && p, gen_t & gen=rng )
{
	// Generate an initial random Rayleigh variable
	T first_result = Rayleigh_rand(sigma);

	// Transform it via Bryan's formula to rein in large values to be less than the max_val
	return (first_result / std::pow(1 + std::pow(first_result / max, p), 1.0 / p));
}

// Returns a Poisson random variable.
template< typename T=int_t >
T Pois_rand( gen_t & gen=rng )
{
	return std::poisson_distribution<T>()(gen);
} // T Pois_rand( gen_t & gen=rng )
template< typename T=int_t, typename T_in=flt_t >
T Pois_rand( T_in && lambda=1., gen_t & gen=rng )
{
	return std::poisson_distribution<T>(std::forward<T_in>(lambda))(gen);
} // T Pois_rand( T_in && lambda=1., gen_t & gen=rng )

#endif // End global function declarations

} // namespace SHE_SIM

#endif
