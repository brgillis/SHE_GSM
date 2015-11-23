/**********************************************************************\
 @file IndPoisson.hpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015 brg

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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDPOISSON_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDPOISSON_HPP_

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/random_functions.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndPoisson: public ParamParam
{
private:

	flt_t _lambda;

public:

	// Constructor and destructor
	IndPoisson( flt_t const & lambda )
	: ParamParam(ParamParam::INDEPENDENT),
	  _lambda(lambda)

	{
	}
	virtual ~IndPoisson() {}

	// Get the name of this
	virtual name_t name() const override { return "poisson"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return Pois_rand(_lambda, gen);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDPOISSON_HPP_