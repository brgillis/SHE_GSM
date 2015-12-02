/**********************************************************************\
 @file common.h
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

#ifndef SHE_SIM_GAL_PARAMS_COMMON_HPP_
#define SHE_SIM_GAL_PARAMS_COMMON_HPP_

#include <Eigen/Core>

#include <memory>
#include <random>
#include <string>
#include <unordered_map>

// General typedefs

namespace SHE_SIM {

typedef short int short_int_t;
typedef int int_t;
typedef long int long_int_t;

typedef double flt_t;

typedef std::string str_t;
typedef str_t name_t;

typedef short_int_t level_t;
typedef std::unique_ptr<level_t> level_ptr_t;
typedef std::unordered_map<name_t,level_ptr_t> generation_level_map_t;

typedef std::ranlux48 gen_t;
typedef gen_t::result_type seed_t;

// Arrays

template< typename T> using array_t = Eigen::Array<T,Eigen::Dynamic,1>;
typedef array_t<flt_t> array_1d_t;
typedef Eigen::Array<flt_t,Eigen::Dynamic,Eigen::Dynamic> array_2d_t;
typedef array_t<array_1d_t> array_1d_array_t;
typedef array_t<array_2d_t> array_2d_array_t;

// Class forward-declarations
class ParamHierarchyLevel;
class ParamGenerator;
class ParamParam;

// Typedefs using these classes

typedef ParamGenerator param_t;
typedef std::unique_ptr<param_t> param_ptr_t;
typedef std::unordered_map<name_t,param_ptr_t> params_t;

typedef ParamParam param_param_t;
typedef std::unique_ptr<param_param_t> param_param_ptr_t;
typedef std::unordered_map<name_t,param_param_ptr_t> param_params_t;

// Typedefs for the type of each object
typedef array_2d_t background_psf_t;

typedef array_2d_array_t binned_observed_flux_distribution_t;
typedef array_2d_array_t binned_psf_t;
typedef array_2d_t core_observed_flux_distribution_t;
typedef array_1d_t core_sed_t;
typedef array_2d_t disk_observed_flux_distribution_t;
typedef array_1d_t disk_sed_t;
typedef array_2d_t observed_flux_distribution_t;
typedef array_2d_t psf_model_t;
typedef array_2d_t pix_galaxy_w_pois_noise_t;

}

#endif // SHE_SIM_GAL_PARAMS_COMMON_HPP_