// ****************************************************************************
/// @file cg.hpp
/// @author Kyle Webster
/// @version 0.2
/// @date 22 Nov 2025
/// @brief Numerics Library - Computer Graphics - @ref cg
/// @details
/// Collection of computer graphics structures and algorithms
// ****************************************************************************
#ifndef CG_HPP
#define CG_HPP
// ** Includes ****************************************************************
#include <variant>
#include <vector>

#include "bra.hpp"
// ****************************************************************************

namespace nl 
{ // ** nl ****************************

/// @namespace cg
/// @brief Computer graphics data structures and algorithms
/// @details
/// Contains:
///  * @ref data_structures
namespace cg 
{ // ** nl::cg ************************

// ****************************************************************************
/// @name data structures

// ****************************************************************************
/// @name spectrums

template <bra::arithmetic T> struct rgb
{
  bra::ℝn<3,T> c;
};
// ** end of spectrums ********************************************************


// ****************************************************************************
/// @name spatial

struct vec 
{
  bra::ℝn<4,float> dir;
};
struct pnt
{
  bra::ℝn<4,float> pos;
};
struct ray 
{
  pnt   p;
  vec   u;
  float t;
};
struct transform 
{
  bra::ℝnxm<3,4,float> M;
  bra::ℝnxm<3,4,float> M_inv;
};

struct bvh {};
// ** end of spatial **********************************************************


// ****************************************************************************
/// @name objects

struct sphere 
{
  transform T;
  float     radius;
};
struct plane 
{
  transform T;
  float     length;
};
struct trimesh 
{
  transform             T;
  std::vector<float>    V;
  std::vector<uint32_t> F;
};
// ** end of objects **********************************************************


// ****************************************************************************
/// @name lights

struct ambientlight 
{
  float irradiance = 0.f;
};
struct pointlight 
{
  float radiant_intensity = 0.f;
  pnt   pos;
};
struct dirlight 
{
  float radiant_intensity = 0.f;
  vec   dir;
};
struct spherelight : sphere 
{
  float radiance;
};
// ** end of lights ***********************************************************


// ****************************************************************************
/// @name materials

struct lambertian 
{

};
struct anisotropic {};
// ** end of materials ********************************************************


using Light = std::variant<ambientlight, pointlight, dirlight, spherelight>;
using Object = std::variant<sphere, plane, trimesh>;
using Material = std::variant<lambertian, anisotropic>;

// ****************************************************************************


// ****************************************************************************
/// @name imaging

template <bra::arithmetic T> struct image 
{
  size_t width, height;
  std::vector<T> data;
};

struct texture {};
// ** end of imaging **********************************************************


// ****************************************************************************
/// @name rendering
struct camera 
{
  transform T;
  float fov;
};

struct scene {};
} // ** end of namespace cg ***********
} // ** end of namespace nl ***********

// ****************************************************************************
#endif // #ifndef CG_HPP