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
#include <string>

#include "json.hpp"

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

// ************************************
/// @name aliases
using materialidx = uint32_t;
using objectidx   = uint32_t;
using meshidx     = uint32_t;

// ****************************************************************************
/// @name data structures

// ****************************************************************************
/// @name spectrums

template <bra::arithmetic T> struct rgb
{
  bra::ℝn<3,T> c;
};
using sRGB   = rgb<uint8_t>;
using linRGB = rgb<float>;
// ** end of spectrums ********************************************************


// ****************************************************************************
/// @name spatial

struct vec 
{
  ℝ4 dir;
  constexpr vec() {dir[3]=0.f;}
  constexpr vec(float (&x)[3]) {for(size_t i=0;i<3;i++)dir[i]=x[i];dir[3]=0.f;}
  constexpr vec(vec const &v)  {for(size_t i=0;i<3;i++)dir[i]=v.dir[i];}
};
struct pnt
{
  ℝ4 pos;
  constexpr pnt() {pos[3]=1.f;}
  constexpr pnt(float (&x)[3]) {for(size_t i=0;i<3;i++)pos[i]=x[i];pos[3]=1.f;}
  constexpr pnt(pnt const &v)  {for(size_t i=0;i<3;i++)pos[i]=v.pos[i];}
};
struct ray 
{
  pnt   p;
  vec   u;
  float t;
};
struct basis
{
  ℝ3 x, y, z;
};
struct transform
{
  bra::ℝnxm<3,4,float> M;
  bra::ℝnxm<3,4,float> M_inv;
};
// ** end of spatial **********************************************************


// ****************************************************************************
/// @name infos

struct hitinfo {};
struct sampleinfo {};
// ****************************************************************************


// ****************************************************************************
/// @name objects

// ************************************
/// @name structures

struct aabb
{
  ℝ3 inf, sup;
};
struct triangle
{
  ℝ3 x0, x1, x2;
};
struct vertex 
{
  vec n;

};
template <typename T> struct bvh {};
struct trimeshdata
{
  bvh<triangle>         bvh;
  std::vector<ℝ3>       V;
  std::vector<uint32_t> F;
};
// ** end of structures ***************

// ************************************
/// @name object instances

struct sphere 
{
  transform   T;
  materialidx mat;
};
struct plane 
{
  transform T;
  materialidx mat;
};
struct trimesh
{
  transform   T;
  materialidx mat;
  meshidx     mesh;
};
// ** end of instances ****************
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
  linRGB Kd;
};
struct anisotropic {};
// ** end of materials ********************************************************


// ****************************************************************************
/// @name imaging

template <bra::arithmetic T> struct image 
{
  size_t width, height;
  std::vector<T> data;
};

template <typename T>
struct texture {};
using valuetex  = texture<float>;
using colourtex = texture<linRGB>;
// ** end of imaging **********************************************************


// ****************************************************************************
/// @name variants

using Light    = std::variant<ambientlight, pointlight, dirlight, spherelight>;
using Object   = std::variant<sphere, plane, trimesh>;
using Mesh     = std::variant<trimeshdata>;
using Material = std::variant<lambertian, anisotropic>;
using Texture  = std::variant<valuetex, colourtex>;
// ****************************************************************************


// ****************************************************************************
/// @name rendering
/// @brief rendering related structures

struct camera 
{
  transform T;
  float fov;
};

struct node
{
  std::string       name;
  std::vector<node> children;
  materialidx       mat;
  objectidx         obj;
};

struct scene 
{
  // main interface components
  camera      cam;
  node        root;
  bvh<Object> obvh;

  // shared storage
  std::vector<Object>   objects;
  std::vector<Mesh>     meshes;
  std::vector<Material> materials;
  std::vector<Light*>   lights;
  std::vector<Texture>  textures;
};
// ** end of data structures **************************************************

/// @namespace intersect
/// @brief intersection code for all objects
namespace intersect
{ // ** nl::cg::intersect *****************************************************

} // ** end of namespace intersect ********************************************


/// @namespace load
/// @brief loading utilities
namespace load
{ // nl::cg::load *************************************************************

/// @todo load .nls file
bool loadNLS(scene &scene, std::string fpath)
{
  
}


/// @todo load gltf node 
bool loadGLTFNode(node *node, scene &scene, std::string fpath) {return false;}
} // ** end of namespace load *************************************************

} // ** end of namespace cg ***********
} // ** end of namespace nl ***********

// ****************************************************************************
#endif // #ifndef CG_HPP