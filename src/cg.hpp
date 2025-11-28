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
  ℝ3x4 M;
  ℝ3x4 M_inv;

  /// @name constructors
  constexpr transform() { M.identity(); M_inv.identity(); }
  constexpr transform(ℝ3 const &x, ℝ3 const &y, ℝ3 const &z, ℝ3 const &p)
  {
    for (int i=0;i<3;i++) {M(i,0)=x[i]; M(i,1)=y[i]; M(i,2)=z[i]; M(i,3)=p[i];}
    /// @todo compute inverse
  }

  /// @name member functions
  static constexpr transform scale(ℝ3 const &x)
  {
    transform T;
    for (int i=0;i<3;i++) { T.M(i,i)=x[i]; T.M_inv(i,i)=1.f/x[i]; }
    return T;
  }
  static constexpr transform rotate(ℝ3 const &u, float degrees)
  { // using Rodrigues' formula R(u,θ) = cosθI+(1-cosθ)*outer(u)+sinθ*skew(u)
    transform T;
    const float θ = radians(degrees);
    const float cosθ = cosf32(θ);
    const float sinθ = sinf32(θ);
    T.M(0,0) = cosθ + (1-cosθ)*u[0]*u[0];
    T.M(0,1) = (1-cosθ)*u[0]*u[1]-sinθ*u[2];
    T.M(0,2) = (1-cosθ)*u[0]*u[2]+sinθ*u[1];
    T.M(1,0) = (1-cosθ)*u[1]*u[0]+sinθ*u[2];
    T.M(1,1) = cosθ + (1-cosθ)*u[1]*u[1];
    T.M(1,2) = (1-cosθ)*u[1]*u[2]-sinθ*u[0];
    T.M(2,0) = (1-cosθ)*u[2]*u[0]-sinθ*u[1];
    T.M(2,1) = (1-cosθ)*u[2]*u[1]+sinθ*u[0];
    T.M(2,2) = cosθ + (1-cosθ)*u[2]*u[2];

    // replace θ with -θ for inverse
    T.M_inv(0,0) = cosθ + (1-cosθ)*u[0]*u[0];
    T.M_inv(0,1) = (1-cosθ)*u[0]*u[1]+sinθ*u[2];
    T.M_inv(0,2) = (1-cosθ)*u[0]*u[2]-sinθ*u[1];
    T.M_inv(1,0) = (1-cosθ)*u[1]*u[0]-sinθ*u[2];
    T.M_inv(1,1) = cosθ + (1-cosθ)*u[1]*u[1];
    T.M_inv(1,2) = (1-cosθ)*u[1]*u[2]+sinθ*u[0];
    T.M_inv(2,0) = (1-cosθ)*u[2]*u[0]+sinθ*u[1];
    T.M_inv(2,1) = (1-cosθ)*u[2]*u[1]-sinθ*u[0];
    T.M_inv(2,2) = cosθ + (1-cosθ)*u[2]*u[2];
    return T;
  }
  static constexpr transform translate(ℝ3 const &x) 
  { 
    transform T; 
    for (int i=0;i<3;i++) { T.M(i,3)=x[i]; T.M_inv(i,3)=-x[i]; } 
    return T;
  }
};

/// @brief compose transform T2 after T1
/// @todo composition operator
constexpr transform operator<<(transform const &T1, transform const &T2)
{

}
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
  uint width, height;
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
using json = nlohmann::json;

template <bra::arithmetic T> void load(T &x, json const &j) {x=j.get<T>();}
void loadℝ3(ℝ3 &x, json const &j) {for(int i=0;i<3;i++)x[i]=j[i].get<float>();}

void loadTransform(transform &T, json const &j)
{
  ℝ3 _scale, _translate, _axis;
  float _deg;
  try { loadℝ3(_scale, j.at("scale")); } catch(...) { _scale=1.f; }
  try { loadℝ3(_translate, j.at("translate")); } catch(...) { _translate=0.f; }
  try 
  { // try to load rotation
    json rot = j.at("rotate"); 
    loadℝ3(_axis, rot.at("axis")); 
    load(_deg, rot.at("degrees")); 
  } catch(...) { _axis={0.f,0.f,1.f}; _deg=0.f; }
  auto scaling = transform::scale(_scale);
  auto rotation = transform::rotate(_axis, _deg);
  auto translation = transform::translate(_translate);
  T = translation*rotation*scaling;
}

void loadCamera(camera &cam, json const &j) 
{
  ℝ3 pos, look_at, up;
  float fov, ar;
  uint width;
  loadℝ3(pos, j.at("pos"));
  loadℝ3(look_at, j.at("look_at"));
  loadℝ3(up, j.at("up"));
  load(fov, j.at("fov"));
  load(width, j.at("width"));
  try { load(ar, j.at("ar")); } catch(...) { ar=1.7778f; }
  cam.fov = fov;
  cam.width = width;
  cam.height = std::ceil(width/ar);
  ℝ3 z = (pos-look_at).normalized();
  ℝ3 x = up^z;
  ℝ3 y = z^x;
  cam.T = transform(x,y,z,pos);
}

void loadLight(Light &light, json const &j)
{

}

void loadObject(Object &obj, json const &j)
{

}

void loadMaterial(Material &mat, json const &j)
{

}

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