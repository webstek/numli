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
#include <fstream>

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

template<typename T>
concept named = requires(T t) {{t.name}->std::convertible_to<std::string>;};
inline std::string getName(const named auto &x) { return x.name; }
template<typename... Ts>
inline std::string getName(const std::variant<Ts...> &v) {
  return std::visit([](auto const &x) { return x.name; }, v);
}

/// @brief named item for other structs to inherit
struct item {std::string name;};

/// @brief list class
/// @warning requires elements to have a public name member
template <typename T>
struct list : std::vector<T>
{
  /// @brief returns index of named element, returns size if not found 
  size_t idxOf(std::string const &name) const
  {
    size_t N = this->size();
    size_t i=0;
    for (; i<N; i++) { if (getName(this->at(i)) == name) break; }
    return i;
  }
  T* find(std::string const &name) 
  { 
    size_t idx = idxOf(name);
    if (idx!=this->size()) { return &this->at(idxOf(name)); }
    return nullptr;
  }
  T const* find(std::string const &name) const 
  {
    size_t idx=idxOf(name);
    if (idx!=this->size()) { return &this->at(idxOf(name)); }
    return nullptr;
  }
};


// ****************************************************************************
/// @name spectrums

template <bra::arithmetic T> struct rgb
{
  bra::ℝn<3,T> c;
  constexpr rgb() {}
  constexpr rgb(T v) {c[0]=v; c[1]=v; c[2]=v;}
  constexpr rgb& operator=(T v) {c[0]=v; c[1]=v; c[2]=v; return *this;}
  constexpr std::string toString()
    {return std::to_string(c[0])+std::to_string(c[1])+std::to_string(c[2]);}
};
using sRGB   = rgb<uint8_t>;
using linRGB = rgb<float>;
// ** end of spectrums ********************************************************


// ****************************************************************************
/// @name spatial

/// @brief computes the inverse transformation matrix of M
constexpr ℝ4x4 inverseTransform(ℝ4x4 const &M)
{
  ℝ3x3 R = {M(0,0),M(0,1),M(0,2),M(1,0),M(1,1),M(1,2),M(2,0),M(2,1),M(2,2)};
  ℝ3x3 R_inv = bra::inverse(R);
  ℝ4x4 M_inv;
  ℝ3 t = -R_inv*bra::column<3>(M,3);
  for (int i=0;i<3;i++) 
  {
    for (int j=0;j<3;j++) { M_inv(i,j)=R_inv(i,j); }
    M_inv(i,3) = t[i];
  }
  M_inv(3,0)=0.f; M_inv(3,1)=0.f; M_inv(3,2)=0.f; M_inv(3,3)=1.f;
  return ℝ4x4(M_inv.elem);
}

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
  ℝ4x4 M;
  ℝ4x4 M_inv;

  /// @name constructors
  constexpr transform() { M.identity(); M_inv.identity(); }
  constexpr transform(ℝ3 const &x, ℝ3 const &y, ℝ3 const &z, ℝ3 const &p)
  {
    for (int i=0;i<3;i++) {M(i,0)=x[i]; M(i,1)=y[i]; M(i,2)=z[i]; M(i,3)=p[i];}
    M(3,0)=0.f; M(3,1)=0.f; M(3,2)=0.f; M(3,3)=1.f;
    M_inv = inverseTransform(M);
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
  transform T;
  T.M = T2.M*T1.M;
  T.M_inv = T1.M_inv*T2.M_inv;
  return T;
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
  bvh<triangle>         _bvh;
  std::vector<ℝ3>       V;
  std::vector<uint32_t> F;
};
// ** end of structures ***************

// ************************************
/// @name object instances
enum class ObjectType {GROUP, SPHERE, PLANE, TRIMESH};
constexpr ObjectType str2obj(std::string s)
{
  if (s=="group")   return ObjectType::GROUP;
  if (s=="sphere")  return ObjectType::SPHERE;
  if (s=="plane")   return ObjectType::PLANE;
  if (s=="trimesh") return ObjectType::TRIMESH;
  throw std::domain_error("No corresponding ObjectType.");
}

struct sphere : item
{
  transform   T;
  materialidx mat;
};
struct plane : item
{
  transform T;
  materialidx mat;
};
struct trimesh : item
{
  transform   T;
  materialidx mat;
  meshidx     mesh;
};
// ** end of instances ****************
// ** end of objects **********************************************************


// ****************************************************************************
/// @name lights
enum class LightType {AMBIENT, POINT, DIR, SPHERE};
constexpr LightType str2light(std::string s)
{
  if (s=="ambient")   return LightType::AMBIENT;
  if (s=="point")     return LightType::POINT;
  if (s=="direction") return LightType::DIR;
  if (s=="sphere")    return LightType::SPHERE;
  throw std::domain_error("No corresponding LightType.");
}

struct ambientlight : item
{
  linRGB irradiance;
};
struct pointlight : item
{
  linRGB radiant_intensity;
  pnt   pos;
};
struct dirlight : item
{
  linRGB radiant_intensity;
  vec   dir;
};
struct spherelight : item
{
  linRGB radiance;
  sphere _sphere;
};
// ** end of lights ***********************************************************


// ****************************************************************************
/// @name materials
enum class MaterialType {NONE, LAMBERTIAN, BLINN, MICROFACET};
constexpr MaterialType str2mat(std::string const &s)
{
  if (s=="lambertian") return MaterialType::LAMBERTIAN;
  if (s=="blinn")      return MaterialType::BLINN;
  if (s=="microfacet") return MaterialType::MICROFACET;
  return MaterialType::NONE;
}

struct lambertian : item
{
  linRGB albedo;
};
struct blinn : item
{
  linRGB Kd, Ks, Kt, Le, reflect, transmit;
  float α, ior;
};
struct microfacet : item
{

};
struct emitter : item
{
  linRGB radiance;
};
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
using Material = std::variant<lambertian, blinn, microfacet, emitter>;
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

struct scene 
{
  // main interface components
  camera      cam;
  bvh<Object> obvh;

  // shared storage
  list<Object>   objects;
  list<Material> materials;
  list<Light>    ideal_lights;
  list<Mesh>     meshes;
  list<Texture>  textures;
  std::vector<Light*> lights;
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

// ************************************
/// @name data loading
template <bra::arithmetic T> inline void load(T &x, json const &j) 
  {x=j.get<T>();}
inline void loadℝ3(ℝ3 &x, json const &j) 
  {for(int i=0;i<3;i++)x[i]=j[i].get<float>();}
inline void loadTransform(transform &T, json const &j)
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
  T = scaling << rotation << translation;
}
// ************************************

// ************************************
/// @brief loads a camera
inline void loadCamera(camera &cam, json const &j) 
{
  ℝ3 pos, look_at, up;
  float fov, ar;
  uint width;
  loadℝ3(pos, j.at("pos"));
  loadℝ3(look_at, j.at("look_at"));
  loadℝ3(up, j.at("up"));
  load(fov, j.at("fov"));
  load(width, j.at("width"));
  try { load(ar, j.at("ar")); } catch(...) { ar=1.7778f; } // default 16:9
  cam.fov = fov;
  cam.width = width;
  cam.height = std::ceil(width/ar);
  ℝ3 z = (pos-look_at).normalized();
  ℝ3 x = up^z;
  ℝ3 y = z^x;
  cam.T = transform(x,y,z,pos);
}
// ************************************

// ************************************
/// @name light loading

inline void loadAmbientLight(ambientlight &l, json const &j) 
  {linRGB irrad; loadℝ3(irrad.c, j.at("irradiance")); l.irradiance = irrad;}
/// @todo
inline void loadPointLight(pointlight &l, json const &j);
/// @todo
inline void loadDirLight(dirlight &l, json const &j);
inline void loadSphereLight(
  spherelight &l, 
  json const &j, 
  list<Material> &mats)
{
  sphere s;
  linRGB radiance;
  loadℝ3(radiance.c, j.at("radiance"));
  loadTransform(s.T, j.at("transform"));
  std::string name = "emitter_"+radiance.toString();
  materialidx mat = mats.idxOf(name);
  if (mat==mats.size()) {emitter m = {name, radiance}; mats.push_back(m);}
  s.mat = mat;
  l.radiance = radiance;
  l._sphere = s;
}
/// @todo point and direction light loading
inline void loadLights(
  scene &scene,
  json const &j) 
{
  size_t n_lights = j.size();
  for (size_t i=0; i<n_lights; i++)
  {
    auto j_light = j[i];
    auto name = j_light.at("name").get<std::string>();
    auto type = str2light(j_light.at("type").get<std::string>());
    switch (type)
    {
    case LightType::AMBIENT:
      {
        ambientlight l; 
        l.name=name; 
        loadAmbientLight(l, j_light); 
        scene.ideal_lights.push_back(l);
        scene.lights.push_back(&scene.ideal_lights.back());
        break;
      }
    case LightType::POINT: break;
    case LightType::DIR: break;
    case LightType::SPHERE:
    {
      spherelight l;
      l.name=name;
      loadSphereLight(l, j_light, scene.materials);
      scene.objects.push_back(l._sphere);
      break;
    }
    }
  }
}
// ************************************

// ************************************
/// @name object loading

/// @todo
inline void loadGroup() {}
inline void loadSphere(sphere &s, json const &j, list<Material> const &mats)
{
  transform T;
  loadTransform(T, j.at("transform"));
  s.T = T;
  s.mat = mats.idxOf(j.at("material").get<std::string>());
}
inline void loadPlane(plane &p, json const &j, list<Material> const &mats) 
{
  transform T;
  loadTransform(T, j.at("transform"));
  p.T = T;
  p.mat = mats.idxOf(j.at("material").get<std::string>());
}
/// @todo
inline void loadTriMesh(trimesh &m, json const &j, list<Material> const &mats);
inline void loadObjects(
  list<Object> &objs, 
  json const &j, 
  list<Material> const &mats)
{
  size_t n_objs = j.size();
  for (size_t i=0; i<n_objs;i++)
  {
    auto j_obj = j[i];
    auto name = j_obj.at("name").get<std::string>();
    auto type = str2obj(j_obj.at("type").get<std::string>());
    switch (type)
    {
    case ObjectType::GROUP:
      /// @bug does not apply transformation of group to children
      {loadObjects(objs, j_obj, mats); break;}
    case ObjectType::SPHERE: 
    {
      sphere s; 
      s.name=name; 
      loadSphere(s, j_obj, mats); 
      objs.push_back(s); 
      break;
    }
    case ObjectType::PLANE:
    {
      plane p; 
      p.name=name; 
      loadPlane(p, j_obj, mats); 
      objs.push_back(p); 
      break;
    }
    case ObjectType::TRIMESH: break;
    }
  }
}
// ************************************

// ************************************
/// @name material loading

inline void loadLambertian(lambertian &m, json const &j)
{
  linRGB albedo;
  loadℝ3(albedo.c, j.at("albedo"));
  m.albedo = albedo;
}
inline void loadBlinn(blinn &m, json const &j)
{
  linRGB Kd, Ks, Kt, Le, reflect, transmit;
  float alpha, ior;
  try {loadℝ3(Kd.c, j.at("Kd"));} catch(...) {Kd.c = 0.f;}
  try {loadℝ3(Ks.c, j.at("Ks"));} catch(...) {Ks.c = 0.f;}
  try {loadℝ3(Kt.c, j.at("Kt"));} catch(...) {Kt.c = 0.f;}
  try {loadℝ3(Le.c, j.at("Le"));} catch(...) {Le.c = 0.f;}
  try {loadℝ3(reflect.c, j.at("reflect"));} catch(...) {reflect.c = 0.f;}
  try {loadℝ3(transmit.c, j.at("transmit"));} catch(...) {transmit.c = 0.f;}
  try {load(alpha, j.at("glossiness"));} catch(...) {alpha = 1024;}
  try {load(ior, j.at("ior"));} catch(...) {ior=1.54;}
  m.Kd = Kd;
  m.Ks = Ks;
  m.Kt = Kt;
  m.Le = Le;
  m.reflect = reflect;
  m.transmit = transmit;
  m.α = alpha;
  m.ior = ior;
}
/// @todo
inline void loadMicrofacet(microfacet &m, json const &j) {(void)m; (void)j;}

inline void loadMaterial(Material &mat, json const &j)
{
  std::visit(
    Overload {
      [&](lambertian &m){ loadLambertian(m, j); },
      [&](blinn &m){ loadBlinn(m, j); },
      [&](microfacet &m){ loadMicrofacet(m, j); },
      [](auto &m){(void)m;}}, // no behaviour fallback
    mat);
}
inline void loadMaterials(list<Material> &mats, json const &j)
{
  size_t n_mats = j.size();
  for (size_t i=0; i<n_mats; i++)
  {
    auto j_mat = j[i];
    auto name = j_mat.at("name").get<std::string>();
    auto type = str2mat(j_mat.at("type").get<std::string>());
    switch (type)
    {
    case MaterialType::NONE: break;
    case MaterialType::LAMBERTIAN:
    {
      lambertian l; 
      l.name=name; 
      loadLambertian(l,j_mat); 
      mats.push_back(l); 
      break;
    }
    case MaterialType::BLINN:
      {blinn b; b.name=name; loadBlinn(b,j_mat); mats.push_back(b); break;}
    case MaterialType::MICROFACET:
    {
      microfacet m; 
      m.name=name; 
      loadMicrofacet(m,j_mat); 
      mats.push_back(m); 
      break;
    }
    } // end of switch
  }
}
// ************************************

/// @todo load gltf node 
inline void loadGLTFNode(scene &scene, std::string fpath);

/// @brief load scene from file
/// @param scene scene to put data in
/// @param fpath path to file to load
/// @return true on successful loading, false otherwise
/// @todo Texture support
inline bool loadNLS(scene &scene, std::string fpath)
{
  std::ifstream file(fpath);
  if (!file.is_open()) 
    {throw std::runtime_error("Could no open file."); return false;}
  json j;
  file >> j;

  loadCamera(scene.cam, j.at("camera"));
  loadMaterials(scene.materials, j.at("materials"));
  loadObjects(scene.objects, j.at("objects"), scene.materials);
  loadLights(scene, j.at("lights"));
  return true;
}



} // ** end of namespace load *************************************************

} // ** end of namespace cg ***********
} // ** end of namespace nl ***********

// ****************************************************************************
#endif // #ifndef CG_HPP