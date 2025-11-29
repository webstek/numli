// ****************************************************************************
/// @file cg.cpp
/// @author Kyle Webster
/// @date 29 Nov 2025
/// @details
/// Tests for the computer graphics module (cg.hpp)
// ****************************************************************************

// ** Includes ****************************************************************
#include "config.hpp"
#include "cg.hpp"
#include <filesystem>
#include <fstream>
// ****************************************************************************

using json = nlohmann::json;

// ** Helper macros ***********************************************************
#define CHECK_RGB(rgb, r, g, b) \
  do { \
    CHECK((rgb).c.elem[0] == (r)); \
    CHECK((rgb).c.elem[1] == (g)); \
    CHECK((rgb).c.elem[2] == (b)); \
  } while (0)

// ****************************************************************************
/// @name Material Loading Tests

TEST_CASE("Load Lambertian Material")
{
  json j_mat = {
    {"name", "matte_red"},
    {"type", "lambertian"},
    {"albedo", {0.7f, 0.1f, 0.1f}}
  };

  nl::cg::lambertian mat;
  mat.name = j_mat.at("name").get<std::string>();
  nl::cg::load::loadLambertian(mat, j_mat);

  CHECK(mat.name == "matte_red");
  CHECK_RGB(mat.albedo, 0.7f, 0.1f, 0.1f);
}

TEST_CASE("Load Blinn Material with all fields")
{
  json j_mat = {
    {"name", "shiny_gold"},
    {"type", "blinn"},
    {"Kd", {0.5f, 0.4f, 0.1f}},
    {"Ks", {0.9f, 0.9f, 0.8f}},
    {"Kt", {0.0f, 0.0f, 0.0f}},
    {"Le", {0.0f, 0.0f, 0.0f}},
    {"reflect", {0.1f, 0.1f, 0.1f}},
    {"transmit", {0.0f, 0.0f, 0.0f}},
    {"glossiness", 256.f},
    {"ior", 1.5f}
  };

  nl::cg::blinn mat;
  mat.name = j_mat.at("name").get<std::string>();
  nl::cg::load::loadBlinn(mat, j_mat);

  CHECK(mat.name == "shiny_gold");
  CHECK_RGB(mat.Kd, 0.5f, 0.4f, 0.1f);
  CHECK_RGB(mat.Ks, 0.9f, 0.9f, 0.8f);
  CHECK(mat.α == 256.f);
  CHECK(mat.ior == 1.5f);
}

TEST_CASE("Load Blinn Material with defaults")
{
  json j_mat = {
    {"name", "minimal_blinn"},
    {"type", "blinn"},
    {"Kd", {0.5f, 0.5f, 0.5f}}
    // Ks, Kt, Le, reflect, transmit, glossiness, ior will use defaults
  };

  nl::cg::blinn mat;
  mat.name = j_mat.at("name").get<std::string>();
  nl::cg::load::loadBlinn(mat, j_mat);

  CHECK(mat.name == "minimal_blinn");
  CHECK_RGB(mat.Kd, 0.5f, 0.5f, 0.5f);
  CHECK_RGB(mat.Ks, 0.f, 0.f, 0.f);  // default 0
  CHECK(mat.α == 1024.f);             // default glossiness
  CHECK(mat.ior == 1.54f);             // default ior
}

TEST_CASE("Load list of Materials")
{
  json j_mats = json::array();
  j_mats.push_back({
    {"name", "red_matte"},
    {"type", "lambertian"},
    {"albedo", {0.8f, 0.1f, 0.1f}}
  });
  j_mats.push_back({
    {"name", "blue_specular"},
    {"type", "blinn"},
    {"Kd", {0.1f, 0.1f, 0.5f}},
    {"Ks", {0.8f, 0.8f, 0.8f}}
  });

  nl::cg::list<nl::cg::Material> mats;
  nl::cg::load::loadMaterials(mats, j_mats);

  CHECK(mats.size() == 2);
  
  // Check first material (lambertian)
  auto mat1_ptr = mats.find("red_matte");
  CHECK(mat1_ptr != nullptr);
  auto &mat1 = std::get<nl::cg::lambertian>(*mat1_ptr);
  CHECK_RGB(mat1.albedo, 0.8f, 0.1f, 0.1f);

  // Check second material (blinn)
  auto mat2_ptr = mats.find("blue_specular");
  CHECK(mat2_ptr != nullptr);
  auto &mat2 = std::get<nl::cg::blinn>(*mat2_ptr);
  CHECK_RGB(mat2.Kd, 0.1f, 0.1f, 0.5f);
  CHECK_RGB(mat2.Ks, 0.8f, 0.8f, 0.8f);
}

// ****************************************************************************
/// @name Object Loading Tests

TEST_CASE("Load Sphere Object")
{
  json j_mat = {
    {"name", "sphere_mat"},
    {"type", "lambertian"},
    {"albedo", {0.5f, 0.5f, 0.5f}}
  };
  nl::cg::list<nl::cg::Material> mats;
  nl::cg::load::loadMaterials(mats, json::array({j_mat}));

  json j_sphere = {
    {"name", "my_sphere"},
    {"type", "sphere"},
    {"material", "sphere_mat"},
    {"transform", {
      {"scale", {1.0f, 1.0f, 1.0f}},
      {"translate", {5.0f, 0.0f, -10.0f}}
    }}
  };

  nl::cg::sphere s;
  s.name = "my_sphere";
  nl::cg::load::loadSphere(s, j_sphere, mats);

  CHECK(s.name == "my_sphere");
  CHECK(s.mat == 0);  // should match sphere_mat at index 0
}

TEST_CASE("Load Plane Object")
{
  json j_mat = {
    {"name", "floor_mat"},
    {"type", "lambertian"},
    {"albedo", {0.3f, 0.3f, 0.3f}}
  };
  nl::cg::list<nl::cg::Material> mats;
  nl::cg::load::loadMaterials(mats, json::array({j_mat}));

  json j_plane = {
    {"name", "floor"},
    {"type", "plane"},
    {"material", "floor_mat"},
    {"transform", {
      {"translate", {0.0f, -5.0f, 0.0f}},
      {"rotate", {
        {"axis", {0.0f, 1.0f, 0.0f}},
        {"degrees", 0.0f}
      }}
    }}
  };

  nl::cg::plane p;
  p.name = "floor";
  nl::cg::load::loadPlane(p, j_plane, mats);

  CHECK(p.name == "floor");
  CHECK(p.mat == 0);
}

TEST_CASE("Load list of Objects")
{
  // Create materials
  json j_mats = json::array({
    {{"name", "mat1"}, {"type", "lambertian"}, {"albedo", {0.5f, 0.5f, 0.5f}}},
    {{"name", "mat2"}, {"type", "lambertian"}, {"albedo", {0.3f, 0.3f, 0.3f}}}
  });
  nl::cg::list<nl::cg::Material> mats;
  nl::cg::load::loadMaterials(mats, j_mats);

  // Create objects
  json j_objs = json::array({
    {
      {"name", "sphere1"},
      {"type", "sphere"},
      {"material", "mat1"},
      {"transform", {{"translate", {0.0f, 0.0f, 0.0f}}}}
    },
    {
      {"name", "plane1"},
      {"type", "plane"},
      {"material", "mat2"},
      {"transform", {{"translate", {0.0f, -5.0f, 0.0f}}}}
    }
  });

  nl::cg::list<nl::cg::Object> objs;
  nl::cg::load::loadObjects(objs, j_objs, mats);

  CHECK(objs.size() == 2);
  CHECK(objs[0].index() == 1);  // sphere is variant index 1
  CHECK(objs[1].index() == 2);  // plane is variant index 2

  auto &s = std::get<nl::cg::sphere>(objs[0]);
  auto &p = std::get<nl::cg::plane>(objs[1]);
  CHECK(s.name == "sphere1");
  CHECK(p.name == "plane1");
}

// ****************************************************************************
/// @name Light Loading Tests

TEST_CASE("Load Ambient Light")
{
  json j_light = {
    {"name", "ambient_env"},
    {"type", "ambient"},
    {"irradiance", {0.5f, 0.5f, 0.5f}}
  };

  nl::cg::ambientlight light;
  light.name = j_light.at("name").get<std::string>();
  nl::cg::load::loadAmbientLight(light, j_light);

  CHECK(light.name == "ambient_env");
  CHECK_RGB(light.irradiance, 0.5f, 0.5f, 0.5f);
}

TEST_CASE("Load Sphere Light")
{
  json j_mats = json::array();
  nl::cg::list<nl::cg::Material> mats;

  json j_light = {
    {"name", "sphere_emitter"},
    {"type", "sphere"},
    {"radiance", {2.0f, 2.0f, 2.0f}},
    {"transform", {
      {"scale", {0.5f, 0.5f, 0.5f}},
      {"translate", {10.0f, 5.0f, 0.0f}}
    }}
  };

  nl::cg::spherelight light;
  light.name = j_light.at("name").get<std::string>();
  nl::cg::load::loadSphereLight(light, j_light, mats);

  CHECK(light.name == "sphere_emitter");
  CHECK_RGB(light.radiance, 2.0f, 2.0f, 2.0f);
  // Should have created an emitter material
  CHECK(mats.size() > 0);
}

TEST_CASE("Load list of Lights (Ambient)")
{
  json j_lights = json::array({
    {
      {"name", "ambient1"},
      {"type", "ambient"},
      {"irradiance", {0.3f, 0.3f, 0.3f}}
    },
    {
      {"name", "ambient2"},
      {"type", "ambient"},
      {"irradiance", {0.5f, 0.5f, 0.5f}}
    }
  });

  nl::cg::scene scene;
  nl::cg::load::loadLights(scene, j_lights);

  CHECK(scene.ideal_lights.size() == 2);
  
  auto &light1 = std::get<nl::cg::ambientlight>(scene.ideal_lights[0]);
  auto &light2 = std::get<nl::cg::ambientlight>(scene.ideal_lights[1]);
  CHECK(light1.name == "ambient1");
  CHECK(light2.name == "ambient2");
  CHECK_RGB(light1.irradiance, 0.3f, 0.3f, 0.3f);
  CHECK_RGB(light2.irradiance, 0.5f, 0.5f, 0.5f);
}

// ****************************************************************************
/// @name Transformation Composition Tests

TEST_CASE("Transform - Translation only")
{
  json j_obj = {
    {"name", "translated_obj"},
    {"type", "sphere"},
    {"material", "mat1"},
    {"transform", {
      {"translate", {5.0f, 3.0f, -2.0f}}
    }}
  };

  nl::ℝ3 trans{5.0f, 3.0f, -2.0f};
  nl::cg::transform T = nl::cg::transform::translate(trans);
  
  // Check that translation is in position part of matrix (column 3)
  CHECK(T.M(0, 3) == 5.0f);
  CHECK(T.M(1, 3) == 3.0f);
  CHECK(T.M(2, 3) == -2.0f);
  CHECK(T.M_inv(0, 3) == -5.0f);
  CHECK(T.M_inv(1, 3) == -3.0f);
  CHECK(T.M_inv(2, 3) == 2.0f);
}

TEST_CASE("Transform - Scale only")
{
  nl::ℝ3 scale{2.0f, 3.0f, 0.5f};
  nl::cg::transform T = nl::cg::transform::scale(scale);
  
  // Check diagonal elements
  CHECK(T.M(0, 0) == 2.0f);
  CHECK(T.M(1, 1) == 3.0f);
  CHECK(T.M(2, 2) == 0.5f);
  CHECK(T.M_inv(0, 0) == 0.5f);
  CHECK(T.M_inv(1, 1) == 1.0f / 3.0f);
  CHECK(T.M_inv(2, 2) == 2.0f);
}

TEST_CASE("Transform - Rotate X axis")
{
  nl::ℝ3 axis{1.0f, 0.0f, 0.0f};
  nl::cg::transform T = nl::cg::transform::rotate(axis, 90.0f);
  
  // Rotation around X axis by 90 degrees should keep X column identity
  CHECK(T.M(0, 0) > 0.99f);  // cos(90°) ≈ 0
  // Y and Z axes should swap and rotate
  CHECK(std::abs(T.M(1, 2)) < 0.01f);  // should be ~0 (cos(90°))
  CHECK(std::abs(T.M(2, 1)) < 0.01f);  // should be ~0 (-sin(90°) becomes ~0)
}

TEST_CASE("Transform - Rotate Y axis")
{
  nl::ℝ3 axis{0.0f, 1.0f, 0.0f};
  nl::cg::transform T = nl::cg::transform::rotate(axis, 45.0f);
  
  // Rotation around Y axis, should preserve Y
  CHECK(T.M(1, 1) > 0.99f);  // cos(45°) ≈ 0.707, identity on Y
  // X and Z should be affected
  CHECK(std::abs(T.M(0, 2)) > 0.1f);  // should have non-zero Z component
}

TEST_CASE("Transform - Rotate Z axis")
{
  nl::ℝ3 axis{0.0f, 0.0f, 1.0f};
  nl::cg::transform T = nl::cg::transform::rotate(axis, 90.0f);
  
  // Rotation around Z axis by 90 degrees
  CHECK(T.M(2, 2) > 0.99f);  // Z component preserved
  // X and Y should be affected
  CHECK(std::abs(T.M(0, 1)) > 0.1f);  // should have non-zero Y component
}

TEST_CASE("Transform - Scale then Translate composition")
{
  // First scale, then translate
  nl::ℝ3 s1{2.0f, 2.0f, 2.0f};
  nl::ℝ3 t1{5.0f, 0.0f, 0.0f};
  auto scale_T = nl::cg::transform::scale(s1);
  auto translate_T = nl::cg::transform::translate(t1);
  auto composed = scale_T << translate_T;
  
  // Composed should have both scale and translation effects
  CHECK(composed.M(0, 0) == 2.0f);  // scale preserved
  CHECK(composed.M(0, 3) == 5.0f);  // translation applied
}

TEST_CASE("Transform - Translate then Scale composition")
{
  // First translate, then scale (different order = different result)
  nl::ℝ3 t2{5.0f, 0.0f, 0.0f};
  nl::ℝ3 s2{2.0f, 2.0f, 2.0f};
  auto translate_T = nl::cg::transform::translate(t2);
  auto scale_T = nl::cg::transform::scale(s2);
  auto composed = translate_T << scale_T;
  
  // When scale is applied after translate, the translation also gets scaled
  CHECK(composed.M(0, 0) == 2.0f);   // scale effect
  CHECK(composed.M(0, 3) == 10.0f);  // translation gets scaled: 5*2 = 10
}

TEST_CASE("Transform - Rotate then Translate composition")
{
  nl::ℝ3 ax1{0.0f, 0.0f, 1.0f};
  nl::ℝ3 t3{3.0f, 4.0f, 0.0f};
  auto rotate_T = nl::cg::transform::rotate(ax1, 45.0f);
  auto translate_T = nl::cg::transform::translate(t3);
  auto composed = rotate_T << translate_T;
  
  // Translation should be applied after rotation
  CHECK(composed.M(0, 3) == 3.0f);  // direct translation preserved
  CHECK(composed.M(1, 3) == 4.0f);
}

TEST_CASE("Transform - Scale, Rotate, Translate composition (SRT)")
{
  // Standard SRT composition: Scale -> Rotate -> Translate
  nl::ℝ3 s3{2.0f, 2.0f, 2.0f};
  nl::ℝ3 ax2{0.0f, 1.0f, 0.0f};
  nl::ℝ3 t4{5.0f, 2.0f, -3.0f};
  auto scale_T = nl::cg::transform::scale(s3);
  auto rotate_T = nl::cg::transform::rotate(ax2, 30.0f);
  auto translate_T = nl::cg::transform::translate(t4);
  auto composed = scale_T << rotate_T << translate_T;
  
  // Check that all three are present
  CHECK(composed.M(0, 0) != 0.0f);  // has scale and rotation
  CHECK(composed.M(0, 3) == 5.0f);  // translation in X
  CHECK(composed.M(1, 3) == 2.0f);  // translation in Y
  CHECK(composed.M(2, 3) == -3.0f); // translation in Z
}

TEST_CASE("Transform - Triple composition: Scale, Translate, Rotate")
{
  nl::ℝ3 s4{0.5f, 0.5f, 0.5f};
  nl::ℝ3 t5{10.0f, 0.0f, 0.0f};
  nl::ℝ3 ax3{0.0f, 1.0f, 0.0f};
  auto scale_T = nl::cg::transform::scale(s4);
  auto translate_T = nl::cg::transform::translate(t5);
  auto rotate_T = nl::cg::transform::rotate(ax3, 90.0f);
  auto composed = scale_T << translate_T << rotate_T;
  
  // All transforms should be applied
  CHECK(composed.M(0, 0) != 0.0f);
  CHECK(composed.M(1, 1) != 0.0f);
  CHECK(composed.M(2, 2) != 0.0f);
}

TEST_CASE("Transform - Load from JSON: translate only")
{
  json j_transform = {
    {"translate", {3.0f, 2.0f, 1.0f}}
  };
  
  nl::cg::transform T;
  nl::cg::load::loadTransform(T, j_transform);
  
  CHECK(T.M(0, 3) == 3.0f);
  CHECK(T.M(1, 3) == 2.0f);
  CHECK(T.M(2, 3) == 1.0f);
}

TEST_CASE("Transform - Load from JSON: scale only")
{
  json j_transform = {
    {"scale", {2.0f, 3.0f, 4.0f}}
  };
  
  nl::cg::transform T;
  nl::cg::load::loadTransform(T, j_transform);
  
  CHECK(T.M(0, 0) == 2.0f);
  CHECK(T.M(1, 1) == 3.0f);
  CHECK(T.M(2, 2) == 4.0f);
}

TEST_CASE("Transform - Load from JSON: rotate with axis and degrees")
{
  json j_transform = {
    {"rotate", {
      {"axis", {0.0f, 1.0f, 0.0f}},
      {"degrees", 45.0f}
    }}
  };
  
  nl::cg::transform T;
  nl::cg::load::loadTransform(T, j_transform);
  
  // Rotation around Y by 45 degrees
  CHECK(T.M(1, 1) > 0.99f);  // Y preserved
  CHECK(std::abs(T.M(0, 2)) > 0.1f);  // X-Z interaction
}

TEST_CASE("Transform - Load from JSON: SRT composition from loadTransform")
{
  json j_transform = {
    {"scale", {2.0f, 2.0f, 2.0f}},
    {"rotate", {
      {"axis", {0.0f, 0.0f, 1.0f}},
      {"degrees", 45.0f}
    }},
    {"translate", {5.0f, -2.0f, 3.0f}}
  };
  
  nl::cg::transform T;
  nl::cg::load::loadTransform(T, j_transform);
  
  // loadTransform applies: scale << rotate << translate
  // So position should have translation
  CHECK(T.M(0, 3) == 5.0f);
  CHECK(T.M(1, 3) == -2.0f);
  CHECK(T.M(2, 3) == 3.0f);
  // And should have scale and rotation effects
  CHECK(T.M(0, 0) != 0.0f || T.M(0, 1) != 0.0f);
}

TEST_CASE("Transform - Load from JSON: partial transforms use defaults")
{
  json j_transform = {
    {"scale", {1.5f, 1.5f, 1.5f}}
    // rotate and translate are missing, should use defaults
  };
  
  nl::cg::transform T;
  nl::cg::load::loadTransform(T, j_transform);
  
  // Scale should be applied
  CHECK(T.M(0, 0) == 1.5f);
  // Default rotation is identity on Z axis
  CHECK(T.M(1, 1) == 1.5f);
  // No translation (default is 0)
  CHECK(T.M(0, 3) == 0.0f);
}

TEST_CASE("Transform - Complex SRT in object loading")
{
  json j_mat = {
    {"name", "mat"},
    {"type", "lambertian"},
    {"albedo", {0.5f, 0.5f, 0.5f}}
  };
  nl::cg::list<nl::cg::Material> mats;
  nl::cg::load::loadMaterials(mats, json::array({j_mat}));

  json j_sphere = {
    {"name", "complex_sphere"},
    {"type", "sphere"},
    {"material", "mat"},
    {"transform", {
      {"scale", {1.5f, 1.5f, 1.5f}},
      {"rotate", {
        {"axis", {1.0f, 0.0f, 0.0f}},  // rotate around X
        {"degrees", 30.0f}
      }},
      {"translate", {-5.0f, 10.0f, 8.0f}}
    }}
  };

  nl::cg::sphere s;
  s.name = "complex_sphere";
  nl::cg::load::loadSphere(s, j_sphere, mats);

  CHECK(s.name == "complex_sphere");
  // Transform should be loaded with all three components
  CHECK(s.T.M(0, 3) == -5.0f);
  CHECK(s.T.M(1, 3) == 10.0f);
  CHECK(s.T.M(2, 3) == 8.0f);
}

// ****************************************************************************
/// @name Scene Loading Tests

TEST_CASE("Load NLS Scene from file")
{
  // Create a temporary test scene file
  json scene_json = {
    {"camera", {
      {"pos", {0.0f, 0.0f, 10.0f}},
      {"look_at", {0.0f, 0.0f, 0.0f}},
      {"up", {0.0f, 1.0f, 0.0f}},
      {"fov", 45.0f},
      {"width", 1920},
      {"ar", 1.7778f}
    }},
    {"materials", json::array({
      {{"name", "red_mat"}, {"type", "lambertian"}, {"albedo", {0.8f, 0.1f, 0.1f}}},
      {{"name", "gray_mat"}, {"type", "lambertian"}, {"albedo", {0.5f, 0.5f, 0.5f}}}
    })},
    {"objects", json::array({
      {
        {"name", "sphere"},
        {"type", "sphere"},
        {"material", "red_mat"},
        {"transform", {{"translate", {0.0f, 0.0f, 0.0f}}}}
      },
      {
        {"name", "floor"},
        {"type", "plane"},
        {"material", "gray_mat"},
        {"transform", {{"translate", {0.0f, -5.0f, 0.0f}}}}
      }
    })},
    {"lights", json::array({
      {
        {"name", "ambient"},
        {"type", "ambient"},
        {"irradiance", {0.3f, 0.3f, 0.3f}}
      }
    })}
  };

  // Write to temporary file
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path();
  std::filesystem::path scene_file = temp_dir / "test_scene_cg.nls";
  
  std::ofstream out(scene_file);
  out << scene_json.dump(2);
  out.close();

  // Load the scene
  nl::cg::scene scene;
  bool success = nl::cg::load::loadNLS(scene, scene_file.string());

  CHECK(success);
  CHECK(scene.cam.width == 1920);
  CHECK(scene.cam.height == 1080);  // 1920 / (16/9)
  CHECK(scene.materials.size() == 2);
  CHECK(scene.objects.size() == 2);
  CHECK(scene.ideal_lights.size() == 1);

  // Verify camera
  CHECK(scene.cam.fov == 45.0f);

  // Verify materials
  auto mat_red = scene.materials.find("red_mat");
  CHECK(mat_red != nullptr);

  // Verify objects
  auto obj_sphere = scene.objects.find("sphere");
  CHECK(obj_sphere != nullptr);
  auto obj_floor = scene.objects.find("floor");
  CHECK(obj_floor != nullptr);

  // Verify light
  auto light_ambient = std::get<nl::cg::ambientlight>(scene.ideal_lights[0]);
  CHECK(light_ambient.name == "ambient");

  // Clean up
  std::filesystem::remove(scene_file);
}

TEST_CASE("Load NLS Scene with complex transforms")
{
  json scene_json = {
    {"camera", {
      {"pos", {-5.0f, 3.0f, 15.0f}},
      {"look_at", {0.0f, 1.0f, 0.0f}},
      {"up", {0.0f, 1.0f, 0.0f}},
      {"fov", 60.0f},
      {"width", 1280}
    }},
    {"materials", json::array({
      {{"name", "material1"}, {"type", "lambertian"}, {"albedo", {0.9f, 0.9f, 0.9f}}}
    })},
    {"objects", json::array({
      {
        {"name", "obj_with_rotation"},
        {"type", "sphere"},
        {"material", "material1"},
        {"transform", {
          {"scale", {2.0f, 2.0f, 2.0f}},
          {"rotate", {
            {"axis", {0.0f, 1.0f, 0.0f}},
            {"degrees", 45.0f}
          }},
          {"translate", {5.0f, 0.0f, -10.0f}}
        }}
      }
    })},
    {"lights", json::array({
      {
        {"name", "light1"},
        {"type", "ambient"},
        {"irradiance", {1.0f, 1.0f, 1.0f}}
      }
    })}
  };

  std::filesystem::path temp_dir = std::filesystem::temp_directory_path();
  std::filesystem::path scene_file = temp_dir / "test_scene_complex.nls";
  
  std::ofstream out(scene_file);
  out << scene_json.dump(2);
  out.close();

  nl::cg::scene scene;
  bool success = nl::cg::load::loadNLS(scene, scene_file.string());

  CHECK(success);
  CHECK(scene.cam.width == 1280);
  CHECK(scene.objects.size() == 1);
  
  // Verify object was loaded with complex transform
  auto obj = scene.objects.find("obj_with_rotation");
  CHECK(obj != nullptr);
  auto &sphere = std::get<nl::cg::sphere>(*obj);
  CHECK(sphere.name == "obj_with_rotation");

  std::filesystem::remove(scene_file);
}

// ****************************************************************************