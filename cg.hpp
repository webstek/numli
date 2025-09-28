// ****************************************************************************
/// @file cg.hpp
/// @author Kyle Webster
/// @version 0.1
/// @date 27 Sep 2025
/// @brief Numerics Library - Computer Graphics - @ref cg
/// @details
/// Collection of computer graphics structures and algorithms
// ****************************************************************************
#include "bra.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
// ****************************************************************************

namespace nl {
namespace cg {
// ** Contents ****************************************************************
  class Mesh;
  Mesh loadPLY(std::string);

// Aliases:
  using vec3 = bra::NdReal<3, double>;      ///< 3d vector
  using vec4 = bra::NdReal<4, double>;      ///< vector in the space PG2
  using point3 = bra::NdReal<3, double>;    ///< 3d point
  using point4 = bra::NdReal<4, double>;    ///< point in the space PG2
  using color = bra::NdReal<3, uint8_t>;    ///< rgb colour (0 to 255)

// Structs:
  /// @brief Representation of a ray, a direction from a point in space.
  typedef struct {
    point3 x;
    vec3 u;
  } Ray;

// Classes:
  /// @brief Triangular mesh
  class Mesh {
  public:
    // Members:
    const uint32_t vertex_count;
    const uint32_t face_count;
    const bra::Mat<double> vertices;            /**< Vertex positions */
    const bra::Mat<uint32_t> vertex_indices;    /**< Face vertex indices */
  };

// Functions:
  /// @brief Produce a Mesh object from a .ply file.
  /// @param path Source file path.
  /// @return mesh Triangular mesh data.
  Mesh loadPLY(std::string path) {
    // prep mesh fields
    uint32_t vert_count = 0;
    uint32_t face_count = 0;
    bra::Mat<double> vertices (3, 1);
    bra::Mat<uint32_t> indices (3, 1);

    // open file stream
    std::ifstream ply;
    ply.open(path);
    if (!ply.is_open()) {
      throw std::runtime_error("failed to open " + path);
    } else {
      // Header is first
      bool is_vertex_line = false;
      bool is_face_line = false;
      std::string line;
      while (std::getline(ply,line)) {
        // make vector of words on the line
        std::vector<std::string> words;
        std::istringstream stream (line);
        std::string word;
        while (stream >> word) words.push_back(word);

        // look for vertex count
        if (words[0] == "element" && words[1] == "vertex") {
          // convert ascii to uint32_t
          vert_count = uint32_t(std::stoul(words[2]));
        }

        // look for face count
        if (words[0] == "element" && words[1] == "face") {
          // convert ascii to uint32_t
          face_count = uint32_t(std::stoul(words[2]));
        }

        if (is_face_line) {
          uint32_t elems[3] = {
            std::stoul(words[1]),
            std::stoul(words[2]),
            std::stoul(words[3])};
          indices.pushCol(elems);
        }

        if (is_vertex_line) {
          double elems[3] = {
            std::stod(words[0]),
            std::stod(words[1]),
            std::stod(words[2])};
          vertices.pushCol(elems);
        }

        // Determine if the next line in the header, vertex data, or triangle
        //   vertex index data.
        if (words[0] == "end_header") is_vertex_line = true;
        if (vertices.cols() == vert_count) {
          is_vertex_line = false;
          is_face_line = true;
        }
      }
      ply.close();
    }

    // make the mesh to return
    Mesh mesh {vert_count, face_count, vertices, indices};
    return mesh;
  }
} // namespace cg
} // namespace nl