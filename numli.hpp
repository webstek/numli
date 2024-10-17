/**
 * Numerics Library - numli
 * A library collecting various useful computational tools.
 */


#ifndef biblio
#define biblio


namespace NL {

namespace cg {
	/**
	 * 3 component vector of doubles
	 */
	class vec3 {
	public:
		double e[3];

		vec3() : e{0,0,0} {};
		vec3(double e0, double e1, double e2) : e{e0, e1, e2} {};
	};

	/**
	 * 3d point type alias of vec3
	 */
	using point3 = vec3;


	/**
	 * 4 component vvector of doubles
	 */
	class vec4 {
	public:
		double e[4];
		vec4() : e{0,0,0,1} {};
		vec4(double e0, double e1, double e2, double e3) : e{e0, e1, e2, e3} {};
	};

	/**
	 * rgb alpha color type alias of vec4
	 */
	using rgba = vec4;
}


}

#endif