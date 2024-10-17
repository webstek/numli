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
}


}

#endif