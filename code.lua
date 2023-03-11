--[[
I want to merge this with the charts
but it requires extra code modules
so for now I'll put all code here
then maybe later I'll move chart-specific code into charts and leave the extras here?
each accepts input vec3 latLonHeight
latLonHeight.x = latitude in degrees
latLonHeight.y = longitude in degrees
latLonHeight.z = height above sealevel, in meters 
--]]
local template = require 'template'
local ModuleSet = require 'modules'

local modules = ModuleSet()

modules:addFromMarkup(template[[
//// MODULE_NAME: M_PI
const float M_PI = <?=math.pi?>;

//// MODULE_NAME: M_SQRT_2
const float M_SQRT_2 = sqrt(2.);

//// MODULE_NAME: rad
//// MODULE_DEPENDS: M_PI

float rad(float d) {
	return d * M_PI / 180.;
}

//// MODULE_NAME: deg
//// MODULE_DEPENDS: M_PI

float deg(float d) {
	return d / M_PI * 180.;
}

//// MODULE_NAME: perp2
vec2 perp2(vec2 a) {
	return vec2(-a.y, a.x);
}

//// MODULE_NAME: isfinite
bool isfinite(float x) {
	return !(isinf(x) || isnan(x));
}

//// MODULE_NAME: WGS84_a

const float WGS84_a = 6378137.;		// equatorial radius

//// MODULE_NAME: chart_sphere
//// MODULE_DEPENDS: M_PI rad WGS84_a

vec3 chart_sphere(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float theta = rad(90. - lat);
	float phi = rad(lon);
	float r = 1. + height / WGS84_a;
	float sinth = sin(theta);
	vec3 pt = r * vec3(
		sinth * cos(phi),
		sinth * sin(phi),
		cos(theta)
	);
	//convert from z-towards-user (3D) to z-up (2D)
	pt.yz = -perp2(pt.yz);	//rotate back so pt is up
	pt.xz = perp2(pt.xz);		//now rotate so prime meridian is along -z instead of +x
	return pt;
}

//// MODULE_DEPENDS: deg

vec3 chartInv_sphere(vec3 pt) {
	// convert from z-up 2D to z-towards-user 3D
	pt.xz = -perp2(pt.xz);
	pt.yz = perp2(pt.yz);
	float r = length(pt);
	float phi = atan(pt.y, pt.x);	//atan2
	float r2 = length(pt.xy);
	float theta = atan(r2 / pt.z);
	float height = (r - 1.) * WGS84_a;
	return vec3(
		deg(theta),
		mod(deg(phi) + 180., 360) - 180.,
		height);
}

//// MODULE_NAME: chart_WGS84
//// MODULE_DEPENDS: M_PI perp2 rad WGS84_a

const float WGS84_b = 6356752.3142;	// polar radius
const float WGS84_esq = 1. - WGS84_b * WGS84_b / (WGS84_a * WGS84_a);
const float WGS84_e = sqrt(WGS84_esq);
const float WGS84_flattening = 1. - WGS84_b / WGS84_a;
const float WGS84_inverseFlattening = 298.257223563;
const float WGS84_eccentricitySquared = (2. * WGS84_inverseFlattening - 1.) / (WGS84_inverseFlattening * WGS84_inverseFlattening);

float WGS84_calc_N(
	float sinTheta
) {
	float denom = sqrt(1. - WGS84_eccentricitySquared * sinTheta * sinTheta);
	return WGS84_a / denom;
}

// TODO instead of making one chart depend on another, put the WGS84 constants in one place

vec3 chart_WGS84(vec3 x) {
	float lat = x.x;
	float lon = x.y;
	float height = x.z;

	float phi = rad(lon);		// spherical φ
	float theta = rad(lat);		// spherical inclination angle (not azumuthal θ)
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);
	
	float N = WGS84_calc_N(sinTheta);
	
	float NPlusH = N + height;
	vec3 y = vec3(
		NPlusH * cosTheta * cos(phi),
		NPlusH * cosTheta * sin(phi),
		(N * (1. - WGS84_eccentricitySquared) + height) * sinTheta
	);
	// at this point we're in meters, matching the geographic-charts code
	// but now I'm going to transform further to match the seismographic-visualization / geo-center-earth code
	y /= WGS84_a;			//convert from meters to normalized coordinates
	y.yz = -perp2(y.yz);	//rotate back so y is up
	y.xz = perp2(y.xz);		//now rotate so prime meridian is along -z instead of +x
	return y;
}

//// MODULE_NAME: chart_cylinder
//// MODULE_DEPENDS: perp2 rad WGS84_a

vec3 chart_cylinder(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float latrad = rad(lat);
	float lonrad = rad(lon);
	float r = WGS84_a + height;
	float x = r * cos(lonrad);
	float y = r * sin(lonrad);
	float z = r * latrad;
	vec3 cartpos = vec3(x, y, z);
	// end of geographic-charts, beginning of vis aligning stuff
	cartpos /= WGS84_a;
	cartpos.yz = -perp2(cartpos.yz);	//rotate back so cartpos is up
	cartpos.xz = perp2(cartpos.xz);		//now rotate so prime meridian is along -z instead of +x
	return cartpos;
}

//// MODULE_NAME: chart_Equirectangular
//// MODULE_DEPENDS: M_PI rad WGS84_a

const float Equirectangular_R = 2. / M_PI;
const float Equirectangular_lambda0 = 0.;
const float Equirectangular_phi0 = 0.;
const float Equirectangular_phi1 = 0.;
const float cos_Equirectangular_phi1 = cos(Equirectangular_phi1);
vec3 chart_Equirectangular(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float latrad = rad(lat);
	float lonrad = rad(lon);
	float x = Equirectangular_R * (lonrad - Equirectangular_lambda0) * cos_Equirectangular_phi1;
	float y = Equirectangular_R * (latrad - Equirectangular_phi0);
	float z = height / WGS84_a;
	return vec3(x,y,z);
}

//// MODULE_NAME: chart_Azimuthal_equidistant 
//// MODULE_DEPENDS: M_PI M_SQRT_2 rad WGS84_a

vec3 chart_Azimuthal_equidistant(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float lonrad = rad(lon);
	float colat = 90. - lat;		//[N,S] => [0,180]
	float azimuthal = colat / 180.;	//[N,S] => [0,1]
	azimuthal *= M_SQRT_2;
	float x = sin(lonrad) * azimuthal;
	float y = -cos(lonrad) * azimuthal;
	float z = height / WGS84_a;
	return vec3(x,y,z);
}
	
//// MODULE_NAME: chart_Lambert_Azimuthal_equal_area
//// MODULE_DEPENDS: M_PI M_SQRT_2 rad
vec3 chart_Lambert_Azimuthal_equal_area(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float colat = 90 - lat;	// spherical friendly
	float polarR = M_SQRT_2 * sin(.5 * rad(colat));
	float lonrad = rad(lon);
	float x = sin(lonrad) * polarR;
	float y = -cos(lonrad) * polarR;
	float z = height / WGS84_a;
	return vec3(x, y, z);
}

//// MODULE_NAME: chart_Mollweide
//// MODULE_DEPENDS: M_PI M_SQRT_2 rad isfinite WGS84_a

const float M_SQRT_8 = sqrt(8.);
const float Mollweide_R = M_PI / 4.;
const float Mollweide_lambda0 = 0.;	// in degrees

vec3 chart_Mollweide(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float lonrad = rad(lon);
	float lambda = lonrad;
	float latrad = rad(lat);
	float phi = latrad;
	float theta;
	if (phi == .5 * M_PI) {
		theta = .5 * M_PI;
	} else {
		theta = phi;
		for (int i = 0; i < 10; ++i) {
			float dtheta = (2. * theta + sin(2. * theta) - M_PI * sin(phi)) / (2. + 2. * cos(theta));
			if (abs(dtheta) < 1e-5) break;
			theta -= dtheta;
		}
	}
	float mollweidex = Mollweide_R * M_SQRT_8 / M_PI * (lambda - Mollweide_lambda0) * cos(theta);
	float mollweidey = Mollweide_R * M_SQRT_2 * sin(theta);
	float mollweidez = height / WGS84_a;
	if (!isfinite(mollweidex)) mollweidex = 0;
	if (!isfinite(mollweidey)) mollweidey = 0;
	if (!isfinite(mollweidez)) mollweidez = 0;
	return vec3(mollweidex, mollweidey, mollweidez);
}

]])

local allChartCode = modules:getCodeAndHeader(
	-- 3D
	'chart_sphere',
	'chart_WGS84',
	'chart_cylinder',
	-- 2D
	'chart_Equirectangular',
	'chart_Azimuthal_equidistant',
	'chart_Lambert_Azimuthal_equal_area',
	'chart_Mollweide'
)

return allChartCode
