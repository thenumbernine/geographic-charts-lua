-- I want to merge this with the charts
-- but it requires extra code modules
-- so for now I'll put all code here
-- then maybe later I'll move chart-specific code into charts and leave the extras here?

local template = require 'template'
local ModuleSet = require 'modules'

local modules = ModuleSet()
modules:addFromMarkup(template[[
//// MODULE_NAME: M_PI
const float M_PI = <?=math.pi?>;

//// MODULE_NAME: rad
//// MODULE_DEPENDS: M_PI
float rad(float d) {
	return d * M_PI / 180.;
}

//// MODULE_NAME: perp2
vec2 perp2(vec2 a) {
	return vec2(-a.y, a.x);
}

//// MODULE_NAME: isfinite
bool isfinite(float x) {
	return !(isinf(x) || isnan(x));
}
]])
	local code_WGS84 = [[
//// MODULE_NAME: chart_WGS84
//// MODULE_DEPENDS: M_PI perp2 rad

const float WGS84_a = 6378137.;		// equatorial radius
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
]]
modules:addFromMarkup(code_WGS84)
	
	local code_cylinder = [[
//// MODULE_NAME: chart_cylinder
//// MODULE_DEPENDS: perp2 rad chart_WGS84

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
]]
modules:addFromMarkup(code_cylinder)

	-- TODO instead of making one chart depend on another, put the WGS84 constants in one place
	local code_Equirectangular = [[
//// MODULE_NAME: chart_Equirectangular
//// MODULE_DEPENDS: M_PI rad chart_WGS84

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
]]
modules:addFromMarkup(code_Equirectangular)

	local code_Azimuthal_equidistant = [[
//// MODULE_NAME: chart_Azimuthal_equidistant 
//// MODULE_DEPENDS: M_PI rad chart_WGS84

vec3 chart_Azimuthal_equidistant(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;
	float latrad = rad(lat);
	float lonrad = rad(lon);
	float azimuthal = M_PI / 2. - latrad;
	float x = -sin(lonrad + M_PI) * azimuthal;
	float y = cos(lonrad + M_PI) * azimuthal;
	float z = height / WGS84_a;
	return vec3(x,y,z);
}
]]
modules:addFromMarkup(code_Azimuthal_equidistant)
	
	local code_Mollweide = [[
//// MODULE_NAME: chart_Mollweide
//// MODULE_DEPENDS: M_PI rad isfinite chart_WGS84

const float M_SQRT_2 = sqrt(2.);
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
]]
modules:addFromMarkup(code_Mollweide)

local allChartCode = modules:getCodeAndHeader(
	'chart_WGS84',
	'chart_cylinder',
	'chart_Equirectangular',
	'chart_Azimuthal_equidistant',
	'chart_Mollweide'
)

return allChartCode
