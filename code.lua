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
<?
local symmath = require 'symmath'
symmath.export.C.numberType = 'float'	-- hmm ... nice to be an arg ...
local charts = require 'geographic-charts'
?>

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

//// MODULE_NAME: sinc

float sinc(float x) {
	if (x == 0.) return 1.;
	return sin(x) / x;
}

//// MODULE_NAME: xformZBackToZUp
//// MODULE_DEPENDS: perp2

vec3 xformZBackToZUp(vec3 pt) {
	//convert from z-towards-user (3D) to z-up (2D)
	pt.yz = -perp2(pt.yz);		//rotate back so pt is up
	pt.xz = perp2(pt.xz);		//now rotate so prime meridian is along -z instead of +x
	return pt;
}

//// MODULE_NAME: xformZUpToZBack
//// MODULE_DEPENDS: perp2

vec3 xformZUpToZBack(vec3 pt) {
	// convert from z-up 2D to z-towards-user 3D
	pt.xz = -perp2(pt.xz);
	pt.yz = perp2(pt.yz);
	return pt;
}

//// MODULE_NAME: isfinite
bool isfinite(float x) {
	return !(isinf(x) || isnan(x));
}

//// MODULE_NAME: WGS84_a

const float WGS84_a = 6378137.;		// equatorial radius

<?=charts.cylinder:getGLSLModule()?>
<?=charts.Equirectangular:getGLSLModule()?>
<?=charts.Mercator:getGLSLModule()?>
<?=charts['Gall-Peters']:getGLSLModule()?>
<?=charts['Lambert cylindrical equal-area']:getGLSLModule()?>
<?=charts['Azimuthal equidistant']:getGLSLModule()?>
<?=charts['Lambert azimuthal equal-area']:getGLSLModule()?>
<?=charts.Sinusoidal:getGLSLModule()?>
<?=charts['Winkel tripel']:getGLSLModule()?>
<?=charts['Kavrayskiy VIII']:getGLSLModule()?>
<?=charts.Wiechel:getGLSLModule()?>
<?=charts.Albers:getGLSLModule()?>
<?=charts.Bonne:getGLSLModule()?>

<?=charts.sphere:getGLSLModule()?>
// still in chart_sphere module ...
//// MODULE_DEPENDS: xformZUpToZBack deg

vec3 chartInv_sphere(vec3 pt) {
	pt = xformZUpToZBack(pt);
	float r = length(pt);
	float lonrad = atan(pt.y, pt.x);
	float latrad = atan(pt.z, length(pt.xy));
	float height = (r - 1.) * WGS84_a;
	return vec3(deg(latrad), deg(lonrad), height);
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

//// MODULE_NAME: chart_Mollweide
//// MODULE_DEPENDS: M_PI M_SQRT_2 rad isfinite WGS84_a

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
	float x = Mollweide_R * 2. * M_SQRT_2 / M_PI * (lambda - Mollweide_lambda0) * cos(theta);
	float y = Mollweide_R * M_SQRT_2 * sin(theta);
	float z = height / WGS84_a;
	if (!isfinite(x)) x = 0;
	if (!isfinite(y)) y = 0;
	if (!isfinite(z)) z = 0;
	return vec3(x, y, z);
}

]])

-- https://en.wikipedia.org/wiki/List_of_map_projections
local allChartCode = modules:getCodeAndHeader(
	-- 3D
	'chart_sphere',
	'chart_WGS84',
	'chart_cylinder',
	-- 2D - rectangular
	'chart_Equirectangular',
	'chart_Mercator',
	-- transverse Mercator ... just puts arctic at mid-top and antarctic at mid-bottom
	'chart_Gall_Peters',
	'chart_Lambert_cylindrical_equal_area',
	-- https://en.wikipedia.org/wiki/Hobo%E2%80%93Dyer_projection
	-- https://en.wikipedia.org/wiki/Behrmann_projection
	-- 2D - rect-ellipse
	'chart_Kavrayskiy_VIII',
	'chart_Winkel_tripel',
	-- https://en.wikipedia.org/wiki/Eckert_IV_projection
	-- https://en.wikipedia.org/wiki/Eckert_VI_projection
	-- https://en.wikipedia.org/wiki/Equal_Earth_projection
	-- 2D - ellipse
	-- https://en.wikipedia.org/wiki/Hammer_projection
	-- requires an integral of a power ... https://en.wikipedia.org/wiki/Tobler_hyperelliptical_projection
	'chart_Mollweide',
	-- 2D - pinched ellipse
	'chart_Sinusoidal',
	-- 2D - circle
	'chart_Azimuthal_equidistant',
	'chart_Lambert_azimuthal_equal_area',
	-- Robinson ... is an interpolation table-based: https://en.wikipedia.org/wiki/Robinson_projection
	-- Stereographic
	'chart_Wiechel',
	-- 2D - conic
	'chart_Albers',
	-- Albers generalizes this: https://en.wikipedia.org/wiki/Lambert_equal-area_conic_projection
	-- https://en.wikipedia.org/wiki/Collignon_projection
	-- https://en.wikipedia.org/wiki/Eckert_II_projection
	'chart_Bonne'
	-- https://en.wikipedia.org/wiki/Bottomley_projection
	-- https://en.wikipedia.org/wiki/Werner_projection
	-- https://en.wikipedia.org/wiki/Strebe_1995_projection
	-- multiple slices
	-- https://en.wikipedia.org/wiki/Boggs_eumorphic_projection
	-- https://en.wikipedia.org/wiki/Goode_homolosine_projection

	-- https://en.wikipedia.org/wiki/Eckert-Greifendorff_projection
)

return allChartCode
