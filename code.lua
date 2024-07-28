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

-- https://en.wikipedia.org/wiki/List_of_map_projections
local table = require 'ext.table'

-- make sure you :build() all charts first
-- I haven't got detection/caching of :build() or I'd just call it in here too
local function codeForCharts(charts)
	local modules = ModuleSet()

	local code = template[[
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
//pt = vec3(pt.x, pt.z, -pt.y);
	pt.xz = perp2(pt.xz);		//now rotate so prime meridian is along -z instead of +x
//pt = vec3(-pt.z, pt.y, pt.x);
//...combined...
//pt = vec3(pt.y, pt.z, pt.x);
//... with inverse ...
//pt = vec3(pt.z, pt.x, pt.y);

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

]]

	local symmath = require 'symmath'
	symmath.export.C.numberType = 'float'	-- hmm ... nice to be an arg ...

	for i,chart in ipairs(charts) do
		code = code .. chart:getGLSLModule() .. '\n'
	end

	modules:addFromMarkup(code)

	return modules:getCodeAndHeader(table.mapi(charts, function(c) return c:getSymbol() end):unpack())
end

return codeForCharts
