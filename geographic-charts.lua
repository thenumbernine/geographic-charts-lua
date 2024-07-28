--[[
I use these often enough that they are going to go into a library soon enough

2D charts will be xy coordinates, z will be for height

3D charts  are z+ = north pole, x+ = prime meridian

each chart has the following:
	.chart(lat, lon, height) = returns the x,y,z for the chart's coordinates.
		lat = [-90,90] in degrees
		lon = [-180,180] in degrees
		height >= 0 in meters
		returns x,y,z in meters
	.basis(lat, lon, height) = returns the ex,ey,ez basis for chart's coordinates ... right-handed system, equal to (d/dlat, d/dlon, -d/dheight) = (north, east, inwards) (to match the WMM basis)


GAAAHHH STANDARDS
physics spherical coordinates: the longitude is φ and the latitude (starting at the north pole and going down) is θ ...
mathematician spherical coordinates: the longitude is θ and the latitude is φ ...
geographic / map charts: the longitude is λ and the latitude is φ ...
so TODO change the calc_* stuff from r_θφ to h_φλ? idk ...

really tempting to just make all charts -- even the 2D ones -- in meters.
and make them in which coordinates?  xyz?
--]]

-- using geographic labels: lat = φ, lon = λ
local table = require 'ext.table'
local class = require 'ext.class'
local math = require 'ext.math'
local vec3d = require 'vec-ffi.vec3d'
local symmath = require 'symmath'
local clnumber = require 'cl.obj.number'

-- input
local latvar = symmath.var'lat'
local lonvar = symmath.var'lon'
--local heightvar = symmath.var'height'
--local latvar = symmath.set.RealInterval(-90, 90, true, true):var'lat'
--local lonvar = symmath.set.RealInterval(-180, 180, true, true):var'lon'
local heightvar = symmath.set.positiveReal:var'height'

-- GLSL-specific, to accept (lat,lon,height) as a vec3
latvar:nameForExporter('C', 'latLonHeight.x')
lonvar:nameForExporter('C', 'latLonHeight.y')
heightvar:nameForExporter('C', 'latLonHeight.z')
-- expressions
local latradval = latvar * symmath.pi / 180
local lonradval = lonvar * symmath.pi / 180

local WGS84_a = 6378137	-- m ... earth equitorial radius

-- tail-call transform list of fields into their values within 't'
local function getfields(t, ...)
	if select('#', ...) == 0 then return end
	local f = ...
	return t[f], getfields(t, select(2, ...))
end

local Chart = class()

function Chart:getCName()
	return (self.name:gsub('[%s%-]', '_'))
end

function Chart:getSymbol()	-- aka GLSL name
	return 'chart_'..self:getCName()
end

-- call this to build numeric fields of each name
-- and to build chart.vars with its fields in-order matching ...
function Chart:buildVars(varkvs)
	self.varnames = table()
	self.vars = {}	-- index by name
	self.varlist = table()	-- index sequence
	for i,kv in ipairs(varkvs or {}) do
		local k, v = next(kv)
		self.varnames[i] = k
		local var = symmath.var(k)
		var:nameForExporter('C', self:getCName()..'_'..k)
		self.vars[k] = var
		self.varlist[i] = var
		self[k] = v
	end
end

-- static method for building code
-- run this Chart:buildFunc()
-- this assigns chart.guiVars, chart.exprOut
-- then it builds chart.exprIn as a concatenation of {latvar,lonvar,heightvar} and whatever is in chart.guiVars
function Chart:buildFunc(exprOut)
	self.exprOut = exprOut
	self.exprIn = table{latvar, lonvar, heightvar}
		:append(self.varlist or {})

	-- TODO for 3D this is in z-back 3D coords?
	-- it doesn't match code.lua's 2D z-up coords
	self.chartFunc = symmath.export.Lua:toFunc{
		input = self.exprIn,
		output = self.exprOut,
	}

	-- for charts with serious problems
	if self.skipBasis then return end

	local y = symmath.Array(table.unpack(self.exprOut))
	self.basisExprs = table{latvar, lonvar, heightvar}:mapi(function(x,i)
		if i == 3 then 	--x == heightvar then
			return (-y):diff(x)()
		else
			return y:diff(x)()
		end
	end)
	if self.normalizeBasisNumerically then
		self.basisFunc = symmath.export.Lua:toFunc{
			input = self.exprIn,
			output = self.basisExprs,
		}
	else
		self.basisNormExprs = self.basisExprs:mapi(function(dy_dx)
			local len = symmath.sqrt(
				dy_dx[1]^2
				+ dy_dx[2]^2
				+ dy_dx[3]^2
			)
			return (dy_dx / len)()
		end)
		self.basisFunc = symmath.export.Lua:toFunc{
			input = self.exprIn,
			output = self.basisNormExprs,
		}
	end
end

-- Be sure to define your `var:nameForExporter('C', ...)` beforehand
-- TODO what about GLSL uniforms?  and providing default values?
function Chart:getGLSLBody()
	return '\t'..symmath.export.C:toCode{
		input = self.exprIn,
		output = self.exprOut,
	}:gsub('\n', '\n\t')
end

function Chart:getGLSLBasisBody()
	if self.skipBasis then
		return table{
			'	vec3 ex = vec3(1., 0., 0.);',
			'	vec3 ey = vec3(0., 1., 0.);',
			'	vec3 ez = vec3(0., 0., 1.);',
		}:concat'\n'
	end

	if self.normalizeBasisNumerically then
		return table{
			'\t'..symmath.export.C:toCode{
				input = self.exprIn,
				-- unwrap these for GLSL ...
				output = {
					self.basisExprs[1][1], self.basisExprs[1][2], self.basisExprs[1][3],
					self.basisExprs[2][1], self.basisExprs[2][2], self.basisExprs[2][3],
					self.basisExprs[3][1], self.basisExprs[3][2], self.basisExprs[3][3],
				},
			}:gsub('\n', '\n\t'),
			'\tvec3 ex = normalize(vec3(out1, out2, out3));',
			'\tvec3 ey = normalize(vec3(out4, out5, out6));',
			'\tvec3 ez = normalize(vec3(out7, out8, out9));',
		}:concat'\n'
	else
		return table{
			'\t'..symmath.export.C:toCode{
				input = self.exprIn,
				output = {
					self.basisNormExprs[1][1], self.basisNormExprs[1][2], self.basisNormExprs[1][3],
					self.basisNormExprs[2][1], self.basisNormExprs[2][2], self.basisNormExprs[2][3],
					self.basisNormExprs[3][1], self.basisNormExprs[3][2], self.basisNormExprs[3][3],
				},
			}:gsub('\n', '\n\t'),
			'\tvec3 ex = vec3(out1, out2, out3);',
			'\tvec3 ey = vec3(out4, out5, out6);',
			'\tvec3 ez = vec3(out7, out8, out9);',
		}:concat'\n'
	end
end

function Chart:getGLSLFunc()
	local escname = self:getCName()
	return self.varnames:mapi(function(k)
		return 'const float '..escname..'_'..k..' = '..clnumber(self[k])..';'
	end):append{
		'vec3 '..self:getSymbol()..'(vec3 latLonHeight) {',
		self:getGLSLBody(),
		'	return vec3(out1, out2, out3);',
		'}',
		'mat3 '..self:getSymbol()..'_basis(vec3 latLonHeight) {',
		self:getGLSLBasisBody(),
		'	return mat3(ex, ey, ez);',
		'}',
	}:concat'\n'
end

-- the '3D' indicates to be sure to insert a `xformZBackToZUp(pt)` last
function Chart:getGLSLFunc3D()
	local escname = self:getCName()
	return self.varnames:mapi(function(k)
		return 'const float '..escname..'_'..k..' = '..clnumber(self[k])..';'
	end):append{
		'vec3 '..self:getSymbol()..'(vec3 latLonHeight) {',
		self:getGLSLBody(),
		'	return xformZBackToZUp(vec3(out1, out2, out3));',
		'}',
		'mat3 '..self:getSymbol()..'_basis(vec3 latLonHeight) {',
		self:getGLSLBasisBody(),
		'	e = mat3(ex, ey, ez);',
		'	e = transpose(e);',
		'	e[0] = xformZBackToZUp(e[0]);',
		'	e[1] = xformZBackToZUp(e[1]);',
		'	e[2] = xformZBackToZUp(e[2]);',
		--'	return mat3(xformZBackToZUp(ex), xformZBackToZUp(ey), xformZBackToZUp(ez));',
		-- right now, for 3d, the chart coords needs to be permuted by the user
		-- so lets force them to permute the basis as well ...
		-- but the downside to not converting is, now you have to convert mid-summing if you are using weighted combinations of graphs ... hmm ...
		--'	return mat3(ex, ey, ez);',
		'}',
	}:concat'\n'
end

function Chart:getGLSLModule()
	return table{
		'//// MODULE_NAME: '..self:getSymbol(),
		'//// MODULE_DEPENDS: M_PI WGS84_a xformZBackToZUp sinc',
		not self.is3D
			and self:getGLSLFunc()
			or self:getGLSLFunc3D(),
	}:concat'\n'
end

function Chart:chart(lat, lon, height)
	return self.chartFunc(lat, lon, height, getfields(self, table.unpack(self.varnames)))
end

-- The vectors are (d/dlat, d/dlon, -d/dheight)
-- Don't forget that GLSL is swapping the z-back for z-up.
-- For a few of these (sphere, cylinder, etc) multiply the output by WGS84_a to put it in meters
function Chart:basis(lat, lon, height)
	if not self.basisFunc then
		local delta = 1e-3
		local x = (vec3d(self:chart(lat+delta, lon, height)) - vec3d(self:chart(lat-delta, lon, height))) * (1 / (2 * delta))
		local y = (vec3d(self:chart(lat, lon+delta, height)) - vec3d(self:chart(lat, lon-delta, height))) * (1 / (2 * delta))
		x = x:normalize()
		y = y:normalize()
		local z = x:cross(y)
		return x,y,z
	else
		local x,y,z = self.basisFunc(lat, lon, height, getfields(self, table.unpack(self.varnames)))
		x = vec3d(table.unpack(x))
		y = vec3d(table.unpack(y))
		z = vec3d(table.unpack(z))
		if self.normalizeBasisNumerically then
			x = x:normalize()
			y = y:normalize()
			z = z:normalize()
		end
		return x,y,z
	end
end

function Chart:updateGUI()
	local ig = require 'imgui'
	if self.varnames then
		for _,n in ipairs(self.varnames) do
			ig.luatableInputFloat(n, self, n)
		end
	end
end


local charts = {
	(function()
		local c = Chart:subclass()

		c.name = 'WGS84'
		c.is3D = true

		-- specific to WGS84:
		c.a = WGS84_a
		c.b = 6356752.3142	-- m ... earth polar radius
		-- ... wait, ellipses use 'a' as the major and 'b' as the minor
		--   which would mean 'a' is the equatorial and 'b' is the polar
		--   did I mix these up?
		--   or did physicists screw up convention again?
		-- symmath vars equivalent
		local avar = symmath.var'a'
		local bvar = symmath.var'b'
		local flatteningVal = 1 - bvar / avar
		local inverseFlatteningVal = 1 / flatteningVal
		local eccentricitySquaredVal = (2 * inverseFlatteningVal - 1) / (inverseFlatteningVal * inverseFlatteningVal)
		local e = math.sqrt(1 - c.b * c.b / (c.a * c.a))
		c.esq = e * e

		local flattening = 1 - c.b / c.a
		c.inverseFlattening = 298.257223563
		c.eccentricitySquared = (2 * c.inverseFlattening - 1) / (c.inverseFlattening * c.inverseFlattening)

		function c:calc_N(sinTheta, equatorialRadius, eccentricitySquared)
			local denom = math.sqrt(1 - eccentricitySquared * sinTheta * sinTheta)
			return equatorialRadius / denom
		end

		function c:calc_dN_dTheta(sinTheta, cosTheta, equatorialRadius, eccentricitySquared)
			local denom = math.sqrt(1 - eccentricitySquared * sinTheta * sinTheta)
			return eccentricitySquared * sinTheta * cosTheta * equatorialRadius / (denom * denom * denom)
		end

		-- |d(x,y,z)/d(h,θ,φ)| for h=0
		function c:dx_dsphere_det_h_eq_0(lat)
			local theta = math.rad(lat)		-- spherical inclination angle (not azumuthal θ)
			local sinTheta = math.sin(theta)
			local cosTheta = math.cos(theta)

			local h = 0
			local N = self:calc_N(sinTheta, self.a, self.eccentricitySquared)
			local dN_dTheta = self:calc_dN_dTheta(sinTheta, cosTheta, self.a, self.eccentricitySquared)
			local cosTheta2 = cosTheta * cosTheta
			return -N * (
				N * cosTheta
				+ self.eccentricitySquared * cosTheta2 * N * cosTheta
				+ self.eccentricitySquared * cosTheta2 * dN_dTheta * sinTheta
			)
		end

		-- lat = [-90,90] in degrees
		-- lon = [-180,180] in degrees
		-- height >= 0 in meters
		-- returns: x,y,z in meters
		function c:chart(lat, lon, height)
			local phi = math.rad(lon)		-- spherical φ
			local theta = math.rad(lat)		-- spherical inclination angle (not azumuthal θ)
			local cosTheta = math.cos(theta)
			local sinTheta = math.sin(theta)

			local N = self:calc_N(sinTheta, self.a, self.eccentricitySquared)

			local NPlusH = N + height
			return
				NPlusH * cosTheta * math.cos(phi),
				NPlusH * cosTheta * math.sin(phi),
				(N * (1 - self.eccentricitySquared) + height) * sinTheta
		end

		-- x,y,z = meters
		-- returns lat (degrees), lon (degrees), height (meters)
		-- lat and lon has the same range as chart()
		--
		-- https://gis.stackexchange.com/questions/28446/computational-most-efficient-way-to-convert-cartesian-to-geodetic-coordinates
		function c:chartInv(x, y, z)
			-- this much is always true
			local phi = math.atan2(y, x);
			local theta			-- spherical inclination angle
			for i=1,10000 do
				-- spherical:
				local r2 = math.sqrt(x*x + y*y);
				local newtheta = math.atan(r2 / z);
				if theta then
					local dtheta = math.abs(newtheta - theta)
					--print(dtheta, z)
					if dtheta < 1e-15 then break end
				end
				theta = newtheta
				x,y,z = self:chart(math.deg(theta), math.deg(phi), 0)
			end

			local NPlusHTimesCosTheta = math.sqrt(x*x + y*y)
			local NPlusH = NPlusHTimesCosTheta / math.cos(theta)
			local height = NPlusH - self:calc_N(math.sin(theta), self.a, self.eccentricitySquared)

			-- lat, lon, height:
			return
				math.deg(theta),
				(math.deg(phi) + 180) % 360 - 180,
				height
		end

		function c:basis(lat, lon, height)
			--return latLonToCartesianTangentSpaceWGS84(lat, lon, height)
			local phi = math.rad(lat)
			local lambda = math.rad(lon)
			local cosLambda = math.cos(lambda)
			local sinLambda = math.sin(lambda)

			local cosPhi = math.cos(phi)
			local sinPhi = math.sin(phi)
			local dphi_cosPhi = -sinPhi
			local dphi_sinPhi = cosPhi

			local rCart = self.a / math.sqrt(1 - self.esq * sinPhi * sinPhi)
			local dphi_rCart = self.a / math.sqrt(1 - self.esq * sinPhi * sinPhi)^3 * self.esq * sinPhi * dphi_sinPhi

			local rCart_over_a = 1 / math.sqrt(1 - self.esq * sinPhi * sinPhi)

			local xp = (rCart + height) * cosPhi
			local dphi_xp = dphi_rCart * cosPhi + (rCart + height) * dphi_cosPhi
			local dheight_xp = cosPhi

			local xp_over_a = (rCart_over_a + height / self.a) * cosPhi

			local zp = (rCart * (1 - self.esq) + height) * sinPhi
			local dphi_zp = (dphi_rCart * (1 - self.esq)) * sinPhi + (rCart * (1 - self.esq) + height) * dphi_sinPhi
			local dheight_zp = sinPhi

			local zp_over_a = (rCart_over_a * (1 - self.esq) + height / self.a) * sinPhi

			local r2D = math.sqrt(xp * xp + zp * zp)
			local dphi_r2D = (xp * dphi_xp + zp * dphi_zp) / r2D
			local dheight_r2D = (xp * dheight_xp + zp * dheight_zp) / r2D

			local r2D_over_a = math.sqrt(xp_over_a * xp_over_a + zp_over_a * zp_over_a)
			local dphi_r2D_over_a = (xp_over_a * dphi_xp + zp_over_a * dphi_zp) / r2D

			local sinPhiSph = zp / r2D
			local dphi_sinPhiSph = (dphi_zp * r2D - zp * dphi_r2D) / (r2D * r2D)
			local dheight_sinPhiSph = (dheight_zp * r2D - zp * dheight_r2D) / (r2D * r2D)

			local cosPhiSph = math.sqrt(1 - sinPhiSph * sinPhiSph)
			--d/du sqrt(1 - x^2) = -x/sqrt(1 - x^2) dx/du
			local dphi_cosPhiSph = -sinPhi / cosPhiSph * dphi_sinPhiSph
			local dheight_cosPhiSph = -sinPhi / cosPhiSph * dheight_sinPhiSph

			--local x = r2D * cosPhiSph / self.a * cosLambda
			--local y = r2D * cosPhiSph / self.a * sinLambda
			--local z = r2D * sinPhiSph / self.a

			local dphi_x = (dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * cosLambda
			local dphi_y = (dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * sinLambda
			local dphi_z = (dphi_r2D_over_a * sinPhiSph + r2D_over_a * dphi_sinPhiSph)

			local dlambda_x = -sinLambda
			local dlambda_y = cosLambda

			local dheight_x = (dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * cosLambda
			local dheight_y = (dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * sinLambda
			local dheight_z = (dheight_r2D * sinPhiSph + r2D * dheight_sinPhiSph)

			return
				vec3d(dphi_x, dphi_y, dphi_z),
				vec3d(dlambda_x, dlambda_y, 0),
				vec3d(-dheight_x, -dheight_y, -dheight_z)
		end

		function c:getGLSLModule()
			return [[

//// MODULE_NAME: chart_WGS84
//// MODULE_DEPENDS: xformZBackToZUp M_PI rad WGS84_a

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

vec3 chart_WGS84(vec3 latLonHeight) {
	float lat = latLonHeight.x;
	float lon = latLonHeight.y;
	float height = latLonHeight.z;

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
	// at this point we're in meters
	// now rotate back
	y = xformZBackToZUp(y);
	return y;
}

mat3 chart_WGS84_basis(vec3 latLonHeight) {
	float phi = rad(latLonHeight.x);
	float lambda = rad(latLonHeight.y);
	float height = latLonHeight.z;

	float cosLambda = cos(lambda);
	float sinLambda = sin(lambda);

	float cosPhi = cos(phi);
	float sinPhi = sin(phi);
	float dphi_cosPhi = -sinPhi;
	float dphi_sinPhi = cosPhi;

	float rCart = WGS84_a / sqrt(1. - WGS84_esq * sinPhi * sinPhi);
	float tmp = sqrt(1. - WGS84_esq * sinPhi * sinPhi);
	float dphi_rCart = WGS84_a / (tmp * tmp * tmp) * WGS84_esq * sinPhi * dphi_sinPhi;

	float rCart_over_a = 1. / sqrt(1. - WGS84_esq * sinPhi * sinPhi);

	float xp = (rCart + height) * cosPhi;
	float dphi_xp = dphi_rCart * cosPhi + (rCart + height) * dphi_cosPhi;
	float dheight_xp = cosPhi;

	float xp_over_a = (rCart_over_a + height / WGS84_a) * cosPhi;

	float zp = (rCart * (1. - WGS84_esq) + height) * sinPhi;
	float dphi_zp = (dphi_rCart * (1. - WGS84_esq)) * sinPhi + (rCart * (1. - WGS84_esq) + height) * dphi_sinPhi;
	float dheight_zp = sinPhi;

	float zp_over_a = (rCart_over_a * (1. - WGS84_esq) + height / WGS84_a) * sinPhi;

	float r2D = sqrt(xp * xp + zp * zp);
	float dphi_r2D = (xp * dphi_xp + zp * dphi_zp) / r2D;
	float dheight_r2D = (xp * dheight_xp + zp * dheight_zp) / r2D;

	float r2D_over_a = sqrt(xp_over_a * xp_over_a + zp_over_a * zp_over_a);
	float dphi_r2D_over_a = (xp_over_a * dphi_xp + zp_over_a * dphi_zp) / r2D;

	float sinPhiSph = zp / r2D;
	float dphi_sinPhiSph = (dphi_zp * r2D - zp * dphi_r2D) / (r2D * r2D);
	float dheight_sinPhiSph = (dheight_zp * r2D - zp * dheight_r2D) / (r2D * r2D);

	float cosPhiSph = sqrt(1. - sinPhiSph * sinPhiSph);
	//d/du sqrt(1 - x^2) = -x/sqrt(1 - x^2) dx/du;
	float dphi_cosPhiSph = -sinPhi / cosPhiSph * dphi_sinPhiSph;
	float dheight_cosPhiSph = -sinPhi / cosPhiSph * dheight_sinPhiSph;

	//float x = r2D * cosPhiSph / WGS84_a * cosLambda;
	//float y = r2D * cosPhiSph / WGS84_a * sinLambda;
	//float z = r2D * sinPhiSph / WGS84_a;

	vec3 dphi = vec3(
		(dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * cosLambda,
		(dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * sinLambda,
		(dphi_r2D_over_a * sinPhiSph + r2D_over_a * dphi_sinPhiSph));

	vec3 dlambda = vec3(
		-sinLambda,
		cosLambda,
		0.);

	vec3 dheight = vec3(
		(dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * cosLambda,
		(dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * sinLambda,
		(dheight_r2D * sinPhiSph + r2D * dheight_sinPhiSph));

	return mat3(dphi, dlambda, -dheight);

}

]]
		end

		return c
	end)(),

	Chart:subclass{
		name = 'sphere',
		is3D = true,
		build = function(c)
			c:buildVars{
				{WGS84_a = WGS84_a},
			}
			local rval = heightvar + c.vars.WGS84_a
			local thetaval = symmath.pi/2 - latradval
			c:buildFunc{
				rval * symmath.sin(thetaval) * symmath.cos(lonradval),
				rval * symmath.sin(thetaval) * symmath.sin(lonradval),
				rval * symmath.cos(thetaval),
			}
			function c:chartInv(x, y, z)
				local r = math.sqrt(x*x + y*y + z*z)
				local phi = math.atan2(y, x)
				local r2 = math.sqrt(x*x + y*y)
				local theta = math.atan(r2 / z)
				local height = r - 1
				return
					math.deg(theta),
					(math.deg(phi) + 180) % 360 - 180,
					height
			end
			function c:getGLSLModule()
				local code = c.super.getGLSLModule(self)
				code = code .. [[

//// MODULE_DEPENDS: xformZUpToZBack deg

vec3 chartInv_sphere(vec3 pt) {
	pt = xformZUpToZBack(pt);
	float r = length(pt);
	float lonrad = atan(pt.y, pt.x);
	float latrad = atan(pt.z, length(pt.xy));
	float height = r - WGS84_a;
	return vec3(deg(latrad), deg(lonrad), height);
}
]]
				return code
			end
		end,
	},

	Chart:subclass{
		name = 'cylinder',
		is3D = true,
		build = function(c)
			c:buildVars{
				{WGS84_a = WGS84_a},
			}
			local rval = heightvar + c.vars.WGS84_a
			c:buildFunc{
				rval * symmath.cos(lonradval),
				rval * symmath.sin(lonradval),
				rval * latradval,
			}
			-- TODO c:chartInv
		end,
	},

	-- https://en.wikipedia.org/wiki/Equirectangular_projection
	Chart:subclass{
		name = 'Equirectangular',
		build = function(c)
			c:buildVars{
				{R = 2/math.pi},
				{lambda0 = 0},
				{phi0 = 0},
				{phi1 = 0},
				{WGS84_a = WGS84_a},
			}
			c:buildFunc{
				c.vars.WGS84_a * c.vars.R * (lonradval - c.vars.lambda0) * symmath.cos(c.vars.phi1),
				c.vars.WGS84_a * c.vars.R * (latradval - c.vars.phi0),
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Mercator_projection
	Chart:subclass{
		name = 'Mercator',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(.5)},
				{Miller = 1},		-- set this to 4/5 to get the Miller projection: https://en.wikipedia.org/wiki/Miller_cylindrical_projection
				{WGS84_a = WGS84_a},
			}
			-- TODO
			-- Mercator is defined so that at phi=-pi/2 and phi=pi/2 the y -> infinity
			-- so ... clamp tan's output?
			-- symmath needs a clamp() and a min() and max() functions ...
			-- until then, just scale by 1-eps
			local eps = 1e-2
			c:buildFunc{
				c.vars.WGS84_a * c.vars.R * lonradval,
				c.vars.WGS84_a * c.vars.R / c.vars.Miller * symmath.log(
					eps + symmath.tan((1-eps) * symmath.pi * (90 + latvar * c.vars.Miller) / 360)
				),
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Gall%E2%80%93Peters_projection
	Chart:subclass{
		name = 'Gall-Peters',
		build = function(c)
			c:buildVars{
				{R = .5},
				{WGS84_a = WGS84_a},
			}
			c:buildFunc{
				c.vars.WGS84_a * c.vars.R * lonradval,
				c.vars.WGS84_a * c.vars.R * 2 * symmath.sin(latradval),
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Lambert_cylindrical_equal-area_projection
	Chart:subclass{
		name = 'Lambert cylindrical equal-area',
		build = function(c)
			c:buildVars{
				{lon0 = 0},
				{WGS84_a = WGS84_a},
			}
			c:buildFunc{
				c.vars.WGS84_a * 1/symmath.sqrt(2) * (lonradval - c.vars.lon0 * symmath.pi / 180),
				c.vars.WGS84_a * 1/symmath.sqrt(2) * symmath.sin(latradval),
				heightvar,
			}
		end,
	},

	--[[
	-- https://en.wikipedia.org/wiki/Hobo%E2%80%93Dyer_projection
	-- ... doesn't have formulas ...
	Chart:subclass{
		name = 'Hobo-Dyer',
		build = function(c)
			c:buildVars
		end,
	},
	--]]

	-- https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
	Chart:subclass{
		name = 'Azimuthal equidistant',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(.5)},
				{WGS84_a = WGS84_a},
			}
			local azimuthalval = c.vars.R * (1 - latvar / 90)
			c:buildFunc{
				c.vars.WGS84_a * symmath.sin(lonradval) * azimuthalval,
				c.vars.WGS84_a * -symmath.cos(lonradval) * azimuthalval,
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
	Chart:subclass{
		name = 'Lambert azimuthal equal-area',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(2)},
				{WGS84_a = WGS84_a},
			}
			local colat = 90 - latvar	-- spherical friendly
			local polarR = c.vars.R * symmath.sin(colat * symmath.pi / 360)
			c:buildFunc{
				c.vars.WGS84_a * symmath.sin(lonradval) * polarR,
				c.vars.WGS84_a * -symmath.cos(lonradval) * polarR,
				heightvar,
			}
		end,
	},

	(function()
		local c = Chart:subclass()

		c.name = 'Mollweide'

		c.R = math.pi / 4
		c.lambda0 = 0	-- in degrees
		function c:updateGUI()
			local ig = require 'imgui'
			ig.luatableInputFloat('R', self, 'R')
			ig.luatableInputFloat('lambda0', self, 'lambda0')
		end

		function c:chart(lat, lon, height)
			local lonrad = math.rad(lon)
			local lambda = lonrad
			local latrad = math.rad(lat)
			local phi = latrad
			local theta
			if phi == .5 * math.pi then
				theta = .5 * math.pi
			else
				theta = phi
				for i=1,10 do
					local dtheta = (2 * theta + math.sin(2 * theta) - math.pi * math.sin(phi)) / (2 + 2 * math.cos(theta))
					if math.abs(dtheta) < 1e-5 then break end
					theta = theta - dtheta
				end
			end
			local mollweidex = WGS84_a * self.R * math.sqrt(8) / math.pi * (lambda - self.lambda0) * math.cos(theta)
			local mollweidey = WGS84_a * self.R * math.sqrt(2) * math.sin(theta)
			local mollweidez = height
			if not math.isfinite(mollweidex) then mollweidex = 0 end
			if not math.isfinite(mollweidey) then mollweidey = 0 end
			if not math.isfinite(mollweidez) then mollweidez = 0 end
			return mollweidex, mollweidey, mollweidez
		end

		function c:basis(lat, lon, height)
			return
				vec3d(0, 1, 0),
				vec3d(1, 0, 0),
				vec3d(0, 0, -1)
		end

		function c:getGLSLModule()
			return [[

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
	float x = WGS84_a * Mollweide_R * 2. * M_SQRT_2 / M_PI * (lambda - Mollweide_lambda0) * cos(theta);
	float y = WGS84_a * Mollweide_R * M_SQRT_2 * sin(theta);
	float z = height;
	if (!isfinite(x)) x = 0.;
	if (!isfinite(y)) y = 0.;
	if (!isfinite(z)) z = 0.;
	return vec3(x, y, z);
}

mat3 chart_Mollweide_basis(vec3 latLonHeight) {
	return mat3(
		vec3(0., 1., 0.),
		vec3(1., 0., 0.),
		vec3(0., 0., -1.));
}

]]
		end

		return c
	end)(),

	-- https://en.wikipedia.org/wiki/Sinusoidal_projection
	Chart:subclass{
		name = 'Sinusoidal',
		build = function(c)
			c:buildVars{
				{lon0 = 0},
				{WGS84_a = WGS84_a},
			}
			c:buildFunc{
				c.vars.WGS84_a * (lonvar - c.vars.lon0) * symmath.pi / 360 * symmath.cos(latradval),
				c.vars.WGS84_a * latradval / 2,
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Winkel_tripel_projection
	Chart:subclass{
		name = 'Winkel tripel',
		normalizeBasisNumerically = true,
		build = function(c)
			c:buildVars{
				{lat1 = 0},
				{WGS84_a = WGS84_a},
			}

			-- TODO move this to symmath?
			-- so sinc has to also be defined in glsl
			local sinc = require 'symmath.Function':subclass()
			sinc.name = 'sinc'
			sinc.realFunc = function(x)
				if x == 0 then return 1 end
				return math.sin(x) / x
			end

			function sinc:evaluateDerivative(deriv, ...)
				local x = table.unpack(self):clone()
				return deriv(x, ...) * (-symmath.sin(x) + x * symmath.cos(x)) / x^2
			end

			local alpha = symmath.acos(symmath.cos(latradval) * symmath.cos(lonradval / 2))
			local sincAlpha = sinc(alpha)
			c:buildFunc{
				c.vars.WGS84_a * (
					lonradval * symmath.cos(c.vars.lat1 * symmath.pi / 180)
					+ 2 * symmath.cos(latradval) * symmath.sin(lonradval / 2) / sincAlpha
				) / 4,
				c.vars.WGS84_a * (latradval + symmath.sin(latradval) / sincAlpha) / 4,
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection
	Chart:subclass{
		name = 'Kavrayskiy VIII',
		normalizeBasisNumerically = true,
		build = function(c)
			c:buildVars{
				{WGS84_a = WGS84_a},
			}
			local latnormval = latvar / 180	--[-1,1]
			c:buildFunc{
				c.vars.WGS84_a * lonradval * symmath.sqrt(symmath.frac(1,3) - latnormval * latnormval),
				c.vars.WGS84_a * latradval,
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Wiechel_projection
	Chart:subclass{
		name = 'Wiechel',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(2)},
				{WGS84_a = WGS84_a},
			}
			local coslon = symmath.cos(lonradval)
			local coslat = symmath.cos(latradval)
			local sinlon = symmath.sin(lonradval)
			local sinlat = symmath.sin(latradval)
			c:buildFunc{
				c.vars.WGS84_a * c.vars.R / 2 * (sinlon * coslat - (1 - sinlat) * coslon),
				c.vars.WGS84_a * c.vars.R / 2 * -(coslon * coslat + (1 - sinlat) * sinlon),
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Albers_projection
	Chart:subclass{
		name = 'Albers',
		normalizeBasisNumerically = true,
		build = function(c)
			c:buildVars{
				{R = 1},
				{lat1 = 15},
				{lat2 = 45},
				{lon0 = 0},
				{lat0 = 0},
				{WGS84_a = WGS84_a},
			}
			local lat1rad = c.vars.lat1 * symmath.pi / 180
			local lat2rad = c.vars.lat2 * symmath.pi / 180
			local lon0rad = c.vars.lon0 * symmath.pi / 180
			local lat0rad = c.vars.lat0 * symmath.pi / 180
			local n = (symmath.sin(lat1rad) + symmath.sin(lat2rad)) / 2
			local C = symmath.cos(lat1rad)^2 + 2 * n * symmath.sin(lat1rad)
			local rho0 = c.vars.R / n * symmath.sqrt(C - 2 * n * symmath.sin(lat0rad))
			local theta = n * (lonradval - lon0rad)
			local rho = c.vars.R / n * symmath.sqrt(C - 2 * n * symmath.sin(latradval))
			c:buildFunc{
				c.vars.WGS84_a * rho * symmath.sin(theta),
				c.vars.WGS84_a * (rho0 - rho * symmath.cos(theta)),
				heightvar,
			}
		end,
	},

	-- https://en.wikipedia.org/wiki/Bonne_projection
	Chart:subclass{
		name = 'Bonne',
		skipBasis = true,
		build = function(c)
			c:buildVars{
				{lon0 = 0},
				{lat1 = 45},
				{WGS84_a = WGS84_a},
			}
			local lon0rad = c.vars.lon0 * symmath.pi / 180
			local lat1rad = c.vars.lat1 * symmath.pi / 180
			local rho = 1 / symmath.tan(lat1rad) + lat1rad - latradval
			local E = (lonradval - lon0rad) * symmath.cos(latradval) / rho
			c:buildFunc{
				c.vars.WGS84_a * rho * symmath.sin(E),
				c.vars.WGS84_a * (1 / symmath.tan(lat1rad) - rho * symmath.cos(E)),
				heightvar,
			}
		end,
		basis = function(c, lat, lon, height)
			return vec3d(0,1,0), vec3d(1,0,0), vec3d(0,0,-1)
		end,
	},

	--[[ TODO:
	-- transverse Mercator ... just puts arctic at mid-top and antarctic at mid-bottom
	-- https://en.wikipedia.org/wiki/Behrmann_projection ... formula not in wikipedia
	-- https://en.wikipedia.org/wiki/Tobler_hyperelliptical_projection ... requires an integral of a power
	-- https://en.wikipedia.org/wiki/Robinson_projection ... Robinson ... is an interpolation table-based
	-- https://en.wikipedia.org/wiki/Lambert_equal-area_conic_projection ... Albers generalizes this
	-- https://en.wikipedia.org/wiki/Eckert_IV_projection ... requires root-finding
	-- https://en.wikipedia.org/wiki/Eckert_VI_projection
	-- https://en.wikipedia.org/wiki/Eckert_II_projection
	-- https://en.wikipedia.org/wiki/Eckert-Greifendorff_projection
	-- https://en.wikipedia.org/wiki/Equal_Earth_projection
	-- https://en.wikipedia.org/wiki/Hammer_projection
	-- Stereographic
	-- https://en.wikipedia.org/wiki/Collignon_projection
	-- https://en.wikipedia.org/wiki/Bottomley_projection
	-- https://en.wikipedia.org/wiki/Werner_projection
	-- https://en.wikipedia.org/wiki/Strebe_1995_projection
	-- multiple slices
	-- https://en.wikipedia.org/wiki/Boggs_eumorphic_projection
	-- https://en.wikipedia.org/wiki/Goode_homolosine_projection

	--]]
}
charts.WGS84_a = WGS84_a
for i=1,#charts do
	local c = charts[i]
	charts[c.name] = c
end
return charts
