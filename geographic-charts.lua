--[[
I use these often enough that they are going to go into a library soon enough

2D charts will be xy coordinates, z will be for height

3D charts  are z+ = north pole, x+ = prime meridian


GAAAHHH STANDARDS
physics spherical coordinates: the longitude is φ and the latitude (starting at the north pole and going down) is θ ...
mathematician spherical coordinates: the longitude is θ and the latitude is φ ...
geographic / map charts: the longitude is λ and the latitude is φ ...
so TODO change the calc_* stuff from r_θφ to h_φλ? idk ...

--]]

-- using geographic labels: lat = φ, lon = λ
local table = require 'ext.table'
local class = require 'ext.class'
local math = require 'ext.math'
local timer = require 'ext.timer'
local vec3d = require 'vec-ffi.vec3d'
local symmath = require 'symmath'

-- input
local latvar = symmath.var'lat'
local lonvar = symmath.var'lon'
local heightvar = symmath.var'height'
-- GLSL-specific, to accept (lat,lon,height) as a vec3
latvar:nameForExporter('C', 'latLonHeight.x')
lonvar:nameForExporter('C', 'latLonHeight.y')
heightvar:nameForExporter('C', 'latLonHeight.z')
-- expressions
local latradval = latvar * symmath.pi / 180
local lonradval = lonvar * symmath.pi / 180

local WGS84_a = 6378137	-- m ... earth equitorial radius
-- TODO global var and remove all WGS84_a's from the var lists (same idea as symmath.pi)
local WGS84_avar = symmath.var'WGS84_a'

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
	self.basisExprs = table{latvar, lonvar, heightvar}:mapi(function(x)
		return y:diff(x)()
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
	return symmath.export.C:toCode{
		input = self.exprIn,
		output = self.exprOut,
	}
end

function Chart:getGLSLFunc()
	local escname = self:getCName()
	return self.varnames:mapi(function(k)
		return 'const float '..escname..'_'..k..' = '..self[k]..';'
	end):append{
		'vec3 chart_'..escname..'(vec3 latLonHeight) {',
		self:getGLSLBody(),
		'	return vec3(out1, out2, out3);',
		'}',
	}:concat'\n'
end

-- the '3D' indicates to be sure to insert a `xformZBackToZUp(pt)` last
function Chart:getGLSLFunc3D()
	local escname = self:getCName()
	return self.varnames:mapi(function(k)
		return 'const float '..escname..'_'..k..' = '..self[k]..';'
	end):append{
		'vec3 chart_'..escname..'(vec3 latLonHeight) {',
		self:getGLSLBody(),
		'	return xformZBackToZUp(vec3(out1, out2, out3));',
		'}',
	}:concat'\n'
end

function Chart:getGLSLModule()
	return table{
		'//// MODULE_NAME: chart_'..self:getCName(),
		'//// MODULE_DEPENDS: M_PI WGS84_a xformZBackToZUp sinc',
		not self.is3D
			and self:getGLSLFunc()
			or self:getGLSLFunc3D(),
	}:concat'\n'
end

function Chart:chart(lat, lon, height)
	return self.chartFunc(lat, lon, height, getfields(self, table.unpack(self.varnames)))
end

-- TODO don't forget that GLSL is swapping the z-back for z-up
-- TODO for a few of these (sphere, cylinder, etc) multiply the output by WGS84_a to put it in meters
function Chart:basis(lat, lon, height)
	local x,y,z = self.basisFunc(lat, lon, height, getfields(self, table.unpack(self.varnames)))
	return vec3d(table.unpack(x)),
			vec3d(table.unpack(y)),
			vec3d(table.unpack(z))
end

function Chart:updateGUI()
	local ig = require 'imgui'
	for _,n in ipairs(self.varnames) do
		ig.luatableInputFloat(n, self, n)
	end
end


local charts = {
	(function()
		local c = class(Chart)

		c.name = 'WGS84'

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

		-- returns x,y,z in meters
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

		return c
	end)(),

	class(Chart, {
		name = 'sphere',
		is3D = true,
		build = function(c)
			c:buildVars()
			local rval = heightvar / WGS84_avar + 1
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
		end,
	}),

	class(Chart, {
		name = 'cylinder',
		is3D = true,
		build = function(c)
			c:buildVars()
			local rval = heightvar / WGS84_avar + 1
			c:buildFunc{
				rval * symmath.cos(lonradval),
				rval * symmath.sin(lonradval),
				rval * latradval,
			}
			-- TODO c:chartInv
		end,
	}),

	-- https://en.wikipedia.org/wiki/Equirectangular_projection
	class(Chart, {
		name = 'Equirectangular',
		build = function(c)
			c:buildVars{
				{R = 2/math.pi},
				{lambda0 = 0},
				{phi0 = 0},
				{phi1 = 0},
			}
			c:buildFunc{
				c.vars.R * (lonradval - c.vars.lambda0) * symmath.cos(c.vars.phi1),
				c.vars.R * (latradval - c.vars.phi0),
				heightvar / WGS84_avar,	-- really tempting to just make all charts -- even the 2D ones -- in meters.
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Mercator_projection
	class(Chart, {
		name = 'Mercator',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(.5)},
			}
			-- TODO
			-- Mercator is defined so that at phi=-pi/2 and phi=pi/2 the y -> infinity
			-- so ... clamp tan's output?
			-- symmath needs a clamp() and a min() and max() functions ...
			-- until then, just scale by 1-eps
			local eps = 1e-2
			c:buildFunc{
				c.vars.R * lonradval,
				c.vars.R * symmath.log(
					eps + symmath.tan((1-eps) * symmath.pi * (90 + latvar) / 360)
				),
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Gall%E2%80%93Peters_projection
	class(Chart, {
		name = 'Gall-Peters',
		build = function(c)
			c:buildVars{
				{R = .5},
			}
			c:buildFunc{
				c.vars.R * lonradval,
				2 * c.vars.R * symmath.sin(latradval),
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Lambert_cylindrical_equal-area_projection
	class(Chart, {
		name = 'Lambert cylindrical equal-area',
		build = function(c)
			c:buildVars{
				{lon0 = 0},
			}
			c:buildFunc{
				1/symmath.sqrt(2) * (lonradval - c.vars.lon0 * symmath.pi / 180),
				1/symmath.sqrt(2) * symmath.sin(latradval),
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
	class(Chart, {
		name = 'Azimuthal equidistant',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(.5)},
			}
			local azimuthalval = c.vars.R * (1 - latvar / 90)
			c:buildFunc{
				symmath.sin(lonradval) * azimuthalval,
				-symmath.cos(lonradval) * azimuthalval,
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
	class(Chart, {
		name = 'Lambert azimuthal equal-area',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(2)},
			}
			local colat = 90 - latvar	-- spherical friendly
			local polarR = c.vars.R * symmath.sin(colat * symmath.pi / 360)
			c:buildFunc{
				symmath.sin(lonradval) * polarR,
				-symmath.cos(lonradval) * polarR,
				heightvar / WGS84_avar,
			}
		end,
	}),

	(function()
		local c = class(Chart)

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
			local mollweidex = self.R * math.sqrt(8) / math.pi * (lambda - self.lambda0) * math.cos(theta)
			local mollweidey = self.R * math.sqrt(2) * math.sin(theta)
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

		return c
	end)(),

	-- https://en.wikipedia.org/wiki/Sinusoidal_projection
	class(Chart, {
		name = 'Sinusoidal',
		build = function(c)
			c:buildVars{
				{lon0 = 0},
			}
			c:buildFunc{
				(lonvar - c.vars.lon0) * symmath.pi / 360 * symmath.cos(latradval),
				latradval / 2,
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Winkel_tripel_projection
	class(Chart, {
		name = 'Winkel tripel',
		normalizeBasisNumerically = true,
		build = function(c)
			c:buildVars{
				{lat1 = 0},
			}

			-- TODO move this to symmath?
			-- so sinc has to also be defined in glsl
			local sinc = class(require 'symmath.Function')
			sinc.name = 'sinc'
			sinc.realFunc = function(x)
				if x == 0 then return 1 end
				return math.sin(x) / x
			end

			local alpha = symmath.acos(symmath.cos(latradval) * symmath.cos(lonradval / 2))
			local sincAlpha = sinc(alpha)
			c:buildFunc{
				(
					lonradval * symmath.cos(c.vars.lat1 * symmath.pi / 180)
					+ 2 * symmath.cos(latradval) * symmath.sin(lonradval / 2) / sincAlpha
				) / 4,
				(latradval + symmath.sin(latradval) / sincAlpha) / 4,
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection
	class(Chart, {
		name = 'Kavrayskiy VIII',
		normalizeBasisNumerically = true,
		build = function(c)
			c:buildVars()
			local latnormval = latvar / 180	--[-1,1]
			c:buildFunc{
				lonradval * symmath.sqrt(symmath.frac(1,3) - latnormval * latnormval),
				latradval,
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Wiechel_projection
	class(Chart, {
		name = 'Wiechel',
		build = function(c)
			c:buildVars{
				{R = math.sqrt(2)},
			}
			local coslon = symmath.cos(lonradval)
			local coslat = symmath.cos(latradval)
			local sinlon = symmath.sin(lonradval)
			local sinlat = symmath.sin(latradval)
			c:buildFunc{
				c.vars.R / 2 * (sinlon * coslat - (1 - sinlat) * coslon),
				-c.vars.R / 2 * (coslon * coslat + (1 - sinlat) * sinlon),
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Albers_projection
	class(Chart, {
		name = 'Albers',
		normalizeBasisNumerically = true,
		build = function(c)
			c:buildVars{
				{R = 1},
				{lat1 = 15},
				{lat2 = 45},
				{lon0 = 0},
				{lat0 = 0},
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
				rho * symmath.sin(theta),
				rho0 - rho * symmath.cos(theta),
				heightvar / WGS84_avar,
			}
		end,
	}),

	-- https://en.wikipedia.org/wiki/Bonne_projection
	class(Chart, {
		name = 'Bonne',
		skipBasis = true,
		build = function(c)
			c:buildVars{
				{lon0 = 0},
				{lat1 = 45},
			}
			local lon0rad = c.vars.lon0 * symmath.pi / 180
			local lat1rad = c.vars.lat1 * symmath.pi / 180
			local rho = 1 / symmath.tan(lat1rad) + lat1rad - latradval
			local E = (lonradval - lon0rad) * symmath.cos(latradval) / rho
			c:buildFunc{
				rho * symmath.sin(E),
				1 / symmath.tan(lat1rad) - rho * symmath.cos(E),
				heightvar / WGS84_avar,
			}
		end,
	}),
}
for i=1,#charts do
	local c = charts[i]
	if c.build then
		timer('building '..c.name, c.build, c)
	end
	charts[c.name] = c
end
return charts
