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
local math = require 'ext.math'
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

local charts = {
	(function()
		local c = {}
		
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

	(function()
		local Rvar = symmath.var'R'
		Rvar:nameForExporter('C', 'WGS84_a')	-- this is the GLSL name
		local rval = heightvar / Rvar + 1
		local thetaval = symmath.pi/2 - latradval
		local xval = rval * symmath.sin(thetaval) * symmath.cos(lonradval)
		local yval = rval * symmath.sin(thetaval) * symmath.sin(lonradval)
		local zval = rval * symmath.cos(thetaval)

		local c = {}
		c.name = 'sphere'
		c.R = WGS84_a
		c.exprIn = {latvar, lonvar, heightvar, Rvar}	-- TODO append guivars ... like radius, subsets of lat lon, etc
		c.exprOut = {xval, yval, zval}

		local f = symmath.export.Lua:toFunc{
			input = c.exprIn,
			output = c.exprOut,
		}

		-- TODO this is in z-back 3D coords
		-- it doesn't match code.lua's 2D z-up coords
		function c:chart(lat, lon, height)
			return f(lat, lon, height, self.R)
		end
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
		function c:basis(lat, lon, height)
			-- TODO diff self.exprOut by latvar, lonvar, heightvar
			-- (make sure exprOut is an symmath.Array for vector operations)
			local theta = math.rad(90 - lat)
			local phi = math.rad(lon)
			local r = 1 + height
			return vec3d(
				math.sin(theta) * math.cos(phi),
				math.sin(theta) * math.sin(phi),
				math.cos(theta)
			), vec3d(
				math.cos(theta) * math.cos(phi),
				math.cos(theta) * math.sin(phi),
				-math.sin(theta)
			), vec3d(
				-math.sin(phi),
				math.cos(phi),
				0
			)
		end
		return c
	end)(),
	
	(function()
		local rval = heightvar + 1
		-- results
		local xval = rval * symmath.cos(lonradval)
		local yval = rval * symmath.sin(lonradval)
		local zval = rval * latradval

		local f = symmath.export.Lua:toFunc{
			input = {
				-- args
				latvar, lonvar, heightvar,
				-- gui vars
			},
			output = {xval, yval, zval},
		}

		local c = {}
		c.name = 'cylinder'
		-- TODO multiply this by WGS84_a to put it in meters
		-- but then update in geo-center-earth and seismograph-visualization
		function c:chart(lat, lon, height)
			return f(lat, lon, height)
		end
		function c:basis(lat, lon, height)
			local lonradval = math.rad(lon)
			local c = math.cos(lonradval)
			local s = math.sin(lonradval)
			return 
				vec3d(0,0,1),
				vec3d(-s, c, 0),
				vec3d(c, s, 0)
		end
		-- TODO c:chartInv
		return c
	end)(),

	(function()
		local c = {}
		c.name = 'Equirectangular'
		
		-- gui vars
		c.R = 2 / math.pi
		c.lambda0 = 0
		c.phi0 = 0
		c.phi1 = 0
		-- symmath vars equivalent
		local Rvar = symmath.var'R'
		local lambda0var = symmath.var'lambda0'
		local phi0var = symmath.var'phi0'
		local phi1var = symmath.var'phi1'

		-- results
		local xval = Rvar * (lonradval - lambda0var) * symmath.cos(phi1var)
		local yval = Rvar * (latradval - phi0var)
		local zval = heightvar / WGS84_a

		local f = symmath.export.Lua:toFunc{
			input = {
				-- args
				latvar, lonvar, heightvar,
				-- gui vars
				Rvar, lambda0var, phi0var, phi1var,
			},
			output = {xval, yval, zval},
		}

		function c:chart(lat, lon, height)
			return f(lat, lon, height, 
				self.R, self.lambda0, self.phi0, self.phi1)
		end
		
		function c:updateGUI()
			local ig = require 'imgui'
			ig.luatableInputFloat('R', self, 'R')
			ig.luatableInputFloat('lambda0', self, 'lambda0')
			ig.luatableInputFloat('phi0', self, 'phi0')
			ig.luatableInputFloat('phi1', self, 'phi1')
		end

		function c:basis(lat, lon, height)
			-- Bx is north, By is east, Bz is down ... smh
			return
				vec3d(0, 1, 0),
				vec3d(1, 0, 0),
				vec3d(0, 0, -1)
		end

		return c
	end)(),

	(function()
		local azimuthalval = symmath.pi / 2 - latradval
		
		-- results
		local xval = symmath.sin(lonradval) * azimuthalval
		local yval = -symmath.cos(lonradval) * azimuthalval
		local zval = heightvar
	
		local f = symmath.export.Lua:toFunc{
			input = {
				-- args
				latvar, lonvar, heightvar,
				-- gui vars
			},
			output = {xval, yval, zval},
		}

		local c = {}
		c.name = 'Azimuthal equidistant'
		
		function c:chart(lat, lon, height)
			return f(lat, lon, height)
		end
	
		function c:basis(lat, lon, height)
			local cosLambda = math.cos(math.rad(lon))
			local sinLambda = math.sin(math.rad(lon))
			return
				vec3d(-cosLambda, -sinLambda, 0),	-- d/dphi
				vec3d(-sinLambda, cosLambda, 0),	-- d/dlambda
				vec3d(0, 0, -1)						-- d/dheight
		end

		return c
	end)(),

	(function()
		local c = {}
		
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
}
for i=1,#charts do
	local chart = charts[i]
	charts[chart.name] = chart
end
return charts
