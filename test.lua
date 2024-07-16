#!/usr/bin/env luajit
local ffi = require 'ffi'
local op = require 'ext.op'
local table = require 'ext.table'
local timer = require 'ext.timer'
local template = require 'template'
local gl = require 'gl'
local GLTex2D = require 'gl.tex2d'
local GLProgram = require 'gl.program'
local glreport = require 'gl.report'
local ig = require 'imgui'
local Image = require 'image'
local vec3d = require 'vec-ffi.vec3d'


local charts = require 'geographic-charts'
local chartCode = require 'geographic-charts.code'
-- TODO validate 'charts' fwd and inverse as well.


local earthtexfn = ... or 'earth-color.png'


local App = require 'imguiapp.withorbit'()

App.title = 'geographic chart demo'
App.viewUseBuiltinMatrixMath = true
App.viewDist = 1.6
App.viewOrthoSize = 2	-- TODO assign in glapp.view

-- TODO combine this, glapp/tests/info.lua, and seismographic-visualization/vis.lua, and a few others, into a consolidated gl.get function (similar to how the CL getters are already defined in the cl lua library)
local function glget(k)
	local int = ffi.new'int[1]'
	gl.glGetIntegerv(assert(gl[k]), int);
	return int[0]
end

-- [[ Winkel tripel uses sinc ... I've got it hooked up to GLSL but not yet to Lua
charts = table.map(charts, function(c,k,t)
	if type(c) == 'table' and c.name == 'Winkel tripel' then return end
	if type(k) == 'number' then return c, #t+1 end
	return c, k
end)
--]]

local chartNames = table.mapi(charts, function(c) return c:getCName() end)

local hasPointSmooth = op.safeindex(gl, 'GL_POINT_SMOOTH')
function App:initGL(...)
	App.super.initGL(self, ...)
	self.view.ortho = true
	self.view.orthoSize = self.viewOrthoSize

	gl.glEnable(gl.GL_DEPTH_TEST)
	if hasPointSmooth then gl.glEnable(gl.GL_POINT_SMOOTH) end	-- not in es
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)

	-- both are too big, max tex size is 16384
	-- and resizing takes too long (and crashes)
	-- so just resize offline
	local image
	timer('loading earth texture', function()
		image = Image(earthtexfn)
	end)
	local maxTextureSize = glget'GL_MAX_TEXTURE_SIZE'
	if image.width > maxTextureSize
	or image.height > maxTextureSize then
		timer('resizing', function()
			image = image:resize(
				math.min(maxTextureSize, image.width),
				math.min(maxTextureSize, image.height)
			)
		end)
	end

glreport'here'
	self.colorTex = GLTex2D{
		image = image,
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		magFilter = gl.GL_LINEAR,
		generateMipmap = true,
	}
glreport'here'
	GLTex2D:unbind()
glreport'here'

	self.globeTexShader = GLProgram{
		version = 'latest',
		header = 'precision highp float;',
		vertexCode = template([[
<?=chartCode?>

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

<? for _,name in ipairs(chartNames) do
?>uniform float weight_<?=name?>;
<? end
?>


in vec3 vertex;
in vec4 color;

out vec4 colorv;

//3D space point that's gonna be used for texturing
// interpolate in 3D so we don't get weird artifats in the texcoord lookup where the texcoords wrap around
out vec3 texcoordptv;

void main() {
	// expect vertex xyz to be lat lon height
	// lat and lon is in degrees
	// height is in meters
	// then generate texcoord etc
	// based on constraints
	//
	// TODO FIXME one of the charts is NaN at the poles ...
	// since this is a linear combination, that's making all charts nan at the poles.
	vec3 pos = 0.
<? for _,name in ipairs(chartNames) do
?>		+ weight_<?=name?> * chart_<?=name?>(vertex)
<? end
?>	;
	
	//from meters to normalized coordinates
	pos /= WGS84_a;

	gl_Position = projectionMatrix * (modelViewMatrix * vec4(pos, 1.));
	colorv = color;

	// ok here ... (lat, lon) to sphere
	texcoordptv = chart_sphere(vertex);
}

]], 	{
			chartCode = chartCode,
			chartNames = chartNames,
		}),
		fragmentCode = template([[
<?=chartCode?>

uniform sampler2D colorTex;

uniform float zeroRoll;
uniform float zeroLon;
uniform float zeroLat;

in vec4 colorv;
in vec3 texcoordptv;

out vec4 fragColor;

vec2 rot2D(vec2 x, float theta) {
	float cth = cos(theta);
	float sth = sin(theta);
	return vec2(
		cth * x.x - sth * x.y,
		sth * x.x + cth * x.y
	);
}

void main() {

	// then rotate sphere by zeroRoll, zeroLon, and maybe a roll too?
	// then back to lat, lon
	// idk what i'm doing
	vec3 pt = texcoordptv;
	pt.xy = rot2D(pt.xy, rad(zeroRoll));
	pt.yz = rot2D(pt.yz, rad(zeroLon));
	pt.xz = rot2D(pt.xz, rad(zeroLat));

	// convert from 3D to lat,lon,heigth
	vec3 invpt = chartInv_sphere(pt);

	// convert from lat,lon in degrees to unit texture coordinates
	float lat = invpt.x;
	float aziFrac = (90. - lat) / 180.;

	float lon = invpt.y;
	float unitLonFrac = (lon + 180.) / 360.;

	vec2 texcoordv = vec2(unitLonFrac, aziFrac);
	fragColor = colorv * texture(colorTex, texcoordv);
}
]],		{
			chartCode = chartCode,
		}),
		uniforms = {
			colorTex = 0,
		},
	}
glreport'here'
	self.globeTexShader:useNone()
glreport'here'
	GLTex2D:unbind()
glreport'here'
end


local weightFields = chartNames:mapi(function(name)
	return 'weight_'..name
end)

local vars = {
	idivs = 100,
	jdivs = 100,
	normalizeWeights = true,
	filterNearest = true,
	zeroRoll = 0,
	zeroLon = 0,
	zeroLat = 0,
	pick = false,
	pickLon = 0,
	pickLat = 0,
}
for _,field in ipairs(weightFields) do
	vars[field] = field == 'weight_Equirectangular' and 1 or 0
end

local function fix3D(v)
	v.y, v.z = v.z, -v.y
	v.x, v.z = -v.z, v.x
end

function App:update()
	gl.glClearColor(0, 0, 0, 1)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	self.globeTexShader:use()
	self.globeTexShader:setUniforms(vars)
	self.globeTexShader:setUniform('modelViewMatrix', self.view.mvMat.ptr)
	self.globeTexShader:setUniform('projectionMatrix', self.view.projMat.ptr)
	self.colorTex:bind()
	if not vars.filterNearest then
		self.colorTex
			:setParameter(gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST)
			:setParameter(gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST)
	else
		self.colorTex
			:setParameter(gl.GL_TEXTURE_MIN_FILTER, gl.GL_LINEAR_MIPMAP_LINEAR)
			:setParameter(gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
	end
	gl.glVertexAttrib4f(self.globeTexShader.attrs.color.loc, 1, 1, 1, 1)
	for j=0,vars.jdivs-1 do
		gl.glBegin(gl.GL_TRIANGLE_STRIP)
		for i=0,vars.idivs do
			local aziFrac = i/vars.idivs
			local azimuthal = aziFrac * math.pi-- azimuthal angle
			local latrad = .5*math.pi - azimuthal	-- latitude
			local lat = math.deg(latrad)

			local unitLonFrac = (j+1)/vars.jdivs
			local lonFrac = unitLonFrac - .5
			local lonrad = lonFrac * 2 * math.pi			-- longitude
			local lon = math.deg(lonrad)
			gl.glVertexAttrib3f(self.globeTexShader.attrs.vertex.loc, lat, lon, 0)

			local unitLonFrac = j/vars.jdivs
			local lonFrac = unitLonFrac - .5
			local lonrad = lonFrac * 2 * math.pi			-- longitude
			local lon = math.deg(lonrad)
			gl.glVertexAttrib3f(self.globeTexShader.attrs.vertex.loc, lat, lon, 0)
		end
		gl.glEnd()
	end
	self.colorTex:unbind()
	self.globeTexShader:useNone()

	do
		local pickHeight = 0
		local colors = table{
			{1,0,0},
			{0,1,0},
			{0,0,1},
		}
		local pos = vec3d(0,0,0)
		local bs = {
			vec3d(0,0,0),
			vec3d(0,0,0),
			vec3d(0,0,0),
		}
		for i=1,#chartNames do
			local c = charts[i]
			local w = vars[weightFields[i]]
			local chartPos = vec3d(c:chart(vars.pickLat, vars.pickLon, pickHeight))
			local cbx, cby, cbz = c:basis(vars.pickLat, vars.pickLon, pickHeight)
			if c.is3D then	-- xformZBackToZUp: same transform as the GLSL 3D -> 2D
				fix3D(chartPos)
				fix3D(cbx)
				fix3D(cby)
				fix3D(cbz)
			end
			pos = pos + chartPos * w
			bs[1] = bs[1] + cbx * w
			bs[2] = bs[2] + cby * w
			bs[3] = bs[3] + cbz * w
		end
		
		for i=1,3 do
			bs[i] = bs[i]:normalize()
		end
		pos = pos + bs[3] * 1e-3	--(1 / charts.WGS84_a)

		local len = .1
		gl.glBegin(gl.GL_LINES)
		for i,b in ipairs(bs) do
			gl.glColor3f(table.unpack(colors[i]))
			gl.glVertex3f(pos.x, pos.y, pos.z)
			gl.glVertex3f(pos.x+b.x*len, pos.y+b.y*len, pos.z+b.z*len)
		end
		gl.glEnd()
	end

	App.super.update(self)
glreport'here'
end

function App:updateGUI()
	ig.igText'view:'
	if ig.igButton'reset view' then
		self.view.ortho = true
		self.view.orthoSize = self.viewOrthoSize
		self.view.angle:set(0,0,0,1)
		self.view.orbit:set(0,0,0)
		self.view.pos:set(0, 0, self.viewDist)
	end
	ig.luatableCheckbox('ortho', self.view, 'ortho')
	ig.igText'resolution:'
	ig.luatableInputInt('idivs', vars, 'idivs')
	ig.luatableInputInt('jdivs', vars, 'jdivs')
	ig.luatableCheckbox('filterNearest', vars, 'filterNearest')
	ig.igText'orientation:'
	ig.luatableSliderFloat('zeroLat', vars, 'zeroLat', -180, 180)
	ig.luatableSliderFloat('zeroLon', vars, 'zeroLon', -180, 180)
	ig.luatableSliderFloat('zeroRoll', vars, 'zeroRoll', -180, 180)
	ig.igText'pick:'
	ig.luatableCheckbox('pick', vars, 'pick')
	ig.luatableInputFloat('pickLon', vars, 'pickLon')
	ig.luatableInputFloat('pickLat', vars, 'pickLat')
	ig.igText'weights:'
	ig.luatableCheckbox('normalize weights', vars, 'normalizeWeights')
	local changed
	for _,field in ipairs(weightFields) do
		if ig.luatableSliderFloat(field, vars, field, 0, 1) then
			changed = field
		end
	end
	if vars.normalizeWeights and changed then
		local restFrac = 1 - vars[changed]
		local totalRest = 0
		for _,field in ipairs(weightFields) do
			if field ~= changed then
				totalRest = totalRest + vars[field]
			end
		end
		for _,field in ipairs(weightFields) do
			if field ~= changed then
				if totalRest == 0 then
					vars[field] = 0
				else
					vars[field] = restFrac * vars[field] / totalRest
				end
			end
		end
	end
end

return App():run()
