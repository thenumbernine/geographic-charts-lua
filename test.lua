#!/usr/bin/env luajit
local ffi = require 'ffi'
local template = require 'template'
local table = require 'ext.table'
local timer = require 'ext.timer'
local gl = require 'gl'
local GLTex2D = require 'gl.tex2d'
local GLProgram = require 'gl.program'
local glreport = require 'gl.report'
local ig = require 'imgui'
local Image = require 'image'

local matrix_ffi = require 'matrix.ffi'
matrix_ffi.real = 'float'	-- default matrix_ffi type

local allChartCode = require 'geographic-charts.code'
-- TODO validate 'charts' fwd and inverse as well.
--local charts = require 'geographic-charts'

local App = require 'imguiapp.withorbit'()

App.title = 'geographic chart demo'
App.viewDist = 1.6
App.viewOrthoSize = 2	-- TODO assign in glapp.view

-- TODO combine this, glapp/tests/info.lua, and seismographic-visualization/vis.lua, and a few others, into a consolidated gl.get function (similar to how the CL getters are already defined in the cl lua library)
local function glget(k)
	local int = ffi.new'int[1]'
	gl.glGetIntegerv(assert(gl[k]), int);
	return int[0]
end

local chartNames = table{
	'sphere',
	'WGS84',
	'cylinder',
	'Equirectangular',
	'Azimuthal_equidistant',
	'Lambert_Azimuthal_equal_area',
	'Mollweide',
}

function App:initGL(...)
	App.super.initGL(self, ...)
	self.view.ortho = true
	self.view.orthoSize = self.viewOrthoSize

	gl.glEnable(gl.GL_DEPTH_TEST)
	gl.glEnable(gl.GL_POINT_SMOOTH)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)

	-- both are too big, max tex size is 16384
	-- and resizing takes too long (and crashes)
	-- so just resize offline
	local image
	timer('loading earth texture', function()
		image = Image'earth-color.png'
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
		minFilter = gl.GL_LINEAR,
		magFilter = gl.GL_LINEAR,
		generateMipmap = true,
	}
glreport'here'
	GLTex2D:unbind()
glreport'here'

	self.modelViewMatrix = matrix_ffi.zeros{4,4}
	self.projectionMatrix = matrix_ffi.zeros{4,4}

	self.globeTexShader = GLProgram{
		vertexCode = template([[
#version 460

<?=allChartCode?>

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

<? for _,name in ipairs(chartNames) do
?>uniform float weight_<?=name?>;
<? end
?>

uniform float zeroLat;
uniform float zeroLon;

in vec3 vertex;
in vec4 color;

out vec4 colorv;
out vec2 texcoordv;

void main() {
	// expect vertex xyz to be lat lon height
	// lat and lon is in degrees
	// height is in meters
	// then generate texcoord etc
	// based on constraints
	vec3 pos = 0.
<? for _,name in ipairs(chartNames) do
?>		+ weight_<?=name?> * chart_<?=name?>(vertex)
<? end
?>	;

	gl_Position = projectionMatrix * (modelViewMatrix * vec4(pos, 1.));
	colorv = color;

	// ok here ... (lat, lon) to sphere
	// then rotate sphere by zeroLat, zeroLon, and maybe a roll too?
	// then back to lat, lon

	float lat = vertex.x + zeroLat;
	float latrad = rad(lat);
	float azimuthal = .5*M_PI - latrad;
	float aziFrac = azimuthal / M_PI;

	float lon = vertex.y + zeroLon;
	float lonrad = rad(lon);
	float lonFrac = lonrad / (2. * M_PI);
	float unitLonFrac = lonFrac + .5;

	texcoordv = vec2(unitLonFrac, aziFrac);
}

]], 	{
			allChartCode = allChartCode,
			chartNames = chartNames,
		}),
		fragmentCode = [[
#version 460
uniform sampler2D colorTex;
in vec4 colorv;
in vec2 texcoordv;
out vec4 fragColor;
void main() {
	fragColor = colorv * texture(colorTex, texcoordv);
}
]],
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



local vars = {
	idivs = 100,
	jdivs = 100,
	normalizeWeights = true,
	zeroLat = 0,
	zeroLon = 0,
}
for _,name in ipairs(chartNames) do
	vars['weight_'..name] = name == 'Equirectangular' and 1 or 0
end

function App:update()
	gl.glClearColor(0, 0, 0, 1)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projectionMatrix.ptr)
	self.globeTexShader:use()
	self.globeTexShader:setUniforms(vars)
	self.globeTexShader:setUniform('modelViewMatrix', self.modelViewMatrix.ptr)
	self.globeTexShader:setUniform('projectionMatrix', self.projectionMatrix.ptr)
	self.colorTex:bind()
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


	App.super.update(self)
glreport'here'
end

local weightFields = chartNames:mapi(function(name)
	return 'weight_'..name
end)

function App:updateGUI()
	if ig.igButton'reset view' then
		self.view.ortho = true
		self.view.orthoSize = self.viewOrthoSize
		self.view.angle:set(0,0,0,1)
		self.view.orbit:set(0,0,0)
		self.view.pos:set(0, 0, self.viewDist)
	end
	ig.luatableInputInt('idivs', vars, 'idivs')
	ig.luatableInputInt('jdivs', vars, 'jdivs')
	ig.luatableSliderFloat('zeroLat', vars, 'zeroLat', -90, 90)
	ig.luatableSliderFloat('zeroLon', vars, 'zeroLon', -180, 180)
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



App():run()
