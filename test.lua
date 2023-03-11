#!/usr/bin/env luajit
local ffi = require 'ffi'
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

local charts = require 'geographic-charts'
local allChartCode = require 'geographic-charts.code'
local wgs84 = charts.WGS84

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
		vertexCode = table{
'#version 460',
allChartCode,
[[
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

uniform float weight_WGS84;
uniform float weight_cylinder;
uniform float weight_Equirectangular;
uniform float weight_Azimuthal_equidistant;
uniform float weight_Mollweide;

in vec3 vertex;
in vec4 color;

out vec4 colorv;
out vec2 texcoordv;

void main() {
	// expect vertex xyz to be lat lon height
	// then generate texcoord etc
	// based on constraints
	vec3 pos = weight_WGS84 * chart_WGS84(vertex)
			+ weight_cylinder * chart_cylinder(vertex)
			+ weight_Equirectangular * chart_Equirectangular(vertex)
			+ weight_Azimuthal_equidistant * chart_Azimuthal_equidistant(vertex)
			+ weight_Mollweide * chart_Mollweide(vertex);

	gl_Position = projectionMatrix * (modelViewMatrix * vec4(pos, 1.));
	colorv = color;

	float lat = vertex.x;
	float latrad = rad(lat);
	float azimuthal = .5*M_PI - latrad;
	float aziFrac = azimuthal / M_PI;

	float lon = vertex.y;
	float lonrad = rad(lon);
	float lonFrac = lonrad / (2. * M_PI);
	float unitLonFrac = lonFrac + .5;

	texcoordv = vec2(unitLonFrac, aziFrac);
}
]]
}:concat'\n',
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



idivs = 100
jdivs = 100
normalizeWeights = true
spheroidCoeff = 0
cylCoeff = 0
equirectCoeff = 1
aziequiCoeff = 0
mollweideCoeff = 0

function App:update()
	gl.glClearColor(0, 0, 0, 1)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projectionMatrix.ptr)
	self.globeTexShader:use()
	self.globeTexShader:setUniforms{
		weight_WGS84 = spheroidCoeff,
		weight_cylinder = cylCoeff,
		weight_Equirectangular = equirectCoeff,
		weight_Azimuthal_equidistant = aziequiCoeff,
		weight_Mollweide = mollweideCoeff,
		modelViewMatrix = self.modelViewMatrix.ptr,
		projectionMatrix = self.projectionMatrix.ptr,
	}
	self.colorTex:bind()
	gl.glVertexAttrib4f(self.globeTexShader.attrs.color.loc, 1, 1, 1, 1)
	for j=0,jdivs-1 do
		gl.glBegin(gl.GL_TRIANGLE_STRIP)
		for i=0,idivs do
			local aziFrac = i/idivs
			local azimuthal = aziFrac * math.pi-- azimuthal angle
			local latrad = .5*math.pi - azimuthal	-- latitude
			local lat = math.deg(latrad)

			local unitLonFrac = (j+1)/jdivs
			local lonFrac = unitLonFrac - .5
			local lonrad = lonFrac * 2 * math.pi			-- longitude
			local lon = math.deg(lonrad)
			gl.glVertexAttrib3f(self.globeTexShader.attrs.vertex.loc, lat, lon, 0)

			local unitLonFrac = j/jdivs
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



local weightFields = {
	'spheroidCoeff',
	'cylCoeff',
	'equirectCoeff',
	'aziequiCoeff',
	'mollweideCoeff',
}

function App:updateGUI()
	if ig.igButton'reset view' then
		self.view.ortho = true
		self.view.orthoSize = self.viewOrthoSize
		self.view.angle:set(0,0,0,1)
		self.view.orbit:set(0,0,0)
		self.view.pos:set(0, 0, self.viewDist)
	end
	ig.luatableInputInt('idivs', _G, 'idivs')
	ig.luatableInputInt('jdivs', _G, 'jdivs')
	ig.luatableCheckbox('normalize weights', _G, 'normalizeWeights')
	local changed
	for _,field in ipairs(weightFields) do
		if ig.luatableSliderFloat(field, _G, field, 0, 1) then
			changed = field
		end
	end
	if normalizeWeights and changed then
		local restFrac = 1 - _G[changed]
		local totalRest = 0
		for _,field in ipairs(weightFields) do
			if field ~= changed then
				totalRest = totalRest + _G[field]
			end
		end
		for _,field in ipairs(weightFields) do
			if field ~= changed then
				if totalRest == 0 then
					_G[field] = 0
				else
					_G[field] = restFrac * _G[field] / totalRest
				end
			end
		end
	end
end



App():run()
