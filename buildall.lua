-- shorthand for when you want to use all charts, and build them all
-- I'm keeping it separate of the charts themselvse because building can take some time
-- [=[ single-threaded
local timer = require 'ext.timer'
local charts = require 'geographic-charts'
for i=1,#charts do
	local c = charts[i]
	if c.build then
		timer('building '..c.name, c.build, c)
	end
end
return charts
--]=]
--[=[ multithreaded
local template = require 'template'
local pool = require 'thread.pool'{
	initcode = template([[
local timer = require 'ext.timer'
local charts = require 'geographic-charts'
]]),
	code = template([[
local chart = charts[tonumber(task)+1]
if c.build then
	timer('building '..c.name, c.build, c)
end
]]),
	donecode = template([[
results = require 'ext.table'{true}:append(charts)
]]),
}

local charts = require 'geographic-charts'

pool:cycle(#charts)

for _,worker in ipairs(pool) do
	for i=1,#charts do
		-- TODO copy across the derived info
		-- but reassign the matching metatables between states
		local dstchart = charts[i]
		local srcchart = worker.thread.lua[[
	return charts[i]
]]
		dstchart.varlist = srcchart.varlist
		dstchart.vars = {}
		dstchart.varnames = table()
		-- TODO it is very well connected with metatables that are local to the thread Lua state.
	end
end
--]=]
