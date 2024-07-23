-- shorthand for when you want to use all charts, and build them all
-- I'm keeping it separate of the charts themselvse because building can take some time
local timer = require 'ext.timer'
local charts = require 'geographic-charts'
for i=1,#charts do
	local c = charts[i]
	if c.build then
		timer('building '..c.name, c.build, c)
	end
end
return charts
