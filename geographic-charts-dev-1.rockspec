package = "geographic-charts"
version = "dev-1"
source = {
	url = "git+https://github.com/thenumbernine/geographic-charts-lua"
}
description = {
	summary = "Geographic Charts",
	detailed = "Geographic Charts",
	homepage = "https://github.com/thenumbernine/geographic-charts-lua",
	license = "MIT",
}
dependencies = {
	"lua >= 5.1"
}
build = {
	type = "builtin",
	modules = {
		["geographic-charts.code"] = "code.lua",
		["geographic-charts"] = "geographic-charts.lua"
	}
}
