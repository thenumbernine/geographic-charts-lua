I am using these often enough that I'm putting them in their own library.

Dependencies:
- `lua-ext` - my lua ext library
- `symmath` - my computer algebra system for differentiation for basis calculation.
- `vec-ffi`
- `imgui` (optional) = if you call a chart's "updateGUI" then it will need this.

Returns a list of charts used for mapping.

Each chart contains:
- name
- chart(lat, lon, height) = returns x, y, z in normalized coordinates.  Multiply by charts.WGS84.a
- basis(lat, lon, height) = returns a set of 3 `vec3d`'s the basis, probably orthonormal.  Lat and lon are in degrees, height is in meters.
- updateGUI() = runs imgui update for modifying variables via gui.
