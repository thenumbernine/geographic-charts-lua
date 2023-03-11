I am using these often enough that I'm putting them in their own library.

Dependencies:
- `lua-ext` - my lua ext library
- `symmath` - my computer algebra system for differentiation for basis calculation.
- `vec-ffi`
- `imgui` (optional) = if you call a chart's "updateGUI" then it will need this.

Returns a list of charts used for mapping.

Each chart contains:
- name
- chart(lat, lon, height) = returns x, y, z in meters
- chartInv(x, y, z) = returns lat, lon, height.  
- basis(lat, lon, height) = returns a set of 3 `vec3d`'s the basis, probably orthonormal.  Lat and lon are in degrees, height is in meters.
- updateGUI() = runs imgui update for modifying variables via gui.
TODO:
- chartGLSL = the code of the chart function, in GLSL, as a module (see modules library).

# Reference

- earth pic is from NASA blue marbel
