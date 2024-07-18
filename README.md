# Geographic Charts

[![Donate via Stripe](https://img.shields.io/badge/Donate-Stripe-green.svg)](https://buy.stripe.com/00gbJZ0OdcNs9zi288)<br>

### [Launch](https://thenumbernine.github.io/glapp/?dir=geographic-charts&file=test.lua)

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

Ok I just went and added a bunch of charts to the GLSL code.
Some of the Lua charts are in symmath, and that means I can procedurally generate them for GLSL.
Would be nice to do this with all charts so I can put the chart in one place, generate Lua and GLSL, use diff to generate the basis, and do things like diff generate the Tissot indicators

# Reference

- earth pic is from NASA blue marble https://visibleearth.nasa.gov/collection/1484/blue-marble
