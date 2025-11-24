# NumLi Scene File Format

`.nls` file are a custom scene format defined in `.json` style.

## Allowed Components

Note: `C` is a child object, `x` is a float, `n` is an integer, `[N]` is an array of `N` numbers, 

| `json`        | Key        | Values                                      |
|---------------|------------|---------------------------------------------|
| Camera        | `camera`   | `pos: [3]`<br> `up: [3]`<br> `look_at: [3]`<br> `fov: x`<br> `ar: x`<br> `width: n`
| Materials     | `material` | array of: `"name": C`, with `C`:<br> `type: "type"`<br> `type_params...`
| material types:| `lambertian`<br> `blinn`<br> `microfacet` | `albedo: [3]`<br> `Kd: [3], Ks: [3], Kt: [3], roughness: x, reflect: [3], Le: [3]`<br> `...`
| Objects       | `objects`  | array of `"name": C`, with `C`:<br> `type: "type"`<br> `type_params...`
| object types: | `group`<br>`sphere`<br>`plane`<br>`trimesh` | `transform: C, children: [objects]`<br> `transform: C, material: "mat"`<br> `transform: C, material: "mat"`<br> `transform: C, material: "mat", source: "file_path"`
| Lights        | `lights`   | array of: `"name": C`, with `C`:<br> `type: "type"`<br> `type_params...`
| light types | `ambient`<br> `point`<br> `direction`<br> `sphere`<br> `plane` | `irradiance: [3]`<br> `pos: [3], radiant_intensity: [3]`<br>`dir: [3], radiant_intensity: [3]`<br> `transform: C, radiance: [3]`<br> `transform: C, radiance: [3]`