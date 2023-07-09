**New Features**
- [x] GitHub actions / CI
- [ ] Unreal Engine 5.2 integration
- [ ] Boost, CGAL, and OpenVDB integration
- [ ] C++ Modules integration
- [ ] Stack frame write-out
- [ ] Installation/setup (Readme file)
- [ ] Dependency list (Readme file)

**Bugs**
- [x] Check for possible double Gamma correction between HDR and blur shaders
- [x] Investigate Texture move constructor/assignment bug
- [x] Fix border shrinking bug during off-screen rendering
- [ ] Fix shadows
- [ ] TeXBox multi-line text issue
- [ ] Remove shader's point/spotlight count reliance.
- [ ] Investigate a possible bug causing asymmetry in point light sources
- [ ] Investigate single triangle and cube render weird issues.
- [ ] Fix parallax occlusion bug

**Minor**
- [ ] Fix segmentation fault when exiting glfw context

**Deprecated**
- [x] Off-screen rendering
- [x] HDR
- [x] Bloom
- [x] Update Texture interface to offer more versatility
- [x] Multi-sampled frame buffer attachments
- [x] Update shadow shaders to account for transparency
- [x] ImGui re-integration
- [x] RenderDoc debugger
- [x] Deprecate String
- [x] Implement post-processing switch
- [x] TeX glyph spacing
- [ ] Resurrect line shader
- [ ] Texture coordinate transformation uniforms in shaders
- [ ] Background texture embedding
- [ ] Add GUI switches for post-processing effects
- [ ] Implement copy constructors for all buffer objects
- [ ] Simple visualisation of lights and cameras
- [ ] Toggles/sliders for all visualiser settings
- [ ] Convert LaTeX textures into SDFs for depth effects
- [ ] ModelGroup loading using Assimp
- [ ] Get GDB debugger working in CLion
- [ ] Implement MSAA frame buffers for directional/point shadow maps
- [ ] Motion blur
- [ ] Deferred rendering
- [ ] Screen-space ambient occlusion (SSAO)
- [ ] PBR pipeline
- [ ] Perlin noise generator
- [ ] Curl noise generator
- [ ] Depth of field
- [ ] Vulkan integration
- [ ] Ray-tracing pipeline
- [ ] Wireframe vs solid models
- [ ] Resurrect spotlights
- [ ] Gnuplot writer
- [ ] Pyplot writer
- [ ] Replace cmake with basil