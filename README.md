MotionFFT
=========


Build instructions
------------------


DGtal:

- get DGtal sources (http://github.com/DGtal-team/DGtal)
- once cloned or unziped, create a build folder in the DGtal source tree
- `cd build`
- `cmake .. -DWITH_QGLVIEWER=true`  (you would need libqglviewer to be installed)
- if no complains, complie with `make`
- let DGtalBuildPath be the build folder of DGtal


MotionFFT tools:

- get the code, create a build folder
- `cd build`
- `cmake .. -DDGtal_DIR=<DGtalBuildPath>`
- `make`

You should have the following tools:

- vol2raw: to convert a VOL file to a (unsigned char) raw file (useful
if you want to import it in paraview)
- sliceViewer: 4-pane visualization tool of a VOL file
- testFFT3D: simple test of FFT and reverse FFT in 3D
- createVol: tool to extrude a shape in time following a linear motion
  (with a simple visualization of the space-time volume), and computes
  its FFT




