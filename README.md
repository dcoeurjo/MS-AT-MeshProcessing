Mumford-Shah Mesh Processing using the Ambrosio-Tortorelli Functional
=====================================================================


Prototype source code. This software intensively uses the DEC source code provided by the SIGGRAPH 2013 Course: https://github.com/dgpdec/course



Building
========

Just use `cmake` to build the project:

```
  mkdir build
  cd build
  cmake -D CMAKE_BUILD_TYPE=Release ..
  make
``` 

The project has several dependencies:


- boost and boost_program_options
- OpenGL with GLUT and GLEW for the viewer
- Suitesparse and Umfpack for the linear algebra
- (optional) OpenMP

(tested on Linux and MacOS, not on Windows)

Usage
=====

Typical usage is
``` 
  ./at-positions-mesh -i helmet_noise.obj
``` 

Then, you can type `n` to process the normal vector field using AT and then `p` to reproject the vertices as described in the paper. If you want to iterate, just add an `x` after the `p` to copy the regularized normals at the previous step as input for the next one. For instance, `npxnpx` performs two steps. Then you can export using a `x` and `w`.

Right clicking on the OpenGL window would give you more options. Note that the viewer is really naive. You may have visual issues if the input normal vector is inverted but the exported mesh should be ok.

Options can be obtained using a `-h` flag:

```
-i [ --input ] arg                    input surface.obj
 --exclusion arg (=0)                  enable feature exclusion in normal
                                       regularization (blue channel)
 --normal-inpainting arg (=0)          enable inpainting in normal
                                       regularization (green channel)
 --normal-preprocess arg (=0)          prefilter normal vector field
 --normal-alpha arg (=0.1)             normal dimensionless alpha
 --normal-lambda arg (=0.05)           normal dimensionless lambda
 -s [ --normal-epsilon-start ] arg (=1)
                                       normal dimensionless start epsilon
 -f [ --normal-epsilon-finish ] arg (=0.125)
                                       normal dimensionless finish epsilon
 -p [ --normal-epsilon-progression ] arg (=3)
                                       normal epsilon progression
 --position-inpainting arg (=0)        enable inpainting in position
                                       regularization (green channel)
 --position-alpha arg (=1)             position dimensionless alpha (w_1 in the paper)
 --position-beta arg (=0.05)            position dimensionless beta (w_2 in the paper)
 -k [ --noise ] arg (=0)               Noise level
```


License
=======


The DEC Package has the following license:


``` 
*
* Copyright 2010 Keenan Crane,Fernando de Goes,Mathieu Desbrun, Peter Schroder.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
* SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
* OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
* ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The views and conclusions contained in the software and documentation are those
* of the author and should not be interpreted as representing official policies,
* either expressed or implied, of any other person or institution.
*
``` 
