/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

#include "Surface.h"

#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>

#include <QOpenGLContext> // must be included before OSG headers

#include <osg/ref_ptr>
#include <osg/Array>
#include <osg/Geometry>
#include <osg/LightModel>
#include <osg/LineWidth>
#include <osg/Point>
#include <osg/PolygonMode>
#include <osg/PrimitiveRestartIndex>
#include <osg/PrimitiveSet>
#include <osg/StateAttribute>
#include <osg/StateSet>

#include "Modeling/MessagesWidget.h"
#include "Util/Helper.h"

std::underlying_type<SurfaceNormalsAverageWeights>::type operator&(const SurfaceNormalsAverageWeights& lhs, const SurfaceNormalsAverageWeights& rhs)
{
  return static_cast<std::underlying_type<SurfaceNormalsAverageWeights>::type>(lhs) & static_cast<std::underlying_type<SurfaceNormalsAverageWeights>::type>(rhs);
}

std::underlying_type<SurfaceNormalsAnimationTypes>::type operator&(const SurfaceNormalsAnimationTypes& lhs, const SurfaceNormalsAnimationTypes& rhs)
{
  return static_cast<std::underlying_type<SurfaceNormalsAnimationTypes>::type>(lhs) & static_cast<std::underlying_type<SurfaceNormalsAnimationTypes>::type>(rhs);
}

template<>
AbstractVisualProperties::Transparency::Type VisualProperties<SurfaceObject>::Transparency::getProperty() const
{
  const SurfaceObject* surface = static_cast<const SurfaceObject*>(mpParent);
  return surface->_transparency.exp;
}

SurfaceObject::SurfaceObject()
    : AbstractVisualizerObjectWithVisualProperties(VisualizerType::surface),
      mClosenessCheckState(SurfaceClosenessCheckState::active),
      mStripsWrappingMethod(SurfaceStripsWrappingMethod::degenerate),
      mNormalsAverageWeights(SurfaceNormalsAverageWeights::bothAreaAndAngle),
      mNormalsAnimationTypes(SurfaceNormalsAnimationTypes::none),
      mPointSize(1.0),
      mLineWidth(1.0),
      mNormalScale(1.0),
      _nu(VisualizerAttribute(0.0)),
      _nv(VisualizerAttribute(0.0)),
      _wireframe(VisualizerAttribute(0.0)),
      _normalized(VisualizerAttribute(0.0)),
      _doublesided(VisualizerAttribute(1.0)),
      _multicolored(VisualizerAttribute(0.0)),
      _transparency(VisualizerAttribute(0.0))
{
}

void SurfaceObject::dumpVisualizerAttributes()
{
  AbstractVisualizerObjectWithVisualProperties::dumpVisualizerAttributes();
  std::cout << "nu " << _nu.getValueString() << std::endl;
  std::cout << "nv " << _nv.getValueString() << std::endl;
  std::cout << "wireframe " << _wireframe.getValueString() << std::endl;
  std::cout << "normalized " << _normalized.getValueString() << std::endl;
  std::cout << "doublesided " << _doublesided.getValueString() << std::endl;
  std::cout << "multicolored " << _multicolored.getValueString() << std::endl;
  std::cout << "transparency " << _transparency.getValueString() << std::endl;
}

void SurfaceObject::fakeTorus(const itype nu, const itype nv,
                              ftype* Vx, ftype* Vy, ftype* Vz,
                              ftype* Nx, ftype* Ny, ftype* Nz,
                              ftype* Cr, ftype* Cg, ftype* Cb) const
{
  constexpr ftype pi = M_PI;

  constexpr ftype R = 1; // Major radius (distance from center of torus to center of tube)
  constexpr ftype r = 0.2; // Minor radius (radius of tube)
  constexpr ftype opening = 0; // Opening angle of torus
  constexpr ftype startAngle = -pi; // Start angle of torus slice
  constexpr ftype  stopAngle = +pi; // Stop  angle of torus slice

  constexpr ftype phi_start = -pi + opening;
  constexpr ftype phi_stop  = +pi - opening;

  const bool normalized = _normalized.exp;
  const bool multicolored = _multicolored.exp;

  ftype alpha;
  ftype beta;

  for (itype u = 0; u < nu; u++) {
    alpha = startAngle + (nu > 1 ? (stopAngle - startAngle) * u / (nu - 1) : 0);
    for (itype v = 0; v < nv; v++) {
      beta = phi_start + (nv > 1 ? (phi_stop - phi_start) * v / (nv - 1) : 0);
      {
        Vz[nv * u + v] = (R + r * std::cos(beta)) * std::cos(alpha);
        Vx[nv * u + v] = (R + r * std::cos(beta)) * std::sin(alpha);
        Vy[nv * u + v] =      r * std::sin(beta);
      }
      if (normalized) {
        Nz[nv * u + v] = std::cos(beta) * std::cos(alpha);
        Nx[nv * u + v] = std::cos(beta) * std::sin(alpha);
        Ny[nv * u + v] = std::sin(beta);
      }
      if (multicolored) {
        Cr[nv * u + v] = 0;
        Cg[nv * u + v] = 1;
        Cb[nv * u + v] = 0;
      }
    }
  }

  if (mClosenessCheckState == SurfaceClosenessCheckState::active) {
    if (phi_stop - pi == phi_start + pi) {
      for (itype u = 0; u < nu; u++) {
        {
          Vx[nv * u + (nv - 1)] = Vx[nv * u + 0];
          Vy[nv * u + (nv - 1)] = Vy[nv * u + 0];
          Vz[nv * u + (nv - 1)] = Vz[nv * u + 0];
        }
        if (normalized) {
          Nx[nv * u + (nv - 1)] = Nx[nv * u + 0];
          Ny[nv * u + (nv - 1)] = Ny[nv * u + 0];
          Nz[nv * u + (nv - 1)] = Nz[nv * u + 0];
        }
      }
    }
    if (stopAngle - pi == startAngle + pi) {
      for (itype v = 0; v < nv; v++) {
        {
          Vx[nv * (nu - 1) + v] = Vx[nv * 0 + v];
          Vy[nv * (nu - 1) + v] = Vy[nv * 0 + v];
          Vz[nv * (nu - 1) + v] = Vz[nv * 0 + v];
        }
        if (normalized) {
          Nx[nv * (nu - 1) + v] = Nx[nv * 0 + v];
          Ny[nv * (nu - 1) + v] = Ny[nv * 0 + v];
          Nz[nv * (nu - 1) + v] = Nz[nv * 0 + v];
        }
      }
    }
  }
}

void SurfaceObject::fakeRectangularBox(const itype nu, const itype nv,
                                       ftype* Vx, ftype* Vy, ftype* Vz,
                                       ftype* Nx, ftype* Ny, ftype* Nz,
                                       ftype* Cr, ftype* Cg, ftype* Cb) const
{
  /* Inspired by https://stackoverflow.com/a/65961094 */
  // Dimensions
  const itype W = nu - 1; // Width
  const itype H = nu - 1; // Height
  const itype L = nu > 1 ? (nv == nu ? 0 : (nv - nu * 3) / 4 - 1) : nv - 1; // Length
  // Conditions
  const itype Wp1 = W + 1;
  const itype Hp1 = H + 1;
  const itype Lp1 = L + 1;
  assert(nu > 0);
  assert(nv > 0);
  if (H == 0) {
    assert(Wp1 == Hp1);
    assert(Hp1 == nu);
    assert(Lp1 == nv);
  } else if (L == 0) {
    assert(nv == nu);
    assert(Wp1 == Hp1);
    assert(Hp1 == nu);
    assert(Lp1 == 1);
  } else {
    assert(nv >= nu * 3 + 4);
    assert(Wp1 == Hp1);
    assert(Hp1 == nu);
    assert(Lp1 * 4 + Wp1 * 2 + Hp1 * 1 == nv);
  }
  // Parameters
  const ftype xscale = 1;
  const ftype yscale = 1;
  const ftype zscale = 1;
  const ftype xmin = -xscale * L / 2;
  const ftype ymin = -yscale * H / 2;
  const ftype zmin = -zscale * W / 2;
  const ftype colors[6][3] = {/*-Z*/{1, 0, 0}, /*-Y*/{0, 0, 1}, /*+Z*/{0, 1, 0}, /*+Y*/{1, 1, 0}, /*-X*/{1, 0, 1}, /*+X*/{0, 1, 1}};
  // Attributes
  const ftype colorR = _color[0].exp;
  const ftype colorG = _color[1].exp;
  const ftype colorB = _color[2].exp;
  const bool normalized = _normalized.exp;
  const bool multicolored = _multicolored.exp;
  // Vertices
  if (H == 0) {
    // Degenerate to a line strip or a single point
    {
      const itype u = 0;
      itype v = 0;
      for (itype l = 0; l <= L; l++, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * l;
          Vy[nv * u + v] = ymin + yscale * 0;
          Vz[nv * u + v] = zmin + zscale * 0;
        }
        if (normalized) {
          Nx[nv * u + v] = 0;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = -1;
        }
        if (multicolored) {
          Cr[nv * u + v] = colorR;
          Cg[nv * u + v] = colorG;
          Cb[nv * u + v] = colorB;
        }
      }
    }
  } else if (L == 0) {
    // Degenerate to a single face
    for (itype h = 0; h <= H; h++) {
      const itype u = h;
      itype v = 0;
      for (itype w = 0; w <= W; w++, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * 0;
          Vy[nv * u + v] = ymin + yscale * h;
          Vz[nv * u + v] = zmin + zscale * w;
        }
        if (normalized) {
          Nx[nv * u + v] = +1;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = 0;
        }
        if (multicolored) {
          Cr[nv * u + v] = colorR;
          Cg[nv * u + v] = colorG;
          Cb[nv * u + v] = colorB;
        }
      }
    }
  } else {
    // Loop over substrips
    for (itype h = 0; h <= H; h++) {
      const itype w = W - h;
      const itype u = h;
      itype v = 0;
      // -Z back face
      for (itype l = 0; l <= L; l++, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * l;
          Vy[nv * u + v] = ymin + yscale * h;
          Vz[nv * u + v] = zmin + zscale * 0;
        }
        if (normalized) {
          Nx[nv * u + v] = 0;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = -1;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[0][0];
          Cg[nv * u + v] = colors[0][1];
          Cb[nv * u + v] = colors[0][2];
        }
      }
      // +X right face above half
      for (itype w = 0; w <= W; w++, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * L;
          Vy[nv * u + v] = ymin + yscale * (H - h > w ? w + h : H);
          Vz[nv * u + v] = zmin + zscale * (W - h > w ? w : W - h);
        }
        if (normalized) {
          Nx[nv * u + v] = +1;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = 0;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[5][0];
          Cg[nv * u + v] = colors[5][1];
          Cb[nv * u + v] = colors[5][2];
        }
      }
      // +Y top face
      for (itype l = L; l >= 0; l--, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * l;
          Vy[nv * u + v] = ymin + yscale * H;
          Vz[nv * u + v] = zmin + zscale * w;
        }
        if (normalized) {
          Nx[nv * u + v] = 0;
          Ny[nv * u + v] = +1;
          Nz[nv * u + v] = 0;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[3][0];
          Cg[nv * u + v] = colors[3][1];
          Cb[nv * u + v] = colors[3][2];
        }
      }
      // -X left face
      for (itype h = H; h >= 0; h--, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * 0;
          Vy[nv * u + v] = ymin + yscale * h;
          Vz[nv * u + v] = zmin + zscale * w;
        }
        if (normalized) {
          Nx[nv * u + v] = -1;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = 0;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[4][0];
          Cg[nv * u + v] = colors[4][1];
          Cb[nv * u + v] = colors[4][2];
        }
      }
      // -Y bottom face
      for (itype l = 0; l <= L; l++, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * l;
          Vy[nv * u + v] = ymin + yscale * 0;
          Vz[nv * u + v] = zmin + zscale * w;
        }
        if (normalized) {
          Nx[nv * u + v] = 0;
          Ny[nv * u + v] = -1;
          Nz[nv * u + v] = 0;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[1][0];
          Cg[nv * u + v] = colors[1][1];
          Cb[nv * u + v] = colors[1][2];
        }
      }
      // +X right face below half
      for (itype w = 0; w <= W; w++, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * L;
          Vy[nv * u + v] = ymin + yscale * (H - h < w ? w + h - H : 0);
          Vz[nv * u + v] = zmin + zscale * (W - h < w ? w : W - h);
        }
        if (normalized) {
          Nx[nv * u + v] = +1;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = 0;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[5][0];
          Cg[nv * u + v] = colors[5][1];
          Cb[nv * u + v] = colors[5][2];
        }
      }
      // +Z front face
      for (itype l = L; l >= 0; l--, v++) {
        {
          Vx[nv * u + v] = xmin + xscale * l;
          Vy[nv * u + v] = ymin + yscale * h;
          Vz[nv * u + v] = zmin + zscale * W;
        }
        if (normalized) {
          Nx[nv * u + v] = 0;
          Ny[nv * u + v] = 0;
          Nz[nv * u + v] = +1;
        }
        if (multicolored) {
          Cr[nv * u + v] = colors[2][0];
          Cg[nv * u + v] = colors[2][1];
          Cb[nv * u + v] = colors[2][2];
        }
      }
    }
  }
}

void SurfaceObject::fakeSphericalArc(const itype nu, const itype nv,
                                     ftype* Vx, ftype* Vy, ftype* Vz,
                                     ftype* Nx, ftype* Ny, ftype* Nz,
                                     ftype* Cr, ftype* Cg, ftype* Cb) const
{
  fakeRectangularBox(nu, nv, Vx, Vy, Vz, Nx, Ny, Nz, Cr, Cg, Cb); // TODO: Implement
}

/**
 * @brief Draw the geometry representing this surface.
 *
 * @details Surface geometry // TODO: Document
 *
 * constexpr itype i = 0; // index of first adjacent facet for top right vertex
 * constexpr itype j = 2; // index of first adjacent facet for top left vertex
 * constexpr itype k = 4; // index of first adjacent facet for bottom right vertex
 * constexpr itype l = 1; // index of second adjacent facet for top left vertex
 * constexpr itype m = 3; // index of second adjacent facet for bottom left vertex
 * constexpr itype n = 5; // index of second adjacent facet for bottom right vertex
 *
 * const ftype area1 = length1 / 2; // surface area = triangle area = half the norm of the cross product
 *
 * const ftype angle11 = std::acos(length11 > 0 ? dot11 / length11 : 0); // corner angle = angle of the corner of the polygon at the vertex
 *
 * @note @c nutnv = nu * nv never overflows because @c nPoints < @c nFacetsTripled, i.e. nu * nv < 3 * 2 * (nu - 1) * (nv - 1),
 * and one/both of them is/are checked for overflow of the same integer type.
 *
 * @note @c itype, @c ftype, @c ltype typedefs are not arbitrary, they must satisfy the following constraints:
 * typedef int itype; // Unsigned integers use modulo arithmetic (for wrapping) which is slower than undefined behavior of under/overflow with signed integers,
 * and anyway the last bit is not needed in practice especially because Modelica Integer type maps to a C signed int type
 * but also it must fit in OSG & OpenGL array types (i.e., conservatively, int)
 * typedef double ftype; // Should be double since Modelica Real primitive type is equivalent to a C double primitive type
 * typedef std::size_t ltype; // Max size is size of address (i.e., pointer size)
 *
 * @note @c Vec2Array, @c Vec3Array, @c Vec4Array are in single precision instead of double precision
 * because @ref osgUtil::LineSegmentIntersector segfaults with @ref osg::Vec3dArray for vertices (it tries to cast it to @ref osg::Vec3Array)
 * https://github.com/openscenegraph/OpenSceneGraph/commit/6e1866ac1857d3466a236a8ad63f91be39e19c71#diff-e5c13d448896b8697db1363ad34ad746a07e66512af001c8d27b087e334281b4L444
 * It still does in latest version, however it might be something else around this cast that leads to a segmentation fault
 * https://github.com/openscenegraph/OpenSceneGraph/blob/master/src/osgUtil/LineSegmentIntersector.cpp#L530
 *
//
List of relevant references:

https://subscription.packtpub.com/book/game-development/9781849512824/4/ch04lvl1sec09/indexing-primitives Book OpenSceneGraph 3.0: Beginner's Guide
https://books.google.ch/books?id=UaRkiFd_yokC Same for free!
http://olmozavala.com/Custom/OpenGL/Tutorials/ProyectoHuracanOpenSceneGraph/Documentos/Documentos_Curso/OSGQSG_Martz.pdf Quick Start Guide for free!
https://weber.itn.liu.se/~karlu20/courses/TNM086-2021/labs/models/openscenegraph_quick_start_guide.pdf Quick Start Guide for free! (3rd revision)

https://github.com/openscenegraph/OpenSceneGraph/blob/34a1d8bc9bba5c415c4ff590b3ea5229fa876ba8/src/osgPlugins/obj/ReaderWriterOBJ.cpp#L528-L771
https://www.khronos.org/opengl/wiki/Vertex_Rendering#Primitive_Restart
https://github.com/openscenegraph/OpenSceneGraph/blob/34a1d8bc9bba5c415c4ff590b3ea5229fa876ba8/include/osg/PrimitiveRestartIndex#L22
https://github.com/openscenegraph/OpenSceneGraph/blob/34a1d8bc9bba5c415c4ff590b3ea5229fa876ba8/src/osg/PrimitiveRestartIndex.cpp#L63
https://github.com/openscenegraph/OpenSceneGraph/tree/master/examples/osgdepthpeeling
https://github.com/openscenegraph/OpenSceneGraph/tree/master/examples/osgoit

https://github.com/openscenegraph/OpenSceneGraph/blob/master/examples/osggeometry/osggeometry.cpp#L139-L141
https://github.com/openscenegraph/OpenSceneGraph/blob/master/src/osg/Shape.cpp#L287

https://learnopengl.com/Getting-started/Hello-Triangle
https://www.glprogramming.com/red/appendixe.html#name2
https://www.glprogramming.com/red/chapter07.html#name2
https://www.glprogramming.com/red/chapter12.html
https://www.glprogramming.com/red/chapter02.html
//

//
See https://www.glprogramming.com/red/chapter07.html#name2 for an example torus drawn with GL_QUAD_STRIP
//

//
See https://www.glprogramming.com/red/chapter02.html#name2 for examples of invalid quadrilaterals
and https://www.glprogramming.com/red/chapter02.html#name8 for some hints to approximate surfaces
"Since OpenGL vertices are always three-dimensional, the points forming the boundary of a particular polygon don't necessarily lie on the same plane in space. (Of course, they do in many cases - if all the z coordinates are zero, for example, or if the polygon is a triangle.) If a polygon's vertices don't lie in the same plane, then after various rotations in space, changes in the viewpoint, and projection onto the display screen, the points might no longer form a simple convex polygon. For example, imagine a four-point quadrilateral where the points are slightly out of plane, and look at it almost edge-on. You can get a nonsimple polygon that resembles a bow tie, as shown in Figure 2-4, which isn't guaranteed to be rendered correctly. This situation isn't all that unusual if you approximate curved surfaces by quadrilaterals made of points lying on the true surface. You can always avoid the problem by using triangles, since any three points always lie on a plane."
=> Do not use GL_QUAD_STRIP but GL_TRIANGLE_STRIP instead (with primitive restart index if necessary)
//

//
See https://www.glprogramming.com/red/appendixe.html#name2 for computing approximate normals
=> Images in MSL documentation seem to just be faceted (otherwise this needs to be done after setting all vertices to have access to neighboring facets' normals)
//

//
See https://www.learnopengles.com/tag/triangle-strips to making use of index buffer objects
=> They suggest using degenerate triangles to connect consecutive row strips, but primitive restart index is an alternative (not available in GLES until 3.1)
//

//
Enum for two methods to connect consecutive row strips: degenerate triangles, primitive restart index
Enum for two methods to compute surface normal vectors: faceted, averaged (average with unit weights or weighted by angle between each face and current vertex BUT this requires the true normal which we do not have since this is what we are trying to approximate, so just average with equal weights
=> Ah! In fact there are some methods out there to do that:
 - https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.534.7649&rep=rep1&type=pdf
 - https://www.bytehazard.com/articles/vertnorm.html
 - https://knowledge.autodesk.com/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2022/ENU/Maya-Modeling/files/GUID-232E99F8-96B4-4870-8BA0-4887C1C8F0F2-htm.html
 - https://github.com/pmnuckels/wnormals
 - https://en.wikipedia.org/wiki/Vertex_normal
)
Note that since we do indexed vertex rendering we have each vertex only once in existence (and pick it multiple times through its index) so we cannot have multiple normals per vertex (which would be necessary to actually set the facet normal to all the vertices of the corresponding face, and each vertex belonging to several faces) but only one normal per vertex, therefore we are forced to average the normal vectors one way or another!
=> By default, unit weights; other weighing possibilities are: area only, angle only, area and angle
//

//
Allocating and passing multidimensional arrays on the heap:
- https://c-faq.com/aryptr/dynmuldimary.html
- https://c-faq.com/aryptr/ary2dfunc3.html
- https://c-faq.com/aryptr/pass2dary.html
- https://c-faq.com/aryptr/fn33.html
- https://stackoverflow.com/a/21944048 as well as https://c-faq.com/aryptr/dynmuldimary.html for contiguous multidimensional arrays (surely better!)
Ragged array (with allocator overhead) vs. Contiguous array (with pointers overhead)
//
 *
 * @return osg::Geometry* The geometry representing this surface.
 */
osg::Geometry* SurfaceObject::drawGeometry() const
{
  osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();

#define SURFACE_DRAW_SPHERICAL_ARC   1
#define SURFACE_DRAW_RECTANGULAR_BOX 0
#define SURFACE_DRAW_TORUS           0

  const char* id = _id.c_str();
#if SURFACE_DRAW_SPHERICAL_ARC
  const itype nu = 4;
  const itype nv = 32;
#elif SURFACE_DRAW_RECTANGULAR_BOX
  const itype nu = 2;
  const itype nv = 14;
#else
  const itype nu = _nu.exp;
  const itype nv = _nv.exp;
#endif
  const bool wireframe = _wireframe.exp;
  const bool normalized = _normalized.exp;
  const bool doublesided = _doublesided.exp;
  const bool multicolored = _multicolored.exp;
  const ftype opacity = 1.0 - _transparency.exp;
  const ftype colorR = _color[0].exp;
  const ftype colorG = _color[1].exp;
  const ftype colorB = _color[2].exp;

  constexpr itype zero  = 0;
  constexpr itype one   = 1;
  constexpr itype two   = 2;
  constexpr itype three = 3;

  if (nu < one || nv < one) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("A dimension is empty for surface \"%1\" "
                                                              "(nu = %2, nv = %3)."))
                                                              .arg(id)
                                                              .arg(nu)
                                                              .arg(nv),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
    return geometry.release();
  }

  constexpr itype ri = 0; // Restart index

  const ftype ps = mPointSize; // Point size
  const ftype lw = mLineWidth; // Line width
  const ftype ns = mNormalScale; // Normal scale

  const bool point = nu == one && nv == one;            // Surface degenerated to a single point
  const bool line = (nu == one || nv == one) && !point; // Surface degenerated to a line strip

  const bool degenerated = point || line; // Surface degenerated to a point or a line
  const bool closed = mClosenessCheckState == SurfaceClosenessCheckState::active; // Surface checked for closeness
  const bool averaged = mNormalsAverageWeights != SurfaceNormalsAverageWeights::none; // Averaged rendering
  const bool faceted = !normalized && !averaged; // Faceted rendering
  const bool objects = !normalized && !degenerated; // Objects allocation

  const bool area  = mNormalsAverageWeights & SurfaceNormalsAverageWeights::area;  // Normals averaged with area
  const bool angle = mNormalsAverageWeights & SurfaceNormalsAverageWeights::angle; // Normals averaged with angle

  const bool vertexes = mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices; // Normals animated for vertices
  const bool facets   = mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets;   // Normals animated for facets

  const bool degenerate = mStripsWrappingMethod == SurfaceStripsWrappingMethod::degenerate; // Degenerate triangles
  const bool restart    = mStripsWrappingMethod == SurfaceStripsWrappingMethod::restart; // Primitive restart index

  const itype o = restart && !faceted && !degenerated ? ri + one : zero; // Index offset

  const itype imax = std::numeric_limits<itype>::max();
  const ltype lmax = std::numeric_limits<ltype>::max();

  const itype num1 = nu - one;
  const itype nvm1 = nv - one;
  const itype num2 = nu - two;

  itype nVertices = 0;
  itype nIndices  = 0;
  itype nPoints   = 0;
  itype nFacets   = 0;

  try {
    if (degenerated || !faceted) {
      if (imax / nu < nv) {
        throw std::overflow_error("[Points counting] Overflow of the number of surface points");
      }
      nPoints = nu * nv;
    }
    if (!degenerated && (faceted || facets)) {
      if (imax / num1 < nvm1 || imax / two < num1 * nvm1) {
        throw std::overflow_error("[Facets counting] Overflow of the number of triangular facets");
      }
      nFacets = two * num1 * nvm1;
    }
    if (degenerated) {
      nVertices = nPoints;
      nIndices  = nPoints;
    } else if (faceted) {
      if (imax / three < nFacets) {
        throw std::overflow_error("[Faceted rendering] Overflow of the number of triangular facets tripled");
      }
      const itype nFacetsTripled = three * nFacets;
      nVertices = nFacetsTripled;
      nIndices  = nFacetsTripled;
    } else {
      if (imax / num1 < nv || imax / two < num1 * nv) {
        throw std::overflow_error("[Normalized or averaged rendering] Overflow of the number of indices in the triangle strip (not counting breaks)");
      }
      const itype nIStrip = two * num1 * nv;
      nVertices = nPoints;
      nIndices  = nIStrip;
    }
    if (vertexes) {
      if (imax / two < nVertices) {
        throw std::overflow_error("[Animation of vertex normals] Overflow of the number of vertices doubled");
      }
      const itype nVerticesDoubled = two * nVertices;
      if (imax - nVerticesDoubled < nVertices) {
        throw std::overflow_error("[Animation of vertex normals] Overflow of the number of vertices");
      }
      nVertices += nVerticesDoubled;
    }
    if (!degenerated && facets) {
      if (imax / two < nFacets) {
        throw std::overflow_error("[Animation of facet normals] Overflow of the number of triangular facets doubled");
      }
      const itype nFacetsDoubled = two * nFacets;
      if (imax - nFacetsDoubled < nVertices) {
        throw std::overflow_error("[Animation of facet normals] Overflow of the number of vertices");
      }
      nVertices += nFacetsDoubled;
    }
    if (!degenerated && !faceted && restart) {
      const itype iOffset = o;
      const itype nMLines = num2;
      if (imax - iOffset < nVertices) {
        throw std::overflow_error("[Primitive restart index] Overflow of the number of vertices");
      }
      if (imax - nMLines < nIndices) {
        throw std::overflow_error("[Primitive restart index] Overflow of the number of indices");
      }
      nVertices += iOffset;
      nIndices  += nMLines;
    }
    if (!degenerated && !faceted && degenerate) {
      const itype nMLines = num2;
      if (imax / two < nMLines) {
        throw std::overflow_error("[Degenerate triangles] Overflow of the number of middle lines (directed along v-dimension) doubled");
      }
      const itype nMLinesDoubled = two * nMLines;
      if (imax - nMLinesDoubled < nIndices) {
        throw std::overflow_error("[Degenerate triangles] Overflow of the number of indices");
      }
      nIndices += nMLinesDoubled;
    }
  } catch (const std::overflow_error& ex) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("Too many vertices for surface \"%1\" "
                                                              "(nu = %2, nv = %3):\n%4."))
                                                              .arg(id)
                                                              .arg(nu)
                                                              .arg(nv)
                                                              .arg(ex.what()),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
    return geometry.release();
  }

  constexpr itype x = 0; // Index of 1st coordinate
  constexpr itype y = 1; // Index of 2nd coordinate
  constexpr itype z = 2; // Index of 3rd coordinate

  constexpr itype nc = 3; // Number of coordinates

  const itype ne = 1 + normalized + multicolored; // Number of user-provided elements
  const itype na = averaged ? 6 : 2; // Number of considered adjacent facets
  const itype nw = area && angle ? 2 : area || angle ? 1 : 0; // Number of area|angle weights
  const itype no = nc + nw; // Number of objects

  {
    try {
      const ltype netnc = ne * nc;
      if (lmax / nu < netnc) {
        throw std::overflow_error("[Array size for elements] Overflow of ne * nc * nu");
      }
      const ltype netnctnu = netnc * nu;
      if (lmax / nv < netnctnu) {
        throw std::overflow_error("[Array size for elements] Overflow of ne * nc * nu * nv");
      }
    } catch (const std::overflow_error& ex) {
      MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                            QString(QObject::tr("Too many vertices for surface \"%1\" "
                                                                "(ne = %2, nc = %3, nu = %4, nv = %5):\n%6."))
                                                                .arg(id)
                                                                .arg(ne)
                                                                .arg(nc)
                                                                .arg(nu)
                                                                .arg(nv)
                                                                .arg(ex.what()),
                                                            Helper::scriptingKind,
                                                            Helper::errorLevel));
      return geometry.release();
    }
  }

  if (objects) {
    try {
      const ltype notna = no * na;
      if (lmax / nu < notna) {
        throw std::overflow_error("[Array size for objects] Overflow of no * na * nu");
      }
      const ltype notnatnu = notna * nu;
      if (lmax / nv < notnatnu) {
        throw std::overflow_error("[Array size for objects] Overflow of no * na * nu * nv");
      }
    } catch (const std::overflow_error& ex) {
      MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                            QString(QObject::tr("Too many vertices for surface \"%1\" "
                                                                "(no = %2, na = %3, nu = %4, nv = %5):\n%6."))
                                                                .arg(id)
                                                                .arg(no)
                                                                .arg(na)
                                                                .arg(nu)
                                                                .arg(nv)
                                                                .arg(ex.what()),
                                                            Helper::scriptingKind,
                                                            Helper::errorLevel));
      return geometry.release();
    }
  }

  const itype nutnv = nu * nv;

  ftype*** E = nullptr;
  ftype*** O = nullptr;

  ftype*** A = nullptr;
  ftype*** W = nullptr;

  {
    const ltype netnc = ne * nc;
    const ltype netnctnutnv = netnc * nutnv;

    ftype*** Ee  = nullptr;
    ftype** Eec  = nullptr;
    ftype* Eecuv = nullptr;

    try {
      Ee  = new ftype**[ne];
      Eec  = new ftype*[netnc];
      Eecuv = new ftype[netnctnutnv]{0};
      for (itype e = 0; e < ne; e++, Eec += nc) {
        Ee[e] = Eec;
        for (itype c = 0; c < nc; c++, Eecuv += nutnv) {
          Ee[e][c] = Eecuv;
        }
      }
    } catch (const std::bad_alloc& ex) {
      MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                            QString(QObject::tr("Not enough memory to allocate elements for surface \"%1\" "
                                                                "(ne = %2, nc = %3, nu = %4, nv = %5):\n%6."))
                                                                .arg(id)
                                                                .arg(ne)
                                                                .arg(nc)
                                                                .arg(nu)
                                                                .arg(nv)
                                                                .arg(ex.what()),
                                                            Helper::scriptingKind,
                                                            Helper::errorLevel));
      delete[] Eecuv;
      delete[] Eec;
      delete[] Ee;
      return geometry.release();
    }

    E = Ee;
  }

  if (objects) {
    const ltype notna = no * na;
    const ltype notnatnutnv = notna * nutnv;

    ftype*** Oo  = nullptr;
    ftype** Ooa  = nullptr;
    ftype* Ooauv = nullptr;

    try {
      Oo  = new ftype**[no];
      Ooa  = new ftype*[notna];
      Ooauv = new ftype[notnatnutnv]{0};
      for (itype o = 0; o < no; o++, Ooa += na) {
        Oo[o] = Ooa;
        for (itype a = 0; a < na; a++, Ooauv += nutnv) {
          Oo[o][a] = Ooauv;
        }
      }
    } catch (const std::bad_alloc& ex) {
      MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                            QString(QObject::tr("Not enough memory to allocate objects for surface \"%1\" "
                                                                "(no = %2, na = %3, nu = %4, nv = %5):\n%6."))
                                                                .arg(id)
                                                                .arg(no)
                                                                .arg(na)
                                                                .arg(nu)
                                                                .arg(nv)
                                                                .arg(ex.what()),
                                                            Helper::scriptingKind,
                                                            Helper::errorLevel));
      delete[] Ooauv;
      delete[] Ooa;
      delete[] Oo;
      delete[] E[0][0];
      delete[] E[0];
      delete[] E;
      return geometry.release();
    }

    O = Oo;

    A = O;
    W = A + nc;
  }

  ftype** V = nullptr;
  ftype** N = nullptr;
  ftype** C = nullptr;

  itype e = 0;
  {
    V = E[e++];
  }
  if (normalized) {
    N = E[e++];
  }
  if (multicolored) {
    C = E[e++];
  }

  ftype* Vx = nullptr;
  ftype* Vy = nullptr;
  ftype* Vz = nullptr;

  ftype* Nx = nullptr;
  ftype* Ny = nullptr;
  ftype* Nz = nullptr;

  ftype* Cr = nullptr;
  ftype* Cg = nullptr;
  ftype* Cb = nullptr;

  {
    Vx = V[x];
    Vy = V[y];
    Vz = V[z];
  }
  if (normalized) {
    Nx = N[x];
    Ny = N[y];
    Nz = N[z];
  }
  if (multicolored) {
    Cr = C[x];
    Cg = C[y];
    Cb = C[z];
  }

  // TODO: Interface with omc instead of drawing fake surfaces
#if SURFACE_DRAW_SPHERICAL_ARC
  fakeSphericalArc  (nu, nv, Vx, Vy, Vz, Nx, Ny, Nz, Cr, Cg, Cb);
#elif SURFACE_DRAW_RECTANGULAR_BOX
  fakeRectangularBox(nu, nv, Vx, Vy, Vz, Nx, Ny, Nz, Cr, Cg, Cb);
#elif SURFACE_DRAW_TORUS
  fakeTorus         (nu, nv, Vx, Vy, Vz, Nx, Ny, Nz, Cr, Cg, Cb);
#endif

  /* Attributes */
  constexpr osg::StateAttribute::GLModeValue mode = osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE | osg::StateAttribute::PROTECTED;
  osg::ref_ptr<osg::StateSet> ss = geometry->getOrCreateStateSet();
  if (point) {
    ss->setAttributeAndModes(new osg::Point(ps), mode);
  }
  if (line || wireframe || vertexes || facets) {
    ss->setAttributeAndModes(new osg::LineWidth(lw), mode);
  }
  if (restart && !faceted && !degenerated) {
    ss->setAttributeAndModes(new osg::PrimitiveRestartIndex(ri), mode);
    ss->setMode(GL_PRIMITIVE_RESTART, mode);
  }
  if (wireframe) {
    ss->setAttributeAndModes(new osg::PolygonMode(osg::PolygonMode::Face::FRONT_AND_BACK, osg::PolygonMode::Mode::LINE), mode);
  }
  if (doublesided) {
    osg::ref_ptr<osg::LightModel> lightModel = new osg::LightModel();
    lightModel->setTwoSided(true);
    ss->setAttributeAndModes(lightModel.get(), mode);
  }

  /* Arrays */
  osg::ref_ptr<Vec3Array> vertices = new Vec3Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<Vec3Array> normals  = new Vec3Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<Vec4Array> colors   = new Vec4Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<Vec2Array> texels   = new Vec2Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<osg::PrimitiveSet> indices;
  osg::ref_ptr<Vec3Array> facetsCenters;
  osg::ref_ptr<Vec3Array> facetsNormals;
  if (facets) {
    facetsCenters = new Vec3Array();
    facetsNormals = new Vec3Array();
  }

  vertices->reserve(nVertices);
  normals ->reserve(nVertices);
  colors  ->reserve(nVertices);
  texels  ->reserve(nVertices);
  if (facets) {
    facetsCenters->reserve(nFacets);
    facetsNormals->reserve(nFacets);
  }

  if (degenerated) {
    const ftype fnutnvm1 = (ftype)(nPoints - one);

    for (itype u = 0; u < nu; u++) {
      for (itype v = 0; v < nv; v++) {
        {
          vertices->push_back(Vec3(Vx[nv * u + v], Vy[nv * u + v], Vz[nv * u + v]));
        }

        if (normalized) {
          normals ->push_back(Vec3(Nx[nv * u + v], Ny[nv * u + v], Nz[nv * u + v]));
        } else {
          normals ->push_back(Vec3(Vx[nv * u + v], Vy[nv * u + v], Vz[nv * u + v]));
        }

        if (multicolored) {
          colors  ->push_back(Vec4(Cr[nv * u + v], Cg[nv * u + v], Cb[nv * u + v], opacity));
        } else {
          colors  ->push_back(Vec4(colorR, colorG, colorB, opacity));
        }

        if (point) {
          texels  ->push_back(Vec2(u, v));
        } else {
          texels  ->push_back(Vec2(u / fnutnvm1, v / fnutnvm1));
        }
      }
    }

    indices = new osg::DrawArrays(point ? osg::PrimitiveSet::POINTS : osg::PrimitiveSet::LINE_STRIP, zero, nIndices);
  } else {
    const ftype fnum1 = (ftype)num1;
    const ftype fnvm1 = (ftype)nvm1;

    if (restart) { // Duplicate a dummy vertex (see OSG commit 353b18b)
      for (itype i = 0; i < o; i++) {
        vertices->push_back(Vec3());
        normals ->push_back(Vec3());
        colors  ->push_back(Vec4());
        texels  ->push_back(Vec2());
      }
    }

    /* Vertices */
    if (!faceted) {
      for (itype u = 0; u < nu; u++) {
        for (itype v = 0; v < nv; v++) {
          vertices->push_back(Vec3(Vx[nv * u + v], Vy[nv * u + v], Vz[nv * u + v]));
        }
      }
    } else {
      for (itype u = 0; u < num1; u++) {
        const itype up0 = u;
        const itype up1 = u + one;
        for (itype v = 0; v < nvm1; v++) {
          const itype vp0 = v;
          const itype vp1 = v + one;
          const ftype Vx00 = Vx[nv * up0 + vp0];
          const ftype Vy00 = Vy[nv * up0 + vp0];
          const ftype Vz00 = Vz[nv * up0 + vp0];
          const ftype Vx01 = Vx[nv * up0 + vp1];
          const ftype Vy01 = Vy[nv * up0 + vp1];
          const ftype Vz01 = Vz[nv * up0 + vp1];
          const ftype Vx10 = Vx[nv * up1 + vp0];
          const ftype Vy10 = Vy[nv * up1 + vp0];
          const ftype Vz10 = Vz[nv * up1 + vp0];
          const ftype Vx11 = Vx[nv * up1 + vp1];
          const ftype Vy11 = Vy[nv * up1 + vp1];
          const ftype Vz11 = Vz[nv * up1 + vp1];
          vertices->push_back(Vec3(Vx00, Vy00, Vz00));
          vertices->push_back(Vec3(Vx10, Vy10, Vz10));
          vertices->push_back(Vec3(Vx01, Vy01, Vz01));
          vertices->push_back(Vec3(Vx10, Vy10, Vz10));
          vertices->push_back(Vec3(Vx11, Vy11, Vz11));
          vertices->push_back(Vec3(Vx01, Vy01, Vz01));
        }
      }
    }

    /* Normals */
    if (normalized) {
      for (itype u = 0; u < nu; u++) {
        for (itype v = 0; v < nv; v++) {
          normals->push_back(Vec3(Nx[nv * u + v], Ny[nv * u + v], Nz[nv * u + v]));
        }
      }
    }
    if (!normalized || facets) {
      constexpr itype i = 0;
      constexpr itype j = 2;
      constexpr itype k = 4;
      constexpr itype l = 1;
      constexpr itype m = 3;
      constexpr itype n = 5;
      for (itype u = 0; u < num1; u++) {
        const itype up0 = u;
        const itype up1 = u + one;
        for (itype v = 0; v < nvm1; v++) {
          const itype vp0 = v;
          const itype vp1 = v + one;
          const ftype Vx00 = Vx[nv * up0 + vp0];
          const ftype Vy00 = Vy[nv * up0 + vp0];
          const ftype Vz00 = Vz[nv * up0 + vp0];
          const ftype Vx01 = Vx[nv * up0 + vp1];
          const ftype Vy01 = Vy[nv * up0 + vp1];
          const ftype Vz01 = Vz[nv * up0 + vp1];
          const ftype Vx10 = Vx[nv * up1 + vp0];
          const ftype Vy10 = Vy[nv * up1 + vp0];
          const ftype Vz10 = Vz[nv * up1 + vp0];
          const ftype Vx11 = Vx[nv * up1 + vp1];
          const ftype Vy11 = Vy[nv * up1 + vp1];
          const ftype Vz11 = Vz[nv * up1 + vp1];
          const ftype du1[nc] = {Vx10 - Vx00, Vy10 - Vy00, Vz10 - Vz00};
          const ftype du2[nc] = {Vx01 - Vx11, Vy01 - Vy11, Vz01 - Vz11};
          const ftype dv1[nc] = {Vx01 - Vx00, Vy01 - Vy00, Vz01 - Vz00};
          const ftype dv2[nc] = {Vx10 - Vx11, Vy10 - Vy11, Vz10 - Vz11};
          const ftype cross1[nc] = {du1[y] * dv1[z] - du1[z] * dv1[y], du1[z] * dv1[x] - du1[x] * dv1[z], du1[x] * dv1[y] - du1[y] * dv1[x]};
          const ftype cross2[nc] = {du2[y] * dv2[z] - du2[z] * dv2[y], du2[z] * dv2[x] - du2[x] * dv2[z], du2[x] * dv2[y] - du2[y] * dv2[x]};
          const ftype length1 = std::sqrt(cross1[x] * cross1[x] + cross1[y] * cross1[y] + cross1[z] * cross1[z]);
          const ftype length2 = std::sqrt(cross2[x] * cross2[x] + cross2[y] * cross2[y] + cross2[z] * cross2[z]);
          const ftype invnorm1 = length1 > 0 ? 1 / length1 : 0;
          const ftype invnorm2 = length2 > 0 ? 1 / length2 : 0;
          const ftype normal1[nc] = {cross1[x] * invnorm1, cross1[y] * invnorm1, cross1[z] * invnorm1};
          const ftype normal2[nc] = {cross2[x] * invnorm2, cross2[y] * invnorm2, cross2[z] * invnorm2};
          if (facets) {
            const ftype center1[nc] = {(Vx00 + Vx10 + Vx01) / 3, (Vy00 + Vy10 + Vy01) / 3, (Vz00 + Vz10 + Vz01) / 3};
            const ftype center2[nc] = {(Vx10 + Vx11 + Vx01) / 3, (Vy10 + Vy11 + Vy01) / 3, (Vz10 + Vz11 + Vz01) / 3};
            facetsCenters->push_back(Vec3(center1[x], center1[y], center1[z]));
            facetsCenters->push_back(Vec3(center2[x], center2[y], center2[z]));
            facetsNormals->push_back(Vec3(normal1[x], normal1[y], normal1[z]));
            facetsNormals->push_back(Vec3(normal2[x], normal2[y], normal2[z]));
          }
          if (!normalized) {
            if (!averaged) {
              for (itype c = 0; c < nc; c++) {
                A[c][i][nv * up0 + vp0] = normal1[c];
                A[c][l][nv * up1 + vp1] = normal2[c];
              }
            } else {
              for (itype c = 0; c < nc; c++) {
                A[c][i][nv * up0 + vp0] = normal1[c];
                A[c][j][nv * up1 + vp0] = normal1[c];
                A[c][k][nv * up0 + vp1] = normal1[c];
                A[c][l][nv * up1 + vp0] = normal2[c];
                A[c][m][nv * up1 + vp1] = normal2[c];
                A[c][n][nv * up0 + vp1] = normal2[c];
              }
              itype w = 0;
              if (area) {
                const ftype area1 = length1 / 2;
                const ftype area2 = length2 / 2;
                W[w][i][nv * up0 + vp0] = area1;
                W[w][j][nv * up1 + vp0] = area1;
                W[w][k][nv * up0 + vp1] = area1;
                W[w][l][nv * up1 + vp0] = area2;
                W[w][m][nv * up1 + vp1] = area2;
                W[w][n][nv * up0 + vp1] = area2;
                w++;
              }
              if (angle) {
                const ftype duv[nc] = {Vx01 - Vx10, Vy01 - Vy10, Vz01 - Vz10};
                const ftype dvu[nc] = {Vx10 - Vx01, Vy10 - Vy01, Vz10 - Vz01};
                const ftype dot11 = dv1[x] * du1[x] + dv1[y] * du1[y] + dv1[z] * du1[z];
                const ftype dot12 = du1[x] * dvu[x] + du1[y] * dvu[y] + du1[z] * dvu[z];
                const ftype dot13 = duv[x] * dv1[x] + duv[y] * dv1[y] + duv[z] * dv1[z];
                const ftype dot21 = dvu[x] * dv2[x] + dvu[y] * dv2[y] + dvu[z] * dv2[z];
                const ftype dot22 = dv2[x] * du2[x] + dv2[y] * du2[y] + dv2[z] * du2[z];
                const ftype dot23 = du2[x] * duv[x] + du2[y] * duv[y] + du2[z] * duv[z];
                const ftype lu1 = std::sqrt(du1[x] * du1[x] + du1[y] * du1[y] + du1[z] * du1[z]);
                const ftype lu2 = std::sqrt(du2[x] * du2[x] + du2[y] * du2[y] + du2[z] * du2[z]);
                const ftype luv = std::sqrt(duv[x] * duv[x] + duv[y] * duv[y] + duv[z] * duv[z]);
                const ftype lv1 = std::sqrt(dv1[x] * dv1[x] + dv1[y] * dv1[y] + dv1[z] * dv1[z]);
                const ftype lv2 = std::sqrt(dv2[x] * dv2[x] + dv2[y] * dv2[y] + dv2[z] * dv2[z]);
                const ftype lvu = std::sqrt(dvu[x] * dvu[x] + dvu[y] * dvu[y] + dvu[z] * dvu[z]);
                const ftype length11 = lv1 * lu1;
                const ftype length12 = lu1 * lvu;
                const ftype length13 = luv * lv1;
                const ftype length21 = lvu * lv2;
                const ftype length22 = lv2 * lu2;
                const ftype length23 = lu2 * luv;
                const ftype angle11 = std::acos(length11 > 0 ? dot11 / length11 : 0);
                const ftype angle12 = std::acos(length12 > 0 ? dot12 / length12 : 0);
                const ftype angle13 = std::acos(length13 > 0 ? dot13 / length13 : 0);
                const ftype angle21 = std::acos(length21 > 0 ? dot21 / length21 : 0);
                const ftype angle22 = std::acos(length22 > 0 ? dot22 / length22 : 0);
                const ftype angle23 = std::acos(length23 > 0 ? dot23 / length23 : 0);
                W[w][i][nv * up0 + vp0] = angle11;
                W[w][j][nv * up1 + vp0] = angle12;
                W[w][k][nv * up0 + vp1] = angle13;
                W[w][l][nv * up1 + vp0] = angle21;
                W[w][m][nv * up1 + vp1] = angle22;
                W[w][n][nv * up0 + vp1] = angle23;
                w++;
              }
            }
          }
        }
      }
      if (!normalized) {
        if (averaged) {
          if (closed) {
            constexpr itype nut0 = 0;
            constexpr itype nvt0 = 0;
            for (itype u = 0; u < nu; u++) {
              if (Vx[nv * u + nvt0] == Vx[nv * u + nvm1] &&
                  Vy[nv * u + nvt0] == Vy[nv * u + nvm1] &&
                  Vz[nv * u + nvt0] == Vz[nv * u + nvm1]) {
                for (itype c = 0; c < nc; c++) {
                  A[c][i][nv * u + nvm1] = A[c][i][nv * u + nvt0];
                  A[c][l][nv * u + nvm1] = A[c][l][nv * u + nvt0];
                  A[c][j][nv * u + nvm1] = A[c][j][nv * u + nvt0];
                  A[c][m][nv * u + nvt0] = A[c][m][nv * u + nvm1];
                  A[c][k][nv * u + nvt0] = A[c][k][nv * u + nvm1];
                  A[c][n][nv * u + nvt0] = A[c][n][nv * u + nvm1];
                }
                for (itype w = 0; w < nw; w++) {
                  W[w][i][nv * u + nvm1] = W[w][i][nv * u + nvt0];
                  W[w][l][nv * u + nvm1] = W[w][l][nv * u + nvt0];
                  W[w][j][nv * u + nvm1] = W[w][j][nv * u + nvt0];
                  W[w][m][nv * u + nvt0] = W[w][m][nv * u + nvm1];
                  W[w][k][nv * u + nvt0] = W[w][k][nv * u + nvm1];
                  W[w][n][nv * u + nvt0] = W[w][n][nv * u + nvm1];
                }
              }
            }
            for (itype v = 0; v < nv; v++) {
              if (Vx[nv * nut0 + v] == Vx[nv * num1 + v] &&
                  Vy[nv * nut0 + v] == Vy[nv * num1 + v] &&
                  Vz[nv * nut0 + v] == Vz[nv * num1 + v]) {
                for (itype c = 0; c < nc; c++) {
                  A[c][i][nv * num1 + v] = A[c][i][nv * nut0 + v];
                  A[c][k][nv * num1 + v] = A[c][k][nv * nut0 + v];
                  A[c][n][nv * num1 + v] = A[c][n][nv * nut0 + v];
                  A[c][l][nv * nut0 + v] = A[c][l][nv * num1 + v];
                  A[c][j][nv * nut0 + v] = A[c][j][nv * num1 + v];
                  A[c][m][nv * nut0 + v] = A[c][m][nv * num1 + v];
                }
                for (itype w = 0; w < nw; w++) {
                  W[w][i][nv * num1 + v] = W[w][i][nv * nut0 + v];
                  W[w][k][nv * num1 + v] = W[w][k][nv * nut0 + v];
                  W[w][n][nv * num1 + v] = W[w][n][nv * nut0 + v];
                  W[w][l][nv * nut0 + v] = W[w][l][nv * num1 + v];
                  W[w][j][nv * nut0 + v] = W[w][j][nv * num1 + v];
                  W[w][m][nv * nut0 + v] = W[w][m][nv * num1 + v];
                }
              }
            }
          }
          for (itype u = 0; u < nu; u++) {
            for (itype v = 0; v < nv; v++) {
              ftype normal[nc] = {0};
              for (itype a = 0; a < na; a++) {
                ftype weight = 1;
                for (itype w = 0; w < nw; w++) {
                  weight *= W[w][a][nv * u + v];
                }
                for (itype c = 0; c < nc; c++) {
                  normal[c] += A[c][a][nv * u + v] * weight;
                }
              }
              const ftype length = std::sqrt(normal[x] * normal[x] + normal[y] * normal[y] + normal[z] * normal[z]);
              const ftype invnorm = length > 0 ? 1 / length : 0;
              for (itype c = 0; c < nc; c++) {
                normal[c] *= invnorm;
              }
              normals->push_back(Vec3(normal[x], normal[y], normal[z]));
            }
          }
        } else {
          for (itype u = 0; u < num1; u++) {
            const itype up0 = u;
            const itype up1 = u + one;
            for (itype v = 0; v < nvm1; v++) {
              const itype vp0 = v;
              const itype vp1 = v + one;
              ftype normal1[nc] = {0};
              ftype normal2[nc] = {0};
              for (itype c = 0; c < nc; c++) {
                normal1[c] += A[c][i][nv * up0 + vp0];
                normal2[c] += A[c][l][nv * up1 + vp1];
              }
              normals->insert(normals->end(), three, Vec3(normal1[x], normal1[y], normal1[z]));
              normals->insert(normals->end(), three, Vec3(normal2[x], normal2[y], normal2[z]));
            }
          }
        }
      }
    }

    /* Colors */
    if (multicolored) {
      if (!faceted) {
        for (itype u = 0; u < nu; u++) {
          for (itype v = 0; v < nv; v++) {
            colors->push_back(Vec4(Cr[nv * u + v], Cg[nv * u + v], Cb[nv * u + v], opacity));
          }
        }
      } else {
        for (itype u = 0; u < num1; u++) {
          const itype up0 = u;
          const itype up1 = u + one;
          for (itype v = 0; v < nvm1; v++) {
            const itype vp0 = v;
            const itype vp1 = v + one;
            const ftype Cr00 = Cr[nv * up0 + vp0];
            const ftype Cg00 = Cg[nv * up0 + vp0];
            const ftype Cb00 = Cb[nv * up0 + vp0];
            const ftype Cr01 = Cr[nv * up0 + vp1];
            const ftype Cg01 = Cg[nv * up0 + vp1];
            const ftype Cb01 = Cb[nv * up0 + vp1];
            const ftype Cr10 = Cr[nv * up1 + vp0];
            const ftype Cg10 = Cg[nv * up1 + vp0];
            const ftype Cb10 = Cb[nv * up1 + vp0];
            const ftype Cr11 = Cr[nv * up1 + vp1];
            const ftype Cg11 = Cg[nv * up1 + vp1];
            const ftype Cb11 = Cb[nv * up1 + vp1];
            colors->push_back(Vec4(Cr00, Cg00, Cb00, opacity));
            colors->push_back(Vec4(Cr10, Cg10, Cb10, opacity));
            colors->push_back(Vec4(Cr01, Cg01, Cb01, opacity));
            colors->push_back(Vec4(Cr10, Cg10, Cb10, opacity));
            colors->push_back(Vec4(Cr11, Cg11, Cb11, opacity));
            colors->push_back(Vec4(Cr01, Cg01, Cb01, opacity));
          }
        }
      }
    } else {
      colors->clear();
      colors->setBinding(osg::Array::BIND_OVERALL);
      colors->push_back(Vec4(colorR, colorG, colorB, opacity));
    }

    /* Texels */
    if (!faceted) {
      for (itype u = 0; u < nu; u++) {
        for (itype v = 0; v < nv; v++) {
          texels->push_back(Vec2(u / fnum1, v / fnvm1));
        }
      }
    } else {
      for (itype u = 0; u < num1; u++) {
        const ftype Tu0 =  u        / fnum1;
        const ftype Tu1 = (u + one) / fnum1;
        for (itype v = 0; v < nvm1; v++) {
          const ftype Tv0 =  v        / fnvm1;
          const ftype Tv1 = (v + one) / fnvm1;
          texels->push_back(Vec2(Tu0, Tv0));
          texels->push_back(Vec2(Tu1, Tv0));
          texels->push_back(Vec2(Tu0, Tv1));
          texels->push_back(Vec2(Tu1, Tv0));
          texels->push_back(Vec2(Tu1, Tv1));
          texels->push_back(Vec2(Tu0, Tv1));
        }
      }
    }

    /* Indices */
    if (!faceted) {
      indices = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLE_STRIP);
      osg::ref_ptr<osg::DrawElements> strip = indices->getDrawElements();
      strip->reserveElements(nIndices);
      const bool num2g0 = num2 > 0;
      const itype num2tnvpo = num2 * nv + o;
      itype up0tnvpopv = o;
      itype up1tnvpo = nv + o;
#define SURFACE_INDICES_V_FIRST()               \
  if (degenerate) {                             \
    strip->addElement(up0tnvpopv);              \
  }
#define SURFACE_INDICES_V_LOOP()                \
  for (; up0tnvpopv < up1tnvpo; up0tnvpopv++) { \
    strip->addElement(up0tnvpopv);              \
    strip->addElement(up0tnvpopv + nv);         \
  }
#define SURFACE_INDICES_V_LAST()                \
  if (degenerate) {                             \
    strip->addElement(up0tnvpopv + nvm1);       \
  } else if (restart) {                         \
    strip->addElement(ri);                      \
  }
      {
        SURFACE_INDICES_V_LOOP();
      }
      if (num2g0) {
        SURFACE_INDICES_V_LAST();
      }
      for (up1tnvpo += nv; up0tnvpopv < num2tnvpo; up1tnvpo += nv) {
        SURFACE_INDICES_V_FIRST();
        SURFACE_INDICES_V_LOOP();
        SURFACE_INDICES_V_LAST();
      }
      if (num2g0) {
        SURFACE_INDICES_V_FIRST();
        SURFACE_INDICES_V_LOOP();
      }
    } else {
      indices = new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, o, nIndices);
    }
  }

  /* Debug */
#define SURFACE_DEBUG_N(o, e, v, n, r, g, b)    \
  const itype l = vertices->size();             \
  const itype c = two * (e - o);                \
  for (itype i = o; i < e; i++) {               \
    const Vec3 vertex = v->at(i);               \
    const Vec3 normal = n->at(i);               \
    const Vec2 texel0 = Vec2(0, 1);             \
    const Vec2 texel1 = Vec2(1, 1);             \
    vertices->push_back(vertex);                \
    vertices->push_back(vertex + normal * ns);  \
    normals ->push_back(normal);                \
    normals ->push_back(normal);                \
    texels  ->push_back(texel0);                \
    texels  ->push_back(texel1);                \
  }                                             \
  if (multicolored) {                           \
    colors->insert(colors->end(), c,            \
        Vec4(r, g, b, opacity));                \
  }                                             \
  osg::ref_ptr<osg::PrimitiveSet> lines =       \
      new osg::DrawArrays(                      \
          osg::PrimitiveSet::LINES, l, c);      \
  geometry->addPrimitiveSet(lines.get());
  if (vertexes) {
    const itype offset = o;
    const itype size = vertices->size();
    SURFACE_DEBUG_N(
        offset, size,
        vertices,
        normals,
        1, 0, 0);
  }
  if (facets) {
    const itype offset = 0;
    const itype size = facetsCenters->size();
    SURFACE_DEBUG_N(
        offset, size,
        facetsCenters,
        facetsNormals,
        0, 0, 1);
  }

  /* Geometry */
  geometry->setUseVertexBufferObjects(true);
  geometry->setVertexArray(vertices.get());
  geometry->setNormalArray(normals.get());
  geometry->setColorArray(colors.get());
  geometry->setTexCoordArray(zero, texels.get());
  geometry->addPrimitiveSet(indices.get());

  if (objects) {
    delete[] O[0][0];
    delete[] O[0];
    delete[] O;
  }
  {
    delete[] E[0][0];
    delete[] E[0];
    delete[] E;
  }

  return geometry.release();
}
