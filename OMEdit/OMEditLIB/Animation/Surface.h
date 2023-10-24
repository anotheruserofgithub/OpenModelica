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

#ifndef SURFACE_H
#define SURFACE_H

#include "AbstractVisualizer.h"

#include <type_traits>

#include <QOpenGLContext> // must be included before OSG headers

#include <osg/Array>
#include <osg/Geometry>

enum class SurfacePolygonsRasterMode {plain = 0b00, wireframe = 0b01, pointcloud = 0b10};
enum class SurfaceClosenessCheckState {inactive = 0b0, active = 0b1};
enum class SurfaceStripsWrappingMethod {restart = 0b0, degenerate = 0b1};
enum class SurfaceNormalsAverageWeights {none = 0b000, equal = 0b001, area = 0b010, angle = 0b100, bothAreaAndAngle = 0b110};
enum class SurfaceNormalsAnimationTypes {none = 0b000, vertices = 0b001, facets = 0b010, bothVerticesAndFacets = 0b011};

std::underlying_type<SurfaceNormalsAverageWeights>::type operator&(const SurfaceNormalsAverageWeights& lhs, const SurfaceNormalsAverageWeights& rhs);
std::underlying_type<SurfaceNormalsAnimationTypes>::type operator&(const SurfaceNormalsAnimationTypes& lhs, const SurfaceNormalsAnimationTypes& rhs);

class SurfaceObject final : public AbstractVisualizerObjectWithVisualProperties<SurfaceObject>
{
public:
  SurfaceObject();
  ~SurfaceObject() = default;
  SurfaceObject(const SurfaceObject&) = default;
  SurfaceObject& operator=(const SurfaceObject&) = default;
  SurfaceObject* asSurface() override final {return this;}
  void dumpVisualizerAttributes() override;
  osg::Geometry* drawGeometry() const;
private:
  typedef int itype;
  typedef double ftype;
  typedef std::size_t ltype;
  typedef osg::Vec2Array Vec2Array;
  typedef osg::Vec3Array Vec3Array;
  typedef osg::Vec4Array Vec4Array;
  typedef Vec2Array::ElementDataType Vec2;
  typedef Vec3Array::ElementDataType Vec3;
  typedef Vec4Array::ElementDataType Vec4;
private: // TODO: Remove
  void fakeTorus         (const itype nu, const itype nv,
                          ftype* Vx, ftype* Vy, ftype* Vz,
                          ftype* Nx, ftype* Ny, ftype* Nz,
                          ftype* Cr, ftype* Cg, ftype* Cb) const;
  void fakeRectangularBox(const itype nu, const itype nv,
                          ftype* Vx, ftype* Vy, ftype* Vz,
                          ftype* Nx, ftype* Ny, ftype* Nz,
                          ftype* Cr, ftype* Cg, ftype* Cb) const;
  void fakeSphericalArc  (const itype nu, const itype nv,
                          ftype* Vx, ftype* Vy, ftype* Vz,
                          ftype* Nx, ftype* Ny, ftype* Nz,
                          ftype* Cr, ftype* Cg, ftype* Cb) const;
private:
  SurfacePolygonsRasterMode mPolygonsRasterMode;
  SurfaceClosenessCheckState mClosenessCheckState;
  SurfaceStripsWrappingMethod mStripsWrappingMethod;
  SurfaceNormalsAverageWeights mNormalsAverageWeights;
  SurfaceNormalsAnimationTypes mNormalsAnimationTypes;
  ftype mPointSize;
  ftype mLineWidth;
  ftype mNormalScale;
  ftype mNormalColorVertex[3];
  ftype mNormalColorFacet[3];
public:
  VisualizerAttribute _nu;
  VisualizerAttribute _nv;
  VisualizerAttribute _wireframe;
  VisualizerAttribute _normalized;
  VisualizerAttribute _doublesided;
  VisualizerAttribute _multicolored;
  VisualizerAttribute _transparency;
};

#endif
