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

std::underlying_type<SurfaceNormalsAverageWeights>::type operator&(const SurfaceNormalsAverageWeights& lhs, const SurfaceNormalsAverageWeights& rhs)
{
  return static_cast<std::underlying_type<SurfaceNormalsAverageWeights>::type>(lhs) & static_cast<std::underlying_type<SurfaceNormalsAverageWeights>::type>(rhs);
}

std::underlying_type<SurfaceNormalsAnimationTypes>::type operator&(const SurfaceNormalsAnimationTypes& lhs, const SurfaceNormalsAnimationTypes& rhs)
{
  return static_cast<std::underlying_type<SurfaceNormalsAnimationTypes>::type>(lhs) & static_cast<std::underlying_type<SurfaceNormalsAnimationTypes>::type>(rhs);
}

SurfaceObject::SurfaceObject()
    : AbstractVisualizerObject(VisualizerType::surface),
      mClosenessCheckState(SurfaceClosenessCheckState::active),
      mStripsWrappingMethod(SurfaceStripsWrappingMethod::degenerate),
      mNormalsAverageWeights(SurfaceNormalsAverageWeights::bothAreaAndAngle),
      mNormalsAnimationTypes(SurfaceNormalsAnimationTypes::bothVerticesAndFacets), // FIXME should be none by default
      _nu(VisualizerAttribute(0.0)),
      _nv(VisualizerAttribute(0.0)),
      _wireframe(VisualizerAttribute(0.0)),
      _normalized(VisualizerAttribute(0.0)),
      _doublesided(VisualizerAttribute(1.0)),
      _multicolored(VisualizerAttribute(0.0)),
      _transparency(VisualizerAttribute(0.0))
{
}

void SurfaceObject::dumpVisualizerAttributes() const
{
  AbstractVisualizerObject::dumpVisualizerAttributes();
  std::cout << "nu " << _nu.getValueString() << std::endl;
  std::cout << "nv " << _nv.getValueString() << std::endl;
  std::cout << "wireframe " << _wireframe.getValueString() << std::endl;
  std::cout << "normalized " << _normalized.getValueString() << std::endl;
  std::cout << "doublesided " << _doublesided.getValueString() << std::endl;
  std::cout << "multicolored " << _multicolored.getValueString() << std::endl;
  std::cout << "transparency " << _transparency.getValueString() << std::endl;
}

#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>

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

#include "Modeling/MessagesWidget.h"
#include "Util/Helper.h"

void SurfaceObject::fakeTorus(const itype nu, const itype nv, ftype** X, ftype** Y, ftype** Z, ftype*** N, ftype*** C) const
{
  constexpr ftype pi = M_PI;

  constexpr ftype R = 1; // Major radius (distance from center of torus to center of tube)
  constexpr ftype r = 0.2; // Minor radius (radius of tube)
  constexpr ftype opening = 0; // Opening angle of torus
  constexpr ftype startAngle = -pi; // Start angle of torus slice
  constexpr ftype  stopAngle = +pi; // Stop  angle of torus slice

  constexpr ftype phi_start = -pi + opening;
  constexpr ftype phi_stop  = +pi - opening;

  ftype alpha;
  ftype beta;

  for (itype u = 0; u < nu; u++) {
    alpha = startAngle + (nu > 1 ? (stopAngle - startAngle) * u / (nu - 1) : 0);
    for (itype v = 0; v < nv; v++) {
      beta = phi_start + (nv > 1 ? (phi_stop - phi_start) * v / (nv - 1) : 0);
      Z[u][v] = (R + r * std::cos(beta)) * std::cos(alpha);
      X[u][v] = (R + r * std::cos(beta)) * std::sin(alpha);
      Y[u][v] =      r * std::sin(beta);
      if (N) {
        N[u][v][2] = std::cos(beta) * std::cos(alpha);
        N[u][v][0] = std::cos(beta) * std::sin(alpha);
        N[u][v][1] = std::sin(beta);
      }
      if (C) {
        C[u][v][0] = 0;
        C[u][v][1] = 1;
        C[u][v][2] = 0;
      }
    }
  }

  if (mClosenessCheckState == SurfaceClosenessCheckState::active) {
    if (phi_stop - pi == phi_start + pi) {
      for (itype u = 0; u < nu; u++) {
        X[u][nv - 1] = X[u][0];
        Y[u][nv - 1] = Y[u][0];
        Z[u][nv - 1] = Z[u][0];
        if (N) {
          N[u][nv - 1][0] = N[u][0][0];
          N[u][nv - 1][1] = N[u][0][1];
          N[u][nv - 1][2] = N[u][0][2];
        }
      }
    }
    if (stopAngle - pi == startAngle + pi) {
      for (itype v = 0; v < nv; v++) {
        X[nu - 1][v] = X[0][v];
        Y[nu - 1][v] = Y[0][v];
        Z[nu - 1][v] = Z[0][v];
        if (N) {
          N[nu - 1][v][0] = N[0][v][0];
          N[nu - 1][v][1] = N[0][v][1];
          N[nu - 1][v][2] = N[0][v][2];
        }
      }
    }
  }
}

osg::Geometry* SurfaceObject::drawGeometry() const
{
  osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();

  const char* id = _id.c_str();
  const itype nu = _nu.exp;
  const itype nv = _nv.exp;
  const bool wireframe = _wireframe.exp;
  const bool normalized = _normalized.exp;
  const bool doublesided = _doublesided.exp;
  const bool multicolored = _multicolored.exp;
  const ftype opacity = 1.0 - _transparency.exp;

  if (nu < 1 || nv < 1) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("A dimension of surface \"%1\" is empty (nu = %2, nv = %3)."))
                                                              .arg(id)
                                                              .arg(nu)
                                                              .arg(nv),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
    return geometry.release();
  }

  const bool degenerated = nu == 1 || nv == 1;
  const bool faceted = !normalized && mNormalsAverageWeights == SurfaceNormalsAverageWeights::none;

  try {
    itype nVerticesBefore = 0;
    itype nIn_dicesBefore = 0;
    itype nLinesElVBefore = 0;
    itype nLinesElFBefore = 0;
    const itype imax = std::numeric_limits<itype>::max();
    const itype nPoints = nu * nv; // Number of surface points
    const itype nMLines = nu - 2; // Number of middle lines directed along v-dimension
    const itype nIStrip = 2 * (nu - 1) * nv; // Number of indices in triangle strip (not counting degenerate triangles)
    const itype nFacets = 2 * (nu - 1) * (nv - 1); // Number of triangular facets
    const itype iOffset = 1; // Index offset
    if (degenerated) {
      if (imax / nu < nv) {
        throw std::overflow_error("Overflow of nPoints (surface degenerated to line or point)");
      }
      const itype capacity = imax - nPoints;
      if (capacity < nVerticesBefore) {
        throw std::overflow_error("Overflow of vertices (surface degenerated to line or point)");
      }
      if (capacity < nIn_dicesBefore) {
        throw std::overflow_error("Overflow of indices (surface degenerated to line or point)");
      }
      nVerticesBefore += nPoints;
      nIn_dicesBefore += nPoints;
    } else if (faceted) {
      if (imax / (nu - 1) < (nv - 1) || imax / 2 < (nu - 1) * (nv - 1)) {
        throw std::overflow_error("Overflow of nFacets (faceted rendering)");
      }
      if (imax / 3 < nFacets) {
        throw std::overflow_error("Overflow of nFacetsTripled (faceted rendering)");
      }
      const itype nFacetsTripled = 3 * nFacets;
      const itype capacity = imax - nFacetsTripled;
      if (capacity < nVerticesBefore) {
        throw std::overflow_error("Overflow of vertices (faceted rendering)");
      }
      if (capacity < nIn_dicesBefore) {
        throw std::overflow_error("Overflow of indices (faceted rendering)");
      }
      nVerticesBefore += nFacetsTripled;
      nIn_dicesBefore += nFacetsTripled;
    } else {
      if (imax / nu < nv) {
        throw std::overflow_error("Overflow of nPoints (averaged rendering)");
      }
      if (imax / (nu - 1) < nv || imax / 2 < (nu - 1) * nv) {
        throw std::overflow_error("Overflow of nIStrip (averaged rendering)");
      }
      const itype capacityV = imax - nPoints;
      const itype capacityI = imax - nIStrip;
      if (capacityV < nVerticesBefore) {
        throw std::overflow_error("Overflow of vertices (averaged rendering)");
      }
      if (capacityI < nIn_dicesBefore) {
        throw std::overflow_error("Overflow of indices (averaged rendering)");
      }
      nVerticesBefore += nPoints;
      nIn_dicesBefore += nIStrip;
    }
    if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices) {
      if (imax / 2 < nVerticesBefore) {
        throw std::overflow_error("Overflow of nVerticesDoubled (animation of vertex normals)");
      }
      const itype nVerticesDoubled = 2 * nVerticesBefore;
      const itype capacity = imax - nVerticesDoubled;
      if (capacity < nVerticesBefore) {
        throw std::overflow_error("Overflow of vertices (animation of vertex normals)");
      }
      if (capacity < nLinesElVBefore) {
        throw std::overflow_error("Overflow of indices of vertex normals (animation of vertex normals)");
      }
      nVerticesBefore += nVerticesDoubled;
      nLinesElVBefore += nVerticesDoubled;
    }
    if (!degenerated && mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
      if (imax / (nu - 1) < (nv - 1) || imax / 2 < (nu - 1) * (nv - 1)) {
        throw std::overflow_error("Overflow of nFacets (animation of facet normals)");
      }
      if (imax / 2 < nFacets) {
        throw std::overflow_error("Overflow of nFacetsDoubled (animation of facet normals)");
      }
      const itype nFacetsDoubled = 2 * nFacets;
      const itype capacity = imax - nFacetsDoubled;
      if (capacity < nVerticesBefore) {
        throw std::overflow_error("Overflow of vertices (animation of facet normals)");
      }
      if (capacity < nLinesElFBefore) {
        throw std::overflow_error("Overflow of indices of facet normals (animation of facet normals)");
      }
      nVerticesBefore += nFacetsDoubled;
      nLinesElFBefore += nFacetsDoubled;
    }
    if (!degenerated && !faceted && mStripsWrappingMethod == SurfaceStripsWrappingMethod::restart) {
      const itype capacityV = imax - iOffset;
      const itype capacityI = imax - nMLines;
      if (capacityV < nVerticesBefore) {
        throw std::overflow_error("Overflow of vertices (primitive restart index)");
      }
      if (capacityI < nIn_dicesBefore) {
        throw std::overflow_error("Overflow of indices (primitive restart index)");
      }
      nVerticesBefore += iOffset;
      nIn_dicesBefore += nMLines;
    }
    if (!degenerated && !faceted && mStripsWrappingMethod == SurfaceStripsWrappingMethod::degenerate) {
      if (imax / 2 < nMLines) {
        throw std::overflow_error("Overflow of nMLinesDoubled (degenerate triangles)");
      }
      const itype nMLinesDoubled = 2 * nMLines;
      const itype capacity = imax - nMLinesDoubled;
      if (capacity < nIn_dicesBefore) {
        throw std::overflow_error("Overflow of indices (degenerate triangles)");
      }
      nIn_dicesBefore += nMLinesDoubled;
    }
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString("nVerticesBefore = %1, nIn_dicesBefore = %2, nLinesElVBefore = %3, nLinesElFBefore = %4")
                                                              .arg(nVerticesBefore)
                                                              .arg(nIn_dicesBefore)
                                                              .arg(nLinesElVBefore)
                                                              .arg(nLinesElFBefore),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
  } catch (const std::overflow_error& ex) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("Too many vertices for surface \"%1\" (nu = %2, nv = %3): %4."))
                                                              .arg(id)
                                                              .arg(nu)
                                                              .arg(nv)
                                                              .arg(ex.what()),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
    return geometry.release();
  }

  itype ne = 0; // TODO number of user-provided elements
  itype nw = 0; // TODO number of adjacent facets/windings
  itype na = 0; // TODO number of area/angle weights

  ne++;
  if (normalized) {
    ne++;
  }
  if (multicolored) {
    ne++;
  }

  if (!degenerated && !normalized) {
    switch (mNormalsAverageWeights) {
      case SurfaceNormalsAverageWeights::none:
        nw = 2;
        break;
      case SurfaceNormalsAverageWeights::equal:
      case SurfaceNormalsAverageWeights::area:
      case SurfaceNormalsAverageWeights::angle:
      case SurfaceNormalsAverageWeights::bothAreaAndAngle:
        nw = 6;
        break;
      default:
        break;
    }

    switch (mNormalsAverageWeights) {
      case SurfaceNormalsAverageWeights::none:
      case SurfaceNormalsAverageWeights::equal:
        na = 0;
        break;
      case SurfaceNormalsAverageWeights::area:
      case SurfaceNormalsAverageWeights::angle:
        na = 1;
        break;
      case SurfaceNormalsAverageWeights::bothAreaAndAngle:
        na = 2;
        break;
      default:
        break;
    }
  }

  const itype no = ne + nw + na; // TODO number of objects

  constexpr itype nc = 3; // TODO number of dimensions/coordinates

  try {
    const std::size_t smax = std::numeric_limits<std::size_t>::max();
    const std::size_t notnc = no * nc;
    const std::size_t notnctnu = no * nc * nu;
    const std::size_t notnctnutnv = no * nc * nu * nv;
    (void)notnctnutnv; // TODO use it or remove it
    if (smax / nc < (std::size_t)no) {
      throw std::overflow_error("Overflow of notnc");
    }
    if (smax / nu < notnc) {
      throw std::overflow_error("Overflow of notnctnu");
    }
    if (smax / nv < notnctnu) {
      throw std::overflow_error("Overflow of notnctnutnv");
    }
  } catch (const std::overflow_error& ex) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("Too many vertices for surface \"%1\" (no = %2, nc = %3, nu = %4, nv = %5): %6."))
                                                              .arg(id)
                                                              .arg(no)
                                                              .arg(nc)
                                                              .arg(nu)
                                                              .arg(nv)
                                                              .arg(ex.what()),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
    return geometry.release();
  }

#if 0 // Ragged array (with allocator overhead)
  ftype**  X = nullptr;
  ftype**  Y = nullptr;
  ftype**  Z = nullptr;
  ftype*** N = nullptr;
  ftype*** C = nullptr;
#else // Contiguous array (with pointers overhead)
  ftype**  Xu   = nullptr;
  ftype**  Yu   = nullptr;
  ftype**  Zu   = nullptr;
  ftype*** Nu   = nullptr;
  ftype*** Cu   = nullptr;
  ftype*   Xuv  = nullptr;
  ftype*   Yuv  = nullptr;
  ftype*   Zuv  = nullptr;
  ftype**  Nuv  = nullptr;
  ftype**  Cuv  = nullptr;
  ftype*   Nuvc = nullptr;
  ftype*   Cuvc = nullptr;
#endif

  try {
#if 0 // Ragged array (with allocator overhead)
    {
      X = new ftype* [nu]();
      Y = new ftype* [nu]();
      Z = new ftype* [nu]();
      for (itype u = 0; u < nu; u++) {
        X[u] = new ftype [nv]();
        Y[u] = new ftype [nv]();
        Z[u] = new ftype [nv]();
      }
    }
    if (normalized) {
      N = new ftype**[nu]();
      for (itype u = 0; u < nu; u++) {
        N[u] = new ftype*[nv]();
        for (itype v = 0; v < nv; v++) {
          N[u][v] = new ftype[nc]();
        }
      }
    }
    if (multicolored) {
      C = new ftype**[nu]();
      for (itype u = 0; u < nu; u++) {
        C[u] = new ftype*[nv]();
        for (itype v = 0; v < nv; v++) {
          C[u][v] = new ftype[nc]();
        }
      }
    }
#else // Contiguous array (with pointers overhead)
    {
      Xu = new ftype* [nu];
      Yu = new ftype* [nu];
      Zu = new ftype* [nu];
      Xuv = new ftype [nu * nv]{0};
      Yuv = new ftype [nu * nv]{0};
      Zuv = new ftype [nu * nv]{0};
    }
    if (normalized) {
      Nu = new ftype**[nu];
      Nuv = new ftype*[nu * nv];
      Nuvc = new ftype[nu * nv * nc]{0};
    }
    if (multicolored) {
      Cu = new ftype**[nu];
      Cuv = new ftype*[nu * nv];
      Cuvc = new ftype[nu * nv * nc]{0};
    }
    {
      for (itype u = 0; u < nu; u++, Xuv += nv, Yuv += nv, Zuv += nv) {
        Xu[u] = Xuv;
        Yu[u] = Yuv;
        Zu[u] = Zuv;
      }
    }
    if (normalized) {
      for (itype u = 0; u < nu; u++, Nuv += nv) {
        Nu[u] = Nuv;
        for (itype v = 0; v < nv; v++, Nuvc += nc) {
          Nu[u][v] = Nuvc;
        }
      }
    }
    if (multicolored) {
      for (itype u = 0; u < nu; u++, Cuv += nv) {
        Cu[u] = Cuv;
        for (itype v = 0; v < nv; v++, Cuvc += nc) {
          Cu[u][v] = Cuvc;
        }
      }
    }
#endif
  } catch (const std::bad_alloc& ex) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("Not enough memory to allocate vertices of surface \"%1\" (nu = %2, nv = %3): %4."))
                                                              .arg(id)
                                                              .arg(nu)
                                                              .arg(nv)
                                                              .arg(ex.what()),
                                                          Helper::scriptingKind,
                                                          Helper::errorLevel));
#if 0 // Ragged array (with allocator overhead)
    if (C) {
      for (itype u = 0; u < nu; u++) {
        if (C[u]) {
          for (itype v = 0; v < nv; v++) {
            if (C[u][v]) {
              delete[] C[u][v];
            }
          }
          delete[] C[u];
        }
      }
      delete[] C;
    }
    if (N) {
      for (itype u = 0; u < nu; u++) {
        if (N[u]) {
          for (itype v = 0; v < nv; v++) {
            if (N[u][v]) {
              delete[] N[u][v];
            }
          }
          delete[] N[u];
        }
      }
      delete[] N;
    }
    if (Z) {
      for (itype u = 0; u < nu; u++) {
        if (Z[u]) {
          delete[] Z[u];
        }
      }
      delete[] Z;
    }
    if (Y) {
      for (itype u = 0; u < nu; u++) {
        if (Y[u]) {
          delete[] Y[u];
        }
      }
      delete[] Y;
    }
    if (X) {
      for (itype u = 0; u < nu; u++) {
        if (X[u]) {
          delete[] X[u];
        }
      }
      delete[] X;
    }
#else // Contiguous array (with pointers overhead)
    delete[] Cuvc;
    delete[] Nuvc;
    delete[] Cuv;
    delete[] Nuv;
    delete[] Zuv;
    delete[] Yuv;
    delete[] Xuv;
    delete[] Cu;
    delete[] Nu;
    delete[] Zu;
    delete[] Yu;
    delete[] Xu;
#endif
    return geometry.release();
  }

#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
  ftype**  X = Xu;
  ftype**  Y = Yu;
  ftype**  Z = Zu;
  ftype*** N = Nu;
  ftype*** C = Cu;
#endif

  fakeTorus(nu, nv, X, Y, Z, N, C);
  // TODO: Fake multicolored surface

  constexpr ftype ps = 1; // Point size // FIXME vendor-specific annotation
  constexpr ftype lw = 1; // Line width // FIXME vendor-specific annotation
  constexpr ftype ns = 0.25; // Normal scale // FIXME vendor-specific annotation

  constexpr itype ri = 0; // Restart index

  constexpr itype x = 0; // Index of 1st coordinate
  constexpr itype y = 1; // Index of 2nd coordinate
  constexpr itype z = 2; // Index of 3rd coordinate

  const bool lines = nu == 1 || nv == 1; // Surface degenerated to a line strip
  const bool point = nu == 1 && nv == 1; // Surface degenerated to a single point

  const bool degenerate = mStripsWrappingMethod == SurfaceStripsWrappingMethod::degenerate; // Degenerate triangles
  const bool restart    = mStripsWrappingMethod == SurfaceStripsWrappingMethod::restart; // Primitive restart index

  const itype o = restart && !lines && (normalized || mNormalsAverageWeights != SurfaceNormalsAverageWeights::none) ? ri + 1 : 0; // Index offset

  /* Attributes */
  constexpr osg::StateAttribute::GLModeValue mode = osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE | osg::StateAttribute::PROTECTED;
  osg::ref_ptr<osg::StateSet> ss = geometry->getOrCreateStateSet();
  if (point) {
    ss->setAttributeAndModes(new osg::Point(ps), mode);
  }
  if ((lines && !point) || wireframe || mNormalsAnimationTypes != SurfaceNormalsAnimationTypes::none) {
    ss->setAttributeAndModes(new osg::LineWidth(lw), mode);
  }
  if (restart && !lines && (normalized || mNormalsAverageWeights != SurfaceNormalsAverageWeights::none)) {
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
  if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
    facetsCenters = new Vec3Array();
    facetsNormals = new Vec3Array();
  }

  if (lines) {
    const itype nutnv = nu * nv;
    const itype nutnvm1 = nutnv - 1;
    const ftype fnutnvm1 = (ftype)nutnvm1;

    const itype nVertices = nutnv;
    const itype nNormals  = nutnv;
    const itype nColors   = nutnv;
    const itype nTexels   = nutnv;
    const itype nIndices  = nutnv;

    vertices->reserve(nVertices);
    normals ->reserve(nNormals );
    colors  ->reserve(nColors  );
    texels  ->reserve(nTexels  );

    for (itype u = 0; u < nu; u++) {
      for (itype v = 0; v < nv; v++) {
        vertices->push_back(Vec3(X[u][v], Y[u][v], Z[u][v]));

        if (normalized) {
          normals->push_back(Vec3(N[u][v][x], N[u][v][y], N[u][v][z]));
        } else {
          normals->push_back(Vec3(X[u][v], Y[u][v], Z[u][v]));
        }

        if (multicolored) {
          colors->push_back(Vec4(C[u][v][x], C[u][v][y], C[u][v][z], opacity));
        } else {
          colors->push_back(Vec4(_color[x].exp, _color[y].exp, _color[z].exp, opacity));
        }

        if (point) {
          texels->push_back(Vec2(u, v));
        } else {
          texels->push_back(Vec2(u / fnutnvm1, v / fnutnvm1));
        }
      }
    }

    indices = new osg::DrawArrays(point ? osg::PrimitiveSet::POINTS : osg::PrimitiveSet::LINE_STRIP, 0, nIndices);
  } else {
    const itype num1 = nu - 1;
    const itype nvm1 = nv - 1;
    const itype num2 = nu - 2;
    const ftype fnum1 = (ftype)num1;
    const ftype fnvm1 = (ftype)nvm1;
    const itype num1tnvm1t2 = (num1 * nvm1) << 1;
    const itype opnutnv = o + nu * nv; // FIXME if not normalized and weighting is none then arrays need to be larger

    const itype nVertices = opnutnv;
    const itype nNormals  = opnutnv;
    const itype nColors   = opnutnv;
    const itype nTexels   = opnutnv;
    const itype nFacets   = num1tnvm1t2;
    const itype nIndices  = opnutnv; // TODO

    vertices->reserve(nVertices);
    normals ->reserve(nNormals );
    colors  ->reserve(nColors  );
    texels  ->reserve(nTexels  );

    if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
      facetsCenters->reserve(nFacets);
      facetsNormals->reserve(nFacets);
    }

    if (restart) { // Duplicate a dummy vertex (see OSG commit 353b18b)
      for (itype i = 0; i < o; i++) {
        vertices->push_back(Vec3());
        normals ->push_back(Vec3());
        colors  ->push_back(Vec4());
        texels  ->push_back(Vec2());
      }
    }

    /* Vertices */
    if (normalized || mNormalsAverageWeights != SurfaceNormalsAverageWeights::none) {
      for (itype u = 0; u < nu; u++) {
        for (itype v = 0; v < nv; v++) {
          vertices->push_back(Vec3(X[u][v], Y[u][v], Z[u][v]));
        }
      }
    } else {
      for (itype u = 0; u < num1; u++) {
        const itype up0 = u;
        const itype up1 = u + 1;
        for (itype v = 0; v < nvm1; v++) {
          const itype vp0 = v;
          const itype vp1 = v + 1;
          const ftype X00 = X[up0][vp0];
          const ftype Y00 = Y[up0][vp0];
          const ftype Z00 = Z[up0][vp0];
          const ftype X01 = X[up0][vp1];
          const ftype Y01 = Y[up0][vp1];
          const ftype Z01 = Z[up0][vp1];
          const ftype X10 = X[up1][vp0];
          const ftype Y10 = Y[up1][vp0];
          const ftype Z10 = Z[up1][vp0];
          const ftype X11 = X[up1][vp1];
          const ftype Y11 = Y[up1][vp1];
          const ftype Z11 = Z[up1][vp1];
          vertices->push_back(Vec3(X00, Y00, Z00));
          vertices->push_back(Vec3(X10, Y10, Z10));
          vertices->push_back(Vec3(X01, Y01, Z01));
          vertices->push_back(Vec3(X10, Y10, Z10));
          vertices->push_back(Vec3(X11, Y11, Z11));
          vertices->push_back(Vec3(X01, Y01, Z01));
        }
      }
    }

    /* Normals */
    if (normalized) {
      for (itype u = 0; u < nu; u++) {
        for (itype v = 0; v < nv; v++) {
          normals->push_back(Vec3(N[u][v][x], N[u][v][y], N[u][v][z]));
        }
      }
    }
    if (!normalized || mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
      itype nw = 6; // TODO number of adjacent facets/windings
      switch (mNormalsAverageWeights) {
        case SurfaceNormalsAverageWeights::none:
          nw = 2;
          break;
        default:
          break;
      }

      ftype**** A = nullptr;

      if (!normalized && nw > 0) {
#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
        ftype**** Au = nullptr;
        ftype*** Auv = nullptr;
        ftype** Auvw = nullptr;
        ftype* Auvwc = nullptr;
#endif

        try {
#if 0 // Ragged array (with allocator overhead)
          A = new ftype***[nu]();
          for (itype u = 0; u < nu; u++) {
            A[u] = new ftype**[nv]();
            for (itype v = 0; v < nv; v++) {
              A[u][v] = new ftype*[nw]();
              for (itype w = 0; w < nw; w++) {
                A[u][v][w] = new ftype[nc]();
              }
            }
          }
#else // Contiguous array (with pointers overhead)
          Au = new ftype***[nu];
          Auv = new ftype**[nu * nv];
          Auvw = new ftype*[nu * nv * nw];
          Auvwc = new ftype[nu * nv * nw * nc]{0};
          for (itype u = 0; u < nu; u++, Auv += nv) {
            Au[u] = Auv;
            for (itype v = 0; v < nv; v++, Auvw += nw) {
              Au[u][v] = Auvw;
              for (itype w = 0; w < nw; w++, Auvwc += nc) {
                Au[u][v][w] = Auvwc;
              }
            }
          }
#endif
        } catch (const std::bad_alloc& ex) {
          MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                                QString(QObject::tr("Not enough memory to allocate adjacent facets normals of surface \"%1\" (nu = %2, nv = %3, nw = %4, nc = %5): %6."))
                                                                    .arg(id)
                                                                    .arg(nu)
                                                                    .arg(nv)
                                                                    .arg(nw)
                                                                    .arg(nc)
                                                                    .arg(ex.what()),
                                                                Helper::scriptingKind,
                                                                Helper::errorLevel));
#if 0 // Ragged array (with allocator overhead)
          if (A) {
            for (itype u = 0; u < nu; u++) {
              if (A[u]) {
                for (itype v = 0; v < nv; v++) {
                  if (A[u][v]) {
                    for (itype w = 0; w < nw; w++) {
                      if (A[u][v][w]) {
                        delete[] A[u][v][w];
                      }
                    }
                    delete[] A[u][v];
                  }
                }
                delete[] A[u];
              }
            }
            delete[] A;
          }
          if (multicolored) {
            for (itype u = 0; u < nu; u++) {
              for (itype v = 0; v < nv; v++) {
                delete[] C[u][v];
              }
              delete[] C[u];
            }
            delete[] C;
          }
          if (normalized) {
            for (itype u = 0; u < nu; u++) {
              for (itype v = 0; v < nv; v++) {
                delete[] N[u][v];
              }
              delete[] N[u];
            }
            delete[] N;
          }
          {
            for (itype u = 0; u < nu; u++) {
              delete[] Z[u];
              delete[] Y[u];
              delete[] X[u];
            }
            delete[] Z;
            delete[] Y;
            delete[] X;
          }
#else // Contiguous array (with pointers overhead)
          delete[] Auvwc;
          delete[] Auvw;
          delete[] Auv;
          delete[] Au;
          if (multicolored) {
            delete[] C[0][0];
            delete[] C[0];
            delete[] C;
          }
          if (normalized) {
            delete[] N[0][0];
            delete[] N[0];
            delete[] N;
          }
          {
            delete[] Z[0];
            delete[] Y[0];
            delete[] X[0];
            delete[] Z;
            delete[] Y;
            delete[] X;
          }
#endif
          return geometry.release();
        }

#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
        A = Au;
#endif
      }

      itype na = 0; // TODO number of area/angle weights
      switch (mNormalsAverageWeights) {
        case SurfaceNormalsAverageWeights::area:
        case SurfaceNormalsAverageWeights::angle:
          na = 1;
          break;
        case SurfaceNormalsAverageWeights::bothAreaAndAngle:
          na = 2;
          break;
        default:
          break;
      }

      ftype**** W = nullptr; // FIXME allocate together with A as a whole contiguous array?

      if (!normalized && na > 0) {
#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
        ftype**** Wu = nullptr;
        ftype*** Wuv = nullptr;
        ftype** Wuvw = nullptr;
        ftype* Wuvwa = nullptr;
#endif

        try {
#if 0 // Ragged array (with allocator overhead)
          W = new ftype***[nu]();
          for (itype u = 0; u < nu; u++) {
            W[u] = new ftype**[nv]();
            for (itype v = 0; v < nv; v++) {
              W[u][v] = new ftype*[nw]();
              for (itype w = 0; w < nw; w++) {
                W[u][v][w] = new ftype[na]();
              }
            }
          }
#else // Contiguous array (with pointers overhead)
          Wu = new ftype***[nu];
          Wuv = new ftype**[nu * nv];
          Wuvw = new ftype*[nu * nv * nw];
          Wuvwa = new ftype[nu * nv * nw * na]{0};
          for (itype u = 0; u < nu; u++, Wuv += nv) {
            Wu[u] = Wuv;
            for (itype v = 0; v < nv; v++, Wuvw += nw) {
              Wu[u][v] = Wuvw;
              for (itype w = 0; w < nw; w++, Wuvwa += na) {
                Wu[u][v][w] = Wuvwa;
              }
            }
          }
#endif
        } catch (const std::bad_alloc& ex) {
          MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                                QString(QObject::tr("Not enough memory to allocate weights for normals of surface \"%1\" (nu = %2, nv = %3, nw = %4, na = %5): %6."))
                                                                    .arg(id)
                                                                    .arg(nu)
                                                                    .arg(nv)
                                                                    .arg(nw)
                                                                    .arg(na)
                                                                    .arg(ex.what()),
                                                                Helper::scriptingKind,
                                                                Helper::errorLevel));
#if 0 // Ragged array (with allocator overhead)
          if (W) {
            for (itype u = 0; u < nu; u++) {
              if (W[u]) {
                for (itype v = 0; v < nv; v++) {
                  if (W[u][v]) {
                    for (itype w = 0; w < nw; w++) {
                      if (W[u][v][w]) {
                        delete[] W[u][v][w];
                      }
                    }
                    delete[] W[u][v];
                  }
                }
                delete[] W[u];
              }
            }
            delete[] W;
          }
          {
            for (itype u = 0; u < nu; u++) {
              for (itype v = 0; v < nv; v++) {
                for (itype w = 0; w < nw; w++) {
                  delete[] A[u][v][w];
                }
                delete[] A[u][v];
              }
              delete[] A[u];
            }
            delete[] A;
          }
          if (multicolored) {
            for (itype u = 0; u < nu; u++) {
              for (itype v = 0; v < nv; v++) {
                delete[] C[u][v];
              }
              delete[] C[u];
            }
            delete[] C;
          }
          if (normalized) {
            for (itype u = 0; u < nu; u++) {
              for (itype v = 0; v < nv; v++) {
                delete[] N[u][v];
              }
              delete[] N[u];
            }
            delete[] N;
          }
          {
            for (itype u = 0; u < nu; u++) {
              delete[] Z[u];
              delete[] Y[u];
              delete[] X[u];
            }
            delete[] Z;
            delete[] Y;
            delete[] X;
          }
#else // Contiguous array (with pointers overhead)
          delete[] Wuvwa;
          delete[] Wuvw;
          delete[] Wuv;
          delete[] Wu;
          {
            delete[] A[0][0][0];
            delete[] A[0][0];
            delete[] A[0];
            delete[] A;
          }
          if (multicolored) {
            delete[] C[0][0];
            delete[] C[0];
            delete[] C;
          }
          if (normalized) {
            delete[] N[0][0];
            delete[] N[0];
            delete[] N;
          }
          {
            delete[] Z[0];
            delete[] Y[0];
            delete[] X[0];
            delete[] Z;
            delete[] Y;
            delete[] X;
          }
#endif
          return geometry.release();
        }

#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
        W = Wu;
#endif
      }

      constexpr itype i = 0; // TODO index of first adjacent facet for top right vertex
      constexpr itype j = 2; // TODO index of first adjacent facet for top left vertex
      constexpr itype k = 4; // TODO index of first adjacent facet for bottom right vertex
      constexpr itype l = 1; // TODO index of second adjacent facet for top left vertex
      constexpr itype m = 3; // TODO index of second adjacent facet for bottom left vertex
      constexpr itype n = 5; // TODO index of second adjacent facet for bottom right vertex

      for (itype u = 0; u < num1; u++) {
        const itype up0 = u;
        const itype up1 = u + 1;
        for (itype v = 0; v < nvm1; v++) {
          const itype vp0 = v;
          const itype vp1 = v + 1;
          const ftype X00 = X[up0][vp0];
          const ftype Y00 = Y[up0][vp0];
          const ftype Z00 = Z[up0][vp0];
          const ftype X01 = X[up0][vp1];
          const ftype Y01 = Y[up0][vp1];
          const ftype Z01 = Z[up0][vp1];
          const ftype X10 = X[up1][vp0];
          const ftype Y10 = Y[up1][vp0];
          const ftype Z10 = Z[up1][vp0];
          const ftype X11 = X[up1][vp1];
          const ftype Y11 = Y[up1][vp1];
          const ftype Z11 = Z[up1][vp1];
          const ftype du1[nc] = {X10 - X00, Y10 - Y00, Z10 - Z00};
          const ftype du2[nc] = {X01 - X11, Y01 - Y11, Z01 - Z11};
          const ftype dv1[nc] = {X01 - X00, Y01 - Y00, Z01 - Z00};
          const ftype dv2[nc] = {X10 - X11, Y10 - Y11, Z10 - Z11};
          const ftype cross1[nc] = {du1[y] * dv1[z] - du1[z] * dv1[y], du1[z] * dv1[x] - du1[x] * dv1[z], du1[x] * dv1[y] - du1[y] * dv1[x]};
          const ftype cross2[nc] = {du2[y] * dv2[z] - du2[z] * dv2[y], du2[z] * dv2[x] - du2[x] * dv2[z], du2[x] * dv2[y] - du2[y] * dv2[x]};
          const ftype length1 = std::sqrt(cross1[x] * cross1[x] + cross1[y] * cross1[y] + cross1[z] * cross1[z]);
          const ftype length2 = std::sqrt(cross2[x] * cross2[x] + cross2[y] * cross2[y] + cross2[z] * cross2[z]);
          const ftype invnorm1 = length1 > 0 ? 1 / length1 : 0;
          const ftype invnorm2 = length2 > 0 ? 1 / length2 : 0;
          const ftype normal1[nc] = {cross1[x] * invnorm1, cross1[y] * invnorm1, cross1[z] * invnorm1};
          const ftype normal2[nc] = {cross2[x] * invnorm2, cross2[y] * invnorm2, cross2[z] * invnorm2};
          if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
            const ftype center1[nc] = {(X00 + X10 + X01) / 3, (Y00 + Y10 + Y01) / 3, (Z00 + Z10 + Z01) / 3};
            const ftype center2[nc] = {(X10 + X11 + X01) / 3, (Y10 + Y11 + Y01) / 3, (Z10 + Z11 + Z01) / 3};
            facetsCenters->push_back(Vec3(center1[x], center1[y], center1[z]));
            facetsCenters->push_back(Vec3(center2[x], center2[y], center2[z]));
            facetsNormals->push_back(Vec3(normal1[x], normal1[y], normal1[z]));
            facetsNormals->push_back(Vec3(normal2[x], normal2[y], normal2[z]));
          }
          if (!normalized) {
            if (mNormalsAverageWeights == SurfaceNormalsAverageWeights::none) {
              ftype** A00 = A[up0][vp0];
              ftype** A11 = A[up1][vp1];
              for (itype c = 0; c < nc; c++) {
                A00[i][c] = normal1[c];
                A11[l][c] = normal2[c];
              }
            } else {
              ftype** A00 = A[up0][vp0];
              ftype** A01 = A[up0][vp1];
              ftype** A10 = A[up1][vp0];
              ftype** A11 = A[up1][vp1];
              for (itype c = 0; c < nc; c++) {
                A00[i][c] = normal1[c];
                A10[j][c] = normal1[c];
                A01[k][c] = normal1[c];
                A10[l][c] = normal2[c];
                A11[m][c] = normal2[c];
                A01[n][c] = normal2[c];
              }
            }
            if (mNormalsAverageWeights & SurfaceNormalsAverageWeights::bothAreaAndAngle) {
              itype a = 0;
              ftype** W00 = W[up0][vp0];
              ftype** W01 = W[up0][vp1];
              ftype** W10 = W[up1][vp0];
              ftype** W11 = W[up1][vp1];
              if (mNormalsAverageWeights & SurfaceNormalsAverageWeights::area) {
                const ftype area1 = length1 / 2; // TODO surface area = triangle area = half the norm of the cross product
                const ftype area2 = length2 / 2; // TODO surface area = triangle area = half the norm of the cross product
                W00[i][a] = area1;
                W10[j][a] = area1;
                W01[k][a] = area1;
                W10[l][a] = area2;
                W11[m][a] = area2;
                W01[n][a] = area2;
                a++;
              }
              if (mNormalsAverageWeights & SurfaceNormalsAverageWeights::angle) {
                const ftype duv[nc] = {X01 - X10, Y01 - Y10, Z01 - Z10};
                const ftype dvu[nc] = {X10 - X01, Y10 - Y01, Z10 - Z01};
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
                const ftype angle11 = std::acos(length11 > 0 ? dot11 / length11 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
                const ftype angle12 = std::acos(length12 > 0 ? dot12 / length12 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
                const ftype angle13 = std::acos(length13 > 0 ? dot13 / length13 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
                const ftype angle21 = std::acos(length21 > 0 ? dot21 / length21 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
                const ftype angle22 = std::acos(length22 > 0 ? dot22 / length22 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
                const ftype angle23 = std::acos(length23 > 0 ? dot23 / length23 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
                W00[i][a] = angle11;
                W10[j][a] = angle12;
                W01[k][a] = angle13;
                W10[l][a] = angle21;
                W11[m][a] = angle22;
                W01[n][a] = angle23;
                a++;
              }
            }
          }
        }
      }

      if (!normalized) {
        if (mNormalsAverageWeights != SurfaceNormalsAverageWeights::none) {
          if (mClosenessCheckState == SurfaceClosenessCheckState::active) {
            constexpr itype nut0 = 0;
            constexpr itype nvt0 = 0;
            for (itype u = 0; u < nu; u++) {
              if (X[u][nvt0] == X[u][nvm1] &&
                  Y[u][nvt0] == Y[u][nvm1] &&
                  Z[u][nvt0] == Z[u][nvm1]) {
                for (itype c = 0; c < nc; c++) {
                  A[u][nvm1][i][c] = A[u][nvt0][i][c];
                  A[u][nvm1][l][c] = A[u][nvt0][l][c];
                  A[u][nvm1][j][c] = A[u][nvt0][j][c];
                  A[u][nvt0][m][c] = A[u][nvm1][m][c];
                  A[u][nvt0][k][c] = A[u][nvm1][k][c];
                  A[u][nvt0][n][c] = A[u][nvm1][n][c];
                }
                for (itype a = 0; a < na; a++) {
                  W[u][nvm1][i][a] = W[u][nvt0][i][a];
                  W[u][nvm1][l][a] = W[u][nvt0][l][a];
                  W[u][nvm1][j][a] = W[u][nvt0][j][a];
                  W[u][nvt0][m][a] = W[u][nvm1][m][a];
                  W[u][nvt0][k][a] = W[u][nvm1][k][a];
                  W[u][nvt0][n][a] = W[u][nvm1][n][a];
                }
              }
            }
            for (itype v = 0; v < nv; v++) {
              if (X[nut0][v] == X[num1][v] &&
                  Y[nut0][v] == Y[num1][v] &&
                  Z[nut0][v] == Z[num1][v]) {
                for (itype c = 0; c < nc; c++) {
                  A[num1][v][i][c] = A[nut0][v][i][c];
                  A[num1][v][k][c] = A[nut0][v][k][c];
                  A[num1][v][n][c] = A[nut0][v][n][c];
                  A[nut0][v][l][c] = A[num1][v][l][c];
                  A[nut0][v][j][c] = A[num1][v][j][c];
                  A[nut0][v][m][c] = A[num1][v][m][c];
                }
                for (itype a = 0; a < na; a++) {
                  W[num1][v][i][a] = W[nut0][v][i][a];
                  W[num1][v][k][a] = W[nut0][v][k][a];
                  W[num1][v][n][a] = W[nut0][v][n][a];
                  W[nut0][v][l][a] = W[num1][v][l][a];
                  W[nut0][v][j][a] = W[num1][v][j][a];
                  W[nut0][v][m][a] = W[num1][v][m][a];
                }
              }
            }
          }

          for (itype u = 0; u < nu; u++) {
            for (itype v = 0; v < nv; v++) {
              ftype normal[nc] = {0};
              for (itype w = 0; w < nw; w++) {
                ftype weight = 1;
                for (itype a = 0; a < na; a++) {
                  weight *= W[u][v][w][a];
                }
                normal[x] += A[u][v][w][x] * weight;
                normal[y] += A[u][v][w][y] * weight;
                normal[z] += A[u][v][w][z] * weight;
              }
              const ftype length = std::sqrt(normal[x] * normal[x] + normal[y] * normal[y] + normal[z] * normal[z]);
              const ftype invnorm = length > 0 ? 1 / length : 0;
              normal[x] *= invnorm;
              normal[y] *= invnorm;
              normal[z] *= invnorm;
              normals->push_back(Vec3(normal[x], normal[y], normal[z]));
            }
          }
        } else {
          for (itype u = 0; u < num1; u++) {
            const itype up0 = u;
            const itype up1 = u + 1;
            for (itype v = 0; v < nvm1; v++) {
              const itype vp0 = v;
              const itype vp1 = v + 1;
              ftype* normal1 = A[up0][vp0][i];
              ftype* normal2 = A[up1][vp1][l];
              normals->insert(normals->end(), 3, Vec3(normal1[x], normal1[y], normal1[z]));
              normals->insert(normals->end(), 3, Vec3(normal2[x], normal2[y], normal2[z]));
            }
          }
        }

        if (na > 0) {
#if 0 // Ragged array (with allocator overhead)
          for (itype u = 0; u < nu; u++) {
            for (itype v = 0; v < nv; v++) {
              for (itype w = 0; w < nw; w++) {
                delete[] W[u][v][w];
              }
              delete[] W[u][v];
            }
            delete[] W[u];
          }
          delete[] W;
#else // Contiguous array (with pointers overhead)
          delete[] W[0][0][0];
          delete[] W[0][0];
          delete[] W[0];
          delete[] W;
#endif
        }
        if (nw > 0) {
#if 0 // Ragged array (with allocator overhead)
          for (itype u = 0; u < nu; u++) {
            for (itype v = 0; v < nv; v++) {
              for (itype w = 0; w < nw; w++) {
                delete[] A[u][v][w];
              }
              delete[] A[u][v];
            }
            delete[] A[u];
          }
          delete[] A;
#else // Contiguous array (with pointers overhead)
          delete[] A[0][0][0];
          delete[] A[0][0];
          delete[] A[0];
          delete[] A;
#endif
        }
      }
    }

    /* Colors */
    if (multicolored) {
      if (normalized || mNormalsAverageWeights != SurfaceNormalsAverageWeights::none) {
        for (itype u = 0; u < nu; u++) {
          for (itype v = 0; v < nv; v++) {
            colors->push_back(Vec4(C[u][v][x], C[u][v][y], C[u][v][z], opacity));
          }
        }
      } else {
        for (itype u = 0; u < num1; u++) {
          const itype up0 = u;
          const itype up1 = u + 1;
          for (itype v = 0; v < nvm1; v++) {
            const itype vp0 = v;
            const itype vp1 = v + 1;
            ftype* C00 = C[up0][vp0];
            ftype* C01 = C[up0][vp1];
            ftype* C10 = C[up1][vp0];
            ftype* C11 = C[up1][vp1];
            colors->push_back(Vec4(C00[x], C00[y], C00[z], opacity));
            colors->push_back(Vec4(C10[x], C10[y], C10[z], opacity));
            colors->push_back(Vec4(C01[x], C01[y], C01[z], opacity));
            colors->push_back(Vec4(C10[x], C10[y], C10[z], opacity));
            colors->push_back(Vec4(C11[x], C11[y], C11[z], opacity));
            colors->push_back(Vec4(C01[x], C01[y], C01[z], opacity));
          }
        }
      }
    } else {
      colors->clear();
      colors->setBinding(osg::Array::BIND_OVERALL);
      colors->push_back(Vec4(_color[x].exp, _color[y].exp, _color[z].exp, opacity));
    }

    /* Texels */
    if (normalized || mNormalsAverageWeights != SurfaceNormalsAverageWeights::none) {
      for (itype u = 0; u < nu; u++) {
        for (itype v = 0; v < nv; v++) {
          texels->push_back(Vec2(u / fnum1, v / fnvm1));
        }
      }
    } else {
      for (itype u = 0; u < num1; u++) {
        const ftype Tu0 =  u      / fnum1;
        const ftype Tu1 = (u + 1) / fnum1;
        for (itype v = 0; v < nvm1; v++) {
          const ftype Tv0 =  v      / fnvm1;
          const ftype Tv1 = (v + 1) / fnvm1;
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
#define SURFACE_FIRST_V()                       \
  if (degenerate) {                             \
    strip->addElement(up0tnvpopv);              \
  }
#define SURFACE_LOOP_V()                        \
  for (; up0tnvpopv < up1tnvpo; up0tnvpopv++) { \
    strip->addElement(up0tnvpopv);              \
    strip->addElement(up0tnvpopv + nv);         \
  }
#define SURFACE_LAST_V()                        \
  if (degenerate) {                             \
    strip->addElement(up0tnvpopv + nvm1);       \
  } else if (restart) {                         \
    strip->addElement(ri);                      \
  }
    if (normalized || mNormalsAverageWeights != SurfaceNormalsAverageWeights::none) {
      indices = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLE_STRIP);
      osg::ref_ptr<osg::DrawElements> strip = indices->getDrawElements();
      strip->reserveElements(nIndices);
      const itype opnum2tnv = o + num2 * nv;
      const bool num2g0 = num2 > 0;
      itype up0tnvpopv = o;
      itype up1tnvpo = nv + o;
      {
        SURFACE_LOOP_V()
      }
      if (num2g0) {
        SURFACE_LAST_V()
      }
      for (up1tnvpo += nv; up0tnvpopv < opnum2tnv; up1tnvpo += nv) {
        SURFACE_FIRST_V()
        SURFACE_LOOP_V()
        SURFACE_LAST_V()
      }
      if (num2g0) {
        SURFACE_FIRST_V()
        SURFACE_LOOP_V()
      }
    } else {
      indices = new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, o, (itype)vertices->size() - o);
    }
  }

  /* Debug */
#define SURFACE_DEBUG_N(o, e, v, n, r, g, b)    \
  const itype l = vertices->size();             \
  const itype c = (e - o) << 1;                 \
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
  if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices) {
    const itype offset = o;
    const itype size = vertices->size();
    SURFACE_DEBUG_N(
        offset, size,
        vertices,
        normals,
        1, 0, 0);
  }
  if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
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
  geometry->setTexCoordArray(0, texels.get());
  geometry->addPrimitiveSet(indices.get());

  const itype idxV = mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices ? 0 : -1;
  const itype idxF = mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets
                   ? mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices ? 1 : 0 : -1;
  const itype idxT = mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets
                   ? mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices ? 2 : 1
                   : mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices ? 1 : 0;
  const itype nLinesElVAfter = idxV < 0 ? 0 : geometry->getPrimitiveSet(idxV)->getNumIndices();
  const itype nLinesElFAfter = idxF < 0 ? 0 : geometry->getPrimitiveSet(idxF)->getNumIndices();
  const itype nIn_dicesAfter = idxT < 0 ? 0 : geometry->getPrimitiveSet(idxT)->getNumIndices();
  const itype nVerticesAfter = geometry->getVertexArray()->getNumElements();
  MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                        QString("nVerticesAfter = %1, nIn_dicesAfter = %2, nLinesElVAfter = %3, nLinesElFAfter = %4")
                                                            .arg(nVerticesAfter)
                                                            .arg(nIn_dicesAfter)
                                                            .arg(nLinesElVAfter)
                                                            .arg(nLinesElFAfter),
                                                        Helper::scriptingKind,
                                                        Helper::errorLevel));

#if 0 // Ragged array (with allocator overhead)
  if (multicolored) {
    for (itype u = 0; u < nu; u++) {
      for (itype v = 0; v < nv; v++) {
        delete[] C[u][v];
      }
      delete[] C[u];
    }
    delete[] C;
  }
  if (normalized) {
    for (itype u = 0; u < nu; u++) {
      for (itype v = 0; v < nv; v++) {
        delete[] N[u][v];
      }
      delete[] N[u];
    }
    delete[] N;
  }
  {
    for (itype u = 0; u < nu; u++) {
      delete[] Z[u];
      delete[] Y[u];
      delete[] X[u];
    }
    delete[] Z;
    delete[] Y;
    delete[] X;
  }
#else // Contiguous array (with pointers overhead)
  if (multicolored) {
    delete[] C[0][0];
    delete[] C[0];
    delete[] C;
  }
  if (normalized) {
    delete[] N[0][0];
    delete[] N[0];
    delete[] N;
  }
  {
    delete[] Z[0];
    delete[] Y[0];
    delete[] X[0];
    delete[] Z;
    delete[] Y;
    delete[] X;
  }
#endif

  return geometry.release();
}

/*
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
*/

/*
See https://www.glprogramming.com/red/chapter07.html#name2 for an example torus drawn with GL_QUAD_STRIP
*/

/*
See https://www.glprogramming.com/red/chapter02.html#name2 for examples of invalid quadrilaterals
and https://www.glprogramming.com/red/chapter02.html#name8 for some hints to approximate surfaces
"Since OpenGL vertices are always three-dimensional, the points forming the boundary of a particular polygon don't necessarily lie on the same plane in space. (Of course, they do in many cases - if all the z coordinates are zero, for example, or if the polygon is a triangle.) If a polygon's vertices don't lie in the same plane, then after various rotations in space, changes in the viewpoint, and projection onto the display screen, the points might no longer form a simple convex polygon. For example, imagine a four-point quadrilateral where the points are slightly out of plane, and look at it almost edge-on. You can get a nonsimple polygon that resembles a bow tie, as shown in Figure 2-4, which isn't guaranteed to be rendered correctly. This situation isn't all that unusual if you approximate curved surfaces by quadrilaterals made of points lying on the true surface. You can always avoid the problem by using triangles, since any three points always lie on a plane."
=> Do not use GL_QUAD_STRIP but GL_TRIANGLE_STRIP instead (with primitive restart index if necessary)
*/

/*
See https://www.glprogramming.com/red/appendixe.html#name2 for computing approximate normals
=> Images in MSL documentation seem to just be faceted (otherwise this needs to be done after setting all vertices to have access to neighboring facets' normals)
*/

/*
See https://www.learnopengles.com/tag/triangle-strips to making use of index buffer objects
=> They suggest using degenerate triangles to connect consecutive row strips, but primitive restart index is an alternative (not available in GLES until 3.1)
*/

/*
Enum for two methods to connect consecutive row strips: degenerate triangles, primitive restart index
Enum for two methods to compute surface normal vectors: faceted, averaged (average with unit weights or weighted by angle between each face and current vertex BUT this requires the true normal which we do not have since this is what we are trying to approximate, so just average with equal weights
=> Ah! In fact there are some methods out there to do that:
 - https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.534.7649&rep=rep1&type=pdf
 - https://www.bytehazard.com/articles/vertnorm.html
 - https://knowledge.autodesk.com/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2022/ENU/Maya-Modeling/files/GUID-232E99F8-96B4-4870-8BA0-4887C1C8F0F2-htm.html
 - https://github.com/pmnuckels/wnormals
 - https://en.wikipedia.org/wiki/Vertex_normal
)
NOTE: Since we do indexed vertex rendering we have each vertex only once in existence (and pick it multiple times through its index) so we cannot have multiple normals per vertex (which would be necessary to actually set the facet normal to all the vertices of the corresponding face, and each vertex belonging to several faces) but only one normal per vertex, therefore we are forced to average the normal vectors one way or another!
=> By default, unit weights; other weighing possibilities are: area only, angle only, area and angle

Allocating and passing multidimensional arrays on the heap:
- https://c-faq.com/aryptr/dynmuldimary.html
- https://c-faq.com/aryptr/ary2dfunc3.html
- https://c-faq.com/aryptr/pass2dary.html
- https://c-faq.com/aryptr/fn33.html
- https://stackoverflow.com/a/21944048 as well as https://c-faq.com/aryptr/dynmuldimary.html for contiguous multidimensional arrays (surely better!)
*/
