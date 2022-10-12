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

#include <type_traits>

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
      mNormalsAnimationTypes(SurfaceNormalsAnimationTypes::bothVerticesAndFacets), // FIXME should be none by default
      _nu(VisualizerAttribute(0.0)),
      _nv(VisualizerAttribute(0.0)),
      _wireframe(VisualizerAttribute(0.0)),
      _normalized(VisualizerAttribute(0.0)),
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
  std::cout << "multicolored " << _multicolored.getValueString() << std::endl;
  std::cout << "transparency " << _transparency.getValueString() << std::endl;
}

#include <cmath>
#include <new>

#include <QOpenGLContext> // must be included before OSG headers

#include <osg/ref_ptr>
#include <osg/Array>
#include <osg/Geometry>
#include <osg/LineWidth>
#include <osg/PolygonMode>
#include <osg/PrimitiveRestartIndex>
#include <osg/PrimitiveSet>
#include <osg/StateAttribute>
#include <osg/StateSet>

#include "Modeling/MessagesWidget.h"
#include "Util/Helper.h"

void SurfaceObject::fakeTorus(const itype nu, const itype nv, ftype** X, ftype** Y, ftype** Z, ftype*** N, ftype*** C) const
{
  (void)N;
  (void)C;

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
    alpha = startAngle + (stopAngle - startAngle) * u / (nu - 1);
    for (itype v = 0; v < nv; v++) {
      beta = phi_start + (phi_stop - phi_start) * v / (nv - 1);
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

osg::Geometry* SurfaceObject::drawGeometry() const
{
  osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();

  constexpr itype nw = 6; // TODO number of adjacent facets/windings
  constexpr itype nc = 3; // TODO number of dimensions/coordinates
  const     itype nu = _nu.exp;
  const     itype nv = _nv.exp;
  const bool wireframe = _wireframe.exp;
  const bool normalized = _normalized.exp;
  const bool multicolored = _multicolored.exp;
  const ftype opacity = 1.0 - _transparency.exp;

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
  } catch (const std::bad_alloc&) {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                          QString(QObject::tr("Not enough memory to allocate vertices of surface \"%1\" (nu = %2, nv = %3)."))
                                                              .arg(_id.c_str())
                                                              .arg(QString::number(nu))
                                                              .arg(QString::number(nv)),
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

  // TODO: Guarantee nu & nv > 1 => Draw line/point otherwise (could technically be a valid degenerate surface)
  // TODO: Guarantee nu * nv < int max minus index offset, or explicitly use a larger type (long long (unsigned) int) if the product is not used as an index

  /* Attributes */ // TODO: #define instead of evaluating conditions every time => instead, let the user choose through vendor-specific annotations
  const bool degenerate = mStripsWrappingMethod == SurfaceStripsWrappingMethod::degenerate; // Degenerate triangles
  const bool restart    = mStripsWrappingMethod == SurfaceStripsWrappingMethod::restart; // Primitive restart index
  constexpr itype r = 0; // Restart index
  const itype o = restart ? r + 1 : 0; // Index offset
  constexpr osg::StateAttribute::GLModeValue mode = osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE | osg::StateAttribute::PROTECTED;
  osg::ref_ptr<osg::StateSet> ss = geometry->getOrCreateStateSet();
  if (restart) {
    ss->setMode(GL_PRIMITIVE_RESTART, mode);
    ss->setAttributeAndModes(new osg::PrimitiveRestartIndex(r), mode);
  }
  if (wireframe) {
    ss->setAttributeAndModes(new osg::PolygonMode(osg::PolygonMode::Face::FRONT_AND_BACK, osg::PolygonMode::Mode::LINE), mode);
  }

  constexpr itype x = 0; // TODO index of 1st coordinate
  constexpr itype y = 1; // TODO index of 2nd coordinate
  constexpr itype z = 2; // TODO index of 3rd coordinate

  const itype num1 = nu - 1;
  const itype nvm1 = nv - 1;
  const itype num2 = nu - 2;
  const itype nutnv = nu * nv;
  const itype opnutnv = o + nutnv;

  const itype nVertices = opnutnv;
  const itype nNormals  = opnutnv;
  const itype nColors   = opnutnv;
  const itype nTexels   = opnutnv;
  const itype nIndices  = opnutnv; // TODO

  osg::ref_ptr<Vec3Array> vertices = new Vec3Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<Vec3Array> normals  = new Vec3Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<Vec4Array> colors   = new Vec4Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<Vec2Array> texels   = new Vec2Array(osg::Array::BIND_PER_VERTEX);
  osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLE_STRIP);

  vertices->reserve(nVertices);
  normals ->reserve(nNormals );
  colors  ->reserve(nColors  );
  texels  ->reserve(nTexels  );
  indices ->reserve(nIndices );

  if (restart) { // Duplicate the first vertex (see OSG commit 353b18b)
    for (itype i = 0; i < o; i++) {
      vertices->push_back(Vec3());
      normals ->push_back(Vec3());
      colors  ->push_back(Vec4());
      texels  ->push_back(Vec2());
    }
  }

  osg::ref_ptr<Vec3Array> facetsCenters = nullptr;
  osg::ref_ptr<Vec3Array> facetsNormals = nullptr;
  if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
    facetsCenters = new Vec3Array();
    facetsNormals = new Vec3Array();
    facetsCenters->reserve(2 * num1 * nvm1); // TODO constexpr
    facetsNormals->reserve(2 * num1 * nvm1); // TODO constexpr
  }

  /* Vertices */
  for (itype u = 0; u < nu; u++) {
    for (itype v = 0; v < nv; v++) {
      vertices->push_back(Vec3(X[u][v], Y[u][v], Z[u][v]));
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
  if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets || !normalized) {
    Vec3::value_type**** A = nullptr;

    if (!normalized) {
#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
      Vec3::value_type**** Au = nullptr;
      Vec3::value_type*** Auv = nullptr;
      Vec3::value_type** Auvw = nullptr;
      Vec3::value_type* Auvwc = nullptr;
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
      } catch (const std::bad_alloc&) {
        MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                              QString(QObject::tr("Not enough memory to allocate adjacent facets normals of surface \"%1\" (nu = %2, nv = %3, nw = %4, nc = %5)."))
                                                                  .arg(_id.c_str())
                                                                  .arg(QString::number(nu))
                                                                  .arg(QString::number(nv))
                                                                  .arg(QString::number(nw))
                                                                  .arg(QString::number(nc)),
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

    Vec3::value_type**** W = nullptr; // FIXME allocate together with A as a whole contiguous array?

    if (!normalized && na > 0) {
#if 0 // Ragged array (with allocator overhead)
#else // Contiguous array (with pointers overhead)
      Vec3::value_type**** Wu = nullptr;
      Vec3::value_type*** Wuv = nullptr;
      Vec3::value_type** Wuvw = nullptr;
      Vec3::value_type* Wuvwa = nullptr;
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
      } catch (const std::bad_alloc&) {
        MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica,
                                                              QString(QObject::tr("Not enough memory to allocate weights for normals of surface \"%1\" (nu = %2, nv = %3, nw = %4, na = %5)."))
                                                                  .arg(_id.c_str())
                                                                  .arg(QString::number(nu))
                                                                  .arg(QString::number(nv))
                                                                  .arg(QString::number(nw))
                                                                  .arg(QString::number(na)),
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
        // TODO: Could avoid creating & destructing so many Vec3 by manually storing and computing things as Vec3::value_type in static arrays
        const Vec3 du1 = Vec3(X10 - X00, Y10 - Y00, Z10 - Z00);
        const Vec3 du2 = Vec3(X01 - X11, Y01 - Y11, Z01 - Z11);
        const Vec3 dv1 = Vec3(X01 - X00, Y01 - Y00, Z01 - Z00);
        const Vec3 dv2 = Vec3(X10 - X11, Y10 - Y11, Z10 - Z11);
        const Vec3 cross1 = du1 ^ dv1;
        const Vec3 cross2 = du2 ^ dv2;
        const Vec3::value_type length1 = cross1.length();
        const Vec3::value_type length2 = cross2.length();
        const Vec3 normal1 = length1 > 0 ? cross1 / length1 : Vec3();
        const Vec3 normal2 = length2 > 0 ? cross2 / length2 : Vec3();
        if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
          const Vec3 center1 = Vec3((X00 + X10 + X01) / 3, (Y00 + Y10 + Y01) / 3, (Z00 + Z10 + Z01) / 3);
          const Vec3 center2 = Vec3((X10 + X11 + X01) / 3, (Y10 + Y11 + Y01) / 3, (Z10 + Z11 + Z01) / 3);
          facetsCenters->push_back(center1);
          facetsCenters->push_back(center2);
          facetsNormals->push_back(normal1);
          facetsNormals->push_back(normal2);
        }
        if (!normalized) {
          {
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
            ftype** W00 = W[up0][vp0];
            ftype** W01 = W[up0][vp1];
            ftype** W10 = W[up1][vp0];
            ftype** W11 = W[up1][vp1];
            itype a = 0;
            if (mNormalsAverageWeights & SurfaceNormalsAverageWeights::area) {
              const Vec3::value_type area1 = length1 / 2; // TODO surface area = triangle area = half the norm of the cross product
              const Vec3::value_type area2 = length2 / 2; // TODO surface area = triangle area = half the norm of the cross product
              W00[i][a] = area1;
              W10[j][a] = area1;
              W01[k][a] = area1;
              W10[l][a] = area2;
              W11[m][a] = area2;
              W01[n][a] = area2;
              a++;
            }
            if (mNormalsAverageWeights & SurfaceNormalsAverageWeights::angle) {
              const Vec3 duv = Vec3(X01 - X10, Y01 - Y10, Z01 - Z10);
              const Vec3 dvu = Vec3(X10 - X01, Y10 - Y01, Z10 - Z01);
              const Vec3::value_type dot11 = dv1 * du1;
              const Vec3::value_type dot12 = du1 * dvu;
              const Vec3::value_type dot13 = duv * dv1;
              const Vec3::value_type dot21 = dvu * dv2;
              const Vec3::value_type dot22 = dv2 * du2;
              const Vec3::value_type dot23 = du2 * duv;
              const Vec3::value_type lu1 = du1.length();
              const Vec3::value_type lu2 = du2.length();
              const Vec3::value_type luv = duv.length();
              const Vec3::value_type lv1 = dv1.length();
              const Vec3::value_type lv2 = dv2.length();
              const Vec3::value_type lvu = dvu.length();
              const Vec3::value_type length11 = lv1 * lu1;
              const Vec3::value_type length12 = lu1 * lvu;
              const Vec3::value_type length13 = luv * lv1;
              const Vec3::value_type length21 = lvu * lv2;
              const Vec3::value_type length22 = lv2 * lu2;
              const Vec3::value_type length23 = lu2 * luv;
              const Vec3::value_type angle11 = std::acos(length11 > 0 ? dot11 / length11 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
              const Vec3::value_type angle12 = std::acos(length12 > 0 ? dot12 / length12 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
              const Vec3::value_type angle13 = std::acos(length13 > 0 ? dot13 / length13 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
              const Vec3::value_type angle21 = std::acos(length21 > 0 ? dot21 / length21 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
              const Vec3::value_type angle22 = std::acos(length22 > 0 ? dot22 / length22 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
              const Vec3::value_type angle23 = std::acos(length23 > 0 ? dot23 / length23 : 0); // TODO corner angle = angle of the corner of the polygon at the vertex
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

      if (mNormalsAverageWeights == SurfaceNormalsAverageWeights::none) {
        // FIXME: Does not seem to work... => Of course because this is by primitive SET (entire triangle strip) not by primitive => Deprecated BIND_PER_PRIMITIVE
        // => The only simple way to go with faceted rendering is to duplicate vertices for each facet in the indices array below, so three normals per triangle
        // normals->setBinding(osg::Array::BIND_PER_PRIMITIVE_SET);
        // for (Vec3& normal : facetsNormals->asVector()) {
        //   normals->push_back(normal);
        // }
      } else {
        for (itype u = 0; u < nu; u++) {
          for (itype v = 0; v < nv; v++) {
            Vec3::value_type normal[nc] = {0};
            for (itype w = 0; w < nw; w++) {
              Vec3::value_type weight = 1;
              for (itype a = 0; a < na; a++) {
                weight *= W[u][v][w][a];
              }
              normal[x] += A[u][v][w][x] * weight;
              normal[y] += A[u][v][w][y] * weight;
              normal[z] += A[u][v][w][z] * weight;
            }
            const Vec3::value_type length = std::sqrt(normal[x] * normal[x] + normal[y] * normal[y] + normal[z] * normal[z]);
            if (length > 0) {
              const Vec3::value_type norm = 1 / length;
              normal[x] *= norm;
              normal[y] *= norm;
              normal[z] *= norm;
            } else {
              normal[x] = 0;
              normal[y] = 0;
              normal[z] = 0;
            }
            normals->push_back(Vec3(normal[x], normal[y], normal[z]));
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
      {
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
    for (itype u = 0; u < nu; u++) {
      for (itype v = 0; v < nv; v++) {
        colors->push_back(Vec4(C[u][v][x], C[u][v][y], C[u][v][z], opacity));
      }
    }
  } else {
    colors->clear();
    colors->push_back(Vec4(_color[x].exp, _color[y].exp, _color[z].exp, opacity));
    colors->setBinding(osg::Array::BIND_OVERALL);
  }

  /* Texels */
  for (itype u = 0; u < nu; u++) {
    for (itype v = 0; v < nv; v++) {
      texels->push_back(Vec2(u / (Vec2::value_type)num1, v / (Vec2::value_type)nvm1));
    }
  }

  /* Indices */
#define SURFACE_FIRST_V()                       \
  if (degenerate) {                             \
    indices->addElement(up0tnvpopv);            \
  }
#define SURFACE_LOOP_V()                        \
  for (; up0tnvpopv < up1tnvpo; up0tnvpopv++) { \
    indices->addElement(up0tnvpopv);            \
    indices->addElement(up0tnvpopv + nv);       \
  }
#define SURFACE_LAST_V()                        \
  if (degenerate) {                             \
    indices->addElement(up0tnvpopv + nvm1);     \
  }                                             \
  if (restart) {                                \
    indices->addElement(r);                     \
  }
#if 0
  const itype opnum1tnv = o + num1 * nv;
  const itype opnum2tnv = o + num2 * nv;
  for (itype up0tnvpopv = o; up0tnvpopv < opnum1tnv;) {
    const bool ug0    = up0tnvpopv > o;
    const bool ulnum2 = up0tnvpopv < opnum2tnv;
    if (ug0) {
      SURFACE_FIRST_V()
    }
    {
      const itype up1tnvpo = up0tnvpopv + nv;
      SURFACE_LOOP_V()
    }
    if (ulnum2) {
      SURFACE_LAST_V()
    }
  }
#else
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
#endif

  /* Debug */
#define SURFACE_DEBUG_N(o, e, v, n, r, g, b)    \
  osg::ref_ptr<osg::DrawElementsUInt> lines =   \
      new osg::DrawElementsUInt(                \
          osg::PrimitiveSet::LINES);            \
  constexpr Vec3::value_type s = 0.25;          \
  itype l = vertices->size();                   \
  for (itype i = o; i < o + e; i++) {           \
    const Vec3 vertex = v->at(i);               \
    const Vec3 normal = n->at(i);               \
    const Vec2 texel0 = Vec2(0, 1);             \
    const Vec2 texel1 = Vec2(1, 1);             \
    vertices->push_back(vertex);                \
    vertices->push_back(vertex + normal * s);   \
    normals ->push_back(normal);                \
    normals ->push_back(normal);                \
    texels  ->push_back(texel0);                \
    texels  ->push_back(texel1);                \
    lines->addElement(l++);                     \
    lines->addElement(l++);                     \
  }                                             \
  if (multicolored) {                           \
    colors->insert(colors->end(), e << 1,       \
        Vec4(r, g, b, opacity));                \
  }                                             \
  geometry->addPrimitiveSet(lines.get());
  if (mNormalsAnimationTypes != SurfaceNormalsAnimationTypes::none) {
    ss->setAttributeAndModes(new osg::LineWidth(5), mode); // TODO constexpr // FIXME line width will change for all lines including wireframe if enabled
    if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::vertices) {
      SURFACE_DEBUG_N(
          o, nutnv,
          vertices,
          normals,
          1, 0, 0);
    }
    if (mNormalsAnimationTypes & SurfaceNormalsAnimationTypes::facets) {
      SURFACE_DEBUG_N(
          (itype)0, (itype)facetsNormals->size(),
          facetsCenters,
          facetsNormals,
          0, 0, 1);
    }
  }

  geometry->setUseVertexBufferObjects(true);
  geometry->setVertexArray(vertices.get());
  geometry->setNormalArray(normals.get());
  geometry->setColorArray(colors.get());
  geometry->setTexCoordArray(0, texels.get());
  geometry->addPrimitiveSet(indices.get());

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
