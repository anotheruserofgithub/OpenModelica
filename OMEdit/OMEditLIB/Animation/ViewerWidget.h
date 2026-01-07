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
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES RECIPIENT'S ACCEPTANCE
 * OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3, ACCORDING TO RECIPIENTS CHOICE.
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
/*
 * @author Adeel Asghar <adeel.asghar@liu.se>
 */

#ifndef VIEWERWIDGET_H
#define VIEWERWIDGET_H

#include <QOpenGLContext> // must be included before OSG headers

#include <osg/ref_ptr>
#include <osgViewer/GraphicsWindow>
#include <osgViewer/CompositeViewer>

#include <OpenThreads/Mutex>

#include <iostream>

#include <QMenu>

#include "AbstractAnimationWindow.h"
#include "AnimationUtil.h"
#include "Util/Helper.h"

/*!
 * \note We need to create two files with same class name since Qt meta object compiler doesn't handle ifdef.
 * OpenGLWidget.h uses QOpenGLWidget and GLWidget.h uses QGLWidget.
 */
#include <QtGlobal>
#if QT_VERSION >= QT_VERSION_CHECK(5, 4, 0)
#include "OpenGLWidget.h"
#else
#include "GLWidget.h"
#endif

class ViewerWidget; // Forward declaration

class GraphicsWindowEmbeddedQt : public osgViewer::GraphicsWindowEmbedded
{
public:
  GraphicsWindowEmbeddedQt(int x, int y, int width, int height, ViewerWidget* pViewerWidget);
  virtual bool realizeImplementation() override {return realized = true;}
  virtual bool isRealizedImplementation() const override {return realized;}
  virtual void requestRedraw() override;
  virtual void requestContinuousUpdate(bool needed = true) override;
  virtual void requestWarpPointer(float x, float y) override;
private:
  bool realized = false;
  ViewerWidget* mpViewerWidget;
  QTimer* mpTimer;
};

/*!
 * This subclassing allows us to forward update requests to the GraphicsWindowEmbeddedQt
 * which then calls the QWidget::update() method of the ViewerWidget only when appropriate.
 */
class View : public osgViewer::View
{
public:
  virtual void requestRedraw() override;
  virtual void requestContinuousUpdate(bool needed = true) override;
  virtual void requestWarpPointer(float x, float y) override;
};

/*!
 * This subclassing allows us to merge update requests coming from multiple View instances.
 * This subclassing also allows us to remove the annoying automatic
 * setting of the CPU affinity to core 0 by osgViewer::ViewerBase,
 * osgViewer::CompositeViewer's base class.
 * Note this: Since OSG 3.6.0 we could call osgViewer::ViewerBase::setUseConfigureAffinity(false)
 * before osgViewer::ViewerBase::realize(), which calls osgViewer::ViewerBase::setUpThreading(),
 * instead of overriding osgViewer::ViewerBase::setUpThreading() with an outdated implementation
 * (see OSG commit 91538d9).
 */
class Viewer : public osgViewer::CompositeViewer
{
public:
  bool requestContinuousUpdate(View* view, bool needed = true);
  virtual void setUpThreading() override;
private:
  QSet<View*> requestContinuousUpdateNeeded;
};

class ViewerWidget : public GLWidget
{
  Q_OBJECT
public:
  ViewerWidget(QWidget *pParent = nullptr, Qt::WindowFlags flags = Qt::WindowFlags());
  View* getSceneView() {return mpSceneView.get();}
  View* getSceneView2() {return mpSceneView2.get();}
  View* getSceneView3() {return mpSceneView3.get();}
  View* getSceneView4() {return mpSceneView4.get();}
  OpenThreads::Mutex* getFrameMutex() {return mpFrameMutex;}
  void frame();
  template<typename T>
  int convertSizeDimensionBack(T dimension);
protected:
  template<typename T>
  int convertSizeDimension(T dimension);
  QPoint convertMousePosition(QMouseEvent *event, bool reverseY = false);
  unsigned int convertMouseButton(QMouseEvent *event);
  QPair<int, int> convertKeyCode(QKeyEvent *event);
  virtual void initializeGL() override;
  virtual void paintGL() override;
  virtual void resizeGL(int width, int height) override;
  virtual void keyPressEvent(QKeyEvent *event) override;
  virtual void keyReleaseEvent(QKeyEvent *event) override;
  virtual void mouseMoveEvent(QMouseEvent *event) override;
  virtual void mousePressEvent(QMouseEvent *event) override;
  virtual void mouseReleaseEvent(QMouseEvent *event) override;
  virtual void mouseDoubleClickEvent(QMouseEvent *event) override;
  virtual void wheelEvent(QWheelEvent *event) override;
  virtual bool event(QEvent *event) override;
  void pickVisualizer(float x, float y);
  void showVisualizerPickContextMenu(const QPoint &pos);
private:
  osg::ref_ptr<GraphicsWindowEmbeddedQt> mpGraphicsWindow;
  osg::ref_ptr<Viewer> mpViewer;
  osg::ref_ptr<Viewer> mpViewer2;
  osg::ref_ptr<View> mpSceneView;
  osg::ref_ptr<View> mpSceneView2;
  osg::ref_ptr<View> mpSceneView3;
  osg::ref_ptr<View> mpSceneView4;
  OpenThreads::Mutex* mpFrameMutex;
  AbstractAnimationWindow* mpAnimationWindow;
  AbstractVisualizerObject* mpSelectedVisualizer;
protected slots:
  void changeVisualizerTransparency();
  void makeVisualizerInvisible();
  void changeVisualizerColor();
  void changeVisualizerSpec();
  void applyCheckerTexture();
  void applyCustomTexture();
  void removeTexture();
  void resetVisualPropertiesForAllVisualizers();
};

#endif
