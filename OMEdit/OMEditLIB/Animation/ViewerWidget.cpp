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

#include <QOpenGLContext> // must be included before OSG headers

#include <osgGA/MultiTouchTrackballManipulator>
#include <osgViewer/ViewerEventHandlers>

#include <QPainter>
#include <QColorDialog>
#include <QInputDialog>
#include <QKeyEvent>
#include <cassert>
#include <QtMath>
#include <QApplication>
#include <QWindow>

#include "ViewerWidget.h"
#include "Modeling/MessagesWidget.h"

/*!
 * \brief GraphicsWindowEmbeddedQt::GraphicsWindowEmbeddedQt
 * \param x
 * \param y
 * \param width
 * \param height
 * \param pViewerWidget
 */
GraphicsWindowEmbeddedQt::GraphicsWindowEmbeddedQt(int x, int y, int width, int height, ViewerWidget* pViewerWidget)
  : GraphicsWindowEmbedded(x, y, width, height)
  , mpViewerWidget(pViewerWidget)
  , mpTimer(new QTimer(mpViewerWidget))
{
  mpTimer->setTimerType(Qt::PreciseTimer);
  QObject::connect(mpTimer, &QTimer::timeout, [=]{OSG_DEBUG << "GraphicsWindowEmbeddedQt::timeout()" << std::endl;mpViewerWidget->update();});
}

/*!
 * \brief GraphicsWindowEmbeddedQt::requestRedraw
 */
void GraphicsWindowEmbeddedQt::requestRedraw()
{OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestRedraw()" << std::endl;
  mpViewerWidget->update();
}

/*!
 * \brief GraphicsWindowEmbeddedQt::requestContinuousUpdate
 * \param needed
 */
void GraphicsWindowEmbeddedQt::requestContinuousUpdate(bool needed)
{OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestContinuousUpdate(" << needed << ")" << std::endl;
  if (needed) {
    mpTimer->start();
  } else {
    mpTimer->stop();
  }
}

/*!
 * \brief GraphicsWindowEmbeddedQt::requestWarpPointer
 * \param x
 * \param y
 */
void GraphicsWindowEmbeddedQt::requestWarpPointer(float x, float y)
{OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestWarpPointer(" << x << ", " << y << ")" << std::endl;
  QPoint pos = QPoint(mpViewerWidget->convertSizeDimensionBack(x), mpViewerWidget->convertSizeDimensionBack(y));
  OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestWarpPointer(): pos = (" << pos.x() << ", " << pos.y() << ")" << std::endl;
  QPoint globalPos = mpViewerWidget->mapToGlobal(pos);
  OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestWarpPointer(): globalPos = (" << globalPos.x() << ", " << globalPos.y() << ")" << std::endl;
  OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestWarpPointer(): QCursor::pos() = (" << QCursor::pos().x() << ", " << QCursor::pos().y() << ")" << std::endl;
  QCursor::setPos(globalPos);
  OSG_DEBUG << "GraphicsWindowEmbeddedQt::requestWarpPointer(): QCursor::pos() = (" << QCursor::pos().x() << ", " << QCursor::pos().y() << ")" << std::endl;
}

/*!
 * \brief View::requestRedraw
 */
void View::requestRedraw()
{OSG_DEBUG << "View::requestRedraw()" << std::endl;
  osgViewer::View::requestRedraw();
  /*osg::ref_ptr<GraphicsWindowEmbeddedQt> gw = dynamic_cast<GraphicsWindowEmbeddedQt*>(_camera->getGraphicsContext());
  if (gw.valid()) {
    gw->requestRedraw();
  }*/
}

/*!
 * \brief View::requestContinuousUpdate
 * \param needed
 */
void View::requestContinuousUpdate(bool needed)
{OSG_DEBUG << "View::requestContinuousUpdate(" << needed << ")" << std::endl;
  osg::ref_ptr<Viewer> viewer = dynamic_cast<Viewer*>(getViewerBase());
  if (viewer.valid()) {
    needed = viewer->requestContinuousUpdate(this, needed);
    OSG_DEBUG << "View::requestContinuousUpdate(): viewer request = " << needed << std::endl;
  }
  osgViewer::View::requestContinuousUpdate(needed);
  osg::ref_ptr<GraphicsWindowEmbeddedQt> gw = dynamic_cast<GraphicsWindowEmbeddedQt*>(_camera->getGraphicsContext());
  if (gw.valid()) {
    needed = false;
    osgViewer::GraphicsWindow::Views views;
    gw->getViews(views);
    for (osg::ref_ptr<osgViewer::View> view : views) {
      if (view.valid() && view->getViewerBase() && view->getViewerBase()->getRequestContinousUpdate()) {
        needed = true;
        break;
      }
    }
    OSG_DEBUG << "View::requestContinuousUpdate(): window request = " << needed << std::endl;
    gw->requestContinuousUpdate(needed);
  }
}

/*!
 * \brief View::requestWarpPointer
 * \param x
 * \param y
 */
void View::requestWarpPointer(float x, float y)
{OSG_DEBUG << "View::requestWarpPointer(" << x << ", " << y << ")" << std::endl;
  osgViewer::View::requestWarpPointer(x, y);
  // Already forwarded to GraphicsWindowEmbeddedQt::requestWarpPointer()
}

/*!
 * \brief Viewer::requestContinuousUpdate
 * \param view
 * \param needed
 * \return
 */
bool Viewer::requestContinuousUpdate(View* view, bool needed)
{OSG_DEBUG << "Viewer::requestContinuousUpdate(" << needed << ")" << std::endl;
  if (needed) {
    requestContinuousUpdateNeeded.insert(view);
  } else {
    requestContinuousUpdateNeeded.remove(view);
  }
  return !requestContinuousUpdateNeeded.isEmpty();
}

/*!
 * \brief Viewer::setUpThreading
 */
void Viewer::setUpThreading()
{
  if (_threadingModel == SingleThreaded) {
    if (_threadsRunning) {
      stopThreading();
    }
  } else {
    if (!_threadsRunning) {
      startThreading();
    }
  }
}

/*!
 * \class ViewerWidget
 * \brief Viewer widget for OpenSceneGraph animations.
 */
/*!
 * \brief ViewerWidget::ViewerWidget
 * \param parent
 * \param flags
 */
ViewerWidget::ViewerWidget(QWidget *parent, Qt::WindowFlags flags)
  : GLWidget(parent, flags)
{
  // Set the number of samples used for multisampling
#if QT_VERSION >= QT_VERSION_CHECK(5, 4, 0)
  QSurfaceFormat format;
  format.setSamples(4);
  setFormat(format);
#else
  QGLFormat format;
  format.setSamples(4);
  setFormat(format);
#endif
  mpGraphicsWindow = new GraphicsWindowEmbeddedQt(x(), y(), width(), height(), this);
  mpViewer = new Viewer();
  mpSceneView = new View();
  mpSceneView2 = new View();
  mpFrameMutex = new OpenThreads::Mutex();
  mpAnimationWindow = qobject_cast<AbstractAnimationWindow*>(parent);
  mpSelectedVisualizer = nullptr;
  // Make sure the event queue has the correct window rectangle size and input range
  mpGraphicsWindow->getEventQueue()->syncWindowRectangleWithGraphicsContext();
  // Add a scene to the viewer
  mpViewer->addView(mpSceneView);
  mpViewer->addView(mpSceneView2);
  // Configure the camera
  osg::ref_ptr<osg::Camera> camera1 = mpSceneView->getCamera();
  camera1->setGraphicsContext(mpGraphicsWindow.get());
  camera1->setClearColor(osg::Vec4(0.95, 0.95, 0.95, 1.0));
  camera1->setViewport(new osg::Viewport(0, 0, width() / 2, height()));
  camera1->setProjectionMatrixAsPerspective(30.0, static_cast<double>(width() / 2) / static_cast<double>(height()), 1.0, 10000.0);
  osg::ref_ptr<osg::Camera> camera2 = mpSceneView2->getCamera();
  camera2->setGraphicsContext(mpGraphicsWindow.get());
  camera2->setClearColor(osg::Vec4(0.95, 0.95, 0.95, 1.0));
  camera2->setViewport(new osg::Viewport(width() / 2, 0, width() / 2, height()));
  camera2->setProjectionMatrixAsPerspective(30.0, static_cast<double>(width() / 2) / static_cast<double>(height()), 1.0, 10000.0);
  // Reverse the mouse wheel zooming
  setCursor(Qt::CrossCursor);
  //setMouseTracking(true);
  osgGA::MultiTouchTrackballManipulator *pMultiTouchTrackballManipulator = new osgGA::MultiTouchTrackballManipulator(osgGA::StandardManipulator::UPDATE_MODEL_SIZE | osgGA::StandardManipulator::COMPUTE_HOME_USING_BBOX | osgGA::StandardManipulator::PROCESS_MOUSE_WHEEL/* | osgGA::StandardManipulator::SET_CENTER_ON_WHEEL_FORWARD_MOVEMENT*/);
  pMultiTouchTrackballManipulator->setWheelZoomFactor(-pMultiTouchTrackballManipulator->getWheelZoomFactor());
  pMultiTouchTrackballManipulator->setAnimationTime(0.0);
  mpSceneView->setCameraManipulator(pMultiTouchTrackballManipulator);
  mpSceneView2->setCameraManipulator(new osgGA::TrackballManipulator);
  // Display OSG statistics by pressing the 's' key (multiple times)
  mpSceneView->addEventHandler(new osgViewer::StatsHandler());
  // Run all OSG traversals in the same thread
  mpViewer->setThreadingModel(osgViewer::CompositeViewer::SingleThreaded);
  // Improve performance for a single-threaded viewer with a single graphics context
  mpViewer->setReleaseContextAtEndOfFrameHint(false);
  // Disable the default setting of Viewer::_done on a QUIT_APPLICATION event
  mpViewer->setQuitEventSetsDone(false);
  // Disable the default setting of Viewer::_done by pressing Escape
  mpViewer->setKeyEventSetsDone(0);
  // Ensure that the scene remains visible
  setMinimumSize(100, 100);
  // This focus policy ensures that the widget will receive keyboard events.
  // The default, Qt::NoFocus, will result in keyboard events that are ignored.
  // This makes the widget grab focus either by tabbing, clicking, or wheeling.
  setFocusPolicy(Qt::WheelFocus);
}

/*!
 * \brief ViewerWidget::convertSizeDimensionBack
 * Converts the size dimension back to account for screen DPI scaling.
 * \param dimension
 * \return
 */
template<typename T>
int ViewerWidget::convertSizeDimensionBack(T dimension)
{
  qreal devicePixelRatio = window() && window()->windowHandle() ? window()->windowHandle()->devicePixelRatio() : qApp->devicePixelRatio();
  return static_cast<int>(static_cast<qreal>(.5) + (static_cast<qreal>(dimension) / devicePixelRatio));
}

/*!
 * \brief ViewerWidget::convertSizeDimension
 * Converts the size dimension to account for screen DPI scaling.
 * \param dimension
 * \return
 */
template<typename T>
int ViewerWidget::convertSizeDimension(T dimension)
{
  qreal devicePixelRatio = window() && window()->windowHandle() ? window()->windowHandle()->devicePixelRatio() : qApp->devicePixelRatio();
  return static_cast<int>(static_cast<qreal>(.5) + (static_cast<qreal>(dimension) * devicePixelRatio));
}

/*!
 * \brief ViewerWidget::convertMousePosition
 * Converts the mouse position to account for screen DPI scaling.
 * \param event
 * \param reverseY
 * \return
 */
QPoint ViewerWidget::convertMousePosition(QMouseEvent *event, bool reverseY)
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
  QPointF position = event->position();
  using T = qreal;
#else
  QPoint position = event->pos();
  using T = int;
#endif
  if (reverseY) {
    position.setY(static_cast<T>(height()) - position.y());
  }
  QPoint point = QPoint(convertSizeDimension(position.x()), convertSizeDimension(position.y()));
  if (reverseY && point.y() > 0) {
    point.ry()--;
  }
  return point;
}

/*!
 * \brief ViewerWidget::convertMouseButton
 * Converts the mouse button to the number expected by OSG.
 * \param event
 * \return
 */
unsigned int ViewerWidget::convertMouseButton(QMouseEvent *event)
{
  // 1 = left mouse button
  // 2 = middle mouse button
  // 3 = right mouse button
  unsigned int button;
  switch (event->button()) {
    case Qt::LeftButton:
      button = 1;
      break;
    case Qt::MiddleButton:
      button = 2;
      break;
    case Qt::RightButton:
      button = 3;
      break;
    default:
      button = 0;
      break;
  }
  return button;
}

/*!
 * \brief ViewerWidget::convertKeyCode
 * Converts the key code to the symbol expected by OSG.
 * \param event
 * \return
 */
QPair<int, int> ViewerWidget::convertKeyCode(QKeyEvent *event)
{
  QString keyString = event->text();
  QByteArray keyByteArray = keyString.toLocal8Bit(); // FIXME see https://wiki.qt.io/Technical_FAQ#How_can_I_convert_a_QString_to_char.2A_and_vice_versa.3F (probably worked because of https://forum.qt.io/post/648180)
  const char* keyData = keyByteArray.constData(); // FIXME see https://doc.qt.io/qt-6/qbytearray.html#data
  int keySymbol = osgGA::GUIEventAdapter::KeySymbol(*keyData);
  int virtualKeySymbol = event->key() == Qt::Key_Control ? osgGA::GUIEventAdapter::KEY_Control_L : 0;
  if (!keySymbol) { // Since OSG 3.6.5 this additional step could be omitted (see OSG commit f4fe1e5)
    keySymbol = virtualKeySymbol;
  }
  return QPair<int, int>(keySymbol, virtualKeySymbol);
}

/*!
 * \brief ViewerWidget::initializeGL
 * Reimplementation of QOpenGLWidget::initializeGL().
 */
void ViewerWidget::initializeGL()
{
  if (!mpViewer->isRealized()) {
    mpViewer->realize();
  }
}

/*!
 * \brief ViewerWidget::paintGL
 * Reimplementation of QOpenGLWidget::paintGL().
 * \note Synchronized frame rendering can lead to a deadlock
 *       in situations where a new paint event is fired while
 *       a frame is currently being rendered, and specifically
 *       MessagesWidget::addPendingMessage() shall be used instead of
 *       MessagesWidget::addGUIMessage() when #mpFrameMutex is locked.
 * \sa ViewerWidget::paintEvent()
 * \sa ViewerWidget::frame()
 */
void ViewerWidget::paintGL()
{
  mpFrameMutex->lock();
  frame();
  mpFrameMutex->unlock();
  MessagesWidget::instance()->showPendingMessages();
}

/*!
 * \brief ViewerWidget::frame
 * Renders the animation frame.
 * \sa ViewerWidget::paintGL()
 */
void ViewerWidget::frame()
{OSG_DEBUG << "ViewerWidget::frame()" << std::endl;
  mpViewer->frame();
}

/*!
 * \brief ViewerWidget::resizeGL
 * Reimplementation of QOpenGLWidget::resizeGL().
 * Resizes the graphics window.
 * \param width
 * \param height
 */
void ViewerWidget::resizeGL(int width, int height)
{
  int x = convertSizeDimension(this->x());
  int y = convertSizeDimension(this->y());
#if QT_VERSION >= QT_VERSION_CHECK(5, 4, 0)
  width = convertSizeDimension(width);
  height = convertSizeDimension(height);
#endif
  mpGraphicsWindow->resized(x, y, width, height);
  mpGraphicsWindow->getEventQueue()->windowResize(x, y, width, height);
}

/*!
 * \brief ViewerWidget::keyPressEvent
 * Passes the QWidget::keyPressEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::keyPressEvent(QKeyEvent *event)
{
  QPair<int, int> key = convertKeyCode(event);
  mpGraphicsWindow->getEventQueue()->keyPress(key.first, key.second);
}

/*!
 * \brief ViewerWidget::keyReleaseEvent
 * Passes the QWidget::keyReleaseEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::keyReleaseEvent(QKeyEvent *event)
{
  QPair<int, int> key = convertKeyCode(event);
  mpGraphicsWindow->getEventQueue()->keyRelease(key.first, key.second);
}

/*!
 * \brief ViewerWidget::mouseMoveEvent
 * Passes the QWidget::mouseMoveEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::mouseMoveEvent(QMouseEvent *event)
{OSG_DEBUG << "ViewerWidget::mouseMoveEvent()" << std::endl;
  QPoint position = convertMousePosition(event);
  mpGraphicsWindow->getEventQueue()->mouseMotion(static_cast<float>(position.x()), static_cast<float>(position.y()));
}

/*!
 * \brief ViewerWidget::mousePressEvent
 * Passes the QWidget::mousePressEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::mousePressEvent(QMouseEvent *event)
{OSG_DEBUG << "ViewerWidget::mousePressEvent()" << std::endl;
  if (event->button() == Qt::RightButton && event->modifiers() == Qt::ShiftModifier) {
    // Qt counts pixels from the upper-left corner and OSG from the lower-left corner (GL viewport), thus pass reverseY = true
    QPoint position = convertMousePosition(event, true);
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    QPoint mousePosition = event->position().toPoint();
#else
    QPoint mousePosition = event->pos();
#endif
    pickVisualizer(static_cast<float>(position.x()), static_cast<float>(position.y()));
    showVisualizerPickContextMenu(mousePosition);
    return;
  }
  unsigned int button = convertMouseButton(event);
  QPoint position = convertMousePosition(event);
  mpGraphicsWindow->getEventQueue()->mouseButtonPress(static_cast<float>(position.x()), static_cast<float>(position.y()), button);
}

/*!
 * \brief ViewerWidget::pickVisualizer
 * Picks the name of the selected visualizer in the osg view
 * \param x - mouse position pixel in x direction in osg system
 * \param y - mouse position pixel in y direction in osg system
 */
void ViewerWidget::pickVisualizer(float x, float y)
{
  mpSelectedVisualizer = nullptr;
  //std::cout<<"pickVisualizer "<<x<<" and "<<y<<std::endl;
  osgUtil::LineSegmentIntersector::Intersections intersections;
  if (mpSceneView->computeIntersections(mpSceneView->getCamera(), osgUtil::Intersector::WINDOW, x, y, intersections)) {
    //take the first intersection with a facette only
    osgUtil::LineSegmentIntersector::Intersections::const_iterator hitr = intersections.cbegin();
    constexpr osg::NodePath::size_type lvl = 2;
    if (hitr->nodePath.size() > lvl && !hitr->nodePath.at(lvl)->getName().empty()) {
      mpSelectedVisualizer = mpAnimationWindow->getVisualization()->getBaseData()->getVisualizerObjectByID(hitr->nodePath.at(lvl)->getName());
      //std::cout<<"Object identified by name "<<mpSelectedVisualizer->_id<<std::endl;
    } else if (hitr->drawable.valid()) {
      mpSelectedVisualizer = mpAnimationWindow->getVisualization()->getBaseData()->getVisualizerObjectByID(hitr->drawable->className());
      //std::cout<<"Object identified by its drawable "<<mpSelectedVisualizer->_id<<std::endl;
    }
  }
}

/*!
 * \brief ViewerWidget::showVisualizerPickContextMenu
 * \param pos
 */
void ViewerWidget::showVisualizerPickContextMenu(const QPoint &pos)
{
  QString name = mpSelectedVisualizer ? QString::fromStdString(mpSelectedVisualizer->_id) : QString();
  //std::cout<<"SHOW CONTEXT "<<name.toStdString()<<" compare "<<QString::compare(name,QString(""))<< std::endl;

  // The context widget
  QMenu contextMenu(tr("Context menu"), this);
  QMenu visualizerMenu(name, &contextMenu);

  visualizerMenu.setIcon(QIcon(":/Resources/icons/animation.svg"));
  QAction action0(QIcon(":/Resources/icons/reset.svg"), tr("Reset Visual Properties"), &contextMenu);
  QAction action1(QIcon(":/Resources/icons/transparency.svg"), tr("Change Transparency"), &visualizerMenu);
  QAction action2(QIcon(":/Resources/icons/invisible.svg"), tr("Make Visualizer Invisible"), &visualizerMenu);
  QAction action3(QIcon(":/Resources/icons/changeColor.svg"), tr("Change Color"), &visualizerMenu);
  QAction action4(QIcon(":/Resources/icons/specularity.svg"), tr("Change Specularity"), &visualizerMenu);
  QAction action5(QIcon(":/Resources/icons/checkered.svg"), tr("Apply Checker Texture"), &visualizerMenu);
  QAction action6(QIcon(":/Resources/icons/texture.svg"), tr("Apply Custom Texture"), &visualizerMenu);
  QAction action7(QIcon(":/Resources/icons/undo.svg"), tr("Remove Texture"), &visualizerMenu);

  connect(&action0, SIGNAL(triggered()), this, SLOT(resetVisualPropertiesForAllVisualizers()));
  connect(&action1, SIGNAL(triggered()), this, SLOT(changeVisualizerTransparency()));
  connect(&action2, SIGNAL(triggered()), this, SLOT(makeVisualizerInvisible()));
  connect(&action3, SIGNAL(triggered()), this, SLOT(changeVisualizerColor()));
  connect(&action4, SIGNAL(triggered()), this, SLOT(changeVisualizerSpec()));
  connect(&action5, SIGNAL(triggered()), this, SLOT(applyCheckerTexture()));
  connect(&action6, SIGNAL(triggered()), this, SLOT(applyCustomTexture()));
  connect(&action7, SIGNAL(triggered()), this, SLOT(removeTexture()));

  // If a visualizer is picked, one can change its properties
  if (mpSelectedVisualizer) {
    action2.setText(tr((std::string("Make ") + mpSelectedVisualizer->getVisualizerType() + " Invisible").c_str()));
    contextMenu.addMenu(&visualizerMenu);
  }

  contextMenu.addAction(&action0);
  visualizerMenu.addAction(&action1);
  visualizerMenu.addAction(&action2);
  visualizerMenu.addSeparator();
  visualizerMenu.addAction(&action3);
  visualizerMenu.addAction(&action4);
  visualizerMenu.addSeparator();
  visualizerMenu.addAction(&action5);
  visualizerMenu.addAction(&action6);
  visualizerMenu.addAction(&action7);

  contextMenu.exec(mapToGlobal(pos));

  mpSelectedVisualizer = nullptr;
}

/*!
 * \brief ViewerWidget::changeVisualizerTransparency
 * opens a number dialog and selects a new transparency for the visualizer
 */
void ViewerWidget::changeVisualizerTransparency()
{
  if (mpSelectedVisualizer) {
    bool ok;
    const int min = 0, max = 100, step = 1; // Unit: [%]
    const int currentTransparency = mpSelectedVisualizer->getVisualProperties()->getTransparency().get() * (max - min) + min;
    const int transparency = QInputDialog::getInt(this, Helper::chooseTransparency, Helper::percentageLabel,
                                                  currentTransparency, min, max, step, &ok);
    if (ok) { // Picked transparency is not OK if the user cancels the dialog
      mpSelectedVisualizer->getVisualProperties()->getTransparency().set((float) (transparency - min) / (max - min));
      mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
    }
  }
}

/*!
 * \brief ViewerWidget::makeVisualizerInvisible
 * suppresses the visualization of this visualizer
 */
void ViewerWidget::makeVisualizerInvisible()
{
  if (mpSelectedVisualizer) {
    mpSelectedVisualizer->getVisualProperties()->getTransparency().set(1.0);
    mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
  }
}

/*!
 * \brief ViewerWidget::changeVisualizerColor
 * opens a color dialog and selects a new color for the visualizer
 */
void ViewerWidget::changeVisualizerColor()
{
  if (mpSelectedVisualizer) {
    const QColor currentColor = mpSelectedVisualizer->getVisualProperties()->getColor().get();
    const QColor color = QColorDialog::getColor(currentColor, this, Helper::chooseColor);
    if (color.isValid()) { // Picked color is invalid if the user cancels the dialog
      mpSelectedVisualizer->getVisualProperties()->getColor().set(color);
      mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
    }
  }
}

/*!
 * \brief ViewerWidget::changeVisualizerSpec
 * opens a number dialog and selects a new specular coefficient for the visualizer
 */
void ViewerWidget::changeVisualizerSpec()
{
  if (mpSelectedVisualizer) {
    bool ok;
    const int min = 0, max = 100, step = 1; // Unit: [%]
    const int currentSpecular = mpSelectedVisualizer->getVisualProperties()->getSpecular().get() * (max - min) + min;
    const int specular = QInputDialog::getInt(this, Helper::chooseSpecularity, Helper::percentageLabel,
                                              currentSpecular, min, max, step, &ok);
    if (ok) { // Picked specular coefficient is not OK if the user cancels the dialog
      mpSelectedVisualizer->getVisualProperties()->getSpecular().set((float) (specular - min) / (max - min));
      mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
    }
  }
}

/*!
 * \brief ViewerWidget::applyCheckerTexture
 * adds a checkered texture to the visualizer
 */
void ViewerWidget::applyCheckerTexture()
{
  if (mpSelectedVisualizer) {
    if (mpSelectedVisualizer->isShape()) {
      ShapeObject* shape = mpSelectedVisualizer->asShape();
      if (isSimpleCADType(shape->_type)) {
        QString msg = tr("Texture feature is not applicable for %1 files.").arg(shape->_type.c_str());
        MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, msg,
                                                              Helper::scriptingKind, Helper::notificationLevel));
        return;
      }
    }
    mpSelectedVisualizer->getVisualProperties()->getTextureImagePath().set(":/Resources/bitmaps/check.png");
    mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
  }
}

/*!
 * \brief ViewerWidget::applyCustomTexture
 * adds a user-defined texture to the visualizer
 */
void ViewerWidget::applyCustomTexture()
{
  if (mpSelectedVisualizer) {
    if (mpSelectedVisualizer->isShape()) {
      ShapeObject* shape = mpSelectedVisualizer->asShape();
      if (isSimpleCADType(shape->_type)) {
        QString msg = tr("Texture feature is not applicable for %1 files.").arg(shape->_type.c_str());
        MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, msg,
                                                              Helper::scriptingKind, Helper::notificationLevel));
        return;
      }
    }
    const QString* currentFileName = nullptr; // File picker starts at the last open directory if any, otherwise the user's home directory
    const QString fileName = StringHandler::getOpenFileName(this, QString("%1 – %2").arg(Helper::applicationName).arg(Helper::chooseFile),
                                                            (QString*) currentFileName, Helper::bitmapFileTypes, nullptr);
    if (!fileName.isEmpty()) { // Picked file name is empty if the user cancels the dialog
      mpSelectedVisualizer->getVisualProperties()->getTextureImagePath().set(fileName.toStdString());
      mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
    }
  }
}

/*!
 * \brief ViewerWidget::removeTexture
 * removes the texture of the visualizer
 */
void ViewerWidget::removeTexture()
{
  if (mpSelectedVisualizer) {
    if (mpSelectedVisualizer->isShape()) {
      ShapeObject* shape = mpSelectedVisualizer->asShape();
      if (isSimpleCADType(shape->_type)) {
        QString msg = tr("Texture feature is not applicable for %1 files.").arg(shape->_type.c_str());
        MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, msg,
                                                              Helper::scriptingKind, Helper::notificationLevel));
        return;
      }
    }
    mpSelectedVisualizer->getVisualProperties()->getTextureImagePath().set("");
    mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(mpSelectedVisualizer);
  }
}

/*!
 * \brief ViewerWidget::resetVisualPropertiesForAllVisualizers
 * sets all visual properties back to default
 */
void ViewerWidget::resetVisualPropertiesForAllVisualizers()
{
  for (AbstractVisualizerObject& visualizer : mpAnimationWindow->getVisualization()->getBaseData()->getVisualizerObjects()) {
    visualizer.getVisualProperties()->resetVisualProperties();
    mpAnimationWindow->getVisualization()->getBaseData()->modifyVisualizer(visualizer);
  }
}

/*!
 * \brief ViewerWidget::mouseReleaseEvent
 * Passes the QWidget::mouseReleaseEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::mouseReleaseEvent(QMouseEvent *event)
{OSG_DEBUG << "ViewerWidget::mouseReleaseEvent()" << std::endl;
  unsigned int button = convertMouseButton(event);
  QPoint position = convertMousePosition(event);
  mpGraphicsWindow->getEventQueue()->mouseButtonRelease(static_cast<float>(position.x()), static_cast<float>(position.y()), button);
}

/*!
 * \brief ViewerWidget::mouseDoubleClickEvent
 * Passes the QWidget::mouseDoubleClickEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::mouseDoubleClickEvent(QMouseEvent *event)
{OSG_DEBUG << "ViewerWidget::mouseDoubleClickEvent()" << std::endl;
  unsigned int button = convertMouseButton(event);
  QPoint position = convertMousePosition(event);
  mpGraphicsWindow->getEventQueue()->mouseDoubleButtonPress(static_cast<float>(position.x()), static_cast<float>(position.y()), button);
}

/*!
 * \brief ViewerWidget::wheelEvent
 * Passes the QWidget::wheelEvent() to Graphics Window.
 * \param event
 */
void ViewerWidget::wheelEvent(QWheelEvent *event)
{OSG_DEBUG << "ViewerWidget::wheelEvent()" << std::endl;
#if QT_VERSION >= QT_VERSION_CHECK(5, 6, 0)
  static QPoint angleDelta = QPoint();
  angleDelta += event->angleDelta();
  QPoint numDegrees = angleDelta / 8;
  QPoint numSteps = numDegrees / 15; // See QWheelEvent documentation
  if (numSteps.isNull()) {
    event->ignore();
    return;
  }
  angleDelta = QPoint();
  bool up = qAbs(numSteps.x()) > qAbs(numSteps.y()) ? numSteps.x() > 0 : numSteps.y() > 0;
#else
  bool up = event->delta() > 0;
#endif
  mpGraphicsWindow->getEventQueue()->mouseScroll(up ? osgGA::GUIEventAdapter::SCROLL_UP : osgGA::GUIEventAdapter::SCROLL_DOWN);
}

/*!
 * \brief ViewerWidget::event
 * Repaints the graphics window on each user interaction.
 * \param event
 * \return
 */
bool ViewerWidget::event(QEvent *event)
{OSG_DEBUG << "ViewerWidget::event(" << event->type() << ")" << std::endl;
  bool handled = GLWidget::event(event);
  // This ensures that the OSG widget is always going to be repainted after the
  // user performed some interaction. Doing this in the event handler ensures
  // that we don't forget about some event and prevents duplicate code.
  switch (event->type()) {
    case QEvent::KeyPress:
    case QEvent::KeyRelease:
    case QEvent::MouseButtonDblClick:
    case QEvent::MouseButtonPress:
    case QEvent::MouseButtonRelease:
    case QEvent::MouseMove:
    case QEvent::Wheel:
      update();
      break;
    default:
      break;
  }OSG_DEBUG << "ViewerWidget::event(return)" << std::endl;
  return handled;
}
