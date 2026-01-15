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
 * @author Volker Waurich <volker.waurich@tu-dresden.de>
 */

#ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <QOpenGLContext> // must be included before OSG headers

#include <osgGA/OrbitManipulator>

#include "AbstractAnimationWindow.h"
#include "Modeling/MessagesWidget.h"
#include "Options/OptionsDialog.h"
#include "ViewerWidget.h"
#include "Visualization.h"
#include "VisualizationMAT.h"
#include "VisualizationCSV.h"
#include "VisualizationFMU.h"

/*!
 * \class AbstractAnimationWindow
 * \brief Abstract animation class defines a QMainWindow for animation.
 */
/*!
 * \brief AbstractAnimationWindow::AbstractAnimationWindow
 * \param pParent
 */
AbstractAnimationWindow::AbstractAnimationWindow(QWidget *pParent)
  : QMainWindow(pParent),
    mPathName(""),
    mFileName(""),
    mpVisualization(nullptr),
    mpViewerWidget(nullptr),
    mpAnimationToolBar(new QToolBar(QString("Animation Toolbar"),this)),
    mpAnimationCameraDockWidget(new QDockWidget(QString("Trackball Manipulator"),this)),
    mpAnimationParameterDockWidget(new QDockWidget(QString("Parameter Settings"),this)),
    mpAnimationChooseFileAction(nullptr),
    mpAnimationInitializeAction(nullptr),
    mpAnimationPlayAction(nullptr),
    mpAnimationPauseAction(nullptr),
    mpAnimationRepeatAction(nullptr),
    mpAnimationSlider(nullptr),
    mpAnimationTimeLabel(nullptr),
    mpTimeTextBox(nullptr),
    mpAnimationSpeedLabel(nullptr),
    mpSpeedComboBox(nullptr),
    mpPerspectiveDropDownBox(nullptr),
    mpRotateCameraLeftAction(nullptr),
    mpRotateCameraRightAction(nullptr),
    mpCenterAtOriginAction(nullptr),
    mpFitInViewAction(nullptr),
    mpCameraControlAction(nullptr),
    mpInteractiveControlAction(nullptr),
    mCameraInitialized(false),
    mSliderRange(1000)
{
  // to distinguish this widget as a subwindow among the plotwindows
  setObjectName(QString("animationWindow"));
  //the viewer widget
  mpViewerWidget = new ViewerWidget(this);
  // we need to set the minimum height so that visualization window is still shown when we cascade windows.
  mpViewerWidget->setMinimumHeight(100);
  // toolbar icon size
  int toolbarIconSize = OptionsDialog::instance()->getGeneralSettingsPage()->getToolbarIconSizeSpinBox()->value();
  mpAnimationToolBar->setIconSize(QSize(toolbarIconSize, toolbarIconSize));
  addToolBar(Qt::TopToolBarArea, mpAnimationToolBar);
  addDockWidget(Qt::RightDockWidgetArea, mpAnimationCameraDockWidget);
  addDockWidget(Qt::RightDockWidgetArea, mpAnimationParameterDockWidget);

  // Viewer layout
  QGridLayout *pGridLayout = new QGridLayout;
  pGridLayout->setContentsMargins(0, 0, 0, 0);
  pGridLayout->addWidget(mpViewerWidget);
  // add the viewer to the frame for boxed rectangle around it.
  QFrame *pCentralWidgetFrame = new QFrame;
  pCentralWidgetFrame->setFrameStyle(QFrame::StyledPanel);
  pCentralWidgetFrame->setLayout(pGridLayout);
  setCentralWidget(pCentralWidgetFrame);
}

/*!
 * \brief AbstractAnimationWindow::openAnimationFile
 * \param fileName
 * \param stashCamera
 */
void AbstractAnimationWindow::openAnimationFile(QString fileName, bool stashCamera)
{
  std::string file = fileName.toStdString();
  if (file.compare("")) {
    std::size_t pos = file.find_last_of("/\\");
    mPathName = file.substr(0, pos + 1);
    mFileName = file.substr(pos + 1, file.length());
    //std::cout<<"file "<<mFileName<<"   path "<<mPathName<<std::endl;
    if (loadVisualization()) {
      // start the widgets
      mpAnimationInitializeAction->setEnabled(true);
      mpAnimationPlayAction->setEnabled(true);
      mpAnimationPauseAction->setEnabled(true);
      mpAnimationRepeatAction->setEnabled(true);
      mpAnimationRepeatAction->setChecked(false);
      mpAnimationSlider->setEnabled(true);
      bool state = mpAnimationSlider->blockSignals(true);
      mpAnimationSlider->setValue(mpVisualization->getTimeManager()->getTimeFraction());
      mpAnimationSlider->blockSignals(state);
      mpTimeTextBox->setEnabled(true);
      mpTimeTextBox->setText(QString::number(mpVisualization->getTimeManager()->getVisTime()));
      mpSpeedComboBox->setEnabled(true);
      mpSpeedComboBox->lineEdit()->setText(QString::number(mpVisualization->getTimeManager()->getSpeedUp()));
      /* Only use isometric view as default for csv file type.
       * Otherwise use side view as default which suits better for Modelica models.
       */
      if (mpVisualization->getVisType() == VisType::CSV) {
        mpPerspectiveDropDownBox->setCurrentIndex(0);
        cameraPositionIsometric();
      } else {
        mpPerspectiveDropDownBox->setCurrentIndex(1);
        cameraPositionSide();
      }

      if(stashCamera && !mCameraInitialized) {         // mCameraInitialized is used to make sure the view is never stashed
        mCameraInitialized = true;      // before the camera is initialized the first time
        stashView();
      }
    }
  }
}

void AbstractAnimationWindow::createActions()
{
  // actions and widgets for the toolbar
  int toolbarIconSize = OptionsDialog::instance()->getGeneralSettingsPage()->getToolbarIconSizeSpinBox()->value();
  // choose file action
  mpAnimationChooseFileAction = new QAction(QIcon(":/Resources/icons/open.svg"), Helper::animationChooseFile, this);
  mpAnimationChooseFileAction->setStatusTip(Helper::animationChooseFileTip);
  connect(mpAnimationChooseFileAction, SIGNAL(triggered()),this, SLOT(chooseAnimationFileSlotFunction()));
  // initialize action
  mpAnimationInitializeAction = new QAction(QIcon(":/Resources/icons/initialize.svg"), Helper::animationInitialize, this);
  mpAnimationInitializeAction->setStatusTip(Helper::animationInitializeTip);
  mpAnimationInitializeAction->setEnabled(false);
  connect(mpAnimationInitializeAction, SIGNAL(triggered()),this, SLOT(initSlotFunction()));
  // animation play action
  mpAnimationPlayAction = new QAction(QIcon(":/Resources/icons/play_animation.svg"), Helper::animationPlay, this);
  mpAnimationPlayAction->setStatusTip(Helper::animationPlayTip);
  mpAnimationPlayAction->setEnabled(false);
  connect(mpAnimationPlayAction, SIGNAL(triggered()),this, SLOT(playSlotFunction()));
  // animation pause action
  mpAnimationPauseAction = new QAction(QIcon(":/Resources/icons/pause.svg"), Helper::animationPause, this);
  mpAnimationPauseAction->setStatusTip(Helper::animationPauseTip);
  mpAnimationPauseAction->setEnabled(false);
  connect(mpAnimationPauseAction, SIGNAL(triggered()),this, SLOT(pauseSlotFunction()));
  // animation repeat action
  mpAnimationRepeatAction = new QAction(QIcon(":/Resources/icons/refresh.svg"), Helper::animationRepeat, this);
  mpAnimationRepeatAction->setStatusTip(Helper::animationRepeatTip);
  mpAnimationRepeatAction->setEnabled(false);
  mpAnimationRepeatAction->setCheckable(true);
  connect(mpAnimationRepeatAction, SIGNAL(triggered(bool)),this, SLOT(repeatSlotFunciton(bool)));
  // animation slide
  mpAnimationSlider = new QSlider(Qt::Horizontal);
  mpAnimationSlider->setMinimum(0);
  mpAnimationSlider->setMaximum(mSliderRange);
  mpAnimationSlider->setSliderPosition(0);
  mpAnimationSlider->setEnabled(false);
  connect(mpAnimationSlider, SIGNAL(valueChanged(int)),this, SLOT(sliderSetTimeSlotFunction(int)));
  // animation time
  QDoubleValidator *pDoubleValidator = new QDoubleValidator(this);
  pDoubleValidator->setBottom(0);
  mpAnimationTimeLabel = new Label;
  mpAnimationTimeLabel->setText(tr("Time [s]:"));
  mpTimeTextBox = new QLineEdit("0.0");
  mpTimeTextBox->setMaximumSize(QSize(toolbarIconSize*2, toolbarIconSize));
  mpTimeTextBox->setEnabled(false);
  mpTimeTextBox->setValidator(pDoubleValidator);
  connect(mpTimeTextBox, SIGNAL(returnPressed()),this, SLOT(jumpToTimeSlotFunction()));
  // animation speed
  mpAnimationSpeedLabel = new Label;
  mpAnimationSpeedLabel->setText(tr("Speed:"));
  mpSpeedComboBox = new QComboBox;
  mpSpeedComboBox->setEditable(true);
  mpSpeedComboBox->addItems(Helper::speedOptions.split(","));
  mpSpeedComboBox->setCurrentIndex(3);
  mpSpeedComboBox->setMaximumSize(QSize(toolbarIconSize*2, toolbarIconSize));
  mpSpeedComboBox->setEnabled(false);
  mpSpeedComboBox->setValidator(pDoubleValidator);
  mpSpeedComboBox->setCompleter(0);
  connect(mpSpeedComboBox, SIGNAL(currentIndexChanged(int)),this, SLOT(setSpeedSlotFunction()));
  connect(mpSpeedComboBox->lineEdit(), SIGNAL(textChanged(QString)),this, SLOT(setSpeedSlotFunction()));
  // perspective drop down
  mpPerspectiveDropDownBox = new QComboBox;
  mpPerspectiveDropDownBox->addItem(QIcon(":/Resources/icons/perspective0.svg"), QString("Isometric"));
  mpPerspectiveDropDownBox->addItem(QIcon(":/Resources/icons/perspective1.svg"), QString("Side"));
  mpPerspectiveDropDownBox->addItem(QIcon(":/Resources/icons/perspective2.svg"), QString("Front"));
  mpPerspectiveDropDownBox->addItem(QIcon(":/Resources/icons/perspective3.svg"), QString("Top"));
  mpPerspectiveDropDownBox->view()->setMinimumWidth(mpPerspectiveDropDownBox->minimumSizeHint().width()); // See https://stackoverflow.com/a/54005207
  connect(mpPerspectiveDropDownBox, SIGNAL(activated(int)), this, SLOT(setPerspective(int)));
  // rotate camera left action
  mpRotateCameraLeftAction = new QAction(QIcon(":/Resources/icons/rotateCameraLeft.svg"), tr("Rotate Left"), this);
  mpRotateCameraLeftAction->setStatusTip(tr("Rotate the scene left"));
  connect(mpRotateCameraLeftAction, SIGNAL(triggered()), this, SLOT(rotateCameraLeft()));
  // rotate camera right action
  mpRotateCameraRightAction = new QAction(QIcon(":/Resources/icons/rotateCameraRight.svg"), tr("Rotate Right"), this);
  mpRotateCameraRightAction->setStatusTip(tr("Rotate the scene right"));
  connect(mpRotateCameraRightAction, SIGNAL(triggered()), this, SLOT(rotateCameraRight()));
  // center at origin action
  mpCenterAtOriginAction = new QAction(QIcon(":/Resources/icons/center-at-origin.svg"), tr("Center at Origin"), this);
  mpCenterAtOriginAction->setStatusTip(tr("Center the scene at the origin"));
  connect(mpCenterAtOriginAction, SIGNAL(triggered()), this, SLOT(centerAtOrigin()));
  // fit in view action
  mpFitInViewAction = new QAction(QIcon(":/Resources/icons/fit-in-view.svg"), tr("Fit in View"), this);
  mpFitInViewAction->setStatusTip(tr("Fit the scene in the view"));
  connect(mpFitInViewAction, SIGNAL(triggered()), this, SLOT(fitInView()));
  // camera control action
  mpCameraControlAction = mpAnimationCameraDockWidget->toggleViewAction();
  mpCameraControlAction->setIcon(QIcon(":/Resources/icons/joystick.svg"));
  mpCameraControlAction->setText(tr("Camera Control"));
  mpCameraControlAction->setStatusTip(tr("Open the camera control panel"));
  //mpAnimationCameraDockWidget->hide();
  // interactive control action
  mpInteractiveControlAction = mpAnimationParameterDockWidget->toggleViewAction();
  mpInteractiveControlAction->setIcon(QIcon(":/Resources/icons/control-panel.svg"));
  mpInteractiveControlAction->setText(tr("Interactive Control"));
  mpInteractiveControlAction->setStatusTip(tr("Open the interactive control panel"));
  mpInteractiveControlAction->setEnabled(false);
  mpAnimationParameterDockWidget->hide();
}

/*!
 * \brief AbstractAnimationWindow::initCameraControlPanel
 */
void AbstractAnimationWindow::initCameraControlPanel()
{
  struct LLabel : public QLabel
  {
    LLabel(const QString &text) : QLabel(text)
    {
      setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    }
  };

  struct RLabel : public QLabel
  {
    RLabel(const QString &text) : QLabel(text)
    {
      setAlignment(Qt::AlignRight | Qt::AlignVCenter);
      setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    }
  };

  struct ADoubleSpinBox : public QDoubleSpinBox
  {
    ADoubleSpinBox() : QDoubleSpinBox()
    {
      setDecimals(3);
      setRange(-DBL_MAX, DBL_MAX);
      setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
      setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
      setMinimumWidth(3 * OptionsDialog::instance()->getGeneralSettingsPage()->getToolbarIconSizeSpinBox()->value());
    }
  };

  struct AScrollArea : public QScrollArea
  {
    QSize sizeHint() const override
    {
      // See https://github.com/qt/qtbase/blob/v6.10.1/src/widgets/widgets/qscrollarea.cpp#L364-L382
      int f = 2 * frameWidth();
      QSize sz(f, f);
      sz += widgetResizable() ? widget()->minimumSizeHint() : widget()->size();
      if (verticalScrollBarPolicy() == Qt::ScrollBarAlwaysOn)
        sz.rwidth() += verticalScrollBar()->sizeHint().width();
      if (horizontalScrollBarPolicy() == Qt::ScrollBarAlwaysOn)
        sz.rheight() += horizontalScrollBar()->sizeHint().height();
      return sz;
    }
  };

  AScrollArea *pScrollArea = new AScrollArea();
  pScrollArea->setWidgetResizable(true);

  QWidget *pWidget = new QWidget();
  pWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

  QVBoxLayout *pVBoxLayout = new QVBoxLayout();
  pVBoxLayout->setSizeConstraint(QLayout::SetMinimumSize);

  // Focal point
  {
    QGroupBox *pFocalPointGroupBox = new QGroupBox(QString("Focal point"));

    QGridLayout *pFocalPointGridLayout = new QGridLayout();

    // Center
    {
      LLabel *pCenterLabel = new LLabel(QString("Center:"));
      pCenterLabel->setToolTip(QString("Vector to the center of rotation in world coordinates."));
      pFocalPointGridLayout->addWidget(pCenterLabel, 0, 0, 1, -1);

      // X axis
      {
        RLabel *pCenterXLabel = new RLabel(QString("x ="));
        pFocalPointGridLayout->addWidget(pCenterXLabel, 1, 0);

        ADoubleSpinBox *pCenterXSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pCenterXSpinBox, 1, 1);
      }

      // Y axis
      {
        RLabel *pCenterYLabel = new RLabel(QString("y ="));
        pFocalPointGridLayout->addWidget(pCenterYLabel, 1, 3);

        ADoubleSpinBox *pCenterYSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pCenterYSpinBox, 1, 4);
      }

      // Z axis
      {
        RLabel *pCenterZLabel = new RLabel(QString("z ="));
        pFocalPointGridLayout->addWidget(pCenterZLabel, 1, 6);

        ADoubleSpinBox *pCenterZSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pCenterZSpinBox, 1, 7);
      }

      // Home center
      {
        QToolButton *pCenterHomeButton = new QToolButton();
        pCenterHomeButton->setAutoRaise(true);
        pCenterHomeButton->setIcon(QIcon(":/Resources/icons/center-at-home-position.svg"));
        pCenterHomeButton->setToolTip(QString("Center the scene at the home position."));
        connect(pCenterHomeButton, SIGNAL(clicked()), this, SLOT(centerAtHomePosition()));
        pFocalPointGridLayout->addWidget(pCenterHomeButton, 1, 9, 1, -1);
      }
    }

    // Rotation
    {
      LLabel *pRotationLabel = new LLabel(QString("Rotation:"));
      pRotationLabel->setToolTip(QString("Quaternion for the orientation with respect to the focal center."));
      pFocalPointGridLayout->addWidget(pRotationLabel, 2, 0, 1, -1);

      // X basis
      {
        RLabel *pRotationXLabel = new RLabel(QString("x ="));
        pFocalPointGridLayout->addWidget(pRotationXLabel, 3, 0);

        ADoubleSpinBox *pRotationXSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pRotationXSpinBox, 3, 1);
      }

      // Y basis
      {
        RLabel *pRotationYLabel = new RLabel(QString("y ="));
        pFocalPointGridLayout->addWidget(pRotationYLabel, 3, 3);

        ADoubleSpinBox *pRotationYSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pRotationYSpinBox, 3, 4);
      }

      // Z basis
      {
        RLabel *pRotationZLabel = new RLabel(QString("z ="));
        pFocalPointGridLayout->addWidget(pRotationZLabel, 3, 6);

        ADoubleSpinBox *pRotationZSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pRotationZSpinBox, 3, 7);
      }

      // W basis
      {
        RLabel *pRotationWLabel = new RLabel(QString("w ="));
        pFocalPointGridLayout->addWidget(pRotationWLabel, 3, 9);

        ADoubleSpinBox *pRotationWSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pRotationWSpinBox, 3, 10);
      }
    }

    // Distance
    {
      LLabel *pDistanceLabel = new LLabel(QString("Distance:"));
      pDistanceLabel->setToolTip(QString("Scalar distance from the focal center to the camera along the line of sight."));
      pFocalPointGridLayout->addWidget(pDistanceLabel, 4, 0, 1, -1);

      // Value
      {
        ADoubleSpinBox *pDistanceSpinBox = new ADoubleSpinBox();
        pFocalPointGridLayout->addWidget(pDistanceSpinBox, 5, 0, 1, -1);
      }
    }

    pFocalPointGroupBox->setLayout(pFocalPointGridLayout);
    pVBoxLayout->addWidget(pFocalPointGroupBox);
  }

  pWidget->setLayout(pVBoxLayout);
  pScrollArea->setWidget(pWidget);
  mpAnimationCameraDockWidget->setWidget(pScrollArea);
}

/*!
 * \brief AbstractAnimationWindow::initInteractiveControlPanel
 */
void AbstractAnimationWindow::initInteractiveControlPanel()
{
  if (getVisualization()) {
    VisualizationFMU* FMUvis = dynamic_cast<VisualizationFMU*>(mpVisualization);
    QWidget *widget = new QWidget(this);
    mSpinBoxVector.clear();
    mStateLabels.clear();

    if (FMUvis->getFMU()->getFMUData()->_stateNames.size()==FMUvis->getFMU()->getFMUData()->_nStates){
    //widgets for the states
    QVBoxLayout *layout = new QVBoxLayout();
    layout->setSpacing(2);

    layout->addWidget(new QLabel(QString("<b>modify state variables</b>"),this));
    for (unsigned int stateIdx = 0; stateIdx < FMUvis->getFMU()->getFMUData()->_nStates; stateIdx++)
    {
      DoubleSpinBoxIndexed* spinBox = new DoubleSpinBoxIndexed(this, stateIdx);
      spinBox->setValue(FMUvis->getFMU()->getFMUData()->_states[stateIdx]);
      spinBox->setMaximum(DBL_MAX);
      spinBox->setMinimum(-DBL_MAX);
      spinBox->setSingleStep(0.1);
      mSpinBoxVector.push_back(spinBox);

      QLabel* stateLabel = new QLabel(QString::number(FMUvis->getFMU()->getFMUData()->_states[stateIdx]), this);
      stateLabel->setMargin(0);
      mStateLabels.push_back(stateLabel);

      QWidget* stateValWidget = new QWidget();
      QHBoxLayout *stateInfoLayout = new QHBoxLayout();
      stateInfoLayout->setContentsMargins(0, 0, 0, 0);
      stateInfoLayout->addWidget(stateLabel);
      stateInfoLayout->addWidget(spinBox);
      stateValWidget->setLayout(stateInfoLayout);

      layout->addSpacing(12);
      layout->addWidget(new QLabel(QString::fromStdString(FMUvis->getFMU()->getFMUData()->_stateNames.at(stateIdx)),this));
      layout->addWidget(stateValWidget);
      connect(mSpinBoxVector.at(stateIdx), SIGNAL(valueChangedFrom(double, int)), this, SLOT(setStateSolveSystem(double, int)));
    }

    //widgets for the inputs
    //layout->addWidget(new QLabel(QString("<b>modify inputs</b>"),this));
    /*
    for (int inputIdx = 0; inputIdx < 0; inputIdx++) {
      std::string name = "some input";
      layout->addWidget(new QLabel(QString::fromStdString("inputName"),this));
      QWidget* valWidget = new QWidget();
      QHBoxLayout *valueLayout = new QHBoxLayout();
      valueLayout->addWidget(new QLabel(QString::fromStdString("min"),this));
      valueLayout->addWidget(new QLabel(QString::fromStdString("value"),this));
      valueLayout->addWidget(new QLabel(QString::fromStdString("max"),this));
      valWidget->setLayout(valueLayout);
      layout->addWidget(valWidget);
      layout->addWidget(new QSlider(Qt::Horizontal, this));
    }
    QPushButton* updateButton = new QPushButton(QIcon(":/Resources/icons/update.svg"), "&update system",this);
    layout->addWidget(updateButton);
    connect(updateButton, SIGNAL(released()), this, SLOT(setStatesSolveSystem()));
    */

    layout->setSizeConstraint(QLayout::SetMinimumSize);
    widget->setLayout(layout);
    QScrollArea *scrollArea = new QScrollArea(this);
    scrollArea->setWidget(widget);
    mpAnimationParameterDockWidget->setWidget(scrollArea);
    } else {
      MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, tr("Information about states could not be determined."),
                                                            Helper::scriptingKind, Helper::errorLevel));
      mpAnimationParameterDockWidget->hide();
    }
  } else {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, tr("Interactive Control needs an FMU ME 2.0"),
                                                          Helper::scriptingKind, Helper::errorLevel));
    mpAnimationParameterDockWidget->hide();
  }
}

/*!
 * \brief AbstractAnimationWindow::updateControlPanelValues
 */
void AbstractAnimationWindow::updateControlPanelValues()
{
  if (getVisualization()) {
    VisualizationFMU* FMUvis = dynamic_cast<VisualizationFMU*>(mpVisualization);
    for (int stateIdx = 0; stateIdx < mSpinBoxVector.size(); stateIdx++) {
      mStateLabels.at(stateIdx)->setText(QString::number(FMUvis->getFMU()->getFMUData()->_states[stateIdx]));
    }
  }
}

/*!
 * \brief AbstractAnimationWindow::setStateSolveSystem
 * \param val The new value
 * \param idx The state that belongs to the value
 */
void AbstractAnimationWindow::setStateSolveSystem(double val, int idx)
{
  if (idx>=0) {
    VisualizationFMU* FMUvis = dynamic_cast<VisualizationFMU*>(mpVisualization);
    if (FMUvis) {
      FMUvis->getFMU()->getFMUData()->_states[idx] = val;
      FMUvis->updateSystem();
      mpViewerWidget->update();
    }
  }
}


/*!
 * \brief AbstractAnimationWindow::clearView
 */
void AbstractAnimationWindow::clearView()
{
  if (mpViewerWidget) {
    mpViewerWidget->getSceneView()->setSceneData(0);
    popView();
    mpViewerWidget->update();
  }
}

void AbstractAnimationWindow::stashView()
{
  if(!mCameraInitialized) return;
  mStashedViewMatrix = mpViewerWidget->getSceneView()->getCameraManipulator()->getMatrix();
}

void AbstractAnimationWindow::popView()
{
  if(!mCameraInitialized) return;
  mpViewerWidget->getSceneView()->getCameraManipulator()->setByMatrix(mStashedViewMatrix);
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::loadVisualization
 * loads the data and the xml scene description
 * \return
 */
bool AbstractAnimationWindow::loadVisualization()
{
  VisType visType = VisType::NONE;
  // Get visualization type.
  if (isFMU(mFileName)) {
    visType = VisType::FMU;
  } else if (isMAT(mFileName)) {
    visType = VisType::MAT;
  } else if (isCSV(mFileName)) {
    visType = VisType::CSV;
  } else {
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, tr("Unknown visualization type."),
                                                          Helper::scriptingKind, Helper::errorLevel));
    return false;
  }
  //load the XML File, build osgTree, get initial values for the visualizers
  bool xmlExists = checkForXMLFile(mFileName, mPathName);
  if (!xmlExists) {
    QString msg = tr("Could not find the visual XML file %1.").arg(QString(assembleXMLFileName(mFileName, mPathName).c_str()));
    MessagesWidget::instance()->addGUIMessage(MessageItem(MessageItem::Modelica, msg, Helper::scriptingKind,
                                                          Helper::errorLevel));
    return false;
  } else {
    //clear previous visualization
    if (mpVisualization) {
      clearView();
      delete mpVisualization;
      mpVisualization = nullptr;
    }
    //init visualization
    if (visType == VisType::MAT) {
      mpVisualization = new VisualizationMAT(mFileName, mPathName);
    } else if (visType == VisType::CSV) {
      mpVisualization = new VisualizationCSV(mFileName, mPathName);
    } else if (visType == VisType::FMU) {
      mpVisualization = new VisualizationFMU(mFileName, mPathName);
    }
    Q_CHECK_PTR(mpVisualization);
    connect(mpVisualization->getTimeManager()->getUpdateSceneTimer(), SIGNAL(timeout()), SLOT(updateScene()));
    mpVisualization->initData();
    mpVisualization->setUpScene();
    mpVisualization->initVisualization();
    //add scene for the chosen visualization
    mpViewerWidget->getSceneView()->setSceneData(mpVisualization->getOMVisScene()->getScene().getRootNode());
    //choose suitable scales for the vector visualizers so that they fit well in the scene
    mpViewerWidget->makeCurrent();
    mpVisualization->getBaseData()->chooseVectorScales(mpViewerWidget->getSceneView(), mpViewerWidget->getFrameMutex(), std::bind(&ViewerWidget::frame, mpViewerWidget));
    mpViewerWidget->doneCurrent();
  }
  //add window title
  setWindowTitle(QString::fromStdString(mFileName));
  //fill trackball manipulator widget
  initCameraControlPanel();
  //open settings dialog for FMU simulation
  if (mpVisualization->getVisType() == VisType::FMU) {
    openFMUSettingsDialog(dynamic_cast<VisualizationFMU*>(mpVisualization));
    mpInteractiveControlAction->setEnabled(true);
    initInteractiveControlPanel();
  }
  else {
    mpInteractiveControlAction->setEnabled(false);
    mpAnimationParameterDockWidget->hide();
  }

  return true;
}

/*!
 * \brief AbstractAnimationWindow::resetCamera
 * resets the camera position
 */
void AbstractAnimationWindow::resetCamera()
{
  mpViewerWidget->getSceneView()->home();
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::cameraPosition
 * sets the camera position with the specified rotation from the focal center, which is reset at the home center,
 * maintaining the current distance to the focal center
 * \param rotation
 */
void AbstractAnimationWindow::cameraPosition(const osg::Quat& rotation)
{
  osg::ref_ptr<osgGA::OrbitManipulator> manipulator = static_cast<osgGA::OrbitManipulator*>(mpViewerWidget->getSceneView()->getCameraManipulator());
  osg::Vec3d eye, center, up;
  manipulator->getHomePosition(eye, center, up);
  manipulator->setCenter(center);
  manipulator->setRotation(rotation);
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::cameraPositionIsometric
 * sets the camera position to isometric view
 */
void AbstractAnimationWindow::cameraPositionIsometric()
{
  cameraPosition(osg::Quat(-std::asin(std::tan(M_PI/6.0)), osg::Vec3d(1, 0, 0), M_PI/4.0, osg::Vec3d(0, 1, 0), 0.0, osg::Vec3d(0, 0, 1)));
}

/*!
 * \brief AbstractAnimationWindow::cameraPositionSide
 * sets the camera position to side view
 */
void AbstractAnimationWindow::cameraPositionSide()
{
  cameraPosition(osg::Quat(0.0, osg::Vec3d(0, 0, 1)));
}

/*!
 * \brief AbstractAnimationWindow::cameraPositionFront
 * sets the camera position to front view
 */
void AbstractAnimationWindow::cameraPositionFront()
{
  cameraPosition(osg::Quat(+M_PI/2.0, osg::Vec3d(0, 1, 0)));
}

/*!
 * \brief AbstractAnimationWindow::cameraPositionTop
 * sets the camera position to top view
 */
void AbstractAnimationWindow::cameraPositionTop()
{
  cameraPosition(osg::Quat(-M_PI/2.0, osg::Vec3d(1, 0, 0)));
}

/*!
 * \brief AbstractAnimationWindow::openFMUSettingsDialog
 * Opens a dialog to set the settings for the FMU visualization
 * \param pVisualizationFMU
 */
void AbstractAnimationWindow::openFMUSettingsDialog(VisualizationFMU* pVisualizationFMU)
{
  FMUSettingsDialog *pFMUSettingsDialog = new FMUSettingsDialog(this, pVisualizationFMU);
  pFMUSettingsDialog->exec();
}

/*!
 * \brief AbstractAnimationWindow::updateSceneFunction
 * updates the visualization objects
 */
void AbstractAnimationWindow::updateScene()
{
  if (mpVisualization) {
    // update scene
    mpVisualization->sceneUpdate();
    mpViewerWidget->update();
    if (!mpVisualization->getTimeManager()->isPaused()) {
      // set time label
      mpTimeTextBox->setText(QString::number(mpVisualization->getTimeManager()->getVisTime()));
      // set time slider
      bool state = mpAnimationSlider->blockSignals(true);
      mpAnimationSlider->setValue(mpVisualization->getTimeManager()->getTimeFraction());
      mpAnimationSlider->blockSignals(state);
    }
    if (mpVisualization->getVisType() == VisType::FMU) {
      // set state labels
      updateControlPanelValues();
    }
  }
}

/*!
 * \brief AbstractAnimationWindow::animationFileSlotFunction
 * opens a file dialog to choose an animation
 */
void AbstractAnimationWindow::chooseAnimationFileSlotFunction()
{
  QString fileName = StringHandler::getOpenFileName(this, QString("%1 – %2").arg(Helper::applicationName).arg(Helper::chooseFile),
                                                    nullptr, Helper::visualizationFileTypes, nullptr);
  if (fileName.isEmpty()) {
    return;
  }
  openAnimationFile(fileName);
}

/*!
 * \brief AbstractAnimationWindow::initSlotFunction
 * slot function for the init button
 */
void AbstractAnimationWindow::initSlotFunction()
{
  mpVisualization->initVisualization();
  mpViewerWidget->update();
  bool state = mpAnimationSlider->blockSignals(true);
  mpAnimationSlider->setValue(mpVisualization->getTimeManager()->getTimeFraction());
  mpAnimationSlider->blockSignals(state);
  mpTimeTextBox->setText(QString::number(mpVisualization->getTimeManager()->getVisTime()));
  if (mpVisualization->getVisType() == VisType::FMU) {
    updateControlPanelValues();
  }
}

/*!
 * \brief AbstractAnimationWindow::playSlotFunction
 * slot function for the play button
 */
void AbstractAnimationWindow::playSlotFunction()
{
  mpVisualization->getTimeManager()->setPause(false);
}

/*!
 * \brief AbstractAnimationWindow::pauseSlotFunction
 * slot function for the pause button
 */
void AbstractAnimationWindow::pauseSlotFunction()
{
  mpVisualization->getTimeManager()->setPause(true);
}

/*!
 * \brief AbstractAnimationWindow::repeatSlotFunciton
 * slot function for the repeat button
 * \param checked
 */
void AbstractAnimationWindow::repeatSlotFunciton(bool checked)
{
  mpVisualization->getTimeManager()->setRepeat(checked);
}

/*!
 * \brief AbstractAnimationWindow::updateSceneTime
 * Updates the scene to the provided point of time
 * \param time The new point of time
 */
void AbstractAnimationWindow::updateSceneTime(double time)
{
  double start = mpVisualization->getTimeManager()->getStartTime();
  double end = mpVisualization->getTimeManager()->getEndTime();
  if (time < start) {
    time = start;
  } else if (time > end) {
    time = end;
  }
  mpVisualization->getTimeManager()->setVisTime(time);
  mpTimeTextBox->setText(QString::number(mpVisualization->getTimeManager()->getVisTime()));
  bool state = mpAnimationSlider->blockSignals(true);
  mpAnimationSlider->setValue(mpVisualization->getTimeManager()->getTimeFraction());
  mpAnimationSlider->blockSignals(state);
  mpVisualization->updateScene(mpVisualization->getTimeManager()->getVisTime());
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::sliderSetTimeSlotFunction
 * slot function for the time slider to jump to the adjusted point of time
 * \param value
 */
void AbstractAnimationWindow::sliderSetTimeSlotFunction(int value)
{
  if (value >= 0) {
    double start = mpVisualization->getTimeManager()->getStartTime();
    double end = mpVisualization->getTimeManager()->getEndTime();
    double time = (end - start) * (value / (double)mSliderRange) + start;
    updateSceneTime(time);
  } else {
    bool state = mpAnimationSlider->blockSignals(true);
    mpAnimationSlider->setValue(mpVisualization->getTimeManager()->getTimeFraction());
    mpAnimationSlider->blockSignals(state);
  }
}

/*!
 * \brief AbstractAnimationWindow::jumpToTimeSlotFunction
 * slot function to jump to the user input point of time
 */
void AbstractAnimationWindow::jumpToTimeSlotFunction()
{
  bool isDouble = false;
  double time = mpTimeTextBox->text().toDouble(&isDouble);
  if (isDouble && time >= 0.0) {
    updateSceneTime(time);
  } else {
    mpTimeTextBox->setText(QString::number(mpVisualization->getTimeManager()->getVisTime()));
  }
}

/*!
 * \brief AbstractAnimationWindow::setSpeedUpSlotFunction
 * slot function to set the user input speed up
 */
void AbstractAnimationWindow::setSpeedSlotFunction()
{
  bool isDouble = false;
  double speed = mpSpeedComboBox->lineEdit()->text().toDouble(&isDouble);
  if (isDouble && speed > 0.0) {
    mpVisualization->getTimeManager()->setSpeedUp(speed);
    mpViewerWidget->update();
  } else {
    mpSpeedComboBox->lineEdit()->setText(QString::number(mpVisualization->getTimeManager()->getSpeedUp()));
  }
}

/*!
 * \brief AbstractAnimationWindow::setPerspective
 * gets the identifier for the chosen perspective and calls the functions
 * \param value
 */
void AbstractAnimationWindow::setPerspective(int value)
{
  switch(value) {
    case 0:
      cameraPositionIsometric();
      break;
    case 1:
      cameraPositionSide();
      break;
    case 2:
      cameraPositionFront();
      break;
    case 3:
      cameraPositionTop();
      break;
  }
}

/*!
 * \brief AbstractAnimationWindow::rotateCamera
 * rotates the camera by the specified angle about the line of sight
 * \param angle
 */
void AbstractAnimationWindow::rotateCamera(double angle)
{
  osg::ref_ptr<osgGA::OrbitManipulator> manipulator = static_cast<osgGA::OrbitManipulator*>(mpViewerWidget->getSceneView()->getCameraManipulator());
  const osg::Quat rotation = osg::Quat(angle, osg::Vec3d(0, 0, 1)) * manipulator->getRotation();
  manipulator->setRotation(rotation);
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::rotateCameraLeft
 * rotates the scene 90 degrees left about the line of sight
 */
void AbstractAnimationWindow::rotateCameraLeft()
{
  rotateCamera(-M_PI/2.0);
}

/*!
 * \brief AbstractAnimationWindow::rotateCameraRight
 * rotates the scene 90 degrees right about the line of sight
 */
void AbstractAnimationWindow::rotateCameraRight()
{
  rotateCamera(+M_PI/2.0);
}

/*!
 * \brief AbstractAnimationWindow::centerAtHomePosition
 * centers the scene at the home position
 */
void AbstractAnimationWindow::centerAtHomePosition()
{
  osg::ref_ptr<osgGA::OrbitManipulator> manipulator = static_cast<osgGA::OrbitManipulator*>(mpViewerWidget->getSceneView()->getCameraManipulator());
  osg::Vec3d eye, center, up;
  manipulator->getHomePosition(eye, center, up);
  manipulator->setCenter(center);
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::centerAtOrigin
 * centers the scene at the origin
 */
void AbstractAnimationWindow::centerAtOrigin()
{
  osg::ref_ptr<osgGA::OrbitManipulator> manipulator = static_cast<osgGA::OrbitManipulator*>(mpViewerWidget->getSceneView()->getCameraManipulator());
  manipulator->setCenter(osg::Vec3d());
  mpViewerWidget->update();
}

/*!
 * \brief AbstractAnimationWindow::fitInView
 * fits the scene in the view
 */
void AbstractAnimationWindow::fitInView()
{
  osg::ref_ptr<osgGA::OrbitManipulator> manipulator = static_cast<osgGA::OrbitManipulator*>(mpViewerWidget->getSceneView()->getCameraManipulator());
  const osg::Quat rotation = manipulator->getRotation();
  mpViewerWidget->getSceneView()->home();
  manipulator->setRotation(rotation);
  mpViewerWidget->update();
}

/*!
 * \brief DoubleSpinBoxIndexed::DoubleSpinBoxIndexed
 * \param pParent
 * \param idx
 */
DoubleSpinBoxIndexed::DoubleSpinBoxIndexed(QWidget *pParent, int idx)
  : QDoubleSpinBox(pParent)
{
  mStateIdx = idx;
  connect( this, SIGNAL(valueChanged(double)),
           this, SLOT(slotForwardValueChanged(double)));
}

/*!
 * \brief DoubleSpinBoxIndexed::slotForwardValueChanged
 * \param val
 */
void DoubleSpinBoxIndexed::slotForwardValueChanged(double val)
{
 emit valueChangedFrom(val, mStateIdx);
}

