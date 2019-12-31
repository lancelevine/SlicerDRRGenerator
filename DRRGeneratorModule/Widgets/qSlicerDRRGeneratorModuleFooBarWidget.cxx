/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerDRRGeneratorModuleFooBarWidget.h"
#include "ui_qSlicerDRRGeneratorModuleFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_DRRGeneratorModule
class qSlicerDRRGeneratorModuleFooBarWidgetPrivate
  : public Ui_qSlicerDRRGeneratorModuleFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerDRRGeneratorModuleFooBarWidget);
protected:
  qSlicerDRRGeneratorModuleFooBarWidget* const q_ptr;

public:
  qSlicerDRRGeneratorModuleFooBarWidgetPrivate(
    qSlicerDRRGeneratorModuleFooBarWidget& object);
  virtual void setupUi(qSlicerDRRGeneratorModuleFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerDRRGeneratorModuleFooBarWidgetPrivate
::qSlicerDRRGeneratorModuleFooBarWidgetPrivate(
  qSlicerDRRGeneratorModuleFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerDRRGeneratorModuleFooBarWidgetPrivate
::setupUi(qSlicerDRRGeneratorModuleFooBarWidget* widget)
{
  this->Ui_qSlicerDRRGeneratorModuleFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerDRRGeneratorModuleFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleFooBarWidget
::qSlicerDRRGeneratorModuleFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerDRRGeneratorModuleFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerDRRGeneratorModuleFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleFooBarWidget
::~qSlicerDRRGeneratorModuleFooBarWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerDRRGeneratorModuleFooBarWidget
::setMRMLScene(vtkMRMLScene* scene)
{
	//this->Superclass::setMRMLScene(scene);
}