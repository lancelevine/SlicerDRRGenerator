/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// DRRGeneratorModule Logic includes
#include <vtkSlicerDRRGeneratorModuleLogic.h>

// DRRGeneratorModule includes
#include "qSlicerDRRGeneratorModuleModule.h"
#include "qSlicerDRRGeneratorModuleModuleWidget.h"

//-----------------------------------------------------------------------------
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
#include <QtPlugin>
Q_EXPORT_PLUGIN2(qSlicerDRRGeneratorModuleModule, qSlicerDRRGeneratorModuleModule);
#endif

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerDRRGeneratorModuleModulePrivate
{
public:
  qSlicerDRRGeneratorModuleModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerDRRGeneratorModuleModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleModulePrivate::qSlicerDRRGeneratorModuleModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerDRRGeneratorModuleModule methods

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleModule::qSlicerDRRGeneratorModuleModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerDRRGeneratorModuleModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleModule::~qSlicerDRRGeneratorModuleModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerDRRGeneratorModuleModule::helpText() const
{
  return "This module can generate a Digitally Reconstructed Radiograph.";
}

//-----------------------------------------------------------------------------
QString qSlicerDRRGeneratorModuleModule::acknowledgementText() const
{
  return "This work was created by a medical student. I apologize for any errors.";
}

//-----------------------------------------------------------------------------
QStringList qSlicerDRRGeneratorModuleModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Lance Levine");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerDRRGeneratorModuleModule::icon() const
{
  return QIcon(":/Icons/DRRGeneratorModule.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerDRRGeneratorModuleModule::categories() const
{
  return QStringList() << "Filtering";
}

//-----------------------------------------------------------------------------
QStringList qSlicerDRRGeneratorModuleModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerDRRGeneratorModuleModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerDRRGeneratorModuleModule
::createWidgetRepresentation()
{
  return new qSlicerDRRGeneratorModuleModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerDRRGeneratorModuleModule::createLogic()
{
  return vtkSlicerDRRGeneratorModuleLogic::New();
}
