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

#ifndef __qSlicerDRRGeneratorModuleModuleWidget_h
#define __qSlicerDRRGeneratorModuleModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerDRRGeneratorModuleModuleExport.h"

#include <vtkMRMLVolumeNode.h>

class qSlicerDRRGeneratorModuleModuleWidgetPrivate;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_DRRGENERATORMODULE_EXPORT qSlicerDRRGeneratorModuleModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerDRRGeneratorModuleModuleWidget(QWidget *parent=0);
  virtual ~qSlicerDRRGeneratorModuleModuleWidget();

public slots:

void setRotationX(double);
void setRotationY(double);
void setRotationZ(double);
void setTranslationX(double);
void setTranslationY(double);
void setTranslationZ(double);
void setCenterX(double);
void setCenterY(double);
void setCenterZ(double);
void setFocalPoint(double);
void setThreshold(double);
void setDRRSizeX(int);
void setDRRSizeY(int);
void setDirection(int);

protected slots :

	void setInputVolume(vtkMRMLNode * volumeNode);
	void setOutputVolume(vtkMRMLNode * volumeNode);
	/*void projectVolume();
	void rotateAndProjectVolume(int x, int y, int z);
	void rotateAndProjectVolume0();
	void rotateAndProjectVolume5();
	void rotateAndProjectVolume15();
	void rotateAndProjectVolume25();
	void projectVolumeDRR3backup(int x, int y, int z);
	void projectVolumeDRR3(int x, int y, int z);
	void projectVolumeDRR0();
	void projectVolumeDRR5();
	void projectVolumeDRR15();
	void projectVolumeDRR25();
	void projectVolumeFilter();
	void cropVolume();
	void rotateVolume();*/
	//void registerVolume3d();
	/*void registerVolumeOld();
	double registerVolume2dSimilarity();
	double registerVolume2dEuler();
	void automaticRegistration();*/
	void GenerateDRR();

protected:
  QScopedPointer<qSlicerDRRGeneratorModuleModuleWidgetPrivate> d_ptr;

  virtual void setup();

  bool debug;

private:
  Q_DECLARE_PRIVATE(qSlicerDRRGeneratorModuleModuleWidget);
  Q_DISABLE_COPY(qSlicerDRRGeneratorModuleModuleWidget);
};

#endif
