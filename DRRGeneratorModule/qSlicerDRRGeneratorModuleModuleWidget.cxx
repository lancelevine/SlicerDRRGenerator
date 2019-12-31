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

// Qt includes
#include <QDebug>

// SlicerQt includes
#include "qSlicerDRRGeneratorModuleModuleWidget.h"
#include "ui_qSlicerDRRGeneratorModuleModuleWidget.h"

#include "qSlicerCoreIOManager.h"

#include "vtkSlicerApplicationLogic.h"
#include "vtkMRMLSegmentationNode.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"

#include "itkMeanProjectionImageFilter.h"
#include "itkProjectionImageFilter.h"


#include "itkMRMLIDImageIO.h"
#include "itkMRMLIDImageIOFactory.h"
#include "itkVTKImageToImageFilter.h"
#include "itkImageToVTKImageFilter.h"

#include "itkVTKImageIO.h"

#include "vtkMatrix4x4.h"
#include "itkFlipImageFilter.h"

#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkJPEGReader.h"

#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkRayCastInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "vtkImageReslice.h"

#include "vtkImageWriter.h"
#include "vtkNrrdReader.h"
#include "itkExtractImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkNrrdImageIO.h"
#include "itkVTKImageIO.h"

#include "vtkImageCast.h"
#include "vtkPNGWriter.h"
#include "vtkImageShiftScale.h"

#include "itkCenteredAffineTransform.h"

#include "itkEuler2DTransform.h"
#include "itkExhaustiveOptimizerv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkImageRegistrationMethodv4.h"

#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkCommand.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkSimilarity2DTransform.h"

#include <vtkTransform.h>
#include "itkSubtractImageFilter.h"

#include <chrono>

class vtkMRMLNode;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerDRRGeneratorModuleModuleWidgetPrivate: public Ui_qSlicerDRRGeneratorModuleModuleWidget
{
public:
  qSlicerDRRGeneratorModuleModuleWidgetPrivate();

  vtkWeakPointer<vtkMRMLVolumeNode> InputVolumeNode;
  vtkWeakPointer<vtkMRMLVolumeNode> OutputVolumeNode;
  double rx;
  double ry;
  double rz;
  double tx;
  double ty;
  double tz;
  double cx;
  double cy;
  double cz;
  double sid;
  double drrthreshold;
  int drrsizex;
  int drrsizey;
  int direction;
};

//-----------------------------------------------------------------------------
// qSlicerDRRGeneratorModuleModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleModuleWidgetPrivate::qSlicerDRRGeneratorModuleModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerDRRGeneratorModuleModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleModuleWidget::qSlicerDRRGeneratorModuleModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerDRRGeneratorModuleModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerDRRGeneratorModuleModuleWidget::~qSlicerDRRGeneratorModuleModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerDRRGeneratorModuleModuleWidget::setup()
{
  Q_D(qSlicerDRRGeneratorModuleModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
 
  connect(d->InputVolumeComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(setInputVolume(vtkMRMLNode*)));

  connect(d->ProjectionComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	  this, SLOT(setOutputVolume(vtkMRMLNode*)));

  connect(d->rotationXSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setRotationX(double)));

  connect(d->rotationYSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setRotationY(double)));

  connect(d->rotationZSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setRotationZ(double)));

  connect(d->translationXSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setTranslationX(double)));

  connect(d->translationYSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setTranslationY(double)));

  connect(d->translationZSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setTranslationZ(double)));

  connect(d->centerXSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setCenterX(double)));

  connect(d->centerYSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setCenterY(double)));

  connect(d->centerZSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setCenterZ(double)));

  connect(d->focalPointSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setFocalPoint(double)));

  connect(d->thresholdSpinBox, SIGNAL(valueChanged(double)),
	  this, SLOT(setThreshold(double)));

  connect(d->sizeXSpinBox, SIGNAL(valueChanged(int)),
	  this, SLOT(setDRRSizeX(int)));

  connect(d->sizeYSpinBox, SIGNAL(valueChanged(int)),
	  this, SLOT(setDRRSizeY(int)));

  connect(d->directionSpinBox, SIGNAL(valueChanged(int)),
	  this, SLOT(setDirection(int)));

  connect(d->GenerateDRRButton, SIGNAL(clicked()),
	  this, SLOT(GenerateDRR()));

  d->rx = 90.00;
  d->ty = 200.00;
  d->sid = 400.00;
  d->drrsizex = 512;
  d->drrsizey = 494;
  d->direction = 2;
}

//-----------------------------------------------------------------------------
/*void qSlicerDRRGeneratorModuleModuleWidget::updateWidgetFromMRML()
{

}*/
//-----------------------------------------------------------------------------
void qSlicerDRRGeneratorModuleModuleWidget::setInputVolume(vtkMRMLNode* volumeNode)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	qDebug() << "function::setInputVolume";

	//qvtkReconnect(d->InputVolumeNode, volumeNode, vtkCommand::ModifiedEvent, this, SLOT(updateVolumeInfo()));
	d->InputVolumeNode = volumeNode;
	//qDebug() << volumeNode->GetID();
}

//-----------------------------------------------------------------------------
void qSlicerDRRGeneratorModuleModuleWidget::setOutputVolume(vtkMRMLNode* volumeNode)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	qDebug() << "function::setOutputVolume";

	//qvtkReconnect(d->InputVolumeNode, volumeNode, vtkCommand::ModifiedEvent, this, SLOT(updateVolumeInfo()));
	d->OutputVolumeNode = volumeNode;
	//qDebug() << volumeNode->GetID();
}

class CommandIterationUpdate : public itk::Command
{
public:
	using Self = CommandIterationUpdate;
	using Superclass = itk::Command;
	using Pointer = itk::SmartPointer<Self>;
	itkNewMacro(Self);
protected:
	CommandIterationUpdate() = default;
public:
	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	using OptimizerPointer = const OptimizerType *;
	void Execute(itk::Object *caller, const itk::EventObject & event) override
	{
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event) override
	{
		auto optimizer = static_cast< OptimizerPointer >(object);
		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};

void qSlicerDRRGeneratorModuleModuleWidget::GenerateDRR()
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);

	using clock = std::chrono::system_clock;
	using ms = std::chrono::milliseconds;

	const auto before = clock::now();

	qDebug() << "function::GenerateDRR()";
	std::cout << "rotation :" << d->rx << " " << d->ry << " " << d->rz << "\n";
	std::cout << "translation :" << d->tx << " " << d->ty << " " << d->tz << "\n";
	std::cout << "center :" << d->cx << " " << d->cy << " " << d->cz << "\n";
	std::cout <<  "size :" << d->drrsizex << " " << d->drrsizey << "\n";
	std::cout << "sid :" << d->sid << "\n";


	const int dim = 3;

	typedef short PType;
	typedef itk::Image< PType, dim > InputImageType;
	typedef itk::Image< PType, 2 > OutputImageType;


	vtkSmartPointer<vtkImageData> inputVTK = vtkSmartPointer<vtkImageData>::New();
	inputVTK = d->InputVolumeNode->GetImageData();

	//vtkNew<vtkMatrix4x4> matr;
	//d->InputVolumeNode->GetIJKToRASMatrix(matr);
	//d->OutputVolumeNode->SetIJKToRASMatrix(matr);

	double bounds29[6];
	inputVTK->GetBounds(bounds29);

	// Rotate about the center of the image
	vtkSmartPointer<vtkTransform> transform29 =
		vtkSmartPointer<vtkTransform>::New();

	// Compute the center of the image
	double center29[3];
	center29[0] = (bounds29[1] + bounds29[0]) / 2.0;
	center29[1] = (bounds29[3] + bounds29[2]) / 2.0;
	center29[2] = (bounds29[5] + bounds29[4]) / 2.0;

	// Rotate about the center
	transform29->Translate(center29[0], center29[1], center29[2]);
	if (d->direction == 0) {
		transform29->RotateWXYZ((180), 1, 0, 0);
	}
	else if (d->direction == 1) {
		transform29->RotateWXYZ((180), 0, 1, 0);
	}
	else {
		transform29->RotateWXYZ((180), 0, 0, 1);
	}
	//transform29->RotateWXYZ((180), 0, 0, 1);
	transform29->Translate(-center29[0], -center29[1], -center29[2]);

	// Reslice does all of the work
	vtkSmartPointer<vtkImageReslice> reslice29 =
		vtkSmartPointer<vtkImageReslice>::New();
	reslice29->SetInputData(inputVTK);
	reslice29->SetResliceTransform(transform29);
	reslice29->SetInterpolationModeToCubic();
	reslice29->SetOutputSpacing(
		inputVTK->GetSpacing()[0],
		inputVTK->GetSpacing()[1],
		inputVTK->GetSpacing()[2]);
	reslice29->SetOutputOrigin(
		inputVTK->GetOrigin()[0],
		inputVTK->GetOrigin()[1],
		inputVTK->GetOrigin()[2]);
	reslice29->SetOutputExtent(inputVTK->GetExtent()); // Use a larger extent than the original image's to prevent clipping
	reslice29->Update();

	// Convert VTK -> ITK
	typedef itk::VTKImageToImageFilter< InputImageType > VTKtoITKFilterType;
	VTKtoITKFilterType::Pointer vtkToItkFilter = VTKtoITKFilterType::New();
	//vtkToItkFilter->SetInput(inputVTK);
	vtkToItkFilter->SetInput(reslice29->GetOutput());
	vtkToItkFilter->Update();
	InputImageType::Pointer image = vtkToItkFilter->GetOutput();

	// Software Guide : BeginLatex
	//
	// Creation of a \code{ResampleImageFilter} enables coordinates for
	// each of the pixels in the DRR image to be generated. These
	// coordinates are used by the \code{RayCastInterpolateImageFunction}
	// to determine the equation of each corresponding ray which is cast
	// through the input volume.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using FilterType = itk::ResampleImageFilter<InputImageType, InputImageType >;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(image);
	filter->SetDefaultPixelValue(0);
	//qDebug() << "function::projectVolume2()::ResampleImageFilter";

	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// An Euler transformation is defined to position the input volume.
	// The \code{ResampleImageFilter} uses this transform to position the
	// output DRR image for the desired view.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using TransformType = itk::CenteredEuler3DTransform< double >;
	TransformType::Pointer transform = TransformType::New();
	transform->SetComputeZYX(true);
	TransformType::OutputVectorType translation;
	translation[0] = d->tx; //left -right Increasing this moves the model right
	translation[1] = d->ty; //THIS WAS ORIGINALLY 0 CHANGED TO 200 Increasing this zooms out further
	translation[2] = d->tz;// tz; up -down Increasing this moves the model downwards
					   // constant for converting degrees into radians
	const double dtr = (std::atan(1.0) * 4.0) / 180.0;
	transform->SetTranslation(translation);
	transform->SetRotation(dtr*d->rx, dtr*d->ry, dtr*d->rz);
	//transform->SetRotation(dtr*0., dtr*0., dtr*0.);
	//transform->SetRotation(dtr*90., dtr*0., dtr*0.); //TODO change 3rd coordinate to 20
													 //transform->SetRotation(dtr*(90.+x), dtr*y, dtr*z);
													 //forward-back, clockwise-ccw, 


	InputImageType::PointType   imOrigin = image->GetOrigin();
	InputImageType::SpacingType imRes = image->GetSpacing();
	using InputImageRegionType = InputImageType::RegionType;
	using InputImageSizeType = InputImageRegionType::SizeType;
	InputImageRegionType imRegion = image->GetBufferedRegion();
	InputImageSizeType   imSize = imRegion.GetSize();
	imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
	imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
	imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;
	TransformType::InputPointType center;
	center[0] = d->cx + imOrigin[0];
	center[1] = d->cy + imOrigin[1];
	center[2] = d->cz + imOrigin[2];
	transform->SetCenter(center);
	/*center[0] = 0 + imOrigin[0];
	center[1] = 0 + imOrigin[1];
	center[2] = 0 + imOrigin[2];
	transform->SetCenter(center);*/
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The \code{RayCastInterpolateImageFunction} is instantiated and passed the transform
	// object. The \code{RayCastInterpolateImageFunction} uses this
	// transform to reposition the x-ray source such that the DRR image
	// and x-ray source move as one around the input volume. This coupling
	// mimics the rigid geometry of the x-ray gantry.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using InterpolatorType =
		itk::RayCastInterpolateImageFunction<InputImageType, double>;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// We can then specify a threshold above which the volume's
	// intensities will be integrated.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet

	interpolator->SetThreshold(d->drrthreshold);

								  // Software Guide : EndCodeSnippet
								  // Software Guide : BeginLatex
								  //
								  // The ray-cast interpolator needs to know the initial position of the
								  // ray source or focal point. In this example we place the input
								  // volume at the origin and halfway between the ray source and the
								  // screen. The distance between the ray source and the screen
								  // is the "source to image distance" \code{sid} and is specified by
								  // the user.
								  //
								  // Software Guide : EndLatex
								  // Software Guide : BeginCodeSnippet
	InterpolatorType::InputPointType focalpoint;
	focalpoint[0] = imOrigin[0];
	focalpoint[1] = imOrigin[1];
	focalpoint[2] = imOrigin[2] - d->sid / 2.;
	interpolator->SetFocalPoint(focalpoint);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// Having initialised the interpolator we pass the object to the
	// resample filter.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	//interpolator->Print(std::cout);
	filter->SetInterpolator(interpolator);
	filter->SetTransform(transform);
	// Software Guide : EndCodeSnippet
	// Software Guide : BeginLatex
	//
	// The size and resolution of the output DRR image is specified via the
	// resample filter.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	// setup the scene
	InputImageType::SizeType   size;
	size[0] = d->drrsizex;//512;//dx;  // number of pixels along X of the 2D DRR image
	size[1] = d->drrsizey;// dy;  // number of pixels along Y of the 2D DRR image
	size[2] = 1;   // only one slice
	filter->SetSize(size);
	InputImageType::SpacingType spacing;
	spacing[0] = 1;// sx;  // pixel spacing along X of the 2D DRR image [mm]
	spacing[1] = 1;// sy;  // pixel spacing along Y of the 2D DRR image [mm]
	spacing[2] = 1; // slice thickness of the 2D DRR image [mm]
	filter->SetOutputSpacing(spacing);
	// Software Guide : EndCodeSnippet

	// Software Guide : BeginLatex
	//
	// In addition the position of the DRR is specified. The default
	// position of the input volume, prior to its transformation is
	// half-way between the ray source and screen and unless specified
	// otherwise the normal from the "screen" to the ray source passes
	// directly through the centre of the DRR.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	double origin[dim];
	origin[0] = imOrigin[0] + 0 - 1.*((double)d->drrsizex - 1.) / 2.;//o2Dx - sx*((double)dx - 1.) / 2.;
	origin[1] = imOrigin[1] + 0 - 1.*((double)d->drrsizey - 1.) / 2.;//o2Dy - sy*((double)dy - 1.) / 2.;
	origin[2] = imOrigin[2] + d->sid / 2.;
	filter->SetOutputOrigin(origin);
	// Software Guide : EndCodeSnippet
	// create writer
	//qDebug() << "function::projectVolume2()::RayCastImageFilter";

	// Software Guide : BeginLatex
	//
	// The output of the resample filter can then be passed to a writer to
	// save the DRR image to a file.
	//
	// Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	using RescaleFilterType = itk::RescaleIntensityImageFilter<
		InputImageType, InputImageType >;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255); //TODO USED TO BE 255
	rescaler->SetInput(filter->GetOutput());
	rescaler->Update();

	typedef itk::ImageFileWriter< InputImageType >  WriterType;
	WriterType::Pointer writer9 = WriterType::New();
	writer9->SetFileName("output_fixed.nrrd");
	writer9->SetInput(rescaler->GetOutput());
	writer9->Update();

	typedef itk::ExtractImageFilter< InputImageType, OutputImageType > SliceType3;
	SliceType3::Pointer slice3 = SliceType3::New();
	slice3->InPlaceOn();
	slice3->SetDirectionCollapseToSubmatrix();
	InputImageType::RegionType inputRegion3 =
		rescaler->GetOutput()->GetLargestPossibleRegion();
	InputImageType::SizeType size3 = inputRegion3.GetSize();
	size3[2] = 0;
	InputImageType::IndexType start3 = inputRegion3.GetIndex();
	const unsigned int sliceNumber3 = 0;
	start3[2] = sliceNumber3;
	InputImageType::RegionType desiredRegion3;
	desiredRegion3.SetSize(size3);
	desiredRegion3.SetIndex(start3);
	slice3->SetExtractionRegion(desiredRegion3);
	//qDebug() << "function::projectVolume2()::" << size3[0] << " " << size3[1] << " ~ " << start3[0] << " " << start3[1];
	slice3->SetInput(rescaler->GetOutput());
	slice3->Update();

	/*WriterType::Pointer writer4 = WriterType::New();
	writer4->SetFileName("output4.jpeg");

	using ExtractFilterType = itk::ExtractImageFilter< OutputImageType, OutputImageType >;
	ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
	extractFilter->SetDirectionCollapseToSubmatrix();

	// set up the extraction region [one slice]
	const OutputImageType * inputImage = rescaler->GetOutput();
	OutputImageType::RegionType inputRegion = inputImage->GetBufferedRegion();
	OutputImageType::SizeType size1 = inputRegion.GetSize();
	size1[2] = 1; // we extract along z direction
	OutputImageType::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = 0;//atoi(argv[3]);
	start[2] = sliceNumber;
	OutputImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	extractFilter->SetExtractionRegion(desiredRegion);
	using PasteFilterType = itk::PasteImageFilter< OutputImageType >;
	PasteFilterType::Pointer pasteFilter = PasteFilterType::New();
	using MedianFilterType = itk::MedianImageFilter< OutputImageType,
	OutputImageType >;
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
	extractFilter->SetInput(inputImage);
	medianFilter->SetInput(extractFilter->GetOutput());
	pasteFilter->SetSourceImage(medianFilter->GetOutput());
	pasteFilter->SetDestinationImage(inputImage);
	pasteFilter->SetDestinationIndex(start);

	OutputImageType::SizeType indexRadius;
	indexRadius[0] = 1; // radius along x
	indexRadius[1] = 1; // radius along y
	indexRadius[2] = 0; // radius along z
	medianFilter->SetRadius(indexRadius);
	medianFilter->UpdateLargestPossibleRegion();
	const OutputImageType * medianImage = medianFilter->GetOutput();
	pasteFilter->SetSourceRegion(medianImage->GetBufferedRegion());
	writer4->SetInput(pasteFilter->GetOutput());
	writer4->Update();*/

	/*OutputImageType::IndexType start1;
	start1[0] = 0;  // first index on X
	start1[1] = 0;  // first index on Y
	start1[2] = 0;  // first index on Z
	OutputImageType::SizeType  size1;
	size1[0] = 490;  // size along X
	size1[1] = 490;  // size along Y
	size1[2] = 1;  // size along Z
	OutputImageType::RegionType region;
	region.SetSize(size1);
	region.SetIndex(start1);

	typedef itk::ExtractImageFilter< OutputImageType, OutputImageType2 > ExtractImageFilterType;
	ExtractImageFilterType::Pointer extract = ExtractImageFilterType::New();
	extract->SetInput(rescaler->GetOutput());
	extract->SetDirectionCollapseToIdentity();
	extract->SetExtractionRegion(region);
	extract->Update();
	using WriterType4 = itk::ImageFileWriter< OutputImageType2 >;
	WriterType4::Pointer writer4 = WriterType4::New();
	writer4->SetFileName("output4.jpeg");
	writer4->SetInput(extract->GetOutput());*/

	using WriterType11 = itk::ImageFileWriter< OutputImageType >;
	WriterType11::Pointer writer11 = WriterType11::New();
	writer11->SetFileName("output_fixed2.nrrd");
	writer11->SetInput(slice3->GetOutput());
	writer11->Update();


	using ITKtoVTKFilterType = itk::ImageToVTKImageFilter< OutputImageType >;
	ITKtoVTKFilterType::Pointer itktovtkfilter = ITKtoVTKFilterType::New();
	itktovtkfilter->SetInput(slice3->GetOutput());
	itktovtkfilter->Update();
	vtkImageData * outputVTK = itktovtkfilter->GetOutput();

	/*vtkSmartPointer<vtkImageReslice> reslice =
	vtkSmartPointer<vtkImageReslice>::New();
	reslice->SetOutputExtent(0, 9, 0, 100, 0, 0);
	reslice->SetInputData(outputVTK);
	reslice->Update();*/

	/*vtkSmartPointer<vtkImageWriter> writer5 =
	vtkSmartPointer<vtkImageWriter>::New();
	writer5->SetInputData(outputVTK);
	writer5->SetFileName("output3.vtk");
	writer5->Write();*/

	vtkSmartPointer<vtkImageCast> castFilter =
		vtkSmartPointer<vtkImageCast>::New();
	castFilter->SetOutputScalarTypeToUnsignedChar();
	castFilter->SetInputData(outputVTK);
	castFilter->Update();

	d->OutputVolumeNode->SetAndObserveImageData(castFilter->GetOutput());

	const auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

	std::cout << "It took " << duration.count() / 1000.0 << "s" << std::endl;
}

void qSlicerDRRGeneratorModuleModuleWidget::setRotationX(double rotationX)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->rx = rotationX;
}

void qSlicerDRRGeneratorModuleModuleWidget::setRotationY(double rotationY)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->ry = rotationY;
}

void qSlicerDRRGeneratorModuleModuleWidget::setRotationZ(double rotationZ)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->rz = rotationZ;
}

void qSlicerDRRGeneratorModuleModuleWidget::setTranslationX(double translationX)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->tx = translationX;
}

void qSlicerDRRGeneratorModuleModuleWidget::setTranslationY(double translationY)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->ty = translationY;
}

void qSlicerDRRGeneratorModuleModuleWidget::setTranslationZ(double translationZ)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->tz = translationZ;
}

void qSlicerDRRGeneratorModuleModuleWidget::setCenterX(double centerX)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->cx = centerX;
}

void qSlicerDRRGeneratorModuleModuleWidget::setCenterY(double centerY)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->cy = centerY;
}

void qSlicerDRRGeneratorModuleModuleWidget::setCenterZ(double centerZ)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->cz = centerZ;
}

void qSlicerDRRGeneratorModuleModuleWidget::setFocalPoint(double focalPoint)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->sid = focalPoint;
}

void qSlicerDRRGeneratorModuleModuleWidget::setThreshold(double threshold)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->drrthreshold = threshold;
}

void qSlicerDRRGeneratorModuleModuleWidget::setDRRSizeX(int sizeX)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->drrsizex = sizeX;
}

void qSlicerDRRGeneratorModuleModuleWidget::setDRRSizeY(int sizeY)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->drrsizey = sizeY;
}

void qSlicerDRRGeneratorModuleModuleWidget::setDirection(int direction)
{
	Q_D(qSlicerDRRGeneratorModuleModuleWidget);
	d->direction = direction;
}
