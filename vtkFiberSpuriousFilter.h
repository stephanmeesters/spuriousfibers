#ifndef vtkFiberSpuriousFilter_h
#define vtkFiberSpuriousFilter_h

#include <iostream>
#include <fstream>

/** Includes - VTK */

#include <vtkPolyDataAlgorithm.h>
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkObjectFactory.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkGenericCell.h>
#include <vtkCallbackCommand.h>
#include <vtkSmartPointer.h>

/** Includes - Custom Files */

#include "SpuriousFiberFilterTypes.h"

#define PI 3.14159265358
#define PISQUARED 9.86960440109
#define PIHALVED = 1.57079632679

/** This class is used to
*/

class vtkFiberSpuriousFilter : public vtkPolyDataAlgorithm
{
	public:

		/** Constructor Call */

		static vtkFiberSpuriousFilter * New();

		/** VTK Macro */

		vtkTypeMacro(vtkFiberSpuriousFilter, vtkPolyDataAlgorithm);

		/** Set a new input volume.
			@param image		New discrete sphere function volume. */

		void SetInputVolume(vtkImageData * image);

		/** Set parameter settings
			@param ps		Parameter settings struct. */

		void SetParameters(ParameterSettings* ps);

        /** Method for curve generation **/

		enum TypeOfCurve
        {
            CURVE_TYPE_GEODESIC = 0,
            CURVE_TYPE_ELASTICA
        };


	protected:

		/** Main entry point of the filter. */

		int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

		/** Constructor. */

		vtkFiberSpuriousFilter();

		/** Destructor. */

		~vtkFiberSpuriousFilter();

		inline double k2(double* x, double* y, double* r, double* v);

		/** Input volume, containing an array defining the spherical directions,
			an array with the radius per direction per voxel, and (optionally)
			a triangles array defining the topology of the glyphs. */

		vtkImageData * inputVolume;
		ParameterSettings* ps;

		double cutoffSquared;

		double orientations[486]; // 162 * x,y,z
        double kernelOutput[14580];

        double kernelDefault;
    
        static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);


}; // class vtkFiberSpuriousFilter


#endif // vtkFiberSpuriousFilter_h
