#ifndef vtkFiberLengthFilter_h
#define vtkFiberLengthFilter_h


/** Includes - STD */

#include <algorithm>
#include <vector>
#include <map>

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
#include <vtkImageCast.h>
#include <vtkObject.h>
#include <vtkMatrix4x4.h>


class vtkFiberLengthFilter : public vtkPolyDataAlgorithm
{
	public:

		/** Constructor Call */

		static vtkFiberLengthFilter * New();

		/** VTK Macro */

		vtkTypeMacro(vtkFiberLengthFilter,vtkPolyDataAlgorithm);

        void SetMaxFiberLength(int b)
        {
            this->fiberLengthMax = b;
        }
    
        void SetMinFiberLength(int b)
        {
            this->fiberLengthMin = b;
        }
    

	protected:

		/** Main entry point of the filter. */

		int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

		/** Constructor. */

		vtkFiberLengthFilter();

		/** Destructor. */

		~vtkFiberLengthFilter();

        /** Fiber length restrictions */
        int fiberLengthMax;
        int fiberLengthMin;

        /** Selected scalar value type */
		int scalarType;
    
    private:
    
        vtkFiberLengthFilter(const vtkFiberLengthFilter&);  // Not implemented.
        void operator=(const vtkFiberLengthFilter&);  // Not implemented.

}; // class vtkFiberLengthFilter


#endif // vtkFiberLengthFilter_h
