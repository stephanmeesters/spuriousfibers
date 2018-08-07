#ifndef vtkFiberSelectAnterior_h
#define vtkFiberSelectAnterior_h


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


class vtkFiberSelectAnterior : public vtkPolyDataAlgorithm
{
	public:

		/** Constructor Call */

		static vtkFiberSelectAnterior * New();

		/** VTK Macro */

		vtkTypeMacro(vtkFiberSelectAnterior,vtkPolyDataAlgorithm);

        void SetFiberTransformationMatrix(vtkMatrix4x4* m)
        {
            this->fiberMatrix = m;
        }

        void SetNumberOfAnteriorFibers(int b)
        {
            this->numberOfAnteriorFibers = b;
        }
    

	protected:

		/** Main entry point of the filter. */

		//virtual void Execute();

		/** Constructor. */

		vtkFiberSelectAnterior();

		/** Destructor. */

		~vtkFiberSelectAnterior();
    
    
        int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    

        /** Fiber transformation matrix */
        vtkMatrix4x4* fiberMatrix;

        /** Number of anterior fibers */
        int numberOfAnteriorFibers;

        /** Selected scalar value type */
		int scalarType;
    
    private:
    
        vtkFiberSelectAnterior(const vtkFiberSelectAnterior&);  // Not implemented.
        void operator=(const vtkFiberSelectAnterior&);  // Not implemented.

}; // class vtkFiberSelectAnterior


#endif // vtkFiberSelectAnterior_h
