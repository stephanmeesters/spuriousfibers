
/** Includes */

#include "vtkFiberSelectAnterior.h"


vtkStandardNewMacro(vtkFiberSelectAnterior);


//-----------------------------[ Constructor ]-----------------------------

vtkFiberSelectAnterior::vtkFiberSelectAnterior()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}


//------------------------------[ Destructor ]-----------------------------

vtkFiberSelectAnterior::~vtkFiberSelectAnterior()
{

}


//-------------------------------[ Execute ]-------------------------------

int vtkFiberSelectAnterior::RequestData(vtkInformation *vtkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector)
{
    // Get the input
    vtkPolyData * input = vtkPolyData::GetData(inputVector[0],0);
    if(!input)
    {
        vtkErrorMacro(<< "Input has not been set.");
        return 0;
    }
    
    // Check if the input contains point data
    vtkPointData * inputPD = input->GetPointData();
    if (!inputPD)
    {
        vtkErrorMacro(<< "Input does not have point data.");
        return 0;
    }
    
    // Get the points of the input
    vtkPoints * inputPoints = input->GetPoints();
    if (!inputPoints)
    {
        vtkErrorMacro(<< "Input does not have points.");
        return 0;
    }
    
    // Get the lines array of the input
    vtkCellArray * inputLines = input->GetLines();
    if (!inputLines)
    {
        vtkErrorMacro(<< "Input does not have lines.");
        return 0;
    }
    
    // Get the output
    vtkPolyData * output = vtkPolyData::GetData(outputVector,0);
    if (!output)
    {
        vtkErrorMacro(<< "Output has not been set.");
        return 0;
    }
    
    // Check if the output contains point data
    vtkPointData * outputPD = output->GetPointData();
    if (!outputPD)
    {
        vtkErrorMacro(<< "Output does not have point data.");
        return 0;
    }
    
    // Get list of scalars
    std::vector<vtkDoubleArray*> outputScalarsList;
    int numberOfScalarTypes = inputPD->GetNumberOfArrays();
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = vtkDoubleArray::New();
        outputScalars->SetName(inputPD->GetArray(i)->GetName());
        outputScalarsList.push_back(outputScalars);
    }
    
    // Create a point set for the output
    vtkPoints * outputPoints = vtkPoints::New();
    output->SetPoints(outputPoints);
    outputPoints->Delete();
    
    // Create a line array for the output
    vtkCellArray * outputLines = vtkCellArray::New();
    output->SetLines(outputLines);
    outputLines->Delete();
    
    // Number of points in the current fiber, and a list of its point IDs
    vtkIdType numberOfPoints;
    vtkIdType * pointList;
    
    // Setup progress bar
    int numberOfCells = inputLines->GetNumberOfCells();
    int progressStep = numberOfCells / 25;
    progressStep += (progressStep == 0) ? 1 : 0;
    this->SetProgressText("Selecting anterior fibers...");
    this->UpdateProgress(0.0);
    
    std::map<double, vtkIdType> fiberMap;
    // Loop through all input fibers and get anterior distance
    for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
    {
        // Update the progress bar
        if ((lineId % progressStep) == 0)
        {
            this->UpdateProgress((double) lineId / (double) numberOfCells);
        }
        
        // Get the data of the current fiber
        vtkCell * currentCell = input->GetCell(lineId);
        int numberOfFiberPoints = currentCell->GetNumberOfPoints();
        
        // Get cell containing fiber
        inputLines->GetNextCell(numberOfPoints, pointList);
        
        // Evaluate if the fiber should be included in the output fibers
        double p[3];
        double maxAnterior = -1e99;
        for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
        {
            // Get the point ID of the current fiber point
            vtkIdType currentPointId = currentCell->GetPointId(pointId);
            
            // Get point coordinates
            inputPoints->GetPoint(currentPointId, p);
            
            // Check if point is within ROI
            double vec[4] = {p[0],p[1],p[2],1};
            fiberMatrix->MultiplyPoint(vec,vec);
            
            if(vec[1] > maxAnterior)
                maxAnterior = vec[1];
        }
        
        // Add the fiber ID and the CM value to the map
        fiberMap.insert( std::pair<double, vtkIdType>(maxAnterior, lineId) );
        //printf("LINE: %f %d \n",maxAnterior,(int)lineId);
    }
    
    int rFiberIndex = 0;
    std::map<double, vtkIdType>::iterator rIter = fiberMap.end();
    
    while (rIter != fiberMap.begin())
    {
        // Decrement the iterator
        rIter--;
        
        // Update the progress bar
        if ((rFiberIndex % progressStep) == 0)
        {
            this->UpdateProgress((double) rFiberIndex / (double) numberOfAnteriorFibers);
        }
        
        // Get the ID of the current fiber
        vtkIdType currentFiberId = rIter->second;
        
        // Get the data of the current fiber
        vtkCell * currentCell = input->GetCell(currentFiberId);
        int numberOfFiberPoints = currentCell->GetNumberOfPoints();
        
        // Create an ID list for the output fiber
        vtkIdList * newFiberList = vtkIdList::New();
        
        // Current point coordinates
        double p[3];
        
        // Loop through all points in the fiber
        for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
        {
            // Get the point ID of the current fiber point
            vtkIdType currentPointId = currentCell->GetPointId(pointId);
            
            // Copy the point coordinates to the output
            inputPoints->GetPoint(currentPointId, p);
            
            vtkIdType newPointId = outputPoints->InsertNextPoint(p);
            newFiberList->InsertNextId(newPointId);
            
            // include old scalar values
            for(int i = 0; i < numberOfScalarTypes; i++)
            {
                // Get the scalar value
                double scalar = inputPD->GetArray(i)->GetTuple1(currentPointId);
                
                // Copy the scalar value to the output
                outputScalarsList.at(i)->InsertNextTuple1(scalar);
            }
            
        }
        
        // Add the new fiber to the output
        outputLines->InsertNextCell(newFiberList);
        
        //printf("%d ---- %d\n",numberOfAnteriorFibers,rFiberIndex);
        
        // Break if we've reached the desired number of fibers
        if (++rFiberIndex == numberOfAnteriorFibers)
            break;
        
    }
    
    // Add back other scalar arrays
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = outputScalarsList.at(i);
        output->GetPointData()->AddArray(outputScalars);
    }
    
    // Finalize the progress bar
    this->UpdateProgress(1.0);
    
    return 1;
}

