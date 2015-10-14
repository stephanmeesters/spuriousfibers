
/** Includes */

#include "vtkFiberLengthFilter.h"


vtkStandardNewMacro(vtkFiberLengthFilter);


//-----------------------------[ Constructor ]-----------------------------

vtkFiberLengthFilter::vtkFiberLengthFilter()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}


//------------------------------[ Destructor ]-----------------------------

vtkFiberLengthFilter::~vtkFiberLengthFilter()
{

}


//-------------------------------[ Execute ]-------------------------------

int vtkFiberLengthFilter::RequestData(vtkInformation *vtkNotUsed(request),
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
    QList<vtkDoubleArray*> outputScalarsList;
    int numberOfScalarTypes = inputPD->GetNumberOfArrays();
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = vtkDoubleArray::New();
        outputScalars->SetName(inputPD->GetArray(i)->GetName());
        outputScalarsList.append(outputScalars);
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
    
    bool applyMin = fiberLengthMin > 0;
    bool applyMax = fiberLengthMax > 0;
    
    QMap<double, vtkIdType> fiberMap;
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
        
        if((applyMax && numberOfFiberPoints > fiberLengthMax) ||
           (applyMin && numberOfFiberPoints < fiberLengthMin))
            continue;
        
        // Evaluate if the fiber should be included in the output fibers
        //double maxAnterior = 0;
        //fiberMap.insert(maxAnterior, lineId);
        
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

