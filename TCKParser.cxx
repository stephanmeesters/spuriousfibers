//
//  TCKParser.cpp
//  SpuriousFibers
//
//  Created by Stephan Meesters on 02/06/15.
//
//

#include "TCKParser.h"
#include <sstream>

vtkPolyData* TCKParser::LoadDataFromFile(std::string filename)
{
    // Create the Qt file handler
    ifstream TCKFile (filename.c_str(), ios::in | ios::binary );

    // Read header
    std::string line;
    std::string end = "END";
    std::string fiberCount = "count:";
    int numfibers;
    while (std::getline(TCKFile, line))
    {
        //printf("%s\n", line.c_str());
        if (line.substr(0,6).find(fiberCount) != std::string::npos)
        {
            numfibers = std::stoi(line.substr(6));
            // printf("number of fibers found: %d\n", numfibers);
        }
        if(line.compare(end) == 0)
        {
            // printf("end found.\n");
            break;
        }
    }

    // printf("pos:%d\n",(int)TCKFile.tellg());

    TCKFile.seekg((int)TCKFile.tellg()-1);

    // Load fiber data
    std::vector< std::vector<float> > fibersList;
    fibersList.resize(numfibers);
    int idx = 0;
    float f1,f2,f3;
    TCKFile.read(reinterpret_cast<char*>(&f1), sizeof(float));
    while(true)
    {
        if(TCKFile.eof())
            break;

        TCKFile.read(reinterpret_cast<char*>(&f1), sizeof(float));
        TCKFile.read(reinterpret_cast<char*>(&f2), sizeof(float));
        TCKFile.read(reinterpret_cast<char*>(&f3), sizeof(float));

        // if(isnan(f1) || isnan(f2) || isnan(f3))
        //     printf("%f %f %f\n", f1,f2,f3);

        if(isnan(f1) && isnan(f2) && isnan(f3))
        {
            idx++;
            // printf("now processing fiber %d\n", idx);
            continue;
        }

        if(isinf(f1) && isinf(f2) && isinf(f3))
        {
            // printf("end of file reached!\n");
            break;
        }

//        std::vector<float> fiber = fibersList.at(idx);

        fibersList[idx].push_back(f1);
        fibersList[idx].push_back(f2);
        fibersList[idx].push_back(f3);
    }

    // close .TCK file
    TCKFile.close();

    
    //
    //  Transformation (not used)
    //
    
    vtkMatrix4x4 * mat = vtkMatrix4x4::New();
    mat->Identity();
    
    vtkTransform* transform = vtkTransform::New();
    transform->SetMatrix(mat);
    //    transform->Translate(-80,-120,-60);
    //    transform->Scale(2,2,2);
    mat = transform->GetMatrix();

    // Create polydata
    vtkPolyData* output = vtkPolyData::New();
    
    // Create a point set for the output
    vtkPoints * outputPoints = vtkPoints::New();
    output->SetPoints(outputPoints);
    outputPoints->Delete();
    
    // Cell array holding the pathways.
    // Each pathway is a single cell (line).
    // Each cell holds the id values to points from the vtkPoints list
    vtkCellArray * outputLines = vtkCellArray::New();
    output->SetLines(outputLines);
    outputLines->Delete();
    
    // Loop over pathways
    int counter = 0;
    int numPathways = fibersList.size();
    for(int i = 0; i < numPathways; i++)
    {
        std::vector<float> fiber = fibersList.at(i);
        
        int numberOfFiberPoints = fiber.size()/3;
        
        // Create a cell representing a fiber
        outputLines->InsertNextCell(numberOfFiberPoints);
        
        // Loop over points in the pathway
        for(int j = 0; j<numberOfFiberPoints; j++)
        {
            outputPoints->InsertNextPoint(fiber[j*3],fiber[j*3+1],fiber[j*3+2]);
            outputLines->InsertCellPoint(counter + j);
        }
        
        counter += numberOfFiberPoints;
    }
    
    // todo: cleanup vars
    
    printf("Fibers loaded succesfully. Nr of fibers=%d\n",(int)fibersList.size());
    
    // return the output PolyData
    return output;
}
