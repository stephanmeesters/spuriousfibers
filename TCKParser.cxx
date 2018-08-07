//
//  TCKParser.cpp
//  SpuriousFibers
//
//  Created by Stephan Meesters on 02/06/15.
//
//

#include "TCKParser.h"

vtkPolyData* TCKParser::LoadDataFromFile(std::string filename)
{
    // Create the Qt file handler
    ifstream TCKFile (filename.c_str(), ios::in | ios::binary );
    
    // Try to open the input file
    if (TCKFile.fail())
    {
        printf("Could not open file %s!\n",filename.c_str());
        return NULL;
    }

    // temp variables
    char c;

    char endBuffer[] = "AAA";
    char end[] = "END";
    for(int i = 0; i<1000; i++)
    {
        TCKFile.read(reinterpret_cast<char*>(&c), sizeof(char));
        
        
        endBuffer[0] = endBuffer[1];
        endBuffer[1] = endBuffer[2];
        endBuffer[2] = c;
        
        if(!strcmp(endBuffer,end))
            break;
    }
    
    double inf=1.0/0.0;
    std::vector< std::vector<float> > fibersList;
    
    int headerPos = TCKFile.tellg();
    int bugfixer = 0;
    int MAXTRIES = 200;
    while(bugfixer < MAXTRIES)
    {
        TCKFile.seekg(headerPos, TCKFile.beg);
        
        for(int i = 0; i<bugfixer; i++)
        {
            TCKFile.read(reinterpret_cast<char*>(&c), sizeof(char));
        }
        
        int j = 0;
        fibersList.clear();
        
        int pos;
        float f;
        while(true)
        {
            pos = TCKFile.tellg();
            TCKFile.read(reinterpret_cast<char*>(&f), sizeof(float));
            if(fabs(f) > 0.01)
                break;
        }
        TCKFile.seekg(pos, TCKFile.beg);
        
        int k = 0;
        bool abort;
        while(true)
        {
            //algo->UpdateProgress((double) j);
            
            std::vector<float> fiber;
            
            float f1,f2,f3;
            k = 0;
            abort = false;
            while(true)
            {
                TCKFile.read(reinterpret_cast<char*>(&f1), sizeof(float));
                TCKFile.read(reinterpret_cast<char*>(&f2), sizeof(float));
                TCKFile.read(reinterpret_cast<char*>(&f3), sizeof(float));
                
                if(f1 != f1 && f2!=f2 && f3!=f3)
                    break;
                
                if(f1 == inf && f2==inf && f3==inf)
                    break;
                
                if(TCKFile.eof())
                    break;
                
                //                f1 += 80;
                //                f1 *= 0.5;
                //
                //                f2 += 120;
                //                f2 *= 0.5;
                //
                //                f3 += 60;
                //                f3 *= 0.5;
                
                fiber.push_back(f1);
                fiber.push_back(f2);
                fiber.push_back(f3);
                
                k++;
                
                if(f1 > 10000 || f2 > 10000 || f3 > 10000)
                {
                    abort = true;
                    break;
                }
            }
            
            if(abort)
            {
                break;
            }
            
            fibersList.push_back(fiber);
            
            j++;
            
            if(TCKFile.eof())
                break;
            
            if(f1 == inf && f2==inf && f3==inf)
                break;
        }
        
        if(abort)
        {
            printf("Bugfix attempt: %d\n",bugfixer);
            bugfixer++;
        }
        else
            break;
        
    }
    
    if(bugfixer == MAXTRIES)
    {
        printf("Error loading %s!\n",filename.c_str());
        return NULL;
    }
    
    //
    //  Transformation
    //
    
    vtkMatrix4x4 * mat = vtkMatrix4x4::New();
    mat->Identity();
    
    vtkTransform* transform = vtkTransform::New();
    transform->SetMatrix(mat);
    //    transform->Translate(-80,-120,-60);
    //    transform->Scale(2,2,2);
    mat = transform->GetMatrix();
    
    // close .TCK file
    TCKFile.close();

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