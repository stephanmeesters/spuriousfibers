//
//  main.cpp
//  SpuriousFibers
//
//  Created by Stephan Meesters on 02/06/15.
//  Copyright (c) 2015 Stephan Meesters. All rights reserved.
//

#include <iostream>

/** Includes - TCLAP */
#include "tclap/CmdLine.h"

/** Includes -- project */
#include "TCKParser.h"
#include "vtkFiberSelectAnterior.h"
#include "SpuriousFiberFilterTypes.h"
#include "vtkFiberSpuriousFilter.h"
#include "vtkFiberLengthFilter.h"

/** Includes - VTK */
#include <vtkPolyData.h>
#include <vtkSplineFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>

/** Includes - OMP */
#ifdef WITH_OPENMP
	#include <omp.h>
#endif

void SpuriousFibers(ParameterSettings* ps)
{
    // load fibers
    printf("Loading fibers from dataset: %s \n",ps->inputFile.c_str());
    
    std::string fn = ps->inputFile;    
    vtkPolyData* polydata = NULL;
    // parse TCK if file if .tck
    if(fn.substr(fn.find_last_of(".") + 1) == "tck")
    {
        polydata = TCKParser::LoadDataFromFile(ps->inputFile);
    }
    
    // parse VTK if file if .vtk
//    else if(infileInfo.suffix() == QString("vtk"))
//    {
//        vtkSmartPointer<vtkPolyDataReader> reader =
//            vtkSmartPointer<vtkPolyDataReader>::New();
//        reader->SetFileName(ps->inputFile.c_str());
//        reader->Update();
//        polydata = static_cast<vtkPolyData*>(reader->GetOutput());
//    }
    
    // unknown format, abort
    else
    {
        printf("Unknown input file extension.\n");
        return;
    }
    
    // check if we actually loaded fibers
    if(polydata == NULL || !polydata || polydata->GetLines() == NULL ||  polydata->GetLines()->GetNumberOfCells() == 0)
    {
        printf("Error loading input fibers.\n");
        return;
    }
    
    //polydata->Print(std::cout);
    
    // Fiber length filter
    vtkFiberLengthFilter* lengthFilter;
    if(ps->applyFilterLength)
    {
        lengthFilter = vtkFiberLengthFilter::New();
        lengthFilter->SetInputData(polydata);
        lengthFilter->SetMaxFiberLength(ps->maxLength);
        lengthFilter->SetMinFiberLength(ps->minLength);
        lengthFilter->Update();
        //lengthFilter->GetOutput()->Print(std::cout);
        printf("Length selection complete. Remaining fibers=%d\n",(int)lengthFilter->GetOutput()->GetLines()->GetNumberOfCells());
    }
    
    // Select anterior fibers
    if((ps->applyFilterLength && ps->numberOfAnteriorFibers >= (int)lengthFilter->GetOutput()->GetLines()->GetNumberOfCells()) ||
       ps->numberOfAnteriorFibers >= (int)polydata->GetLines()->GetNumberOfCells())
        ps->applySelectAnterior = false; // select all fibers, so don't bother applying the filter
    vtkFiberSelectAnterior* anteriorFilter;
    if(ps->applySelectAnterior)
    {
        anteriorFilter = vtkFiberSelectAnterior::New();
        if(ps->applyFilterLength)
            anteriorFilter->SetInputData(lengthFilter->GetOutput());
        else
            anteriorFilter->SetInputData(polydata);
        anteriorFilter->SetNumberOfAnteriorFibers(ps->numberOfAnteriorFibers);
        
        // Get the transformation matrix
        vtkMatrix4x4* transformationMatrix = vtkMatrix4x4::New();
        transformationMatrix->Identity();
        anteriorFilter->SetFiberTransformationMatrix(transformationMatrix);
        anteriorFilter->Update();
        printf("Anterior fiber selection complete (#=%d).\n",ps->numberOfAnteriorFibers);
        
        //anteriorFilter->GetOutput()->Print(std::cout);
    }
    
    // Spline sample
    vtkSplineFilter* splineFilter;
    if(ps->applySubsampling)
    {
        splineFilter = vtkSplineFilter::New();
        if(ps->applySelectAnterior)
            splineFilter->SetInputData(anteriorFilter->GetOutput());
        else if(ps->applyFilterLength)
            splineFilter->SetInputData(lengthFilter->GetOutput());
        else
            splineFilter->SetInputData(polydata);
        splineFilter->SetSubdivideToLength();
        splineFilter->SetLength(ps->samplingStep);
        splineFilter->SetGenerateTCoordsToOff();
        splineFilter->Update();
        printf("Fiber subsampling complete.\n");
        
        //splineFilter->GetOutput()->Print(std::cout);
    }
    
    // check if there are any fibers remaining, otherwise abort
    
    
    // Spurious fiber filter
    vtkFiberSpuriousFilter* scoringFilter = vtkFiberSpuriousFilter::New();
    if(ps->applySubsampling)
        scoringFilter->SetInputData(splineFilter->GetOutput());
    else if(ps->applySelectAnterior)
        scoringFilter->SetInputData(anteriorFilter->GetOutput());
    else if(ps->applyFilterLength)
        scoringFilter->SetInputData(lengthFilter->GetOutput());
    else
        scoringFilter->SetInputData(polydata);
    scoringFilter->SetParameters(ps);
    scoringFilter->Update();
    printf("Scoring complete.\n");
    
    // Write the results
    vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
    writer->SetFileName(ps->outputFile.c_str());
    writer->SetInputData(scoringFilter->GetOutput());
    writer->SetFileTypeToASCII();
    writer->Write();
    writer->Delete();
    printf("Results written to %s.\n",ps->outputFile.c_str());
}

int main(int argc, const char * argv[])
{
    try
    {
        // Command line arguments setup
        TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
        TCLAP::ValueArg<std::string> path_tck("",
                                              "fibers",
                                              "Path to the .tck or .vtk file",
                                              true,
                                              "",
                                              "Path to .tck file");
        TCLAP::ValueArg<std::string> path_output("",
                                              "output",
                                              "Path to the output .vtk file",
                                              true,
                                              "",
                                              "Path to output .vtk file");
        TCLAP::ValueArg<float> param_d33("",
                                              "d33",
                                              "Value of D33 parameter. Default value: 1.0",
                                              false,
                                              1.0,
                                              "Value of D33 parameter");
        TCLAP::ValueArg<float> param_d44("",
                                         "d44",
                                         "Value of D44 parameter. Default value: 0.04",
                                         false,
                                         0.04,
                                         "Value of D44 parameter");
        TCLAP::ValueArg<float> param_t("",
                                         "t",
                                         "Value of diffusion time parameter. Default value: 1.4",
                                         false,
                                         1.4,
                                         "Value of diffusion time parameter");
        TCLAP::ValueArg<bool> param_bothdirs("",
                                       "apply-both-dirs",
                                       "Calculate kernel in both forward and backward directions. Default value: false",
                                       false,
                                       false,
                                       "Calculate kernel in both forward and backward directions");
        TCLAP::ValueArg<float> param_subsampling("",
                                             "subsample-step",
                                             "Step to use for fiber subsampling (use 0 to disable). Default value: 1.0",
                                             false,
                                             1.0,
                                             "Step to use for fiber subsampling (use 0 to disable)");
        TCLAP::ValueArg<int> param_anteriorfibers("",
                                                "num-anterior",
                                                "Number of anterior fibers to select (use 0 to disable). Default value: 0",
                                                false,
                                                0,
                                                "Number of anterior fibers to select (use 0 to disable)");
        TCLAP::ValueArg<int> param_maxlength("",
                                                  "max-length",
                                                  "Maximum length of fiber (use 0 to disable). Default value: 0",
                                                  false,
                                                  0,
                                                  "Maximum length of fiber (use 0 to disable)");
        TCLAP::ValueArg<int> param_minlength("",
                                             "min-length",
                                             "Minimum length of fiber (use 0 to disable). Default value: 0",
                                             false,
                                             0,
                                             "Minimum length of fiber (use 0 to disable)");
        TCLAP::ValueArg<float> param_mindist("",
                                                 "min-dist",
                                                 "Kernel distance cutoff. Default value: 15",
                                                 false,
                                                 15,
                                                 "Kernel distance cutoff");
		TCLAP::ValueArg<float> param_nthreads("",
                                                 "nthreads",
                                                 "Number of threads for multicore processing: 1",
                                                 false,
                                                 1,
                                                 "Number of threads");
        TCLAP::SwitchArg param_verbose("","verbose","Verbose mode",false);
        
        cmd.add( param_minlength);
        cmd.add( param_maxlength );
        cmd.add( param_verbose );
        cmd.add( param_mindist );
        cmd.add( param_anteriorfibers);
        cmd.add( param_subsampling);
        cmd.add( param_bothdirs );
        cmd.add( param_t );
        cmd.add( param_d44 );
        cmd.add( param_d33 );
        cmd.add( path_output );
        cmd.add( path_tck );
		cmd.add( param_nthreads );
        
        // Parse the args.
        cmd.parse( argc, argv );
        
        // create parameter settings struct
        ParameterSettings* ps = new ParameterSettings;
        ps->applyKernelInBothDirs = param_bothdirs.getValue();
        ps->numberOfAnteriorFibers = param_anteriorfibers.getValue();
        ps->applySelectAnterior = ps->numberOfAnteriorFibers > 0;
        ps->samplingStep = param_subsampling.getValue();
        ps->applySubsampling = ps->samplingStep > 0;
        ps->D33 = param_d33.getValue();
        ps->D44 = param_d44.getValue();
        ps->t = param_t.getValue();
        ps->verbose = param_verbose.getValue();
        ps->minDist = param_mindist.getValue();
        ps->inputFile = path_tck.getValue();
        ps->outputFile = path_output.getValue();
        ps->cutoff = 3.5;
        ps->maxLength = param_maxlength.getValue();
        ps->minLength = param_minlength.getValue();
        ps->applyFilterLength = ps->maxLength > 0 || ps->minLength > 0;

		// OpenMP settings
#ifdef WITH_OPENMP
 		omp_set_dynamic(0);     // Explicitly disable dynamic teams
    	omp_set_num_threads(param_nthreads.getValue()); // Use 4 threads for all consecutive parallel regions
#endif
        
        // Run spurious fiber filter
        SpuriousFibers(ps);
        
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    
    return 0;
}












