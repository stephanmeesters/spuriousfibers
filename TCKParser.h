//
//  TCKParser.h
//  SpuriousFibers
//
//  Created by Stephan Meesters on 02/06/15.
//
//

#ifndef __SpuriousFibers__TCKParser__
#define __SpuriousFibers__TCKParser__

/** Includes - VTK */
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkAlgorithm.h>

/** Includes - C++ */
#include <iostream>
#include <fstream>
#include <assert.h>

/** Data structs **/

struct StatHeader
{
    int i1;             /* is this statistic good as a luminance encoding for paths? (not used) */
    int i2;             /* is this statistic stored per point? */
    int i3;             /* does it show up in the stat panel by default? (not used) */
    char c1[255];       /* aggregate name */
    char c2[255];       /* local name */
    int i4;             /* unique ID */
};

struct AlgoHeader
{
    char c1[255];
    char c2[255];
    int i1;
};

struct Pathway
{
    int headerSize;
    int numPoints;
    int algoInt;
    int seedPointIndex;
    double* pathStats;
    double* points;         // {{x_1,y_1,z_1},{x_2,y_2,z_2},...,{x_n,y_n,z_n}} -- n = number of points
    double* pointStats;     // {{s1_1,s1_2,...s1_n},{s2_1,s2_2,...s2_n},...,{sm_1,sm_2,...sm_n}} -- m = number of point stats
};

class TCKParser
{
public:
    static vtkPolyData* LoadDataFromFile(std::string filename);
};


#endif /* defined(__SpuriousFibers__TCKParser__) */
