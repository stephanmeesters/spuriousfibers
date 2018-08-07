#ifndef bmia_SpuriousFiberFilterTypes_h
#define bmia_SpuriousFiberFilterTypes_h

/** Holding parameter settings */
typedef struct
{
    std::string outputFiberDataName;

    double D33;
    double D44;
    double t;
    bool applyKernelInBothDirs;
    double minDist;
    bool requireRecompute;
    double cutoff;
    int maxLength;
    int minLength;

    // saved results
    double avgScoreTotal;
    double* fiberMinScores;
    double* fiberScores;
    int* fiberStartIds;

    bool applySubsampling;
    double samplingStep;

    bool applySelectAnterior;
    int numberOfAnteriorFibers;
    
    bool applyFilterLength;
    
    bool verbose;
    
    std::string inputFile;
    std::string outputFile;

} ParameterSettings;

#endif  // bmia_SpuriousFiberFilterTypes_h
