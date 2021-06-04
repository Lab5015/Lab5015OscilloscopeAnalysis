#ifndef WFMREADER_H
#define WFMREADER_H

#include "stdint.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>




class WFMReader
{
public:
  
  //---ctor
  WFMReader(const std::string& binFileName,
            const bool& verbosity = false);
  
  //---dtor
  ~WFMReader() {};
  
  //---methods
  void CloseFile() { if(binFile) binFile->close(); };
  int GetNFrames() { return nFrames; };
  std::pair<std::vector<double>,std::vector<double> > GetFrame(const int& iFrame, const double& ped = 0., const double& polarity = 1.);
  int GetNPoints() { return impDim1Size; };
  double GetVUnit() { return expDim1Scale; };
  double GetTUnit() { return impDim1Scale/1E-09; };
  double GetTMin()  { return impDim1Offset/1.E-09; };
  double GetTMax()  { return ((impDim1Size-1)*impDim1Scale+impDim1Offset)/1.E-09; };
  
  //---members
  std::ifstream* binFile;
  
  std::string binFileName;
  
  char* buf_staticHeader;
  char* buf_waveformHeader;
  char* buf_expDim1;
  char* buf_expDim2;
  char* buf_impDim1;
  char* buf_impDim2;
  char* buf_timeBaseInfo1;
  char* buf_timeBaseInfo2;
  char* buf_updateSpec;
  char* buf_curveInfo;  
  
  std::string version;
  int nBytesLeft;
  int bytesPerPoint;
  int byteOffset;
  int nFrames;
  
  int wfmSetType;
  int wfmCount;
  int reqFFrames;
  int acqFFrames;
  
  double expDim1Scale;
  double expDim1Offset;
  double expDim1Resolution;
  int expDim1Size;
  int expDim1DataSize;
  std::string expDim1Units;
  
  double expDim2Scale;
  double expDim2Offset;
  double expDim2Resolution;
  int expDim2Size;
  int expDim2DataSize;
  std::string expDim2Units;
  
  double impDim1Scale;
  double impDim1Offset;
  double impDim1Resolution;
  int impDim1Size;
  std::string impDim1Units;
  
  double impDim2Scale;
  double impDim2Offset;
  double impDim2Resolution;
  int impDim2Size;
  std::string impDim2Units;
  
  int realPointOffset;
  double TTOffset;
  double fracSec;
  unsigned int gmtSec;
};

#endif



WFMReader::WFMReader(const std::string& binFileName,
                     const bool& verbosity):
binFileName(binFileName)
{
  binFile = new std::ifstream(binFileName.c_str(), std::ios::in | std::ios::binary);
  
  // if the file was successfully open
  if( (*binFile) )
  {
    //-------------
    // static header
    if( verbosity ) std::cout << "------------------ static header ------------------" << std::endl;
    
    buf_staticHeader = new char[78];
    binFile->read(buf_staticHeader,78);
    
    
    char cversion[8];
    for(int ii = 0; ii < 8; ++ii)
      cversion[ii] = buf_staticHeader[2+ii];
    version = std::string(cversion);
    if( verbosity ) std::cout << "version: " << version << std::endl;
    
    nBytesLeft = (*(int*)&buf_staticHeader[11]);
    if( verbosity ) std::cout << "nBytesLeft: " << nBytesLeft << " (" << nBytesLeft/1024./1024. << " MB)" << std::endl;
    
    bytesPerPoint = (int)((*(char*)&buf_staticHeader[15]));
    if( verbosity ) std::cout << "bytesPerPoint: " << bytesPerPoint << std::endl;

    byteOffset = (*(int*)&buf_staticHeader[16]);
    if( verbosity ) std::cout << "byteOffset: " << byteOffset << std::endl;
    
    nFrames = (*(int*)&buf_staticHeader[72]);
    if( verbosity ) std::cout << "nFrames: " << nFrames << std::endl;
    
    
    //----------------
    // waveform header
    if( verbosity ) std::cout << "----------------- waveform header -----------------" << std::endl;
    
    buf_waveformHeader = new char[90];
    binFile->read(buf_waveformHeader,90);
    
    wfmSetType = (*(int*)&buf_waveformHeader[0]);
    if( verbosity ) std::cout << "waveformSetType (0=single, 1=FastFrame): " << wfmSetType << std::endl;
    
    wfmCount = (*(int*)&buf_waveformHeader[4]);
    if( verbosity ) std::cout << "wfmCount = " << wfmCount << std::endl;
    
    reqFFrames = (*(int*)&buf_waveformHeader[68]);
    if( verbosity ) std::cout << "reqFFrames = " << reqFFrames << std::endl;
    
    acqFFrames = (*(int*)&buf_waveformHeader[72]);
    if( verbosity ) std::cout << "acqFFrames = " << acqFFrames << std::endl;
    
    
    //---------------------
    // explicit dimension 1
    if( verbosity ) std::cout << "------------------- exp. dim. 1 -------------------" << std::endl;
    
    buf_expDim1 = new char[160];
    binFile->read(buf_expDim1,160);
    
    expDim1Scale =  (*(double*)&buf_expDim1[0]);
    if( verbosity ) std::cout << "expDim1Scale: "  << expDim1Scale << " ( " << expDim1Scale*256./10.*1000. << " mV/div )" << std::endl;
    
    expDim1Offset =  (*(double*)&buf_expDim1[8]);
    if( verbosity ) std::cout << "expDim1Offset: " << expDim1Offset << std::endl;
    
    expDim1Resolution =  (*(double*)&buf_expDim1[56]);
    if( verbosity ) std::cout << "expDim1Resolution: " << expDim1Resolution << std::endl;
    
    expDim1Size =  (*(int*)&buf_expDim1[16]);
    if( verbosity ) std::cout << "expDim1Size: "   << expDim1Size << std::endl;
    
    expDim1DataSize =  (*(int*)&buf_expDim1[72]);
    if( verbosity ) std::cout << "expDim1DataSize: " << expDim1DataSize << std::endl;
    
    char cexpDim1Units[20];
    for(int ii = 0; ii < 20; ++ii)
      cexpDim1Units[ii] = buf_expDim1[20+ii];
    expDim1Units = std::string(cexpDim1Units);
    if( verbosity ) std::cout << "expDim1Units: \"" << expDim1Units << "\"" << std::endl;
    
    
    //---------------------
    // explicit dimension 2
    if( verbosity ) std::cout << "------------------- exp. dim. 2 -------------------" << std::endl;
    
    buf_expDim2 = new char[160];
    binFile->read(buf_expDim2,160);
    
    expDim2Scale =  (*(double*)&buf_expDim2[0]);
    if( verbosity ) std::cout << "expDim2Scale: "  << expDim2Scale <<  std::endl;
    
    expDim2Offset =  (*(double*)&buf_expDim2[8]);
    if( verbosity ) std::cout << "expDim2Offset: " << expDim2Offset << std::endl;
    
    expDim2Resolution =  (*(double*)&buf_expDim2[56]);
    if( verbosity ) std::cout << "expDim2Resolution: " << expDim2Resolution << std::endl;
    
    expDim2Size =  (*(int*)&buf_expDim2[16]);
    if( verbosity ) std::cout << "expDim2Size: "   << expDim2Size << std::endl;
    
    expDim2DataSize =  (*(int*)&buf_expDim2[72]);
    if( verbosity ) std::cout << "expDim2DataSize: " << expDim2DataSize << std::endl;
    
    char cexpDim2Units[20];
    for(int ii = 0; ii < 20; ++ii)
      cexpDim2Units[ii] = buf_expDim2[20+ii];
    expDim2Units = std::string(cexpDim2Units);
    if( verbosity ) std::cout << "expDim2Units: \"" << expDim2Units << "\"" << std::endl;
    
    
    //---------------------
    // implicit dimension 1
    if( verbosity ) std::cout << "------------------- imp. dim. 1 -------------------" << std::endl;
    
    buf_impDim1 = new char[136];
    binFile->read(buf_impDim1,136);
    
    impDim1Scale =  (*(double*)&buf_impDim1[0]);
    if( verbosity ) std::cout << "impDim1Scale: "  << impDim1Scale <<  std::endl;
    
    impDim1Offset =  (*(double*)&buf_impDim1[8]);
    if( verbosity ) std::cout << "impDim1Offset: " << impDim1Offset << std::endl;
    
    impDim1Resolution =  (*(double*)&buf_impDim1[56]);
    if( verbosity ) std::cout << "impDim1Resolution: " << impDim1Resolution << std::endl;
    
    impDim1Size =  (*(int*)&buf_impDim1[16]);
    if( verbosity ) std::cout << "impDim1Size: "   << impDim1Size << std::endl;     
    
    
    char cimpDim1Units[20];
    for(int ii = 0; ii < 20; ++ii)
      cimpDim1Units[ii] = buf_impDim1[20+ii];
    impDim1Units = std::string(cimpDim1Units);
    if( verbosity ) std::cout << "impDim1Units: \"" << impDim1Units << "\"" << std::endl;
    
    
    //---------------------
    // implicit dimension 2
    if( verbosity ) std::cout << "------------------- imp. dim. 2 -------------------" << std::endl;
    
    buf_impDim2 = new char[136];
    binFile->read(buf_impDim2,136);
    
    impDim2Scale =  (*(double*)&buf_impDim2[0]);
    if( verbosity ) std::cout << "impDim2Scale: "  << impDim2Scale <<  std::endl;
    
    impDim2Offset =  (*(double*)&buf_impDim2[8]);
    if( verbosity ) std::cout << "impDim2Offset: " << impDim2Offset << std::endl;
    
    impDim2Resolution =  (*(double*)&buf_impDim2[56]);
    if( verbosity ) std::cout << "impDim2Resolution: " << impDim2Resolution << std::endl;
    
    impDim2Size =  (*(int*)&buf_impDim2[16]);
    if( verbosity ) std::cout << "impDim2Size: "   << impDim2Size << std::endl;
    
    char cimpDim2Units[20];
    for(int ii = 0; ii < 20; ++ii)
      cimpDim2Units[ii] = buf_impDim2[20+ii];
    impDim2Units = std::string(cimpDim2Units);
    if( verbosity ) std::cout << "impDim2Units: \"" << impDim2Units << "\"" << std::endl;
    
    
    //----------------------
    // time base information
    if( verbosity ) std::cout << "------------------ timeBaseInfo ------------------" << std::endl;
    
    buf_timeBaseInfo1 = new char[12];
    binFile->read(buf_timeBaseInfo1,12);
    
    buf_timeBaseInfo2 = new char[12];
    binFile->read(buf_timeBaseInfo2,12);
    
    /*
    //---------------------
    // update specification
    if( verbosity ) std::cout << "------------------- updateSpec ------------------" << std::endl;
    
    buf_updateSpec = new char[24];
    binFile->read(buf_updateSpec,24);
    
    realPointOffset = (*(int*)&buf_updateSpec[0]);
    if( verbosity ) std::cout << "realPointOffset: " << realPointOffset << std::endl;
    
    TTOffset = (*(double*)&buf_updateSpec[4]);
    if( verbosity ) std::cout << "TTOffset: " << TTOffset << std::endl;
    
    fracSec = (*(double*)&buf_updateSpec[12]);
    if( verbosity ) std::cout << "fracSec: " << fracSec << std::endl;
    
    gmtSec =  (*(unsigned int*)&buf_updateSpec[20]);
    if( verbosity ) std::cout << "gmtSec: " << gmtSec << std::endl;
    
    
    //------------------
    // curve information
    if( verbosity ) std::cout << "-------------------- curveInfo ------------------" << std::endl;
    
    buf_curveInfo = new char[30];
    binFile->read(buf_curveInfo,30);
    */
    
    if( verbosity ) std::cout << "---------------------- end  ----------------------" << std::endl;
    binFile->close();
  }
  else
  {
    std::cerr << "WFMReader::Error: failed to open binary file " << binFileName << std::endl;
  }
}




std::pair<std::vector<double>,std::vector<double> >  WFMReader::GetFrame(const int& iFrame, const double& ped, const double& polarity)
{
  std::pair<std::vector<double>,std::vector<double> > ret;
  
  
  if( (iFrame > nFrames) || iFrame < 0 )
  {
    std::cerr << "WFMReader::GetFrame::Error: frame " << iFrame << " does not exist" << std::endl;  
    return ret;
  }
  
  /*
  int bytesToSkip1 = 838 + iFrame*(54);
  if((iFrame-1)%4 == 0) bytesToSkip1 += 1;
  
  if( !binFile || !binFile->is_open() )
  {
    binFile = new std::ifstream(binFileName.c_str(), std::ios::in | std::ios::binary);
    binFile -> ignore(bytesToSkip1);
  }
  
  char* buf_updateSpec = new char[24];
  binFile->read(buf_updateSpec,24);
  
  realPointOffset = (*(int*)&buf_updateSpec[0]);
  std::cout << "realPointOffset: " << realPointOffset << std::endl;
  
  TTOffset = (*(double*)&buf_updateSpec[4]);
  std::cout << "TTOffset: " << TTOffset << std::endl;

  fracSec = (*(double*)&buf_updateSpec[12]);
  std::cout << "fracSec: " << fracSec << std::endl;
  
  gmtSec =  (*(unsigned int*)&buf_updateSpec[20]);
  std::cout << "gmtSec: " << gmtSec << std::endl;
  
  
  std::bitset<192> b;

  for(int i = 0; i < 24; ++i)
  {
    uint8_t cur = (*(uint8_t*)&buf_updateSpec[i]);
    int offset = i * 8;

    for(int bit = 0; bit < 8; ++bit)
    {
      b[offset] = cur & 1;
      ++offset;   // Move to next bit in b
      cur >>= 1;  // Move to next bit in array
    }
  }

  std::cout << b << std::endl;
  
  binFile -> close();
  */
  
  int bytesToSkip = byteOffset + impDim1Size*bytesPerPoint*iFrame;
  
  if( !binFile || !binFile->is_open() )
  {
    binFile = new std::ifstream(binFileName.c_str(), std::ios::in | std::ios::binary);
    binFile -> ignore(bytesToSkip);
  }
  
  char buf_curve[impDim1Size*bytesPerPoint];
  binFile->read(buf_curve,sizeof(buf_curve));
  
  int8_t val_curve;
  for(int ii = 0; ii < impDim1Size; ++ii)
  {
    val_curve = (*(int8_t*)&buf_curve[ii*bytesPerPoint]);
    ret.first.push_back( (ii*impDim1Scale+impDim1Offset)/1.E-09 );
    ret.second.push_back( (double(val_curve)*expDim1Scale+expDim1Offset-ped)*polarity );
    //std::cout << "ii: " << ii << "   t: " <<  (ii*impDim1Scale+impDim1Offset) << "   V: " << int(val_curve) << std::endl;
  }
  
  if( iFrame == (nFrames-1) )
  {
    binFile -> close();
  }
  
  return ret;
}
