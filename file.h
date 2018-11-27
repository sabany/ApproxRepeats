/*
 *  Purpose: Provides classes for interacting with files.
 *
 *  Modified by Sarah Banyassady, August 2015
 *
 *   Added classes and functions:
 *   - seedData
 *   - tableEntry
 *   - seedFileReadInfo
 *   - compare_query
 *   - getStartnEndOfSequence
 *   - removeSubsets
 *
 *   Modified classes and functions:
 *   - memExt
 *   - printMemonTerminal
 *   - removeDuplicates
 *
 * ================================================================= *
 *  file.h : Header file with supporting class definitions           *
 *                                                                   *
 *  E-MEM: An efficient (MUMmer-like) tool to retrieve Maximum Exact *
 *         Matches using hashing based algorithm                     *
 *                                                                   *
 *  Copyright (c) 2014, Nilesh Khiste                                *
 *  All rights reserved                                              *
 *                                                                   * 
 *  This program is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public License as   *
 *  published by the Free Software Foundation, either version 3 of   *
 *  the License, or (at your option) any later version.              *
 *                                                                   *
 *  This program is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 *  GNU General Public License for more details.                     *
 *                                                                   *
 *  You should have received a copy of the GNU General Public        *
 *  License along with this program.                                 *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#ifndef __FILE_H_INCLUDED__
#define __FILE_H_INCLUDED__

#include <unistd.h>

using namespace std;
typedef uint32_t intType;       //type of integer for sequence positions

#define CHARS2BITS(x)       (2*(x))     //to convert char position to bit position
#define DATATYPE_WIDTH      64      //number of bits
#define NUM_TMP_FILES       100     //number of temporary files
#define NUM_TMP_FILES_REV   50
#define NEGINFINITY         -2000000000

class commonData {
  public:
    static int32_t numThreads;
    static int32_t minMemLen;
    static int32_t maxDistance;
    static int32_t matchScore;
    static int32_t misMatScore;
    static int32_t indelScore;
    static int32_t xDrop;
    static int32_t minScore;
    static int32_t stepSize;
    static int32_t d;
    static int32_t ignoreN;
    static int32_t fColOutput;
    static int32_t lenInHeader;
    static int32_t relQueryPos;
    static double  kBlast;
    static double  lambdaBlast;
    static double  expectationValue;
};

int32_t commonData::numThreads  = 1;
int32_t commonData::minMemLen   = 0;
int32_t commonData::maxDistance = 5;
int32_t commonData::matchScore  = 2;
int32_t commonData::misMatScore= -2;
int32_t commonData::indelScore = -3;
int32_t commonData::xDrop       = 5;
int32_t commonData::minScore = 0;
int32_t commonData::stepSize = 1;
int32_t commonData::d        = 1;
int32_t commonData::ignoreN  = 0;
int32_t commonData::fColOutput  = 0;
int32_t commonData::lenInHeader = 0;
int32_t commonData::relQueryPos = 0;
double  commonData::kBlast       = -1.0;
double  commonData::lambdaBlast  = -1.0;
double  commonData::expectationValue = 10.00;

class seedData {
public:
    uint32_t length;
    uint32_t weight;        //number of 1's in seed
    uint32_t binSeedSize;
    uint64_t * binSeed;
    
    seedData() {
        weight=0;
        length=0;
        binSeedSize=0;
        binSeed=NULL;
    }
    
    ~seedData() {
        if(binSeed)
        delete [] binSeed;
    }
    
    void allocBinArray(uint32_t length) {
        binSeedSize = ceil(static_cast<float>(length)/32);
        try{
            binSeed=new uint64_t[binSeedSize];
        }catch(std::bad_alloc& ba){
            std::cerr << "seedData _ bad_alloc caught: " << ba.what() <<endl;
        }
    }
    
    void swap(seedData * seed) {
        
        uint32_t tempW, tempL, tempSize;
        uint64_t * tempBin;
        
        tempW = this->weight;
        this->weight = seed->weight;
        seed->weight = tempW;
        
        tempL = this->length;
        this->length = seed->length;
        seed->length = tempL;
        
        tempSize = this->binSeedSize;
        this->binSeedSize = seed->binSeedSize;
        seed->binSeedSize = tempSize;
        
        tempBin = this->binSeed;
        this->binSeed = seed->binSeed;
        seed->binSeed = tempBin;
    }
};

/* Entry of dynammic table which is used for extension */
class tableEntry {
public:
    int32_t score;
    int32_t editDist;
    
    tableEntry()
    {
        score = NEGINFINITY;
        editDist = 0;
    };
    
    tableEntry(int32_t s, int32_t d)
    {
        score = s;
        editDist = d;
    };
};

class seedFileReadInfo {
    fstream file;
    string filename;
    uint32_t numSeeds;
    
    void sortSeedLength() {
        
        seedData tempSeed;
        uint32_t j;
        for(uint32_t i=1; i<numSeeds ; i++) {
            j = i;
            while (j>0 && seeds[j-1].length > seeds[j].length) {
                seeds[j-1].swap(&seeds[j]);
                j--;
            }
        }
    }
    
    void clearFileFlag() {
        file.clear();
        file.seekg(0, ios::beg);
    }
    
public:
    seedData * seeds;
    uint32_t minLen;    //length of the shortest seed
    uint32_t maxLen;    //length of the longest seed
    
    seedFileReadInfo() {
        numSeeds=0;
        seeds=NULL;
        minLen=NEGINFINITY;
        maxLen=0;
    }
    
    void destroy() {
        numSeeds=0;
        if(seeds)
        delete [] seeds;
    }
    
    uint32_t getNumSeeds() {
        return numSeeds;
    }
    
    void openFile(string s) {
        file.open(s, ios::in);
        if(!file.is_open()) {
            cout << "ERROR: unable to open "<< s << " file" << endl;
            exit(EXIT_FAILURE);
        }
        filename=s;
    }
    
    void setNumSeeds()  {
        string line;
        while ( getline(file,line).good() )
            numSeeds++;
        
        if(numSeeds == 0){
            cout << "ERROR: "<< filename << " file is empty" << endl;
            exit(EXIT_SUCCESS);
        }
        try{
            seeds = new seedData[numSeeds];
        }catch(std::bad_alloc& ba){
            std::cerr << "seedFileReadInfo _ bad_alloc caught: " << ba.what() <<endl;
        }
        
        clearFileFlag();
        for (uint32_t i=0; getline(file, line).good() ;i++)
            seeds[i].allocBinArray((uint32_t)line.length());
    }
    
    void readSeeds()  {
        clearFileFlag();
        string line;
        uint32_t binSeedPos;
        
        for (uint32_t i=0; getline(file, line).good() ;i++)
        {
            binSeedPos=0;
            seeds[i].binSeed[binSeedPos] = 0;
            for( std::string::iterator it=line.begin(); it!=line.end(); ++it)
            {
                switch(*it)
                {
                    case '1':
                        seeds[i].binSeed[binSeedPos] <<= 2;
                        seeds[i].binSeed[binSeedPos] |= 3;
                        seeds[i].weight++;
                        break;
                    case '*':
                        seeds[i].binSeed[binSeedPos] <<= 2;
                }
                seeds[i].length++;
                
                if ((seeds[i].length%32)==0){
                    binSeedPos++;
                    seeds[i].binSeed[binSeedPos] = 0;
                }
            }
            
            if(seeds[i].weight > 32){
                cout << "ERROR: seed weights cannot be greater than 32" << endl;
                exit(EXIT_SUCCESS);
            }
            
            if(minLen > seeds[i].length)  minLen=seeds[i].length;
            if(maxLen < seeds[i].length)  maxLen=seeds[i].length;
        }
        file.close();
        sortSeedLength();
    }

};


class seqData {
  public:
    uint64_t start;
    uint64_t end;
    char seq[32];
    seqData() 
    {
        start=0;
        end=0;
        memset(&seq, 0, 32);
    };

    bool operator ()(const seqData &obj1, const seqData &obj2)
    {
      return (obj2.start>obj1.end?true:false);
    }
};

class mapObject {
    public:
      uint64_t left;
      uint64_t right;
      mapObject() {
          left=0;
          right=0;
      }

      mapObject(uint64_t x, uint64_t y) {
          left=x;
          right=y;
      }

      bool operator()(const uint64_t &x, const mapObject &y) 
      {
          return x < y.left;
      }
};

class seqFileReadInfo {
      fstream file;
      string filename;
      uint64_t size;     //size of the sequence/subsequence being processed
      string strTmp, strName;
      uint64_t binReadSize;
      uint64_t binReadsLocation;
      uint64_t currPos;
      uint64_t numSequences;
   
      void processTmpString(uint64_t &sz, uint64_t &blockNCount)
      {
          string line = strTmp;
          strTmp.clear();
          totalBases=0;
          binReadsLocation=0;
          processInput(line, sz, blockNCount);
      }

      /*
       * Converts a character sequence into an array of integers.
       * Input: character string
       * Output: array of integers, total number of bases
       */
      void processInput(string &str, uint64_t &sz, uint64_t &blockNCount)
      {
          int chooseLetter=0;
          uint64_t k=0;

          if (!totalBases) {
              for (k=0; k<binReadSize; ++k)
                 binReads[k]=0;
          }

          /* Processing the sequences by encoding the base pairs into 2 bits */
          for ( std::string::iterator it=str.begin(); it!=str.end(); ++it)
          {
              if (totalBases == sz){ //sz=size+minSize
                  strTmp += *it;
                  continue;
              }else if (totalBases >= size) {  //to make borders overlap in case commonData::d != 1 or commonData::numThreads != 1
                  strTmp += *it;
              }
              switch(*it)
              {
                  case 'A':
                  case 'a':
                      binReads[binReadsLocation] <<= 2;
                      if (commonData::ignoreN && blockNCount){
                         blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                         blockNCount=0;
                      }
                      numOfAs++;
                      break;
                  case 'C':
                  case 'c':
                      binReads[binReadsLocation] <<= 2;
                      binReads[binReadsLocation] |= 1;
                      if (commonData::ignoreN && blockNCount){
                         blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                         blockNCount=0;
                      }
                      numOfCs++;
                      break;
                  case 'G':
                  case 'g':
                      binReads[binReadsLocation] <<= 2;
                      binReads[binReadsLocation] |= 2;
                      if (commonData::ignoreN && blockNCount){
                         blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                         blockNCount=0;
                      }
                      numOfGs++;
                      break;
                  case 'T':
                  case 't':
                      binReads[binReadsLocation] <<= 2;
                      binReads[binReadsLocation] |= 3;
                      if (commonData::ignoreN && blockNCount){
                         blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                         blockNCount=0;
                      }
                      numOfTs++;
                      break;
                  default:
                      if(!blockNCount)
                          blockNCount=totalBases+1;
                      chooseLetter = rand() % 4;
                      if (chooseLetter == 0)
                          binReads[binReadsLocation] <<= 2;
                      else if (chooseLetter == 1)
                      {
                          binReads[binReadsLocation] <<= 2;
                          binReads[binReadsLocation] |= 1;
                      }
                      else if (chooseLetter == 2)
                      {
                          binReads[binReadsLocation] <<= 2;
                          binReads[binReadsLocation] |= 2;
                      }
                      else
                      {
                         binReads[binReadsLocation] <<= 2;
                         binReads[binReadsLocation] |= 3;
                      }
              }
              totalBases++;
              if ((totalBases%32)==0){
                  binReadsLocation++;
              }
          }
      }

    public:
      uint64_t *binReads;   //an array of integers containing binay representation of the sequence chunks
      intType totalBases;   //total number of bases in the sequence
      intType numOfAs;
      intType numOfCs;
      intType numOfGs;
      intType numOfTs;
      std::vector <mapObject> blockOfNs;
      
      seqFileReadInfo() {
          size=0;
          currPos=0;
          binReadSize=0;
          binReadsLocation=0;
          numSequences=0;
          totalBases=0;
          numOfAs=0;
          numOfCs=0;
          numOfGs=0;
          numOfTs=0;
          binReads=NULL;
      }

      uint64_t &getNumSequences() {
          return numSequences;
      }
    
      string getFileName() {
          return filename;
      }

      void openFile(string s) {
          file.open(s, ios::in);
          if(!file.is_open()) {
              cout << "ERROR: unable to open "<< s << " file" << endl;
              exit(EXIT_FAILURE);
          }
          filename = s;
      }

      void setReverseFile() {
          char buffer[50];
          memset(buffer,0,50);
          sprintf(buffer, "./%d_tmp/revComp", getpid());
          file.close();
          openFile(buffer);
      }
    
      void destroy() {
          currPos=0;
          binReadSize=0;
          binReadsLocation=0;
          totalBases=0;
          numOfAs=0;
          numOfCs=0;
          numOfGs=0;
          numOfTs=0;
          strName.clear();
          strTmp.clear();
          clearMapForNs();
          if(binReads)
          delete [] binReads;
      }

      void clearFileFlag() {
          file.clear();
          file.seekg(0, ios::beg);
      } 
    
      uint64_t allocBinArray()
      {
          size = size/commonData::d;
          binReadSize = floor((size+commonData::minMemLen/2+commonData::d)/32+4);
          try{
              binReads = new uint64_t[binReadSize];
          }catch(std::bad_alloc& ba){
              std::cerr << "seqFileReadInfo _ bad_alloc caught: " << ba.what() <<endl;
          }
          return size;
      }
  
      void clearMapForNs() {
          blockOfNs.clear();
      } 

      void clearTmpString()
      {
          strTmp.clear();
          strName.clear();
          clearMapForNs();
      }
      
      void getKmerLeftnRightBoundForNs(intType currKmerPos, mapObject &bounds)
      {
          uint64_t right=0;
          /*
           * Since we do all computation with bits, all our
           * positions are even. Here I return 1 (odd position),
           * an indication of no Ns towards left   
           */

          if (!blockOfNs.size()){
              bounds.left=0x1;
              bounds.right=CHARS2BITS(totalBases-1);
              return;
          }

          vector<mapObject>::iterator it;
          it=upper_bound(blockOfNs.begin(), blockOfNs.end(), currKmerPos, mapObject()); 
          /* No N block beyond this point */
          if (it == blockOfNs.end())
              right = CHARS2BITS(totalBases-1);
          else
              right = (*it).left-2;

          /* This function never gets a position which is N */
          if (!currKmerPos || it==blockOfNs.begin()){
              bounds.left=0x1;
              bounds.right=right;
              return;
          }

          --it;

          bounds.left=(*it).right+2;
          bounds.right=right;
          return;
      }

    
      bool checkKmerForNs(intType currKmerPos, vector<mapObject>::iterator &it, uint32_t seedLength, uint32_t & double_check)
      {
          if (!blockOfNs.size())
              return false;

          while(it != blockOfNs.end())
          {
              if ((*it).left>currKmerPos)
                  break;
              else
                  ++it;
          }    
    
          /* No N block beyond this point
           * end() returns an iterator referring to the past-the-end element,
           * the theoretical element that would follow the last element
           */
          if (it == blockOfNs.end()){
              --it;
              /* Current position within N block */
              if (((*it).left <=currKmerPos) && (currKmerPos <= (*it).right)){
                  ++it;
                  return true;
              }else{
                  ++it;
                  return false;
              }
          }
 
          if ((*it).left > (currKmerPos+seedLength-2)){
              if (it != blockOfNs.begin()){
                  --it;
                  if ((*it).right < currKmerPos){
                      ++it;
                      double_check = 1; //the same currKmerPos with a longer seed needs to be checked again
                      return false;
                  }else {
                      ++it;
                      return true;
                  }
              }else {
                  double_check = 1;
                  return false;
              }
          }else {
              return true;
          }
      }


      void setCurrPos() {
          currPos+=size;
      } 

      uint64_t getCurrPos() {
          return currPos;
      } 

      void resetCurrPos() {
          currPos=0;
      } 

      bool readChunks()
      {
          string line;
          uint64_t blockNCount=0;
          int minSize = commonData::minMemLen/2-1;
          uint64_t sz=size+minSize;
          /* Process anything remaining from the last iteration */
          processTmpString(sz, blockNCount);
          
          while(getline( file, line ).good() ){
              if(line[0] == '>' || (totalBases == sz)){
                  if( !strName.empty()){    //process what we read from the last entry
                     if(line[0] != '>') {
                          processInput(line, sz, blockNCount);
                      }
                      if (totalBases == sz) {
                          if ((totalBases%32)!=0){
                              intType offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                              binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                              binReadsLocation++;
                              binReads[binReadsLocation]=0;
                          }
                          if (commonData::ignoreN && blockNCount){
                              blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                              blockNCount=0;
                          }
                          return true;
                      }
                  }
                  if( !line.empty() ){
                      strName = line.substr(1);
                  }
              }else if( !strName.empty() ){
                  processInput(line, sz, blockNCount);
              }
          }

          if( !strName.empty() ){   //process what we read from the last entry
              if ((totalBases%32)!=0){
                  intType offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                  binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                  binReadsLocation++;
                  binReads[binReadsLocation]=0;
              }
              if (commonData::ignoreN && blockNCount){
                  blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                  blockNCount=0;
              }
              if (!strTmp.size())
                  strName.clear();
              return true;
          }
          return false;
      }

      void flipCharacter(char &in, char &out)
      {
          switch(in)
          {
              case 'A':
              case 'a':
                  out='T';
                  break;
              case 'C':
              case 'c':
                  out='G';
                  break;
              case 'G':
              case 'g':
                  out='C';
                  break;
              case 'T':
              case 't':
                  out='A';
                  break;
              default:
                  out=in;
          }
      }

      void flipNswap(string &content)
      {
          string::iterator itBeg = content.begin();
          string::iterator itEnd = --content.end();
          char beg=0, end=0;
          uint64_t d=0;
          while ((d=distance(itBeg,itEnd))) 
          {
              flipCharacter(*itBeg, end);
              flipCharacter(*itEnd, beg);
              (*itEnd)=end;
              (*itBeg)=beg;
              ++itBeg;
              --itEnd;
              if(d==1)
                  break;
          }
          if (!d)
              flipCharacter(*itEnd, *itEnd);
      }

      void writeReverseComplementString(string &name, string &content, fstream &file) {
          file << ">" << name << endl;
          flipNswap(content);
          file << content ;
      }

      /* Writes the sequence/pos mapping of the reference and query in
       * the first two TmpFiles. It also sets numSequences and size 
       */
      bool generateSeqPosMap(vector<seqData> &seqInfo, uint32_t revComplement) {
          uint64_t i=0,j=0;
          seqData s;
          string line, name,content;
          fstream revFile;

          if (revComplement) {
              char buffer[50];
              memset(buffer,0,50);
              sprintf(buffer, "./%d_tmp/revComp", getpid());
              revFile.open(buffer, ios::out);
              if (!revFile.is_open())
              {
                  cout << "ERROR: unable to open temporary reverse complement file" << endl;
                  return false;
              }
          }

          clearFileFlag();
          while(getline(file, line).good()){

              if(line[0] == '>'){
                  if(!name.empty()) {
                      s.start=CHARS2BITS(j);
                      s.end=CHARS2BITS(i-1);
                      strncpy(s.seq,name.c_str(),30);
                      s.seq[31]='\0';
                      numSequences++;
                      seqInfo.push_back(s);
                      j=i;
                      if (revComplement){
                          writeReverseComplementString(name, content, revFile);
                          content.clear();
                      }
                      name.clear();
                  }
                  if(!line.empty())
                      name=line.substr(1);
              } else if( !name.empty() ) {
                    i+=line.length();
                    size+=line.length();
                    if (revComplement) {
                        content += "\n";
                        content += line;
                    }
              }
          }
          if(size == 0) {
              cout << "ERROR: "<< filename << " file is empty" << endl;
              return false;
          }
          
          if( !name.empty() ) {
              s.start=CHARS2BITS(j);
              s.end=CHARS2BITS(i-1);
              strncpy(s.seq,name.c_str(),30);
              s.seq[31]='\0';
              numSequences++;
              seqInfo.push_back(s);
              if (revComplement){
                  writeReverseComplementString(name, content, revFile);
                  content.clear();
              }
              name.clear();
          }
          if (revComplement)
              revFile.close();
          return true;
      }
};

class MemExt {
public:
    intType lR;
    intType lQ;
    uint16_t Rsize;
    uint16_t Qsize;
    int32_t score;
    MemExt() {
    }
    
    ~MemExt() {
    }
    
    MemExt(intType lr, intType lq, uint16_t rs, uint16_t qs, int32_t s){
        lR=lr;
        lQ=lq;
        Rsize=rs;
        Qsize=qs;
        score=s;
    }
    
    void copy(MemExt source){
        this->lR = source.lR;
        this->Rsize = source.Rsize;
        this->lQ = source.lQ;
        this->Qsize = source.Qsize;
        this->score=source.score;
    }
    
    bool operator () (const MemExt &obj1, const MemExt &obj2)
    {
        if (obj1.lQ<obj2.lQ)
            return true;
        else if (obj1.lQ>obj2.lQ)
            return false;
        else{
            if (obj1.lR<obj2.lR)
                return true;
            else if (obj1.lR>obj2.lR)
                return false;
            else{     //obj1.lQ = obj2.lQ & obj1.lR = obj2.lR
                if (obj1.Qsize>obj2.Qsize)
                    return true;
                else if (obj1.Qsize<obj2.Qsize)
                    return false;
                else{
                    if(obj1.Rsize>obj2.Rsize)
                        return true;
                    return false;
                }
            }
        }
    }
};


class tmpFilesInfo {
    fstream *TmpFiles;
    uint64_t numMemsInFile;
    uint64_t numRevMemsInFile;
    vector <MemExt> MemInBorders;
    vector<seqData> refSeqInfo;
    vector<seqData> querySeqInfo;

    /* Check if the match straddles borders of the sequences returns true (in case commonData::d != 1 or commonData::numThreads != 1) */
    bool checkMEMExt(intType &lr, intType &rr, intType &lq, intType &rq, seqFileReadInfo &QueryFile, seqFileReadInfo &RefFile) {
      if ((!lq && QueryFile.getCurrPos()) || rq == CHARS2BITS(QueryFile.totalBases-1)) {
         return true;
      }else if((!lr && RefFile.getCurrPos()) || rr == CHARS2BITS(RefFile.totalBases-1)) {
         return true;
      } 
      return false;
    }
    
    void writeToFile(intType lQ, intType rQ, intType lR, intType rR, int32_t s, uint32_t &revComplement) {
        MemExt m;
        m.lQ=lQ;
        m.lR=lR;
        m.Qsize=static_cast<uint16_t>(rQ-lQ+2);
        m.Rsize=static_cast<uint16_t>(rR-lR+2);
        m.score=s;
        if (IS_MATCH_BOTH_DEF(revComplement))
            TmpFiles[(m.lQ/numRevMemsInFile)+NUM_TMP_FILES].write((char *)&m, sizeof(MemExt));
        else
            TmpFiles[m.lQ/numMemsInFile].write((char *)&m, sizeof(MemExt));

    }
    
    void writeToVector(intType lQ, intType rQ, intType lR, intType rR, int32_t s) {
        MemExt m;
        m.lQ=lQ;
        m.lR=lR;
        m.Qsize=static_cast<uint16_t>(rQ-lQ+2);
        m.Rsize=static_cast<uint16_t>(rR-lR+2);
        m.score=s;
        MemInBorders.push_back(m);
    }

  public:
 
    tmpFilesInfo(int numFiles) {
        TmpFiles = new fstream[numFiles];
    }

    ~tmpFilesInfo() {
        delete [] TmpFiles;
    }
   
    void setNumMemsInFileVec(uint64_t size) {
        numMemsInFile =  (2*(size+commonData::minMemLen/2)*commonData::d+commonData::d)/NUM_TMP_FILES;
        numRevMemsInFile = (2*(size+commonData::minMemLen/2)*commonData::d+commonData::d)/NUM_TMP_FILES_REV;
    }
   
    static bool compare_reference (const MemExt &obj1, const MemExt &obj2)
    {
      return (obj1.lR>=obj2.lR?false:true);
    }
    
    static bool compare_query (const MemExt &obj1, const MemExt &obj2)
    {
        if (obj1.lQ<obj2.lQ)
            return true;
        else if (obj1.lQ>obj2.lQ)
            return false;
        else{
            if (obj1.lR<obj2.lR)
                return true;
            else
                return false;
        }
    }

    static bool myUnique(const MemExt &obj1, const MemExt &obj2)
    {
        if((obj1.lQ<=obj2.lQ) && (obj1.lR<=obj2.lR) && ((obj1.lQ+obj1.Qsize) >= (obj2.lQ+obj2.Qsize)) && ((obj1.lR+obj1.Rsize) >= (obj2.lR+obj2.Rsize)))
            return true;
        else
            return false;
    }

    void openFiles(ios_base::openmode mode, int32_t numFiles, uint32_t &revComplement) {
        char buffer[50];
        memset(buffer,0,50);
        static int flag=0;
        sprintf(buffer, "./%d_tmp", getpid());
        if (!flag) {
            if(mkdir(buffer, S_IRWXU|S_IRGRP|S_IXGRP))
            {
                cout << "ERROR: unable to open temporary directory" << endl;
                exit(EXIT_FAILURE);
            }
            flag=1;
        }
 
        for (int32_t i=0; i<numFiles; i++) {
            /* Temporary file to hold the mems */
            sprintf(buffer, "./%d_tmp/%d", getpid(),i);
            TmpFiles[i].open(buffer, mode);
            if (!TmpFiles[i].is_open())
            {
                cout << "ERROR: unable to open temporary file" << endl;
                EXIT_AND_REMOVE(i, revComplement)
            }
        }
    }
    
    void closeFiles(int32_t numFiles) {
        for (int32_t i=0; i<numFiles; i++)
           TmpFiles[i].close();
    }
    
    vector<seqData>& getRefInfoVec() {
        return refSeqInfo;
    }
    
    vector<seqData>& getQueryInfoVec() {
        return querySeqInfo;
    }
    
    void writeMemInTmpFilesVecs(intType &lRef, intType &rRef, intType &lQue, intType &rQue, int32_t s, seqFileReadInfo &QueryFile, seqFileReadInfo &RefFile, uint32_t &revComplement) {
       MemExt m;
       intType currPosQ = CHARS2BITS(QueryFile.getCurrPos());
       intType currPosR = CHARS2BITS(RefFile.getCurrPos());
        
       if (!(commonData::d==1 && commonData::numThreads==1) && checkMEMExt(lRef, rRef, lQue, rQue, QueryFile, RefFile)){
           #pragma omp critical
           writeToVector(currPosQ+lQue, currPosQ+rQue, currPosR+lRef, currPosR+rRef, s);
       }
       else{
            #pragma omp critical
            writeToFile(currPosQ+lQue, currPosQ+rQue, currPosR+lRef, currPosR+rRef, s, revComplement);
       }
    }

    void printQueryHeader(vector<seqData>::iterator &itQ, uint32_t &revComplement)
    {

        if (revComplement & 0x1){
            if (commonData::lenInHeader)
                cout << "> " << strtok((*itQ).seq," ") << " Reverse" << " Len = " << ((*itQ).end-(*itQ).start+2)/2 << endl;
            else
                cout << "> " << strtok((*itQ).seq," ") << " Reverse" << endl;
        }else{
            if (commonData::lenInHeader)
                cout << "> " << strtok((*itQ).seq," ") << " Len = " << ((*itQ).end-(*itQ).start+2)/2 << endl;
            else
                cout << "> " << strtok((*itQ).seq," ") << endl;
        }
    }
    
    
    void getStartnEndOfSequence(intType refPos, intType quePos, seqData &refSeq, seqData &queSeq)
    {
        if (refSeqInfo.size() == 0 || querySeqInfo.size() == 0)
            return;
        
        static vector<seqData>::iterator itQ = querySeqInfo.begin();
        while(quePos >= (*itQ).end)
            ++itQ;
        queSeq = *itQ;
    
        vector<seqData>::iterator itR;
        seqData s;
        s.start=refPos;
        itR = lower_bound(refSeqInfo.begin(), refSeqInfo.end(), s, seqData());
        refSeq = (*itR);
        return;
    }

    void printMemOnTerminal(MemExt m, uint32_t &revComplement) {
        intType &lRef = m.lR;
        intType rRef = (m.lR + m.Rsize - 2);
        intType &lQue = m.lQ;
        intType rQue = (m.lQ + m.Qsize - 2);
        static int flag=0;
        vector<seqData>::iterator itR;
        static vector<seqData>::iterator itQ=querySeqInfo.begin();
        seqData s;
        
        /* Print the remianing query sequences - if any */
        if (!lRef && !lQue && !rRef && !rQue) {
            /* No matches found - simply return */
            if (!flag){
                printQueryHeader(itQ, revComplement);
            }
            while(itQ != --querySeqInfo.end()){
                ++itQ;
                printQueryHeader(itQ, revComplement);
            }
            itQ=querySeqInfo.begin();
            flag=0;
            return;
        }
        
        s.start=lRef;
        s.end=rRef;

        /* Process relative position for reference sequence */
        itR = lower_bound(refSeqInfo.begin(), refSeqInfo.end(), s, seqData());
        
        if ((*itR).start <= lRef && (*itR).end >= rRef){
            // MEM within acutal sequence
            // s------e--s------e--s------e
            // s--|--|e
            lRef?lRef-=((*itR).start):lRef;
            rRef-=((*itR).start);
        }
       
        /* Print the first query sequence */
        if (!flag){
            printQueryHeader(itQ, revComplement);
            flag=1;
        }
        /* Process relative position for query sequence */
        while(lQue >= (*itQ).end){
            ++itQ;
            printQueryHeader(itQ, revComplement);
        }
        if ((*itQ).start <= lQue && (*itQ).end >= rQue){
            // MEM within acutal sequence
            // s------e--s------e--s------e
            // s--|--|e
            lQue?lQue-=((*itQ).start):lQue;
            rQue-=((*itQ).start);
        }

        if (refSeqInfo.size() == 1 && !commonData::fColOutput) {
            if ((revComplement & 0x1) && commonData::relQueryPos)
                    
                cout << setw(7) << " " << "(" << std::left << setw(10) << (lRef+2)/2 << ", " << std::left << setw(10) << (rRef+2)/2 << ")   "
                << " (" << std::left << setw(10) << ((*itQ).end-(*itQ).start-lQue+2)/2 << ", " << std::left << setw(10)
                << ((*itQ).end-(*itQ).start-rQue+2)/2 << ")\t" << "score: " << setw(10) << m.score << " length: " << m.Rsize/2 << ", "
                << m.Qsize/2 << endl;
            else
                cout << setw(7) << " " << "(" << std::left << setw(10) << (lRef+2)/2 << ", " << std::left << setw(10) << (rRef+2)/2 << ")   "
                << " (" << std::left << setw(10) << (lQue+2)/2 << ", " << std::left << setw(10) << (rQue+2)/2 << ")\t" << "score: "<< setw(10)
                << m.score << " length:" << m.Rsize/2 << ", " << m.Qsize/2 << endl;
        }else{
            if ((revComplement & 0x1) && commonData::relQueryPos)
                cout << " " << setw(15) << std::left <<strtok((*itR).seq, " ") << " (" << std::left << setw(10) << (lRef+2)/2 << ", " << std::left
                << setw(10) << (rRef+2)/2 << ")   " << " (" << std::left << setw(10) << ((*itQ).end-(*itQ).start-lQue+2)/2 << ", " << std::left
                << setw(10) << ((*itQ).end-(*itQ).start-rQue+2)/2 << ")\t" << "score: " << setw(10) << m.score << " length:" << m.Rsize/2
                << ", " << m.Qsize/2 << endl;
            else
                cout << " " << setw(15) << std::left <<strtok((*itR).seq, " ") << " (" << std::left << setw(10) << (lRef+2)/2 << ", " << std::left
                << setw(10) << (rRef+2)/2 << ")   "<< " (" << std::left << setw(10) << (lQue+2)/2 << ", " << std::left << setw(10) << (rQue+2)/2
                << ")\t" << "score: " << setw(10) << m.score << " length:" << m.Rsize/2 << ", " << m.Qsize/2 << endl;
        }
    }
    
    /* This is the exact same function as "mergeMemExtVector" of e-mem */
    void mergeMemInBorders (uint32_t &revComplement) {
        int flag=0;
        MemExt m;
        if (commonData::d==1 && commonData::numThreads==1)
            return;
        
        if (MemInBorders.size() > 1) {
            do {
                flag=0;
                sort(MemInBorders.begin(), MemInBorders.end(), compare_query);
                for (vector<MemExt>::iterator it=MemInBorders.begin(); it != --MemInBorders.end(); ++it) {
                    vector<MemExt>::iterator dup = it;
                    ++dup;
                    for (; dup != MemInBorders.end(); ++dup) {
                        if((*dup).lQ + static_cast<intType>(commonData::minMemLen-2) > (*it).lQ+static_cast<intType>((*it).Qsize) )
                            break;
                        if((*dup).lQ + static_cast<intType>(commonData::minMemLen-2) == (*it).lQ+static_cast<intType>((*it).Qsize) ) {
                            if((*dup).lR + static_cast<intType>(commonData::minMemLen-2) == (*it).lR+static_cast<intType>((*it).Rsize) ) {
                                flag=1;
                                (*it).Qsize=(*dup).Qsize;
                                (*it).Rsize=(*dup).Rsize;
                                MemInBorders.erase(dup);
                                break;
                            }
                        }
                    }
                    if (flag)
                        break;
                }
                
                sort(MemInBorders.begin(), MemInBorders.end(), compare_reference);
                for (vector<MemExt>::iterator it=MemInBorders.begin(); it != --MemInBorders.end(); ++it) {
                    vector<MemExt>::iterator dup = it;
                    ++dup;
                    for (; dup != MemInBorders.end(); ++dup) {
                        if((*dup).lR + static_cast<intType>(commonData::minMemLen-2) > (*it).lR+static_cast<intType>((*it).Rsize) )
                            break;
                        if((*dup).lR + static_cast<intType>(commonData::minMemLen-2) == (*it).lR+static_cast<intType>((*it).Rsize) ) {
                            if((*dup).lQ + static_cast<intType>(commonData::minMemLen-2) == (*it).lQ+static_cast<intType>((*it).Qsize) ) {
                                flag=1;
                                (*it).Qsize=(*dup).Qsize;
                                (*it).Rsize=(*dup).Rsize;
                                MemInBorders.erase(dup);
                                break;
                            }
                        }
                    }
                    if (flag)
                        break;
                }
            } while (flag);
        }
        
        for (vector<MemExt>::iterator it=MemInBorders.begin(); it != MemInBorders.end(); ++it) {
             writeToFile((*it).lQ, (*it).lQ+(*it).Qsize-2, (*it).lR, (*it).lR+(*it).Rsize-2, (*it).score, revComplement);
        }
        if(!MemInBorders.empty())
            MemInBorders.clear();
    }
    
    /* if IS_MATCH_BOTH_DEF:
     *   forFile contains forward repeats
     *   revFile contains reverse repeats
     *   This function prints out the content of these two files.
     */
    void outputInMummerFormat(fstream &forFile, fstream &revFile) {
        string line, last_line;
        fstream *filePtr;
        static int first=1;
        char buffer[50];

        forFile.clear();
        revFile.clear();   
        forFile.seekg(ios::beg);   
        revFile.seekg(ios::beg);

        filePtr = &forFile;
        if(getline((*filePtr), line).good()) 
            cout << line << endl;

        while(getline((*filePtr), line).good()) {
            if(line[0] == '>'){
                if (last_line.size())
                    cout << last_line << endl;
                last_line = line;
                if ((*filePtr) == forFile) {
                    filePtr = &revFile;
                    if (first) {
                        if(getline((*filePtr), line).good()) 
                            cout << line << endl;
                        first=0;
                    }
                }else
                    filePtr = &forFile;
                continue;
            }
            cout << line << endl;
        }
        forFile.close();

        cout << last_line << endl;
        filePtr = &revFile;
        while(getline((*filePtr), line).good())
            cout << line << endl;
        revFile.close();
        
        memset(buffer,0,50);
        sprintf(buffer, "./%d_tmp/%d", getpid(), NUM_TMP_FILES+NUM_TMP_FILES_REV);
        remove(buffer);
        memset(buffer,0,50);
        sprintf(buffer, "./%d_tmp/%d", getpid(), NUM_TMP_FILES+NUM_TMP_FILES_REV+1);
        remove(buffer);
    }
    
    
    void removeSubsets(vector<MemExt> * vec){
        vector<MemExt>::iterator it=vec->begin();
        MemExt curre = vec->front();
        
        for(++it; it!=vec->end() ; ++it)
            if ( ((*it).lQ+(*it).Qsize <= curre.lQ+curre.Qsize) && ((*it).lR >= curre.lR) && ((*it).lR+(*it).Rsize <= curre.lR+curre.Rsize) ) {
                vec->erase(it);
                it--;
            }
    }
    
    void removeDuplicates(uint32_t revComplement, uint32_t &options) {
        streambuf *coutbuf=std::cout.rdbuf();
        fstream outFile;
        int numFiles=0;
        MemExt mem, curreM(0,0,2,2,0);
        vector<MemExt> * finalVec, * MemExtVector;
        vector<MemExt>::iterator last;
        char buffer[50];
        
        if(IS_MATCH_BOTH_DEF(revComplement))
            numFiles=NUM_TMP_FILES+NUM_TMP_FILES_REV;
        else
            numFiles=NUM_TMP_FILES;

        openFiles(ios::in|ios::binary, numFiles, revComplement);
        
        if (IS_MATCH_BOTH_DEF(revComplement)) {
            memset(buffer,0,50);
            sprintf(buffer, "./%d_tmp/%d", getpid(), numFiles);
            TmpFiles[numFiles].open(buffer, ios::out|ios::in|ios::trunc);
            memset(buffer,0,50);
            sprintf(buffer, "./%d_tmp/%d", getpid(), numFiles+1);
            TmpFiles[numFiles+1].open(buffer, ios::out|ios::in|ios::trunc);
        }
        
        /* Redirect std::cout to the Output file */
        if (IS_OUT_FILE_DEF(options)) {
            outFile.open("Output",ios::out);
            std::cout.rdbuf(outFile.rdbuf());
        }
        
        /* Indication that reverse complement is being processed */
        if (IS_MATCH_REV_DEF(revComplement))
            revComplement|=0x1;
        
        /* Redirect std::cout to a file */
        if (IS_MATCH_BOTH_DEF(revComplement))
            std::cout.rdbuf(TmpFiles[numFiles].rdbuf());
        
        finalVec = new vector<MemExt>;
        for (int32_t i=0; i<NUM_TMP_FILES; i++){
            
            MemExtVector=new vector<MemExt>;
            
            memset(buffer,0,50);
            sprintf(buffer, "./%d_tmp/%d", getpid(),i);
            
            while(!TmpFiles[i].read((char *)&mem, sizeof (MemExt)).eof()) {
                MemExtVector->push_back(mem);
            }
            sort(MemExtVector->begin(), MemExtVector->end(), MemExt());
            last=unique(MemExtVector->begin(), MemExtVector->end(), myUnique);
            MemExtVector->resize(distance(MemExtVector->begin(),last));
            
            TmpFiles[i].close();
            remove(buffer);
                
            for (vector<MemExt>::iterator it=MemExtVector->begin(); it!=MemExtVector->end(); ++it) {
                
                while(finalVec->size() >= 2 && (*it).lQ+commonData::minMemLen > (curreM.lQ+curreM.Qsize)) {
                    
                    removeSubsets(finalVec);
                    printMemOnTerminal(finalVec->front(), revComplement);
                    numOfOutputRepeats++;
                    
                    finalVec->erase(finalVec->begin());
                    if(finalVec->empty())
                        break;
                    curreM = finalVec->front();
                }
                mem.copy(*it);
                if (finalVec->empty())
                    curreM.copy(mem);
                finalVec->push_back(mem);
            }
            if(!MemExtVector->empty())
                MemExtVector->clear();
            delete MemExtVector;
        }
        while( finalVec->size() >= 2)  {
            
            removeSubsets(finalVec);
            printMemOnTerminal(finalVec->front(), revComplement);
            numOfOutputRepeats++;
            finalVec->erase(finalVec->begin());
        }
        if(!finalVec->empty()) {
                    
            printMemOnTerminal(finalVec->front(), revComplement);
            numOfOutputRepeats++;
            finalVec->pop_back();
        }
        if(!finalVec->empty())
            finalVec->clear();
        delete finalVec;
    

        if (IS_MATCH_BOTH_DEF(revComplement)) {
            
            /* Output any unsued query sequence */
            mem.lR=mem.lQ=0;
            mem.Rsize=mem.Qsize=2;
            printMemOnTerminal(mem, revComplement);
            
            /* Processing reverse complement files now */
            revComplement|=0x1;
            /* Redirect output to reverse complement file */
            std::cout.rdbuf(TmpFiles[numFiles+1].rdbuf());
        }
        
        finalVec = new vector<MemExt>;
        for (int32_t i=NUM_TMP_FILES; i<numFiles; i++){
            
            MemExtVector=new vector<MemExt>;
            
            memset(buffer,0,50);
            sprintf(buffer, "./%d_tmp/%d", getpid(),i);
            
            while(!TmpFiles[i].read((char *)&mem, sizeof (MemExt)).eof()) {
                MemExtVector->push_back(mem);
            }
            sort(MemExtVector->begin(), MemExtVector->end(), MemExt());
            last=unique(MemExtVector->begin(), MemExtVector->end(), myUnique);
            MemExtVector->resize(distance(MemExtVector->begin(),last));
            
            TmpFiles[i].close();
            remove(buffer);
                
            for (vector<MemExt>::iterator it=MemExtVector->begin(); it!=MemExtVector->end(); ++it) {
                
                while(finalVec->size() >= 2 && (*it).lQ+commonData::minMemLen > (curreM.lQ+curreM.Qsize)) {
                
                    removeSubsets(finalVec);
                    printMemOnTerminal(finalVec->front(), revComplement);
                    numOfOutputRepeats++;
                
                    finalVec->erase(finalVec->begin());
                    if(finalVec->empty())
                        break;
                    curreM = finalVec->front();
                }
                mem.copy(*it);
                if (finalVec->empty())
                    curreM.copy(mem);
                finalVec->push_back(mem);
            }
            if(!MemExtVector->empty())
                MemExtVector->clear();
            delete MemExtVector;
        }
        while( finalVec->size() >= 2)  {
            
            removeSubsets(finalVec);
            printMemOnTerminal(finalVec->front(), revComplement);
            numOfOutputRepeats++;
            finalVec->erase(finalVec->begin());
        }
        if(!finalVec->empty()) {
            
            printMemOnTerminal(finalVec->front(), revComplement);
            numOfOutputRepeats++;
            finalVec->pop_back();
        }
        if(!finalVec->empty())
            finalVec->clear();
        delete finalVec;
        
        /* Output any unsued query sequence */
        mem.lR=mem.lQ=0;
        mem.Rsize=mem.Qsize=2;
        printMemOnTerminal(/*refSeqInfo, querySeqInfo,*/ mem, revComplement);
        
        /* Restore std::cout */
        if (IS_MATCH_BOTH_DEF(revComplement)){
            if(IS_OUT_FILE_DEF(options))
                std::cout.rdbuf(outFile.rdbuf());
            else
                std::cout.rdbuf(coutbuf);
            outputInMummerFormat(TmpFiles[numFiles], TmpFiles[numFiles+1]);
        }
        
        if(IS_OUT_FILE_DEF(options)) {
            std::cout.rdbuf(coutbuf);
            outFile.close();
        }
        
        if(revComplement) {
            memset(buffer,0,50);
            sprintf(buffer, "./%d_tmp/revComp", getpid());
            remove(buffer);
        }
        memset(buffer,0,50);
        sprintf(buffer, "./%d_tmp", getpid());
        remove(buffer);
    }
};

#endif