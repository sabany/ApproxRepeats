/*  
 *  Purpose: Retrieve pairs of approximate repeats between two DNA
 *           sequences using multiple spaced seeds.
 *
 *  Modified by Sarah Banyassady, August 2015
 *
 *   Added functions:
 *   - gaplessExtension
 *   - gappedExtension
 *
 *   Modified functions:
 *   - helpReportMem
 *   - buildRefHash
 *   - reportMem
 *   - processQuery
 *   - processReference
 *   - checkCommandLineOptions
 *   - print_help_msg
 *   - main
 *
 * ================================================================= *
 *  e-mem.cpp : Main program                                         *
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

#include <iostream>
#include <fstream>
#include <cstdint>         
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <iterator>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>

#include "e-mem.h"
#include "file.h"
#include "qlist.h"
#include "stats.h"

using namespace std;

/*
 * Function builds a hash table for a reference sequence.
 * Input: empty refHash
 * Output: populated refHash
 */
void buildRefHash(rightKnode** &refHash_r, leftKnode** &refHash_l, intType** &occ, seqFileReadInfo &RefFile, seedFileReadInfo &SeedFile)
{
    intType j=0, offset=0;
    intType currKmerPos=0;
    vector<mapObject>::iterator it;
    it = upper_bound(RefFile.blockOfNs.begin(), RefFile.blockOfNs.end(), currKmerPos, mapObject());
    uint32_t seedL, double_check = 0;
    
    uint32_t minSeedLen = CHARS2BITS(SeedFile.minLen);
    uint32_t maxSeedLen = CHARS2BITS(SeedFile.maxLen);
    intType totalBits = CHARS2BITS(RefFile.totalBases-1);
    uint32_t binSeedSize = ceil(maxSeedLen/64.0);
    uint64_t * currKmers = new uint64_t[binSeedSize];
    uint64_t * finalKmers = new uint64_t[binSeedSize];
    int32_t stepCounter = 0;
    
    cout<<"Indexing "<<RefFile.getFileName()<<"..."<<endl;
    for (intType currKmerPos=0; currKmerPos<=totalBits; currKmerPos+=2)
    {
        if (currKmerPos + minSeedLen - 2 > totalBits)
            break;
        
        while (currKmerPos + maxSeedLen - 2 > totalBits){
            maxSeedLen -= 2;
            binSeedSize = ceil(maxSeedLen/64.0);
        }

        double_check = 0;
        if(RefFile.checkKmerForNs(currKmerPos, it, minSeedLen, double_check)){

            stepCounter++;
            if (stepCounter == commonData::stepSize)  stepCounter = 0;
            continue;
        }

        offset = currKmerPos%DATATYPE_WIDTH;
        j=currKmerPos/DATATYPE_WIDTH;
        
        for(uint32_t k=j; k<binSeedSize+j ; k++) {
            
            currKmers[k-j] = (RefFile.binReads[k] << offset);
            
            if(offset)    currKmers[k-j] |= (RefFile.binReads[k+1] >> (DATATYPE_WIDTH-offset));
        }
        
         /* Add kmer to the hash table */
        for(uint32_t i=0; i<SeedFile.getNumSeeds() ; i++) {
            
            if( static_cast<int32_t>(i%commonData::stepSize) == stepCounter ) {
            
                seedL = CHARS2BITS(SeedFile.seeds[i].length);
            
                if ((currKmerPos + seedL - 2) > totalBits)
                    break;
            
                if(double_check) {
                    double_check = 0;
                    if(RefFile.checkKmerForNs(currKmerPos, it, seedL, double_check))
                        break;
                }
            
                for(uint32_t k=0; k<binSeedSize ; k++)
                    finalKmers[k]=currKmers[k];

                refHash_r[i]->addKmerNode(finalKmers, (currKmerPos+2)/2, RefFile, SeedFile.seeds[i], refHash_l[i], occ[i]);
            }
        }
        stepCounter++;
        if (stepCounter == commonData::stepSize)  stepCounter = 0;
    }
    cout<<"Done!"<<endl;
    
    if(currKmers)     delete [] currKmers;
    if(finalKmers)    delete [] finalKmers;
}

/*
 * Perfroms gapless extension on a hit pair
 * Output: true if the extended pair scores at least "minScore",
 *         and false otherwise.
 */
bool gaplessExtension(intType currRPos, intType currQPos, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, uint32_t seedL, mapObject &RefBound, mapObject &QueryBound)
{
    intType lRef=currRPos, lQue=currQPos;
    intType offsetR=0,offsetQ=0;
    intType rRef=currRPos + seedL - 2, rQue=currQPos + seedL - 2;
    uint64_t currR=0, currQ=0, currRtemp=0, currQtemp=0;
    int32_t i=0,j=0, score=0, extScore=0, maxScore=0;
    
    /* Evaluate the seed score */
    offsetR = (currRPos)%DATATYPE_WIDTH;
    i = (currRPos)/DATATYPE_WIDTH;
    offsetQ = (currQPos)%DATATYPE_WIDTH;
    j = (currQPos)/DATATYPE_WIDTH;
    currR = RefFile.binReads[i];
    currQ = QueryFile.binReads[j];
    for(uint32_t k=0; k<seedL; k+=2)
    {
        currRtemp = (currR << offsetR);
        currQtemp = (currQ << offsetQ);
        if((currRtemp & seed_mask[31]) == (currQtemp & seed_mask[31]))
            score += commonData::matchScore;
        else
            score += commonData::misMatScore;
                
        currRPos += 2;
        offsetR += 2;
        if(offsetR == DATATYPE_WIDTH){
            offsetR = 0;
            currR = RefFile.binReads[++i];
        }
        
        currQPos += 2;
        offsetQ += 2;
        if(offsetQ == DATATYPE_WIDTH){
            offsetQ = 0;
            currQ = QueryFile.binReads[++j];
        }
    }
    
    /* Extend to the right */
    extScore = maxScore = 0;
    while (extScore > (maxScore - commonData::xDrop) && (rRef < RefBound.right) && (rQue < QueryBound.right))
    {
        currRtemp = (currR << offsetR);
        currQtemp = (currQ << offsetQ);
        if((currRtemp & seed_mask[31]) == (currQtemp & seed_mask[31]))
            extScore += commonData::matchScore;
        else
            extScore += commonData::misMatScore;
                
        if(extScore > maxScore) maxScore = extScore;
        
        rRef += 2;
        rQue += 2;
        if ((rRef >= RefBound.right) || (rQue >= QueryBound.right))
            break;
        
        currRPos += 2;
        offsetR += 2;
        if(offsetR == DATATYPE_WIDTH){
            offsetR = 0;
            currR = RefFile.binReads[++i];
        }
        
        currQPos += 2;
        offsetQ += 2;
        if(offsetQ == DATATYPE_WIDTH){
            offsetQ = 0;
            currQ = QueryFile.binReads[++j];
        }
    }
    score += extScore;
    
    /* Extend to the left */
    if(lRef && lQue)
    {
        currRPos = lRef;
        currQPos = lQue;
        offsetR = (currRPos)%DATATYPE_WIDTH;
        i = (currRPos)/DATATYPE_WIDTH;
        offsetQ = (currQPos)%DATATYPE_WIDTH;
        j = (currQPos)/DATATYPE_WIDTH;
        currR = RefFile.binReads[offsetR?i:--i];
        currQ = QueryFile.binReads[offsetQ?j:--j];
    
        if(!offsetR) offsetR = DATATYPE_WIDTH;
        if(!offsetQ) offsetQ = DATATYPE_WIDTH;
            
        extScore = maxScore = 0;
        while (extScore > (maxScore - commonData::xDrop) && lRef && lQue && (QueryBound.left < lQue) && (RefBound.left < lRef))

        {
            currRtemp = currR >> (DATATYPE_WIDTH-offsetR);
            currQtemp = currQ >> (DATATYPE_WIDTH-offsetQ);
            if((currRtemp & seed_mask[0]) == (currQtemp & seed_mask[0]))
                extScore += commonData::matchScore;
            else
                extScore += commonData::misMatScore;
        
            if(extScore > maxScore) maxScore = extScore;
        
            lRef -= 2;
            lQue -= 2;
            if ((lRef <= 0) || (lQue <= 0) || (lRef <= RefBound.left) || (lQue <= QueryBound.left))
                break;
        
            currRPos -= 2;
            offsetR -= 2;
            if(offsetR == 0){
                offsetR = DATATYPE_WIDTH;
                currR = RefFile.binReads[--i];
            }
                
            currQPos -= 2;
            offsetQ -= 2;
            if(offsetQ == 0){
                offsetQ = DATATYPE_WIDTH;
                currQ = QueryFile.binReads[--j];
            }
        }
        score += extScore;
    }

    if(score < commonData::minScore)
        return false;
    else
        return true;
}

/*
 * Perfroms gapped extension on a hit pair
 * Output: false if the extended pair is shorter than "minMemLen" or its
 *         edit distance is larger than "maxDistance", and true otherwise.
 */
bool gappedExtension(intType &lRef, intType &lQue, intType &rRef, intType &rQue, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, uint32_t seedL, mapObject &RefBound, mapObject &QueryBound, int32_t &finalScore)
{
    intType midRef = (lRef + seedL - 1) & 0xFFFFFFFFFFFFFFFE,
            midQue = (lQue + seedL - 1) & 0xFFFFFFFFFFFFFFFE;  //& 0xFFFFFFFFFFFFFFFE makes the number even
    
    intType offsetR=0,offsetQ=0, currRPos, currQPos;
    uint64_t currR=0, currQ=0;
    
    vector<tableEntry> prevAntiDiag, antiDiag, nextAntiDiag;
    vector<tableEntry>::iterator it;
    vector<tableEntry>::reverse_iterator Rit;
    
    int32_t score, localMaxScore, maxScore, dist=0, editDist;
    int32_t i,j, k, L,U,preL,preU, itN, iR, iQ;
    
    /* Extend to the right */
    int32_t rRefLen = (RefBound.right - midRef + 2)/2, rQueLen = (QueryBound.right - midQue + 2)/2;
    score=0, localMaxScore=0, maxScore=0, finalScore = 0;
    k=0, L=0, U=0, preL=0, preU=0;
    antiDiag.push_back(tableEntry(0,0));
    
    do{
        k++;
        nextAntiDiag.reserve(k+1);
        for(i=0; i<L ;i++) nextAntiDiag.push_back(tableEntry());
        
        for(i=L; i<=U+1 ; i++){
            
            j=k-i;
            
            score = NEGINFINITY;
            if (preL<i && i<=preU+1 && j>0){
                
                currRPos = midRef + CHARS2BITS(i) - 2;
                currQPos = midQue + CHARS2BITS(j) - 2;
                
                offsetR = (currRPos)%DATATYPE_WIDTH;
                iR = (currRPos)/DATATYPE_WIDTH;
                offsetQ = (currQPos)%DATATYPE_WIDTH;
                iQ = (currQPos)/DATATYPE_WIDTH;
                
                currR = (RefFile.binReads[iR] << offsetR);
                currQ = (QueryFile.binReads[iQ] << offsetQ);
                
                if((currR & seed_mask[31]) == (currQ & seed_mask[31]))
                {
                    score = prevAntiDiag.at(i-1).score + commonData::matchScore;
                    dist = prevAntiDiag.at(i-1).editDist;
                }
                else
                {
                    score = prevAntiDiag.at(i-1).score + commonData::misMatScore;
                    dist = prevAntiDiag.at(i-1).editDist + 1;
                }
            }
            if (i<=U && score < (antiDiag.at(i).score + commonData::indelScore)){
                score = antiDiag.at(i).score + commonData::indelScore;
                dist = antiDiag.at(i).editDist + 1;
            }
            if (L<i && score < (antiDiag.at(i-1).score + commonData::indelScore)){
                score = antiDiag.at(i-1).score + commonData::indelScore;
                dist = antiDiag.at(i-1).editDist + 1;
            }
            
            if (score > localMaxScore) localMaxScore = score;
            if (score < maxScore - commonData::xDrop) score = NEGINFINITY;
            
            nextAntiDiag.push_back(tableEntry(score,dist));
        }
        
        preL = L, itN=0;
        for(it = nextAntiDiag.begin(); it != nextAntiDiag.end(); ++it){
            if((*it).score > NEGINFINITY)
                break;
            itN++;
        }
        L = itN;
        
        preU = U, itN=nextAntiDiag.size()-1;
        for(Rit = nextAntiDiag.rbegin(); Rit != nextAntiDiag.rend() ; ++Rit){
            if((*Rit).score > NEGINFINITY)
                break;
            itN--;
        }
        U = itN;
        
        if( k+1-rQueLen > L) L = k+1-rQueLen;
        if( rRefLen-1 < U ) U = rRefLen-1;
        
        maxScore = localMaxScore;
        prevAntiDiag = antiDiag;
        antiDiag = nextAntiDiag;
        nextAntiDiag.clear();
        
    }while(L<=U+1);
    
    /* Entries  of the last antiDiagonal are all -INFINITY */
    if (U == -1){
        antiDiag = prevAntiDiag;
        k--;
    }
    
    itN = 0;
    it = antiDiag.begin();
    maxScore = (*it).score, dist = (*it).editDist;
    for(; it != antiDiag.end(); ++it, ++itN)
    {
        if((*it).score > maxScore){
            maxScore = (*it).score;
            dist = (*it).editDist;
            i = itN;
        }
    }
    j=k-i;
    editDist = dist;
    finalScore = maxScore;

    rRef = midRef + CHARS2BITS(i) - 2;
    rQue = midQue + CHARS2BITS(j) - 2;
    
    if (editDist > commonData::maxDistance)
        return false;
    
    if ((rRef+2) < static_cast<uint32_t>(commonData::minMemLen))
        return false;
    
    if ((rQue+2) < static_cast<uint32_t>(commonData::minMemLen))
        return false;
    
    /* Extend to the left */
    int32_t lRefLen = (midRef - RefBound.left)/2, lQueLen = (midQue - QueryBound.left)/2;
    score=0, localMaxScore=0, maxScore=0;
    k=0, L=0, U=0, preL=0, preU=0;
    
    prevAntiDiag.clear();
    antiDiag.clear();
    nextAntiDiag.clear();
    antiDiag.push_back(tableEntry(0,0));
    
    do{
        k++;
        nextAntiDiag.reserve(k+1);
        for(i=0; i<L ;i++) nextAntiDiag.push_back(tableEntry());
        
        for(i=L; i<=U+1 ; i++){
            
            j=k-i;
            
            score = NEGINFINITY;
            if (preL<i && i<=preU+1 && j>0){
                
                currRPos = midRef - CHARS2BITS(i);
                currQPos = midQue - CHARS2BITS(j);
                
                offsetR = (currRPos)%DATATYPE_WIDTH;
                iR = (currRPos)/DATATYPE_WIDTH;
                offsetQ = (currQPos)%DATATYPE_WIDTH;
                iQ = (currQPos)/DATATYPE_WIDTH;
                
                currR = (RefFile.binReads[iR] << offsetR);
                currQ = (QueryFile.binReads[iQ] << offsetQ);
                
                if((currR & seed_mask[31]) == (currQ & seed_mask[31]))
                {
                    score = prevAntiDiag.at(i-1).score + commonData::matchScore;
                    dist = prevAntiDiag.at(i-1).editDist;
                }
                else
                {
                    score = prevAntiDiag.at(i-1).score + commonData::misMatScore;
                    dist = prevAntiDiag.at(i-1).editDist + 1;
                }
            }
            if (i<=U && score < (antiDiag.at(i).score + commonData::indelScore)){
                score = antiDiag.at(i).score + commonData::indelScore;
                dist = antiDiag.at(i).editDist + 1;
            }
            if (L<i && score < (antiDiag.at(i-1).score + commonData::indelScore)){
                score = antiDiag.at(i-1).score + commonData::indelScore;
                dist = antiDiag.at(i-1).editDist + 1;
            }
            
            if (score > localMaxScore) localMaxScore = score;
            if (score < maxScore - commonData::xDrop) score = NEGINFINITY;
            
            nextAntiDiag.push_back(tableEntry(score,dist));
        }
        
        preL = L, itN=0;
        for(it = nextAntiDiag.begin(); it != nextAntiDiag.end(); ++it){
            if((*it).score > NEGINFINITY)
                break;
            itN++;
        }
        L = itN;
        
        preU = U, itN=nextAntiDiag.size()-1;
        for(Rit = nextAntiDiag.rbegin(); Rit != nextAntiDiag.rend() ; ++Rit){
            if((*Rit).score > NEGINFINITY)
                break;
            itN--;
        }
        U = itN;
        
        if( k+1-lQueLen > L) L = k+1-lQueLen;
        if( lRefLen-1 < U ) U = lRefLen-1;
        
        maxScore = localMaxScore;
        prevAntiDiag = antiDiag;
        antiDiag = nextAntiDiag;
        nextAntiDiag.clear();
        
    }while(L<=U+1);
    
    /* Entries  of the last antiDiagonal are all -INFINITY */
    if (U == -1){
        antiDiag = prevAntiDiag;
        k--;
    }
    
    itN = 0;
    it = antiDiag.begin();
    maxScore = (*it).score, dist = (*it).editDist;
    for(; it != antiDiag.end(); ++it, ++itN)
    {
        if((*it).score > maxScore){
            maxScore = (*it).score;
            dist = (*it).editDist;
            i = itN;
        }
    }
    j=k-i;
    editDist += dist;
    finalScore += maxScore;
    
    lRef = midRef - CHARS2BITS(i);
    lQue = midQue - CHARS2BITS(j);
    
    if (editDist > commonData::maxDistance)
        return false;
    
    if ((rQue-lQue+2) < static_cast<uint32_t>(commonData::minMemLen))
        return false;
    
    if ((rRef-lRef+2) < static_cast<uint32_t>(commonData::minMemLen))
        return false;
    return true;
}

/*
 * Calls the gapless and the gapped extension processes, then checks 
 * whether the extended pair is completely contained in a previously
 * detected repeat.
 */
void helpReportMem(intType currRPos, intType currQPos, queryList* &currQueryMEMs, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, uint32_t seedL, tmpFilesInfo &arrayTmp, mapObject &RefNpos, mapObject &QueryNpos, uint32_t &revComplement)
{
    /*
     * lRef and lQue are local variables for left extension of
     * reference and query sequence respectively. rRef and rQue
     * are their right counterparts.
     */
    intType lRef = currRPos, lQue = currQPos;
    intType rRef = currRPos + seedL - 2, rQue = currQPos + seedL - 2;
    seqData refSeq, queSeq;
    mapObject RefBound, QueryBound;
    int32_t finalScore;

    arrayTmp.getStartnEndOfSequence(currRPos, currQPos, refSeq, queSeq);

    if (!(((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && rQue<=QueryNpos.right))
        QueryFile.getKmerLeftnRightBoundForNs(lQue, QueryNpos);
    
    if (queSeq.start >= ((QueryNpos.left==0x1)?0:QueryNpos.left))
          QueryBound.left = queSeq.start;
    else  QueryBound.left = QueryNpos.left;
    
    if (queSeq.end <= QueryNpos.right)
          QueryBound.right = queSeq.end;
    else  QueryBound.right = QueryNpos.right;
    
    if ((QueryBound.right-QueryBound.left+2) < static_cast<uint32_t>(commonData::minMemLen))
        return;

    if (!(((RefNpos.left==0x1)?true:RefNpos.left<=lRef) && rRef<=RefNpos.right))
        RefFile.getKmerLeftnRightBoundForNs(lRef, RefNpos);
    
    if (refSeq.start >= ((RefNpos.left==0x1)?0:RefNpos.left))
          RefBound.left = refSeq.start;
    else  RefBound.left = RefNpos.left;
    
    if (refSeq.end <= RefNpos.right)
          RefBound.right = refSeq.end;
    else  RefBound.right = RefNpos.right;
    
    if ((RefBound.right-RefBound.left+2) < static_cast<uint32_t>(commonData::minMemLen))
        return;
    
    if (!gaplessExtension(currRPos, currQPos, RefFile, QueryFile, seedL, RefBound, QueryBound))
        return;
    
    if (!gappedExtension(lRef, lQue, rRef, rQue, RefFile, QueryFile, (seedL/2), RefBound, QueryBound, finalScore))
        return;

    if(!(currQueryMEMs->checkRedundantMEM(&currQueryMEMs, lRef, rRef, lQue, rQue)))
    {
        arrayTmp.writeMemInTmpFilesVecs(lRef, rRef, lQue, rQue, finalScore, QueryFile, RefFile, revComplement);
        currQueryMEMs->ListAdd(&currQueryMEMs, lQue, rQue,lRef, rRef);
    }
}

void reportMEM(rightKnode** &refHash_r, leftKnode** &refHash_l, intType** &occ, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, seedFileReadInfo &SeedFile, tmpFilesInfo &arrayTmp, uint32_t &revComplement)
{
  uint32_t minSeedLen = CHARS2BITS(SeedFile.minLen);
  uint32_t maxSeedLen = CHARS2BITS(SeedFile.maxLen);
  intType totalQBits =  CHARS2BITS(QueryFile.totalBases-1);
    
  #pragma omp parallel num_threads(commonData::numThreads) firstprivate(maxSeedLen)
  {
      queryList *currQueryMEMs = NULL;
      intType j=0, offset=0, occTail, RefPos;
      mapObject QueryNpos, RefNpos;
      vector<mapObject>::iterator it;
      uint32_t seedL;
      it = upper_bound(QueryFile.blockOfNs.begin(), QueryFile.blockOfNs.end(), 0, mapObject());
      
      uint32_t binSeedSize = ceil(maxSeedLen/64.0);
      uint64_t * currKmers = new uint64_t[binSeedSize];
      uint64_t * finalKmers = new uint64_t[binSeedSize];
      
      int32_t stepCounter = 0;
      #pragma omp for
      for (intType currKmerPos=0; currKmerPos<=totalQBits; currKmerPos+=2)
      {
          if ((currKmerPos + minSeedLen - 2) > totalQBits)
              continue;
          
          while((currKmerPos + maxSeedLen - 2) > totalQBits){
              maxSeedLen -= 2;
              binSeedSize = ceil(maxSeedLen/64.0);
          }
        
          uint32_t double_check = 0;
          if(QueryFile.checkKmerForNs(currKmerPos, it, minSeedLen, double_check)) {
              /* Do not process this Kmer, Ns in it */
              stepCounter++;
              if (stepCounter == commonData::stepSize)  stepCounter = 0;
              continue;
          }

          j=currKmerPos/DATATYPE_WIDTH;
          offset = currKmerPos%DATATYPE_WIDTH;
          
          for(uint32_t k=j; k<binSeedSize+j ; k++)
          {
              currKmers[k-j] = (QueryFile.binReads[k] << offset);
              
              if(offset)   currKmers[k-j] |= (QueryFile.binReads[k+1] >> (DATATYPE_WIDTH-offset));
          }

          /* Find the K-mer in the refHash */
          for(uint32_t i=0; i<SeedFile.getNumSeeds() ; i++)
          {
              if( static_cast<int32_t>(i%commonData::stepSize) == stepCounter ) {
                  
                  seedL = CHARS2BITS(SeedFile.seeds[i].length);
                  
                  if ((currKmerPos + seedL - 2) > totalQBits)
                      break;
                  
                  if(double_check) {
                      double_check = 0;
                      if(QueryFile.checkKmerForNs(currKmerPos, it, seedL, double_check))
                          break;
                  }
                  
                  for(uint32_t k=0; k<binSeedSize ; k++)
                      finalKmers[k]=currKmers[k];
                  
                  /* Find the subSequence of query in the hash table */
                  occTail = 0;
                  if(refHash_r[i]->findKmer(finalKmers, occTail, SeedFile.seeds[i], RefFile, refHash_l[i]))
                      while(occTail != 0)
                      {
                          RefPos = CHARS2BITS(occTail-1);
                          helpReportMem(RefPos, currKmerPos, currQueryMEMs, RefFile, QueryFile, seedL, arrayTmp, RefNpos, QueryNpos, revComplement);
                          occTail = occ[i][occTail];
                      }
              }
          }
          stepCounter++;
          if (stepCounter == commonData::stepSize)  stepCounter = 0;
      }
      currQueryMEMs->ListFree(&currQueryMEMs);
      if(currKmers)      delete [] currKmers;
      if(finalKmers)     delete [] finalKmers;
  }
}

void processQuery(rightKnode** &refHash_r, leftKnode** &refHash_l, intType** &occ, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, seedFileReadInfo &SeedFile, tmpFilesInfo &arrayTmp, uint32_t &revComplement, int32_t numTmpFiles)
{
    if (IS_MATCH_REV_DEF(revComplement) || IS_MATCH_BOTH_DEF(revComplement))
        cout<<"Scanning Reverse Sequence..."<<endl;
    else
        cout<<"Scanning "<<QueryFile.getFileName()<<"..."<<endl;
    
    QueryFile.clearFileFlag();
    QueryFile.resetCurrPos();
    for (int32_t i=0; i<commonData::d; i++) {
        if(QueryFile.readChunks()){
            
            /* statistics */
            fillNbLetters(QueryFile, 1);
            if (!setMinScore(RefFile.totalBases, QueryFile.totalBases))
                EXIT_AND_REMOVE(numTmpFiles, revComplement)
            
            reportMEM(refHash_r, refHash_l, occ, RefFile, QueryFile, SeedFile, arrayTmp, revComplement);
            QueryFile.setCurrPos();
            QueryFile.clearMapForNs();
        }
        else
            break;
    }
    QueryFile.clearTmpString();
    cout<<"Done!"<<endl;
}

void processReference(seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, seedFileReadInfo &SeedFile, tmpFilesInfo &arrayTmp, uint32_t &revComplement, int32_t numTmpFiles, uint32_t lastCall)
{
    uint32_t n;
    static rightKnode ** refHash_r;
    static leftKnode ** refHash_l;
    static intType ** occ;
    
    if(!IS_MATCH_BOTH_DEF(revComplement)){
    
        intType numberOfKmers = RefFile.totalBases - SeedFile.minLen + 1;

        /* Set the size of the hash table to the numberofKmers/numberofSeeds */
        for (n=0; n<450; ++n)
        {
            if (hashTableSize[n] > (numberOfKmers/SeedFile.getNumSeeds()))
            {
                rightKnode::currHashTabSize = (intType)hashTableSize[n];
                break;
            }
        }
    
        /* Set the initial size of second-level hash tables */
        uint32_t powOf2L = 1, powOf2U = 2;
        while((2 * SeedFile.getNumSeeds() - powOf2L) > powOf2U)
        {
            powOf2L = powOf2U;
            powOf2U <<=1 ;
        }
        rightKnode::initialSubTabSize = powOf2L;


        /* Create the refHash for K-mers */
        try{
            refHash_r = new rightKnode*[SeedFile.getNumSeeds()];
        }catch(std::bad_alloc& ba){
            std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
        }
        try{
            refHash_l = new leftKnode*[SeedFile.getNumSeeds()];
        }catch(std::bad_alloc& ba){
            std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
        }
    
        MemoryForHashTable += SeedFile.getNumSeeds() * (sizeof(rightKnode*)+sizeof(leftKnode*));
    
        for(n=0; n<SeedFile.getNumSeeds() ; n++)
        {
            try{
                refHash_r[n] = new rightKnode[rightKnode::currHashTabSize];
            }catch(std::bad_alloc& ba){
                std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
            }
        
            try{
                refHash_l[n] = new leftKnode[rightKnode::currHashTabSize];
            }catch(std::bad_alloc& ba){
                std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
            }
            MemoryForHashTable += rightKnode::currHashTabSize * (sizeof(rightKnode)+sizeof(leftKnode));
        }
    
        /* Allocate memory for lists of occurrences */
        try{
            occ=new intType*[SeedFile.getNumSeeds()];
        }catch(std::bad_alloc& ba){
            std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
        }
    
        MemoryForHashTable += SeedFile.getNumSeeds() * sizeof(intType*);
    
        for(n=0; n<SeedFile.getNumSeeds() ; n++)
        {
            try{
                occ[n]=new intType[numberOfKmers+1];
            }catch(std::bad_alloc& ba){
                std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
            }
            memset(occ[n],0,(numberOfKmers+1)*sizeof(intType));
        
            MemoryForHashTable += (numberOfKmers+1)*sizeof(intType);
        }

        buildRefHash(refHash_r, refHash_l, occ, RefFile, SeedFile);
        
        cout << "The Hash Tables Size: " << MemoryForHashTable/1024 << " Kbyte" << endl;
        MemoryForHashTable = 0;
    }
    processQuery(refHash_r, refHash_l, occ, RefFile, QueryFile, SeedFile, arrayTmp, revComplement, numTmpFiles);
    
    if(lastCall)
    {
        for(n=0; n<SeedFile.getNumSeeds(); n++)
        {
            if(refHash_r[n]) delete [] refHash_r[n];
            if(refHash_l[n]) delete [] refHash_l[n];
            if(occ[n])       delete [] occ[n];
        }
        if(refHash_r)  delete [] refHash_r;
        if(refHash_l)  delete [] refHash_l;
        if(occ)        delete [] occ;
    }
    return;
}

bool is_numeric(const string &str)
{
    return all_of(str.begin(), str.end(), ::isdigit);
}

void checkCommandLineOptions(uint32_t &options)
{
    if (!IS_REF_FILE_DEF(options)){
        cout << "ERROR: reference file must be passed!" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (!IS_QUERY_FILE_DEF(options)){
        cout << "ERROR: query file must be passed!" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (!IS_SEED_FILE_DEF(options)){
        cout << "ERROR: seed file must be passed!" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (IS_SPLIT_SIZE_DEF(options)){
        if (commonData::d <= 0){
            cout << "ERROR: -d cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_NUM_THREADS_DEF(options)){
        if (commonData::numThreads <= 0){
            cout << "ERROR: -t cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_LENGTH_DEF(options)){
        if (commonData::minMemLen <= 2){
            cout << "ERROR: -l cannot be less than or equal to one!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_STEP_SIZE_DEF(options)){
        if (commonData::stepSize <= 0){
            cout << "ERROR: -k cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_MAX_DISTANCE_DEF(options)){
        if (commonData::maxDistance < 0){
            cout << "ERROR: -D cannot be less than zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_XDROP_DEF(options)){
        if (commonData::xDrop < 0){
            cout << "ERROR: Xdrop cannot be less than zero!" << endl;
            exit(EXIT_FAILURE);
        }
        else if (commonData::xDrop > 100){
            cout << "ERROR: Xdrop cannot be greater than 100" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_EVALUE_DEF(options)){
        if (commonData::expectationValue < 1e-100){
            cout << "ERROR: E-value cannot be less than 1e-100!" << endl;
            exit(EXIT_FAILURE);
        }
        else if (commonData::expectationValue > 1e+10){
            cout << "ERROR: E-value cannot be greater than 1e+10" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (IS_MATCH_REV_DEF(options) && IS_MATCH_BOTH_DEF(options)) {
        cout << "ERROR: option -b and option -r exclude each other!" << endl;
        exit(EXIT_FAILURE);
    }
  
    if(IS_RELREV_QUEPOS_DEF(options)) {
        if (!IS_MATCH_REV_DEF(options) && !IS_MATCH_BOTH_DEF(options)) {
            cout << "ERROR: option -c requires either option -r or - b" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void print_help_msg()
{
    cout << endl;
    cout << "ApproxRepeats finds and outputs the position, length, and score" << endl;
    cout << "of approximate repeats between <reference-file> and <query-file>" << endl;
    cout << endl;
    cout << "Usage: ./ApproxRepeats [options]  <reference-file>  <query-file>  <seed-file>" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-n" << "Match only the characters a, c, g, or t." << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "They can be in upper or in lower case." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-l length" << "Set the minimum length of a repeat." << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "The default length is the length of the shortest seed." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-s scr1 scr2 scr3" << "Set the scoring matrix for [ Match, Mismatch, Indel ]." << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "The default is [ 2, -2, -3 ]." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-D edit-distance" << "Set the maximum edit distance between two repeats." << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "The default distance is 5." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-x xdrop" << "Set the Xdrop threshold score. The default is 5." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-e evalue" << "Set the E-value threshold. The default is 10." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-b" << "Compute forward and reverse complement repeats." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-r" << "Only compute reverse complement repeats." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-c" << "Report the query-position of a reverse complement repeat" << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "relative to the original query sequence." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-o" << "Print the output in the 'Output' file." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-f" << "Force printing the reference sequence name beside each pair" << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "regardless of the number of reference sequence input." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-L" << "Show the length of the query sequences on the header line." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-k step-size" << "Set the step size for indexing/scanning the two sequences." << endl;
    cout << setw(2) << " " << std::left << setw(20) << " " << "The default step is 1." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-d split-size" << "Set the split size. The default value is 1." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-t threads" << "Set the number of threads. The default is 1 thread." << endl;
    cout << setw(2) << " " << std::left << setw(20) << "-h" << "Show possible options." << endl;
}

int main (int argc, char *argv[])
{    
    int32_t i=0, n=1, numTmpFiles;
    uint32_t options=0, revComplement=0, lastCall=0;
    seqFileReadInfo RefFile, QueryFile;
    seedFileReadInfo SeedFile;
    
    /* Check arguments */
    if (argc==1 || argc==2 || argc==3){
       print_help_msg();
       exit(EXIT_SUCCESS);
    }
    
    while(argv[n]) {
        if(boost::equals(argv[n],"-l")){
            if (IS_LENGTH_DEF(options)) {
                cout << "ERROR: Length argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_LENGTH(options);
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -l option!" << endl;
                exit(EXIT_FAILURE);
            }
            commonData::minMemLen = 2*std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-t")){
            if (IS_NUM_THREADS_DEF(options)) {
                cout << "ERROR: Number of threads argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -t option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_NUM_THREADS(options);
            commonData::numThreads = std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-D")){
            if (IS_MAX_DISTANCE_DEF(options)) {
                cout << "ERROR: Maximum distance argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -D option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_MAX_DISTANCE(options);
            commonData::maxDistance = std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-s")){
            if (IS_SCORE_MAT_DEF(options)) {
                cout << "ERROR: Scoring matrix passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] ||  !argv[n+2] || !argv[n+3]){
                cout << "ERROR: -s option must be followed by 3 integers !" << endl;
                exit(EXIT_FAILURE);
            }
            if ((!is_numeric(argv[n+1]) && argv[n+1][0] != '-') ||
                (!is_numeric(argv[n+2]) && argv[n+2][0] != '-') ||
                (!is_numeric(argv[n+3]) && argv[n+3][0] != '-')){
                cout << "ERROR: Invalid value for -s option!" << endl;
                exit(EXIT_FAILURE);
            }
            commonData::matchScore = std::stoi(argv[n+1]);
            commonData::misMatScore = std::stoi(argv[n+2]);
            commonData::indelScore = std::stoi(argv[n+3]);
            SET_SCORE_MAT(options);
            n+=4;
        }else if (boost::equals(argv[n],"-x")){
            if (IS_XDROP_DEF(options)) {
                cout << "ERROR: Xdrop argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -x option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_XDROP(options);
            commonData::xDrop = std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-e")){
            if (IS_EVALUE_DEF(options)) {
                cout << "ERROR: E-value argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !strtod(argv[n+1], NULL)){
                cout << "ERROR: Invalid value for -e option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_EVALUE(options);
            commonData::expectationValue = strtod(argv[n+1], NULL);
            n+=2;
        }else if (argv[n][0] != '-'){
            /* These are files */
            if (!IS_REF_FILE_DEF(options)) {
                /* Open referencead file provided by the user */
                SET_REF_FILE(options);
                RefFile.openFile(argv[n]);
                n+=1;
                continue;
            }
            if (!IS_QUERY_FILE_DEF(options)) {
                /* Open query file provided by the user */
                SET_QUERY_FILE(options);
                QueryFile.openFile(argv[n]);
                n+=1;
                continue;
            }
            if (!IS_SEED_FILE_DEF(options)) {
                /* Open seed file provided by the user */
                SET_SEED_FILE(options);
                SeedFile.openFile(argv[n]);
                n+=1;
                continue;
            }
            cout << "ERROR: More input files than expected!" << endl;
            exit(EXIT_FAILURE);
        }else if (boost::equals(argv[n],"-r")){
            if (IS_MATCH_REV_DEF(options)) {
                cout << "ERROR: Reverse repeat argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_MATCH_REV(options);
            n+=1;
        }else if (boost::equals(argv[n],"-b")){
            if (IS_MATCH_BOTH_DEF(options)) {
                cout << "ERROR: option -b passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_MATCH_BOTH(options);
            n+=1;
        }else if (boost::equals(argv[n],"-n")){
            if (IS_IGNORE_N_DEF(options)) {
                cout << "ERROR: Ignore N's argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_IGNORE_N(options);
            commonData::ignoreN = 1;
            n+=1;
        }else if (boost::equals(argv[n],"-c")){
            if (IS_RELREV_QUEPOS_DEF(options)) {
                cout << "ERROR: option -c passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_RELREV_QUEPOS(options);
            n+=1;
            commonData::relQueryPos = 1;
        }else if (boost::equals(argv[n],"-f")){
            if (IS_FCOL_OUTPUT_DEF(options)) {
                cout << "ERROR: option -f passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FCOL_OUTPUT(options);
            n+=1;
            commonData::fColOutput = 1;
        }else if (boost::equals(argv[n],"-o")){
            if (IS_OUT_FILE_DEF(options)) {
                cout << "ERROR: option -o passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_OUT_FILE(options);
            n+=1;
        }else if (boost::equals(argv[n],"-L")){
            if (IS_LEN_IN_HEADER_DEF(options)) {
                cout << "ERROR: option -L passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_LEN_IN_HEADER(options);
            n+=1;
            commonData::lenInHeader = 1;
        }else if (boost::equals(argv[n],"-k")){
            if (IS_STEP_SIZE_DEF(options)) {
                cout << "ERROR: Step size argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -k option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_STEP_SIZE(options);
            commonData::stepSize = 2*std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-d")){
            if (IS_SPLIT_SIZE_DEF(options)) {
                cout << "ERROR: Split size argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -d option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_SPLIT_SIZE(options);
            commonData::d = std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-mum")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-mumcand")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-mumreference")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-maxmatch")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-h")){
            print_help_msg();
            exit(EXIT_SUCCESS);
        }else {
            cout << "ERROR: Invalid option." << endl << flush;
            print_help_msg();
            exit( EXIT_FAILURE );
        }
    }

    checkCommandLineOptions(options);
    
    SeedFile.setNumSeeds();
    SeedFile.readSeeds();
    
    /* statistics */
    fillSubstitutionMatrix();
    if (!IS_LENGTH_DEF(options))
        commonData::minMemLen = CHARS2BITS(SeedFile.minLen);
    numTmpFiles = IS_MATCH_BOTH_DEF(options)?(NUM_TMP_FILES+NUM_TMP_FILES_REV+2):NUM_TMP_FILES;
    
    tmpFilesInfo arrayTmp(numTmpFiles);
    arrayTmp.openFiles(ios::out|ios::binary, numTmpFiles, revComplement);
    
    if ( !RefFile.generateSeqPosMap(arrayTmp.getRefInfoVec(), 0) )
        EXIT_AND_REMOVE(numTmpFiles, revComplement)
        
    if ( !QueryFile.generateSeqPosMap(arrayTmp.getQueryInfoVec(), (IS_MATCH_REV_DEF(options) || IS_MATCH_BOTH_DEF(options))) )
        EXIT_AND_REMOVE(numTmpFiles, revComplement)
    
    /* Only reverse complement repeats */
    if (IS_MATCH_REV_DEF(options)){
        QueryFile.setReverseFile();
        SET_MATCH_REV(revComplement);   //revComplement = MATCH_REV
    }
    /* Only reverse complement repeats | Only forward repeats */
    if(!IS_MATCH_BOTH_DEF(options))
        lastCall = 1;
    
    arrayTmp.setNumMemsInFileVec(QueryFile.allocBinArray());
    RefFile.allocBinArray();
    RefFile.clearFileFlag();
    
    while (true)
    {
        for (i=0; i<commonData::d; i++) {
            if(RefFile.readChunks()){
                /* statistics */
                fillNbLetters(RefFile, 0);
                processReference(RefFile, QueryFile, SeedFile, arrayTmp, revComplement, numTmpFiles, lastCall);
                RefFile.setCurrPos();
                RefFile.clearMapForNs();
            }
            else
                break;
        }
        
        /*
         * Process MemInBorders list
         */
        
        arrayTmp.mergeMemInBorders(revComplement);

        if (revComplement)
            break;
        if (IS_MATCH_BOTH_DEF(options)){
            SET_MATCH_BOTH(revComplement);  //revComplement = MATCH_BOTH
            lastCall = 1;
            RefFile.clearFileFlag();
            RefFile.resetCurrPos();
            RefFile.totalBases=0;
            QueryFile.setReverseFile();
            QueryFile.totalBases=0;
        }
        else
            break;
    }
    
    /*
     * Free up the allocated arrays
     */
    arrayTmp.closeFiles(numTmpFiles);
    RefFile.destroy();
    QueryFile.destroy();
    SeedFile.destroy();
    arrayTmp.removeDuplicates(revComplement, options);
    
    if (IS_OUT_FILE_DEF(options))
        cout << "Output file written successfully!" << endl;

    cout<< "\nThe number of detected repeats: " << numOfOutputRepeats << endl;
    return 0;
}