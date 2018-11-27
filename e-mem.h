/*  
 *  Purpose: Provides classes for creating a two-level hash table.
 *
 *  Modified by Sarah Banyassady, August 2015
 *
 *   Added classes and functions:
 *   - leftKnode
 *   - getHashIndex
 *   - get2ndLevelHashIndex
 *   - retrieveKey
 *   - reSizeSubTable
 *   - reHashSubTable
 *
 *   Modified functions:
 *   - getHashKey
 *   - findKmer
 *   - addKmerNode
 *
 *   Moved global variables and macros defintions to "global_var.h".
 *
 * ================================================================= *
 *  e-mem.h : Header file for main program                           *
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

#ifndef __EMEM_H_INCLUDED__
#define __EMEM_H_INCLUDED__

#include "global_var.h"
#include "file.h"

/* The left part of the main hash table storing the size counteres */
class leftKnode {
public:
    uint8_t offset;
    
    leftKnode(uint8_t a = 0){
        offset = a;
    }
};

/* The right part of the main hash table storing the pointers to subTables */
class rightKnode {
    intType * subTab_tail;
    
    uint64_t getHashKey(uint64_t * currKmer, seedData &seed)
    {
        uint64_t partialKey, key=0;
        uint32_t numShifts=0, currW=0;
        
        uint32_t llast = seed.length - (seed.binSeedSize-1)*32;     //size of the last chunck of seed
        int32_t binPos = seed.binSeedSize-1;
        currKmer[binPos] >>= DATATYPE_WIDTH-CHARS2BITS(llast);
        
        for(uint32_t i=0; i<llast ; i++)
        {
            if ((seed.binSeed[binPos] & seed_mask[i]) == 0) numShifts++;
            
            else  key |= ((currKmer[binPos] & seed_mask[i]) >> CHARS2BITS(numShifts));
        }
        currW = llast - numShifts;
        
        for(binPos=seed.binSeedSize-2; binPos>=0 ; binPos--)
        {
            partialKey = numShifts = 0;
            for(uint32_t i=0; i<32 ; i++)
            {
                if ((seed.binSeed[binPos] & seed_mask[i]) == 0) numShifts++;
                
                else partialKey |= ((currKmer[binPos] & seed_mask[i]) >> CHARS2BITS(numShifts));
            }
            key |= (partialKey << CHARS2BITS(currW));
            currW += (32 - numShifts);
        }
        return key;
    }
    
    intType getHashIndex(uint64_t key)
    {
        intType index = key % currHashTabSize;
        return index;
    }
    
    uint32_t get2ndLevelHashIndex(uint64_t key, uint32_t & match_found, seedData &seed, seqFileReadInfo &RefFile, uint8_t offset)
    {
        uint32_t subTabSize = initialSubTabSize<<offset;
        uint32_t index = key % subTabSize;
        uint32_t step = (subTabSize>1)?((key % (subTabSize>>1))|1):0;
        uint32_t count=1;
        while (this->subTab_tail[index] != 0)
        {
            if(count > subTabSize){
                match_found = 2;
                return 0;
            }
            if (key == this->retrieveKey(this->subTab_tail[index], seed, RefFile))
            {
                match_found = 1;
                return index;
            }
            else
            {
                index = (index + (count*step)) % subTabSize;
                count++;
            }
        }
        return index;   //no collision
    }
    
    uint64_t retrieveKey(intType pos, seedData &seed, seqFileReadInfo &RefFile)
    {
        pos = CHARS2BITS(pos)-2;
        
        intType j = pos / DATATYPE_WIDTH;
        intType offset = pos % DATATYPE_WIDTH;
        
        uint64_t * Kmer = new uint64_t[seed.binSeedSize];
        for(uint32_t i=0; i<seed.binSeedSize ; i++, j++)
        {
            Kmer[i] = (RefFile.binReads[j]<<offset);
            
            if (offset)    Kmer[i] |= (RefFile.binReads[j+1] >> (DATATYPE_WIDTH-offset));
        }
        
        uint64_t key = this->getHashKey(Kmer, seed);
        delete [] Kmer;
        return key;
    }
    
    void reSizeSubTable(intType index,seedData &seed,seqFileReadInfo &RefFile, leftKnode * &refHash_l)
    {
        uint32_t newSize, oldSize = initialSubTabSize<<refHash_l[index].offset;
        
        /* Copy current information */
        intType * tempP = NULL;
        uint64_t * tempKeys = NULL;
        try{
            tempP = new intType[oldSize];
        }catch(std::bad_alloc& ba){
            std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
        }
        try{
            tempKeys = new uint64_t[oldSize];
        }catch(std::bad_alloc& ba){
            std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
        }
        
        for(uint32_t i=0; i < oldSize ; i++)
        {
            tempP[i] = this->subTab_tail[i];
            if (tempP[i])
                tempKeys[i] = this->retrieveKey(this->subTab_tail[i], seed, RefFile);
        }
        
        do{
            refHash_l[index].offset++;
            newSize = initialSubTabSize<<refHash_l[index].offset;
            
            /* Allocate new space */
            delete [] this->subTab_tail;
            try{
                this->subTab_tail = new intType[newSize];
            }catch(std::bad_alloc& ba){
                std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
            }
            memset(this->subTab_tail, 0, newSize*sizeof(intType));
            
            /* Re-hash subtable */
        }while(!reHashSubTable(oldSize, refHash_l[index].offset, tempP, tempKeys, seed, RefFile));
        
        MemoryForHashTable += (newSize-oldSize) * ( sizeof(intType) );
        
        delete [] tempP;
        delete [] tempKeys;
        return;
    }
    
    bool reHashSubTable(uint32_t oldSize, uint8_t offset, intType * tempP, uint64_t * tempKeys, seedData &seed, seqFileReadInfo &RefFile)
    {
        uint32_t index2, match;
        for(uint32_t i=0; i<oldSize ; i++)
        {
            if(tempP[i])
            {
                match = 0;
                index2 = this->get2ndLevelHashIndex(tempKeys[i], match, seed, RefFile, offset);
                if (match == 0)
                    this->subTab_tail[index2] = tempP[i];
                else if (match == 2)
                {
                    cout<<"Resized again!!"<<endl;
                    return false;
                }
            }
        }
        return true;
    }

public:
    static intType currHashTabSize;
    static uint32_t initialSubTabSize;

    bool findKmer(uint64_t * currKmer, intType &tail, seedData &seed, seqFileReadInfo &RefFile, leftKnode * &refHash_l)
    {
        uint64_t key = this->getHashKey(currKmer, seed);
        intType index = this->getHashIndex(key);
        uint32_t match=0;
        
        if (this[index].subTab_tail == NULL)
            return false;
        
        uint32_t index2 = this[index].get2ndLevelHashIndex(key, match, seed, RefFile, refHash_l[index].offset);
        if(match == 1)
        {
            tail = this[index].subTab_tail[index2];
            return true;
        }
        else
            return false;
    }
    
    void addKmerNode(uint64_t * currKmer, intType currKmerPos, seqFileReadInfo &RefFile, seedData &seed, leftKnode * &refHash_l, intType * &occ)
    {
        uint64_t key = this->getHashKey(currKmer, seed);
        intType index = this->getHashIndex(key);
        uint32_t index2, match = 0;

        if(this[index].subTab_tail == NULL)
        {
            try{
                this[index].subTab_tail = new intType[initialSubTabSize];
            }catch(std::bad_alloc& ba){
                std::cerr << MemoryForHashTable << " bad_alloc caught: " << ba.what() <<endl;
            }
            memset(this[index].subTab_tail, 0, initialSubTabSize*sizeof(intType));

            MemoryForHashTable += initialSubTabSize * ( sizeof(intType) );
            
            index2 = this[index].get2ndLevelHashIndex(key, match, seed, RefFile, refHash_l[index].offset);
            this[index].subTab_tail[index2]= currKmerPos;
        }
        else
        {
            match = 0;
            index2=this[index].get2ndLevelHashIndex(key, match, seed, RefFile, refHash_l[index].offset);

            /* The slot is empty */
            if(match == 0)
                this[index].subTab_tail[index2]= currKmerPos;

            /* A match is found: must be hashed into the same slot */
            else if (match == 1)
            {
                occ[currKmerPos] = this[index].subTab_tail[index2];
                this[index].subTab_tail[index2] = currKmerPos;
            }
            
            else    //if(match == 2)
            {
                this[index].reSizeSubTable(index, seed, RefFile, refHash_l);
                
                match = 0;
                index2 = this[index].get2ndLevelHashIndex(key, match, seed, RefFile, refHash_l[index].offset);
                if(match == 0)
                    this[index].subTab_tail[index2 ]= currKmerPos;
                else
                    cout<<"A hit was missed!"<<endl;
            }
        }
        
    }
    
    rightKnode(intType * a = NULL)  {
        subTab_tail=a;
    }
    
    ~rightKnode()  {
        if(subTab_tail)
            delete [] subTab_tail;
    }
};

intType rightKnode::currHashTabSize = 0;
uint32_t rightKnode::initialSubTabSize = 0;

#endif