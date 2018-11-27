/*  
 *  Purpose: Provides classes for storing list of approximate repeats
 *           between query and reference sequences in memory.
 *
 *  Modified by Sarah Banyassady, August 2015
 *
 *  Modified classes:
 *   - MEMS_LIST:
 *      Modified to not store the key value of a match which  can  be
 *      used to retrive the position of its corresponding exact match
 *      in the reference sequence.  Instead, directly store the  left
 *      and right positions of the corresponding approximate match in
 *      the reference sequence.
 *   - queryList:
 *      Modified to comply with the changes in the MEMS_LIST class.
 *
 * ================================================================= *
 *  qlist.h : Header file with supporting class definitions          *
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

#ifndef __QLIST_H_INCLUDED__
#define __QLIST_H_INCLUDED__

class qneryList; 
 
class MEMS_LIST {
    friend class queryList;
    intType rLeft;
    intType rRight;
    MEMS_LIST *next;
};

class queryList { 
    intType left;
    intType right;
    class queryList* next;
    class MEMS_LIST* mems;

    queryList *queryList_Alloc (intType l=0, intType r=0, intType rl=0, intType rr=0)
    {
        queryList *q = new queryList;
        q->left = l;
        q->right = r;
        q->next = NULL;
        q->mems = new MEMS_LIST;
        q->mems->rLeft = rl;
        q->mems->rRight = rr;
        q->mems->next = NULL;
        return q;
    }
  
 public:
    void ListFree(queryList** listRef)
    {
        if (!*listRef)
            return;
        while(*listRef)
        {
            /* Free all MEMs node found using this Kmer */
            MEMS_LIST *remMems = (*listRef)->mems, *remMemNext=NULL;
            while(remMems) {
                remMemNext = remMems->next;
                delete remMems;
                remMems = remMemNext;
            }
            queryList *tmp = (*listRef)->next;
            delete *listRef;
            *listRef=tmp;
        }
    }

    void ListAdd(queryList** listRef, intType l, intType r,intType rl, intType rr)
    {
        queryList *prev_node=*listRef;
        queryList *node=*listRef;

        if (*listRef == NULL) {
            *listRef = queryList_Alloc (l, r, rl, rr);
            return;
        }

        while(node)
        {
            if (node->right > r) {
                queryList *p = queryList_Alloc(l, r, rl, rr);
                p->next = node;
                if (node == this)
                   *listRef = p;
                else
                   prev_node->next = p;
                return;
            }else if (node->right == r) {
                if(node->left == l){
                    /* Add MEM_LIST item */
                    MEMS_LIST *mems=new MEMS_LIST;
                    mems->rLeft = rl;
                    mems->rRight = rr;

                    MEMS_LIST *prev_mem = node->mems;
                    MEMS_LIST *curr_mem = node->mems;
                    while(curr_mem)
                    {
                        if(curr_mem->rRight >= rr){
                            mems->next = curr_mem;
                            if (curr_mem == prev_mem)
                                node->mems = mems;
                            else
                                prev_mem->next = mems;
                            return;
                        }
                        prev_mem = curr_mem;
                        curr_mem = curr_mem->next;
                        if(!curr_mem) { //end of mems list
                            mems->next = NULL;
                            prev_mem->next = mems;
                            return;
                        }
                    }
                }
            }

            prev_node = node;
            node=node->next;
            if (!node) {    //end of list
                queryList *p = queryList_Alloc(l, r, rl, rr);
                prev_node->next = p;
                return;
            }
         }
         return;
    }

    bool checkRedundantMEM(queryList ** listRef, intType refLeft, intType refRight, intType QueryLeft, intType QueryRight)
    {
        /* Find MEM positions from currQueryMEMs list */
        queryList *p = *listRef; 
        while (p) {
            if (QueryLeft >= p->left && QueryRight <= p->right) {
                MEMS_LIST *mems = p->mems;
                while (mems) {
                    if(refLeft >= mems->rLeft && refRight <= mems->rRight)  //found
                        return true;
                    else
                        mems = mems->next;
                }
                
            }else if ((QueryLeft+commonData::minMemLen-2) > p->right) {
                /* Free all MEMs node found using this Kmer */
                MEMS_LIST *remMems = p->mems, *remMemNext=NULL;
                p->mems = NULL;
                while(remMems) {
                    remMemNext = remMems->next;
                    delete remMems;
                    remMems = remMemNext;
                }
                *listRef = p->next;
                delete p;
                p = *listRef;
                continue;
            }
            p=p->next;
        }
        return false;
    }
};

#endif