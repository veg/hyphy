/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (spond@ucsd.edu)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@ucsd.edu)
 
 Module Developers:
 Lance Hepler (nlhepler@gmail.com)
 Martin Smith (martin.audacis@gmail.com)
 
 Significant contributions from:
 Spencer V Muse (muse@stat.ncsu.edu)
 Simon DW Frost (sdf22@cam.ac.uk)
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
*/


 #define MIN(a,b) ((a)>(b)?(b):(a))
__kernel void LeafKernel(  __global float* node_cache,                 // argument 0
                           __global const float* model,                // argument 1
                           __global const float* nodRes_cache,         // argument 2
                           __global const long* nodFlag_cache,         // argument 3
                           long sites,                                 // argument 4
                           long characters,                            // argument 5
                           long childNodeIndex,                        // argument 6
                           long parentNodeIndex,                       // argument 7
                           long roundCharacters,                       // argument 8
                           int intTagState,                            // argument 9
                           int nodeID,                                 // argument 10
                           __global int* scalings,                     // argument 11
                           float scalar,                               // argument 12
                           float uFlowThresh                           // argument 13
                           )
{
   int gx = get_global_id(0); // pchar
   if (gx > characters) return;
   int gy = get_global_id(1); // site
   if (gy > sites) return;
   long parentCharacterIndex = parentNodeIndex*sites*roundCharacters + gy*roundCharacters + gx;
   float privateParentScratch = 1.0f;
   int scale = 0;
   if (intTagState == 1)
   {
       privateParentScratch = node_cache[parentCharacterIndex];
       scale = scalings[parentCharacterIndex];
   }
   long siteState = nodFlag_cache[childNodeIndex*sites + gy];
   privateParentScratch *= model[nodeID*roundCharacters*roundCharacters + siteState*roundCharacters + gx];
   if (gy < sites && gx < characters)
   {
       node_cache[parentCharacterIndex] = privateParentScratch;
       scalings[parentCharacterIndex] = scale;
   }
}
__kernel void AmbigKernel(     __global float* node_cache,                 // argument 0
                               __global const float* model,                // argument 1
                               __global const float* nodRes_cache,         // argument 2
                               __global const long* nodFlag_cache,         // argument 3
                               long sites,                                 // argument 4
                               long characters,                            // argument 5
                               long childNodeIndex,                        // argument 6
                               long parentNodeIndex,                       // argument 7
                               long roundCharacters,                       // argument 8
                               int intTagState,                            // argument 9
                               int nodeID,                                 // argument 10
                               __global int* scalings,                     // argument 11
                               float scalar,                               // argument 12
                               float uFlowThresh                           // argument 13
                               )
{
   // thread index
   int tx = get_local_id(0);
   int ty = get_local_id(1);
   // global index
   int gx = get_global_id(0);
   int gy = get_global_id(1);
   long parentCharacterIndex = parentNodeIndex*sites*roundCharacters + gy*roundCharacters + gx;
   float privateParentScratch = 1.0f;
   int scale = 0;
   if (intTagState == 1 && gy < sites && gx < characters)
   {
       privateParentScratch = node_cache[parentCharacterIndex];
       scale = scalings[parentCharacterIndex];
   }
   float sum = 0.f;
   float childSum = 0.f;
   int scaleScratch = 0;
   __local float childScratch[BLOCK_SIZE][BLOCK_SIZE];
   __local float modelScratch[BLOCK_SIZE][BLOCK_SIZE];
   int siteState = nodFlag_cache[childNodeIndex*sites + gy];
   int ambig = 0;
   if (siteState < 0)
   {
       ambig = 1;
       siteState = -siteState-1;
   }
   int cChar = 0;
   for (int charBlock = 0; charBlock < 64/BLOCK_SIZE; charBlock++)
   {
       if (ambig && gy < sites && gx < characters)
           childScratch[ty][tx] =  nodRes_cache[siteState*characters + (charBlock*BLOCK_SIZE) + tx];
           //childScratch[ty][tx] =
            //   nodRes_cache[siteState*characters + (charBlock*BLOCK_SIZE) + tx];
       else if (gy < sites && gx < characters)
       {
           if (charBlock*BLOCK_SIZE + tx == siteState)
               childScratch[ty][tx] = 1;
           else
               childScratch[ty][tx] = 0;
       }
       else
               childScratch[ty][tx] = 0;
       modelScratch[ty][tx] = model[nodeID*roundCharacters*roundCharacters + roundCharacters*((charBlock*BLOCK_SIZE)+ty) + gx];
       barrier(CLK_LOCAL_MEM_FENCE);
       for (int myChar = 0; myChar < MIN(BLOCK_SIZE, (characters-cChar)); myChar++)
       {
           sum += childScratch[ty][myChar] * modelScratch[myChar][tx];
           childSum += childScratch[ty][myChar];
       }
       barrier(CLK_LOCAL_MEM_FENCE);
       cChar += BLOCK_SIZE;
   }
   while (childSum < 1 && childSum != 0)
   {
       childSum *= scalar;
       sum *= scalar;
       scaleScratch++;
   }
   scale += scaleScratch;
   privateParentScratch *= sum;
   if (gy < sites && gx < characters)
   {
       scalings    [parentCharacterIndex]  = scale;
       node_cache  [parentCharacterIndex]  = privateParentScratch;
   }
   /*
   */
}
__kernel void InternalKernel(  __global float* node_cache,                 // argument 0
                               __global const float* model,                // argument 1
                               __global const float* nodRes_cache,         // argument 2
                               long sites,                                 // argument 3
                               long characters,                            // argument 4
                               long childNodeIndex,                        // argument 5
                               long parentNodeIndex,                       // argument 6
                               long roundCharacters,                       // argument 7
                               int intTagState,                            // argument 8
                               int nodeID,                                 // argument 9
                               __global float* root_cache,                 // argument 10
                               __global int* scalings,                     // argument 11
                               float scalar,                               // argument 12
                               float uFlowThresh,                          // argument 13
                               __global int* root_scalings                 // argument 10
                               )
{
   // thread index
   short tx = get_local_id(0);   //local pchar
   short ty = get_local_id(1);   //local site
   // global index
   short gx = get_global_id(0);
   int gy = get_global_id(1);
   long parentCharacterIndex = parentNodeIndex*sites*roundCharacters + gy*roundCharacters + gx;
   float privateParentScratch = 1.0f;
   short scale = 0;
   if (intTagState == 1 && gy < sites && gx < characters)
   {
       privateParentScratch = node_cache[parentCharacterIndex];
       scale = scalings[parentCharacterIndex];
   }
   float sum = 0.f;
   float childSum = 0.f;
   int scaleScratch = scalings[childNodeIndex*sites*roundCharacters + gy*roundCharacters + gx];
   __local float  childScratch[BLOCK_SIZE][BLOCK_SIZE];
   __local float  modelScratch[BLOCK_SIZE][BLOCK_SIZE];
   short cChar = 0;
   for (int charBlock = 0; charBlock < 64/BLOCK_SIZE; charBlock++)
   {
       childScratch[ty][tx] =
           node_cache[childNodeIndex*sites*roundCharacters + roundCharacters*gy + (charBlock*BLOCK_SIZE) + tx];
       modelScratch[ty][tx] = model[nodeID*roundCharacters*roundCharacters + roundCharacters*((charBlock*BLOCK_SIZE)+ty) + gx];
       barrier(CLK_LOCAL_MEM_FENCE);
       for (int myChar = 0; myChar < MIN(BLOCK_SIZE, (characters-cChar)); myChar++)
       {
           sum += childScratch[ty][myChar] * modelScratch[myChar][tx];
           childSum += childScratch[ty][myChar];
       }
       barrier(CLK_LOCAL_MEM_FENCE);
       cChar += BLOCK_SIZE;
   }
   while (childSum < 1 && childSum != 0)
   {
       childSum *= scalar;
       sum *= scalar;
       scaleScratch++;
   }
   scale += scaleScratch;
   privateParentScratch *= sum;
   if (gy < sites && gx < characters)
   {
       scalings     [parentCharacterIndex]  = scale;
       root_scalings[gy*roundCharacters+gx] = scale;
       node_cache   [parentCharacterIndex]  = privateParentScratch;
       root_cache   [gy*roundCharacters+gx] = privateParentScratch;
   }
}
__kernel void ResultKernel (   __global int* freq_cache,                   // argument 0
                               __global float* prob_cache,                 // argument 1
                               __global fpoint* result_cache,              // argument 2
                               __global float* root_cache,                 // argument 3
                               __global int* root_scalings,                // argument 4
                               long sites,                                 // argument 5
                               long roundCharacters,                       // argument 6
                               float scalar,                               // argument 7
                               long characters                             // argument 8
                           )
{
   // shrink the work group to sites, rather than sites x characters
   #ifdef __GPUResults__
   int site = get_global_id(0);
   int localSite = get_local_id(0);
   __local fpoint resultScratch[BLOCK_SIZE*BLOCK_SIZE];
   resultScratch[localSite] = 0.0;
   while (site < sites)
   {
       result_cache[site] = 0.0;
       fpoint acc = 0.0;
       int scale = root_scalings[site*roundCharacters];
       for (int rChar = 0; rChar < characters; rChar++)
       {
           acc += root_cache[site*roundCharacters + rChar] * prob_cache[rChar];
       }
       //resultScratch[localSite] += (native_log(acc)-scale*native_log(scalar)) * freq_cache[site];
       resultScratch[localSite] += (log(acc)-scale*log(scalar)) * freq_cache[site];
       //result_cache[site] += (log(acc)-scale*log(scalar)) * freq_cache[site];
       site += get_global_size(0);
   }
   barrier(CLK_LOCAL_MEM_FENCE);
   for (int offset = get_local_size(0)/2; offset > 0; offset >>= 1)
   {
       if (localSite < offset)
       {
           fpoint other = resultScratch[localSite + offset];
           fpoint mine  = resultScratch[localSite];
           resultScratch[localSite] = mine+other;
       }
       barrier(CLK_LOCAL_MEM_FENCE);
   }
   // TODO: this would probably be faster if I saved them further apart to reduce bank conflicts
   if (localSite == 0) result_cache[get_group_id(0)] = resultScratch[0];
   #else
   if (get_global_id(0) != 0) return;
   int site = get_global_id(1);
   result_cache[site] = 0.0;
   if (get_group_id(1) >= get_local_size(0)*get_local_size(1)) return;
   while (site < sites)
   {
       float acc = 0.0;
       int scale = root_scalings[site*roundCharacters];
       for (int rChar = 0; rChar < characters; rChar++)
       {
           acc += root_cache[site*roundCharacters + rChar] * prob_cache[rChar];
       }
       result_cache[site] += (native_log(acc)-scale*native_log(scalar)) * freq_cache[site];
       site += get_local_size(0)*get_local_size(1);
   }
   barrier(CLK_LOCAL_MEM_FENCE);
   #endif
   /*
   */
}
__kernel void ReductionKernel ( __global double* result_cache               // argument 1
                           )
{
   int groupNum = get_local_id(0);
   __local double resultScratch[BLOCK_SIZE*BLOCK_SIZE];
   resultScratch[groupNum] = result_cache[groupNum];
   barrier(CLK_LOCAL_MEM_FENCE);
   for (int offset = get_local_size(0)/2; offset > 0; offset >>= 1)
   {
       if (groupNum < offset)
       {
           double other = resultScratch[groupNum + offset];
           double mine  = resultScratch[groupNum];
           double sum  = mine + other;
				//if (offset > 1)										
           	resultScratch[groupNum]  = sum;
				//else										
			//{
           //	resultScratch[2]  = mine;
           //	resultScratch[3]  = other;
           //	resultScratch[4]  = (sum);
       	//}
       }
       barrier(CLK_LOCAL_MEM_FENCE);
   }
   result_cache[groupNum] = resultScratch[groupNum];
   //
   /*
   if (get_group_id(0) != 0) return;
   if (get_group_id(1) != 0) return;
   int groupDim = get_local_size(0)*get_local_size(1);
   int groupNum = get_local_id(1)*get_local_size(0)+get_local_id(0);
   __local float resultScratch[BLOCK_SIZE*BLOCK_SIZE];
   resultScratch[groupNum] = result_cache[groupNum];
   //result_cache[groupNum] = 0.f;
   barrier(CLK_LOCAL_MEM_FENCE);
   for (int offset = groupDim/2; offset > 0; offset >>= 1)
   {
   //int offset = groupDim/2;
       if (groupNum < offset)
       {
           float other = resultScratch[groupNum + offset];
           float mine  = resultScratch[groupNum];
				if (offset > 1)	
           	resultScratch[groupNum]  = mine + other;
				else		
			{
           	resultScratch[2]  = mine;
           	resultScratch[3]  = other;
       	}
           //resultScratch[groupNum + offset]  = 0.f;
       }
       barrier(CLK_LOCAL_MEM_FENCE);
   }
   //if (groupNum == 0) result_cache[0] = resultScratch[0];
   result_cache[groupNum] = resultScratch[groupNum];
   */
}
