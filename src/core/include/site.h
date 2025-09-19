/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
   Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
   Art FY Poon    (apoon42@uwo.ca)
   Steven Weaver (sweaver@temple.edu)

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

#ifndef _GENSITE_
#define _GENSITE_
//#pragma once

#include "avllist.h"
#include "avllistx.h"
#include "avllistxl.h"
#include "list.h"
#include "parser.h"
#include "simplelist.h"
#include "stdlib.h"
#include "translation_table.h"


/**
 * @brief A struct to store the state of a data set file
 */
struct FileState {

  /**
   * @brief The translation table for the file
   */
  _TranslationTable *translationTable;
  /**
   * @brief The current species being read
   */
  long curSpecies,
       /**
        * @brief The total number of species read
        */
       totalSpeciesRead,
       /**
        * @brief The total number of sites read
        */
       totalSitesRead,
       /**
        * @brief The total number of species expected
        */
       totalSpeciesExpected,
       /**
        * @brief The total number of sites expected
        */
       totalSitesExpected,
       /**
        * @brief The current site being read
        */
       curSite,
       /**
        * @brief The maximum string length
        */
       maxStringLength,
       /**
        * @brief The position in the source
        */
       pInSrc,
       /**
        * @brief The current line in the file
        */
       currentFileLine;
  /**
   * @brief Whether the file is accepting commands
   */
  bool acceptingCommands,
       /**
        * @brief Whether all species have been defined
        */
       allSpeciesDefined,
       /**
        * @brief Whether the file is interleaved
        */
       interleaved,
       /**
        * @brief Whether to auto-detect the file type
        */
       autoDetect,
       /**
        * @brief Whether to skip in NEXUS
        */
       isSkippingInNEXUS;
  /**
   * @brief The file type
   */
  int fileType,
      /**
       * @brief The base length
       */
      baseLength;
  /**
   * @brief The repeat character
   */
  char repeat,
       /**
        * @brief The skip character
        */
       skip;

  /**
   * @brief The source string
   */
  _String *theSource,
          /**
           * @brief The namespace
           */
          *theNamespace;
    
  /**
   * @brief The line buffer
   */
  _StringBuffer lineBuffer;
  /**
   * @brief The raw lines format
   */
  _SimpleList rawLinesFormat;

};
//_________________________________________________________

/**
 * @brief A class to represent a site in a sequence
 */
class _Site : public _StringBuffer {

public:

  /**
   * @brief Construct a new _Site object
   */
  _Site(void);

  /**
   * @brief Construct a new _Site object
   *
   * @param length The length of the site
   * @param ref_no The reference number
   */
  _Site(unsigned long length, long ref_no);

  /**
   * @brief Construct a new _Site object
   *
   * @param s The string to construct from
   */
  _Site(_String const &);

  /**
   * @brief Construct a new _Site object
   *
   * @param c The character to construct from
   */
  _Site(char);

  /**
   * @brief Construct a new _Site object
   *
   * @param l The long to construct from
   */
  _Site(long);

  /**
   * @brief Destroy the _Site object
   */
  virtual ~_Site(void);

  /**
   * @brief Complete the site
   */
  void Complete(void);

  /**
   * @brief Get the reference number
   *
   * @return long The reference number
   */
  long GetRefNo     (void) const { return refNo < 0L ? -refNo - 2L : refNo - 2L; }
  /**
   * @brief Check if the site is complete
   *
   * @return true if the site is complete, false otherwise
   */
  bool IsComplete   (void) const    { return refNo < 0L; }
  /**
   * @brief Set the reference number
   *
   * @param r The reference number
   */
  void SetRefNo     (long r) { refNo = -r - 2L; }
    
  /**
   * @brief The new operator for the _Site class
   *
   * @param size The size to allocate
   * @return void* A pointer to the allocated memory
   */
  void *      operator new       (size_t size);
  /**
   * @brief The delete operator for the _Site class
   *
   * @param p A pointer to the memory to deallocate
   */
  void        operator delete    (void * p);


private:
  long refNo; // if this site contains a reference to another one
  // if refNo is negative, then shows whether the definition of this datatype
  // has been completed
};

#endif
