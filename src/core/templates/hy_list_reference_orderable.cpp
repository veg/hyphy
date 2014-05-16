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

/*
 Constructors
*/

template<typename PAYLOAD>
_hyListReference<PAYLOAD>::_hyListReference () : _hyListOrderable<PAYLOAD*> () {
}

template<typename PAYLOAD>
_hyListReference<PAYLOAD>::_hyListReference (const unsigned long items) : _hyListOrderable<PAYLOAD*> (items) {
  
}

template<typename PAYLOAD>
_hyListReference<PAYLOAD>::_hyListReference (const PAYLOAD& value) : _hyListOrderable<PAYLOAD*> (&value) {
  value.AddAReference ();
}

template<typename PAYLOAD>
_hyListReference<PAYLOAD>::_hyListReference (const _hyListReference<PAYLOAD> &l, const long from, const long to) : _hyListOrderable<PAYLOAD*> (l,from,to) {
  for (unsigned long k = 0; k < this->lLength; k++) {
    this->lData[k]->AddAReference();
      // 20140113: CHECK THAT THIS DOESN'T HAVE TO BE A makeDynamic CALL
  }
}

  // Data constructor (variable number of long constants)
template<typename PAYLOAD>
_hyListReference<PAYLOAD>::_hyListReference(const unsigned long number, const PAYLOAD* items[]) : 
_hyListOrderable<PAYLOAD> (number, items) {
}

template<typename PAYLOAD>
_hyListReference<PAYLOAD>::~_hyListReference (void) {
  for (unsigned long k = 0; k < this->lLength; k++) {
    if (this->lData[k]->CanFreeMe()) {
      delete this->lData[k];
    } else {
      this->lData[k]->RemoveAReference();
    }
  }
}

/*
 Storage/reference managament functions
 */

template<typename PAYLOAD>
void _hyListReference<PAYLOAD>::Clear (bool deallocate_memory) {
  for (unsigned long item = 0UL; item < this->lLength; item++) {
    if (this->lData[item]->CanFreeMe()) {
      delete this->lData[item];
    } else {
      this->lData[item]->RemoveAReference();
    }
  }
  this->_hyListOrderable<PAYLOAD*>::Clear(deallocate_memory);
}


/*
 Accessor/storage functions
 */

template<typename PAYLOAD>
void _hyListReference<PAYLOAD>::AppendNewInstance (PAYLOAD* the_ref) {
  this->append (the_ref);
}

/*
 Operator overloads
*/

  //Assignment operator
template<typename PAYLOAD>
_hyListReference<PAYLOAD> const _hyListReference<PAYLOAD>::operator=(const _hyListReference<PAYLOAD>& l)
{
   Clear();
   Clone (&l);
   return *this;
}

  //Append a list
template<typename PAYLOAD>
_hyListReference<PAYLOAD> const _hyListReference<PAYLOAD>::operator&(const _hyListReference<PAYLOAD>& l)
{
  _hyListReference <PAYLOAD> res (this->llength + l.lLength);
  
  res << (*this);
  res << l;
  return res;
}

  //Append a list
template<typename PAYLOAD>
_hyListReference<PAYLOAD> const _hyListReference<PAYLOAD>::operator&(const PAYLOAD* l)
{
  _hyListReference <PAYLOAD> res (this->llength + 1UL);
  
  res << (*this);
  (*res) && l;
  return res;
}

  //append an element
template<typename PAYLOAD>
void _hyListReference<PAYLOAD>::operator<<(PAYLOAD* item)
{
  this->append (item);
  item->AddAReference();
}

  //append the dynamic copy of item
template<typename PAYLOAD>
void _hyListReference<PAYLOAD>::operator&&(const PAYLOAD* item)
{
   this->append ((PAYLOAD*) item->makeDynamic());
}


  //append all elements in the list
template<typename PAYLOAD>
void _hyListReference<PAYLOAD>::operator<<(const _hyListReference<PAYLOAD>& l)
{
  this->RequestSpace (this->lLength + l.lLength);
  for (unsigned long k = 0UL; k < l.lLength; k++) {
    (*this) << l.lData[k];
  }
}

/*
 Element Comparison Functions
*/

template<typename PAYLOAD>
bool _hyListReference<PAYLOAD>::ItemEqualToValue(unsigned long index, const PAYLOAD* & value) const
{
  return this->lData[index]->Equal (value);
}

template<typename PAYLOAD>
long _hyListReference<PAYLOAD>::Compare(const long i, const long j) const {
  
  return this->lData[i]->Compare (this->lData[j]);
}

