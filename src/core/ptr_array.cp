/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2006  
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)
						 
Some of the Original Code by William A Casey.

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


//------------------------------------------------------------------------

template <class array_data> void ptr_array<array_data >::add(array_data in){

		length++;
        if (length > 1)
		{
			 array_data *temp = new array_data[length];
			 for (long i=0;i < length - 1; i++)
				temp[i] = data[i];
			 delete [] data;
			 data = temp;
			 data[length-1] = in;
		}
        else
		{
			 data = new array_data[1];
			 data[0] = in;
		}
}

//------------------------------------------------------------------------

template <class array_data>	void ptr_array<array_data >::prepend(array_data in){
	
		length++;
        if (length > 1)
		{
			 array_data *temp = new array_data[length];
			 for (long i=1;i < (length); i++)
				temp[i] = data[i-1];
			 delete [] data;
			 data = temp;
			 data[0] = in;
		}
        else
		{
			 data = new array_data[1];
			 data[0] = in;
		}
}

//------------------------------------------------------------------------

template <class array_data>	void ptr_array<array_data>::delete_entry(int index){

  array_data *temp;
  if (length > 0){
	  length--;
      if (length) {
          temp = new array_data [length];
          for (long i=0; i < index-1 ; i++) temp[i] = data[i];
          for (long k=index-1; k < length; k++) temp[k] = data[k+1];
      } else {
        temp = nil;
      } 
      delete [] data;
      data = temp;
  }
}
///____________________________________________________________
