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

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/


//------------------------------------------------------------------------

template <class array_data> void ptr_array<array_data >::add(array_data in){

        array_data *temp;

		length++;
        if (length > 1)
		{
			 temp = new array_data[length];
			 for (long i=0;i < length - 1; i++)
				temp[i] = data[i];
			 delete data;
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

        array_data *temp;
  
		length++;
        if (length > 1)
		{
			 temp = new array_data[length];
			 for (long i=1;i < (length); i++)
				temp[i] = data[i-1];
			 delete data;
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
      temp = new array_data [length];
      for (long i=0; i < index-1 ; i++) temp[i] = data[i];
      for (long k=index-1; k < length; k++) temp[k] = data[k+1];
      delete data;
      data = temp;
  }
}
///____________________________________________________________
