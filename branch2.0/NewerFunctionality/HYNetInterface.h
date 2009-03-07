/*

	This is a file which declares NN interfacing with HBL model structures.
	
	September 2003, Sergei L Kosakovsky Pond

*/

#ifndef	__HYNNINTERFACE__
#define	__HYNNINTERFACE__

#include "hy_strings.h"
#include "batchlan.h"
#include "matrix.h"
#include "SerangNet.h"

/*---------------------------------------------------------*/

//class	 _MatrixNN: public _SimpleList 
//{

//	Net * 
//}

/*---------------------------------------------------------*/


void	 	 TrainModelNN 	(_String*, _String*);
void		 LoadModelNN  	(_String*, _String*);
_Matrix*	 ComputeModel	(long);

/*---------------------------------------------------------*/

void		 NNMatrixSampler(long, _Matrix&, _SimpleList&, _SimpleList&, _Matrix*, _List&, _List&);

extern	 _String	  		ModelTrainNNFlag,
					  		ModelLoadNNFlag;
					  	

#endif