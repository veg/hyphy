/*--------------------------------------------------------*/

function determineCoordinateTicks (x1,x2) 
{
	_range 	   = x2-x1;
	_log10Tick = Log(_range)/Log(10) /* round to the next smallest integer */
				 $1;
				 
	_log10Tick = Exp(_log10Tick*Log(10));
	if (_range/_log10Tick < 4)
	{
		_log10Tick = _log10Tick / (((4*_log10Tick/_range)+0.5)$1);
	}
	return _log10Tick;
}


/*--------------------------------------------------------*/

function boost (value, direction, factor)
{
	if (direction)
	{
		if (value > 0)
		{
			return value*factor;
		}
		else
		{
			return value/factor;
		}
	}
	
	if (value > 0)
	{
		return value/factor;
	}
	else
	{
		return value*factor;
	}
} 


/*--------------------------------------------------------*/

ExecuteAFile  ( "PostScript.bf");


function ScatterPlot		 (xy&, 			/* Nx2 matrix with x,y,value points to plot */
							  xyranges, 	/* 2x2 matrix {{x_min, x_max}{y_min, y_max} 
							  				   will be adjusted to cover the data in xy if needed*/
							  fontFace, 	/* font to use for plotting */
							  plotDim, 		/* 1x3 matrix {{width, height,font_size}} of the plot in points 
							  				   if embedLabels is 1, then this must be a 1x4 matrix; last entry is the font size 
							  				   for the embedded labels */
							  colors, 		/* Nx3 matrix of RGB colors to plot each point with */
							  shapes, 		/* Nx1 matrix of shapes to plot for each point */
							  labels,  		/* 1x3 matrix of strings: plot-label, x-axis label, y-axis label*/
							  pointLabels	/* Nx1 matrix of strings with labels for every point */,
							  embedLabels   /* whether or not to plot points (0) or text labels (1) */,
							  centroid		/* 2x1 point of the centroid */,
							  hull			/* Kx1 list of points (indices in 'xy' which must be traversed 
							  				   counter clock-wise to obtain the convex hull */,
							  doWrappers    /* should PS prefix and suffix be included */
							  )
{
	
	
	psDensityPlot = ""; psDensityPlot*1024;
	
	plotHeight = Max (100, plotDim[1]);
	plotWidth  = Max (100, plotDim[0]);
	
	plotOriginX = 4.5*plotDim[2];
	plotOriginY = 3.5*plotDim[2];
	
	xMin		= xyranges[0][0];
	xMax		= xyranges[0][1];
	yMin		= xyranges[1][0];
	yMax		= xyranges[1][1];
	
	plotSpanX   = plotWidth + 5*plotDim[2];
	plotSpanY	= plotHeight + 4*plotDim[2];
	
	if (doWrappers)
	{
		psDensityPlot * _HYPSPageHeader (plotSpanX,plotSpanY, "Density Plot");
		psDensityPlot * "\n";
		psDensityPlot * _HYPSTextCommands(0);
	}
	
	psDensityPlot * _HYPSSetFont (fontFace, plotDim[2]);
	psDensityPlot * "\n";
	
	psDensityPlot * "\n 1 setlinewidth 1 setlinecap 0 setlinejoin 0 0 0 setrgbcolor";
	psDensityPlot * ("\n " + plotOriginX + " " + plotOriginY + " " + plotWidth + " " + plotHeight + " rectstroke\n");
	
	/* adjust data ranges if necessary */
	
	_x = Rows (xy);

	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		xMin = Min(xMin,xy[_dataPoint][0]);
		xMax = Max(xMax,xy[_dataPoint][0]);
		yMin = Min(yMin,xy[_dataPoint][1]);
		yMax = Max(yMax,xy[_dataPoint][1]);
	}
	
	xMin = boost(xMin,0,1.1);
	yMin = boost(yMin,0,1.1);
	xMax = boost(xMax,1,1.1);
	yMax = boost(yMax,1,1.1);
	
	diff = (yMax-yMin) - (xMax-xMin);
	if (diff > 0)
	{
		xMin = xMin - diff/2;
		xMax = xMax + diff/2;
	}
	else
	{
		yMin = yMin + diff/2;
		yMax = yMax - diff/2;
	}
	
	px = plotWidth /(xMax-xMin);
	py = plotHeight/(yMax-yMin);
	
	
	_hullSize = Abs(hull);
	if (_hullSize>1)
	{
		psDensityPlot * ("0.5 0.5 0.5 setrgbcolor\n[3] 0 setdash\n");
		for (_dataPoint = 1; _dataPoint < _hullSize; _dataPoint = _dataPoint + 1)
		{
			psDensityPlot * ("newpath " + (plotOriginX+(xy[hull[_dataPoint-1]][0]-xMin)*px) + " " 
										+ (plotOriginY+(xy[hull[_dataPoint-1]][1]-yMin)*py) + " moveto "
										+ (plotOriginX+(xy[hull[_dataPoint]][0]-xMin)*px) + " " 
										+ (plotOriginY+(xy[hull[_dataPoint]][1]-yMin)*py) 
										+ " lineto stroke\n");
		
		}
		psDensityPlot * ("0 0 0 setrgbcolor\n[] 0 setdash\n");
	}

	if (embedLabels)
	{			
		psDensityPlot * _HYPSSetFont ("Helvetica", plotDim[3]);
		
		for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
		{
			psDensityPlot * (""+ colors[_dataPoint][0] + " " + colors[_dataPoint][1] + " " + colors[_dataPoint][2] + " setrgbcolor\n");

			aLabel = pointLabels[_dataPoint];
			aLabel = aLabel ^ {{"\\(","\("}};
			aLabel = aLabel ^ {{"\\)","\)"}};
			psDensityPlot * ("newpath " + (plotOriginX+(xy[_dataPoint][0]-xMin)*px) + " " 
										+ (plotOriginY+(xy[_dataPoint][1]-yMin)*py) + " moveto (" 
										+ aLabel + ") show\n");
			
			
		
		}

		psDensityPlot * _HYPSSetFont ("Times-Roman", plotDim[2]);
	}
	else
	{
		for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
		{
			psDensityPlot * (""+ colors[_dataPoint][0] + " " + colors[_dataPoint][1] + " " + colors[_dataPoint][2] + " setrgbcolor\n");
			myX_coord = plotOriginX+(xy[_dataPoint][0]-xMin)*px;
			myY_coord = plotOriginY+(xy[_dataPoint][1]-yMin)*py;
			psDensityPlot * ("newpath " + (myX_coord) + " " 
										+ (myY_coord) + " " 
										+ "3 0 360 arc fill\n");
										
		}	
	}
	

	if (Rows(centroid))
	{
		psDensityPlot * ("\n 0.75 0.75 0.75 setrgbcolor " + (plotOriginX+px*(centroid[0]-xMin)) + " " + (plotOriginX+py*(centroid[1]-yMin)) + " 8 0 360 arc fill\n");
		psDensityPlot * ("\n 0.25 0.25 0.25 setrgbcolor " + (plotOriginX+px*(centroid[0]-xMin)) + " " + (plotOriginX+py*(centroid[1]-yMin)) + " 8 0 360 arc stroke\n");
	}

	xscaler = determineCoordinateTicks (xMin,xMax);
	_x	= ((xMin/xscaler)$1)*xscaler;
	psDensityPlot * ("0 0 0 setrgbcolor\n");
	while (_x < xMax)
	{
		xStep = (plotOriginX + px*(_x-xMin));
		psDensityPlot * ("" +  xStep + " " + (2.5*plotDim[2]) + " (" + Format(_x,4,2) + ") centertext\n");  
		psDensityPlot * ("" +  xStep + " " + (plotOriginY+0.25*plotDim[2]) + " moveto 0 "
							+ (-0.25*plotDim[2]) +" rlineto stroke\n");  
		_x = _x + xscaler;
	}
	
	
	yscaler = determineCoordinateTicks (yMin,yMax);
	_y	= ((yMin/yscaler)$1)*yscaler;
	while (_y < yMax)
	{
		yStep = (plotOriginY + py*(_y-yMin));
		psDensityPlot * ("" +  (4*plotDim[2]) + " " + yStep + " (" + Format(_y,4,2) + ") righttext\n");  
		psDensityPlot * ("" +  plotOriginX    + " " + yStep + " moveto "+(0.25*plotDim[2]) +" 0 rlineto stroke\n");  
		_y = _y + yscaler;
	}

	psDensityPlot * ("" + (plotOriginX+plotWidth/2) + " " + (0.5*plotDim[2]) +" (" + labels[1] + ") centertext\n");
	psDensityPlot * ("" + (plotOriginY+plotHeight/2) + " " + (-1.5*plotDim[2]) +" ("+ labels[2] + ") vcentertext\n");	
	psDensityPlot * 0;
	
	if (doWrappers)
	{
		psDensityPlot * "\nshowpage\n";
	}
	return psDensityPlot;
}

/*--------------------------------------------------------------------------------------------------------------------------*/

function StackedBarPlot		 (xy&, 			/* x axis followed by K columns of y values*/
							  xyranges, 	/* 2x2 matrix {{x_min, x_max}{y_min, y_max} 
							  				   will be adjusted to cover the data in xy if needed*/
							  fontFace,		/* use this font */
							  plotDim, 		/* 1x3 matrix {{width, height,font_size}} of the plot in points 
							  			    */
							  colors, 		/* Kx3 matrix of RGB colors to plot each point with */
							  labels,  		/* 3x1 matrix of strings: plot-label, x-axis label, y-axis label*/
							  dataLabels,   /* Kx1 matrix of strings to label the observations with */
							  doWrappers,   /* should PS prefix and suffix be included */
							  lastLabelSP   /* use the last column to scatter plot over stacked bars */
							  )
{
	
	
	psDensityPlot = ""; psDensityPlot*1024;
	
	plotHeight = Max (100, plotDim[1]);
	plotWidth  = Max (100, plotDim[0]);
	
	plotOriginX = 4.5*plotDim[2];
	plotOriginY = 3.5*plotDim[2];
	
	xMin		= xyranges[0][0];
	xMax		= xyranges[0][1];
	yMin		= xyranges[1][0];
	yMax		= xyranges[1][1];
	
	
	_yColumns     		 = Columns (xy)-1-(lastLabelSP>0);

	legendWidth		   = 0;
	for (_dataPoint = 0; _dataPoint < _yColumns; _dataPoint = _dataPoint + 1)
	{
		px = _HYPSGetStringWidth (dataLabels[_dataPoint]) * plotDim[2];
		if (px > legendWitdh)
		{
			legendWitdh = px;
		}
	}
	
	if (legendWitdh)
	{
		legendWitdh = legendWitdh + 2*plotDim[2];
	}

	plotSpanX	= plotWidth;
	plotSpanY   = plotHeight;
	plotWidth   = plotWidth  - 5*plotDim[2]-legendWitdh;
	plotHeight	= plotHeight - 4*plotDim[2];
	
	if (doWrappers)
	{
		psDensityPlot * _HYPSPageHeader (plotSpanX,plotSpanY, "Stacked Bar Plot");
		psDensityPlot * "\n";
		psDensityPlot * _HYPSTextCommands(0);
	}
	
	psDensityPlot * _HYPSSetFont (fontFace, plotDim[2]);
	psDensityPlot * "\n";
	
	psDensityPlot * "\n 1 setlinewidth 1 setlinecap 0 setlinejoin 0 0 0 setrgbcolor";
	psDensityPlot * ("\n " + plotOriginX + " " + plotOriginY + " " + plotWidth + " " + plotHeight + " rectstroke\n");
	
	
	_x 		      		 = Rows (xy);
	_yTotalHeight 		 = {1,_x};
	
	barWidth			 = 1e100;
	_yIterator			 = (xy[-1][0])%0;
	for (_dataPoint = 1; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		barWidth = Min (barWidth,_yIterator[_dataPoint] - _yIterator[_dataPoint-1]);
	}

	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		xMin = Min(xMin,xy[_dataPoint][0]);
		xMax = Max(xMax,xy[_dataPoint][0]);
		for (_yIterator = 1; _yIterator <= _yColumns; _yIterator = _yIterator + 1)
		{	
			_yTotalHeight[_dataPoint] = _yTotalHeight[_dataPoint] + xy[_dataPoint][_yIterator];
		}
		yMin = Min(yMin,_yTotalHeight[_dataPoint]);
		yMax = Max(yMax,_yTotalHeight[_dataPoint]);
		if (lastLabelSP)
		{
			yMin = Min(yMin,xy[_dataPoint][_yIterator]);
			yMax = Max(yMax,xy[_dataPoint][_yIterator]);	
		}
	}
		
	/*diff = (yMax-yMin) - (xMax-xMin);
	if (diff > 0)
	{
		xMin = xMin - diff/2;
		xMax = xMax + diff/2;
	}
	else
	{
		yMin = yMin + diff/2;
		yMax = yMax - diff/2;
	}*/
	
	plotWidth 	= plotWidth  - 2;
	plotHeight  = plotHeight - 2;
	px 			= plotWidth /(xMax - xMin + barWidth);
	py 			= plotHeight/(yMax - yMin + 2);
	barWidth 	= barWidth * px;
	xShift		= 1;	
	
	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		myX_coord = plotOriginX+(xy[_dataPoint][0]-xMin)*px+xShift;
		myY_coord = plotOriginY;
		
		for (_yIterator = 0; _yIterator < _yColumns; _yIterator = _yIterator + 1)
		{	
			_rectHeight = xy[_dataPoint][_yIterator+1]*py;
			if (_rectHeight > 0)
			{
				psDensityPlot * (""+ colors[_yIterator][0] + " " + colors[_yIterator][1] + " " + colors[_yIterator][2] + " setrgbcolor\n");
				psDensityPlot * (""+ myX_coord + " " + myY_coord + " " + barWidth*1.5 + " " + _rectHeight + " rectfill\n");
				myY_coord = myY_coord + _rectHeight;
			}
		}									
	}	
	
	if (lastLabelSP)
	{
		for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
		{
			myX_coord = plotOriginX+(xy[_dataPoint][0]-xMin)*px+xShift;
			myY_coord = plotOriginY + xy[_dataPoint][_yIterator+1]*py;
			psDensityPlot * (""+ colors[_yIterator][0] + " " + colors[_yIterator][1] + " " + colors[_yIterator][2] + " setrgbcolor\n");
			psDensityPlot * ("newpath " + (myX_coord) + " " 
										+ (myY_coord) + " " 
										+ "1 0 360 arc fill\n");
		}
	}

	plotWidth 	= plotWidth  + 2;
	plotHeight  = plotHeight + 2;


	xscaler = determineCoordinateTicks (xMin,xMax);
	_x	= ((xMin/xscaler)$1)*xscaler;
	psDensityPlot * ("0 0 0 setrgbcolor\n");
	plottedZero = (_x == 0);
	
	while (_x < xMax)
	{
		xStep = (plotOriginX + px*(_x-xMin));
		psDensityPlot * ("" +  xStep + " " + (2.5*plotDim[2]) + " (" + Format(_x,0,0) + ") centertext\n");  
		psDensityPlot * ("" +  xStep + " " + (plotOriginY+0.25*plotDim[2]) + " moveto 0 "
							+ (-0.25*plotDim[2]) +" rlineto stroke\n");  
		_x = _x + xscaler;
	}
	
	
	yscaler = determineCoordinateTicks (yMin,yMax);
	_y	= ((yMin/yscaler)$1)*yscaler;
	if (plottedZero && (_y == 0))
	{
		_y = yscaler;
	}
	
	while (_y < yMax)
	{
		yStep = (plotOriginY + py*(_y-yMin));
		psDensityPlot * ("" +  (4*plotDim[2]) + " " + yStep + " (" + Format(_y,0,0) + ") righttext\n");  
		psDensityPlot * ("" +  plotOriginX    + " " + yStep + " moveto "+(0.25*plotDim[2]) +" 0 rlineto stroke\n");  
		_y = _y + yscaler;
	}

	psDensityPlot * ("" + (plotOriginX+plotWidth/2) + " " + (0.5*plotDim[2]) +" (" + labels[1] + ") centertext\n");
	psDensityPlot * ("" + (plotOriginY+plotHeight/2) + " " + (-1.5*plotDim[2]) +" ("+ labels[2] + ") vcentertext\n");	

	if (legendWitdh)
	{
		yLoc = plotOriginY + plotHeight - 1.5*plotDim[2];
		xLoc = plotOriginX +  _x * px + 0.5*plotDim[2];
		
		for (_segment = 0; _segment < _yColumns; _segment = _segment + 1)
		{
			psDensityPlot * ("" + colors[_segment][0] + " " + colors[_segment][1] + " " + colors[_segment][2] + " setrgbcolor\n");
			psDensityPlot * ("newpath " + xLoc + " " 
									+ yLoc + " " 
									+ plotDim[2] + " " 
									+ plotDim[2] + " " 
									+ " rectfill\n");
			psDensityPlot * ("0 0 0 setrgbcolor\n");
			psDensityPlot * ("newpath " + xLoc + " " 
									+ yLoc + " " 
									+ plotDim[2] + " " 
									+ plotDim[2] + " " 
									+ " rectstroke\n");
									
			psDensityPlot * ("newpath " + (xLoc + plotDim[2]*1.5) + " " + (yLoc+plotDim[2]/6) + " moveto (" + dataLabels[_segment] + ") show stroke\n");

			yLoc = yLoc - plotDim[2] * 1.5;
			
		}
	}

	psDensityPlot * 0;
	
	if (doWrappers)
	{
		psDensityPlot * "\nshowpage\n";
	}
	return psDensityPlot;
}

/*--------------------------------------------------------------------------------------------------------------------------*/

function CircleGraphPlot	  (edgeWeights&, 			/* a KxK matrix of graph edge weights */
							  fontFace,		/* use this font */
							  plotDim, 		/* 1x4 matrix {{width, height,font_size,maxlinewidth}} of the plot in points */
							  dataLabels,   /* Kx1 matrix of strings to label the nodes with */
							  colors,		/* Kx3 matrix of colors for the nodes  */
							  doWrappers    /* should PS prefix and suffix be included */
							  )
{
	
	
	psCircleGraphPlot = ""; 
	psCircleGraphPlot*1024;
	
	plotHeight = Max (100, plotDim[1]);
	plotWidth  = Max (100, plotDim[0]);
	
	plotOriginX = plotDim[2];
	plotOriginY = plotDim[2];
	
	plotSpanX	= plotWidth;
	plotSpanY   = plotHeight;
	plotWidth   = plotWidth  - 2*plotDim[2];
	plotHeight	= plotHeight - 2*plotDim[2];
	
	if (doWrappers)
	{
		psCircleGraphPlot * _HYPSPageHeader (plotSpanX,plotSpanY, "Circular graph plot");
		psCircleGraphPlot * "\n";
		psCircleGraphPlot * _HYPSTextCommands(0);
	}
	
	psCircleGraphPlot * _HYPSSetFont (fontFace, plotDim[2]);
	psCircleGraphPlot * "\n";
	
	psCircleGraphPlot * "\n 1 setlinewidth 1 setlinecap 0 setlinejoin 0 0 0 setrgbcolor";
	psCircleGraphPlot * ("\n " + plotOriginX + " " + plotOriginY + " " + plotWidth + " " + plotHeight + " rectstroke\n");
	
	_x 		      		 = Rows (edgeWeights);
	_yTotalHeight 		 = {1,_x};
	maxlinewidth		 = Max(plotDim[3],1);
	
	_nodeCenters		 = {_x,2};
	_anglePerNode		 = 2*3.1415926/_x;
	_radialCenter		 = {{plotWidth__/2+plotOriginX__,plotHeight__/2+plotOriginY__}
							{plotWidth__/2,plotHeight__/2}
						   };
	
	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		myX_coord = _radialCenter[0][0] + _radialCenter[1][0]*Cos(_anglePerNode*_dataPoint);
		myY_coord = _radialCenter[0][1] + _radialCenter[1][1]*Sin(_anglePerNode*_dataPoint);
		_nodeCenters[_dataPoint][0] = myX_coord;
		_nodeCenters[_dataPoint][1] = myY_coord;
	}

	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		for (_dataPoint2 = _dataPoint+1; _dataPoint2 < _x; _dataPoint2 = _dataPoint2 + 1)
		{
			myEW = edgeWeights[_dataPoint][_dataPoint2];
			if (myEW > 0.95)
			{
				psCircleGraphPlot * ("newpath " + _nodeCenters[_dataPoint][0] + " "
												+ _nodeCenters[_dataPoint][1] + " moveto " 
												+ myEW * maxlinewidth + " setlinewidth " 
												+ _nodeCenters[_dataPoint2][0] + " "
												+ _nodeCenters[_dataPoint2][1] + " lineto stroke\n");
			}								
		}
	}
		
	/*
	plotWidth 	= plotWidth  - 2;
	plotHeight  = plotHeight - 2;
	px 			= plotWidth /(xMax - xMin + barWidth);
	py 			= plotHeight/(yMax - yMin + 2);
	barWidth 	= barWidth * px;
	xShift		= 1;	
	
	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		myX_coord = plotOriginX+(xy[_dataPoint][0]-xMin)*px+xShift;
		myY_coord = plotOriginY;
		
		for (_yIterator = 0; _yIterator < _yColumns; _yIterator = _yIterator + 1)
		{	
			_rectHeight = xy[_dataPoint][_yIterator+1]*py;
			if (_rectHeight > 0)
			{
				psCircleGraphPlot * (""+ colors[_yIterator][0] + " " + colors[_yIterator][1] + " " + colors[_yIterator][2] + " setrgbcolor\n");
				psCircleGraphPlot * (""+ myX_coord + " " + myY_coord + " " + barWidth*1.5 + " " + _rectHeight + " rectfill\n");
				myY_coord = myY_coord + _rectHeight;
			}
		}									
	}	
	
	if (lastLabelSP)
	{
		for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
		{
			myX_coord = plotOriginX+(xy[_dataPoint][0]-xMin)*px+xShift;
			myY_coord = plotOriginY + xy[_dataPoint][_yIterator+1]*py;
			psCircleGraphPlot * (""+ colors[_yIterator][0] + " " + colors[_yIterator][1] + " " + colors[_yIterator][2] + " setrgbcolor\n");
			psCircleGraphPlot * ("newpath " + (myX_coord) + " " 
										+ (myY_coord) + " " 
										+ "1 0 360 arc fill\n");
		}
	}

	plotWidth 	= plotWidth  + 2;
	plotHeight  = plotHeight + 2;


	xscaler = determineCoordinateTicks (xMin,xMax);
	_x	= ((xMin/xscaler)$1)*xscaler;
	psCircleGraphPlot * ("0 0 0 setrgbcolor\n");
	plottedZero = (_x == 0);
	
	while (_x < xMax)
	{
		xStep = (plotOriginX + px*(_x-xMin));
		psCircleGraphPlot * ("" +  xStep + " " + (2.5*plotDim[2]) + " (" + Format(_x,0,0) + ") centertext\n");  
		psCircleGraphPlot * ("" +  xStep + " " + (plotOriginY+0.25*plotDim[2]) + " moveto 0 "
							+ (-0.25*plotDim[2]) +" rlineto stroke\n");  
		_x = _x + xscaler;
	}
	
	
	yscaler = determineCoordinateTicks (yMin,yMax);
	_y	= ((yMin/yscaler)$1)*yscaler;
	if (plottedZero && (_y == 0))
	{
		_y = yscaler;
	}
	
	while (_y < yMax)
	{
		yStep = (plotOriginY + py*(_y-yMin));
		psCircleGraphPlot * ("" +  (4*plotDim[2]) + " " + yStep + " (" + Format(_y,0,0) + ") righttext\n");  
		psCircleGraphPlot * ("" +  plotOriginX    + " " + yStep + " moveto "+(0.25*plotDim[2]) +" 0 rlineto stroke\n");  
		_y = _y + yscaler;
	}

	psCircleGraphPlot * ("" + (plotOriginX+plotWidth/2) + " " + (0.5*plotDim[2]) +" (" + labels[1] + ") centertext\n");
	psCircleGraphPlot * ("" + (plotOriginY+plotHeight/2) + " " + (-1.5*plotDim[2]) +" ("+ labels[2] + ") vcentertext\n");	

	if (legendWitdh)
	{
		yLoc = plotOriginY + plotHeight - 1.5*plotDim[2];
		xLoc = plotOriginX +  _x * px + 0.5*plotDim[2];
		
		for (_segment = 0; _segment < _yColumns; _segment = _segment + 1)
		{
			psCircleGraphPlot * ("" + colors[_segment][0] + " " + colors[_segment][1] + " " + colors[_segment][2] + " setrgbcolor\n");
			psCircleGraphPlot * ("newpath " + xLoc + " " 
									+ yLoc + " " 
									+ plotDim[2] + " " 
									+ plotDim[2] + " " 
									+ " rectfill\n");
			psCircleGraphPlot * ("0 0 0 setrgbcolor\n");
			psCircleGraphPlot * ("newpath " + xLoc + " " 
									+ yLoc + " " 
									+ plotDim[2] + " " 
									+ plotDim[2] + " " 
									+ " rectstroke\n");
									
			psCircleGraphPlot * ("newpath " + (xLoc + plotDim[2]*1.5) + " " + (yLoc+plotDim[2]/6) + " moveto (" + dataLabels[_segment] + ") show stroke\n");

			yLoc = yLoc - plotDim[2] * 1.5;
			
		}
	}*/

	psCircleGraphPlot * 0;
	
	if (doWrappers)
	{
		psCircleGraphPlot * "\nshowpage\n";
	}
	return psCircleGraphPlot;
}

/*--------------------------------------------------------------------------------------------------------------------------*/


function generateDensityPlot (data_matrix&, /* Nx3 matrix with x,y,value points to plot */
							  xyranges, 	/* 2x3 matrix {{x_min, x_max, steps}{y_min, y_max, steps} */
							  zrange,   	/* 1x2 matrix {{z_min, z_max}}, used for color scaling 
							  			   	z_min and z_max will be automatically adjusted to include
							  			   	at least the range of values from data_matrix*/
							  plotDim, 		/* 1x3 matrix {{width, height,font_size}} of the plot in points */
							  colors, 		/* 2x3 matrix {{R_base,G_base,B_base}{R_max,G_max,R_max
							  			 		the colors are linearly interpolated from base (min intensity)*/
							  labels,  		/* 1x3 matrix of strings: plot-label, x-axis label, y-axis label*/
							  circles	    /* Nx3 matrix (N could be 0) with coordinates/radii of circles to place on the map */
							  )
{
	
	psDensityPlot = ""; psDensityPlot*1024;
	
	plotHeight = Max (100, plotDim[1]);
	plotWidth  = Max (100, plotDim[0]);
	
	plotOriginX = 4.5*plotDim[2];
	plotOriginY = 3.5*plotDim[2];
	
	xBoxes		= xyranges[0][2];
	yBoxes		= xyranges[1][2];
	xStep		= (xyranges[0][1]-xyranges[0][0])/xBoxes;
	yStep		= (xyranges[1][1]-xyranges[1][0])/yBoxes;
	px			= plotWidth/xBoxes;
	py			= plotHeight/yBoxes;
	zMin		= 1e100;
	zMax		= -1e100;

	psDensityPlot * _HYPSPageHeader (plotWidth + 5*plotDim[2], plotHeight + 4*plotDim[2], "Density Plot");
	psDensityPlot * "\n";
	doPercentage  = Columns (plotDim) > 3;
	if (doPercentage)
	{
		secondarySize  = Min(Min (px,py)/2,plotDim[3]);
		secondaryShift = secondarySize/2;
		
		psDensityPlot * _HYPSSetFont ("Times-Roman", secondarySize);
	}
	else
	{
		psDensityPlot * _HYPSSetFont ("Times-Roman", plotDim[2]);
	}
	psDensityPlot * "\n";
	psDensityPlot * _HYPSTextCommands(0);
	
	psDensityPlot * "\n 1 setlinewidth 1 setlinecap 0 setlinejoin 0 0 0 setrgbcolor";
	psDensityPlot * ("\n " + plotOriginX + " " + plotOriginY + " " + plotWidth + " " + plotHeight + " rectstroke\n");
	
	zValues		= {xBoxes,yBoxes};
	
	/* compute and condense the data matrix */
	
	
	_x = Rows (data_matrix);
	
	totalSum = 0;
	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		xSquare = Max(0,Min(((data_matrix[_dataPoint][0]-xyranges[0][0])/xStep)$1,xBoxes-1));
		ySquare = Max(0,Min(((data_matrix[_dataPoint][1]-xyranges[1][0])/yStep)$1,yBoxes-1));
		zValues [xSquare][ySquare] = zValues [xSquare][ySquare] + data_matrix[_dataPoint][2];
		zMin = Min(zMin,zValues [xSquare][ySquare]);
		zMax = Max(zMax,zValues [xSquare][ySquare]);
		totalSum = totalSum + data_matrix[_dataPoint][2];
	}
	
	zMin = Min(zMin,zrange[0]);
	zMax = Max(zMax,zrange[1])-zMin;
	
	for (_x = 0; _x < xBoxes; _x = _x+1)
	{
		for (_y = 0; _y < yBoxes; _y = _y+1)
		{
			meColor = Sqrt((zValues [_x][_y]-zMin)/zMax);
			meColor2 = colors[0][-1]*(1-meColor) + colors[1][-1]*meColor;
			psDensityPlot * (""+ meColor2[0] + " " + meColor2[1] + " " + meColor2[2] + " setrgbcolor\n");
			psDensityPlot * ("" + (plotOriginX+_x*px) + " " + (plotOriginY+_y*py) + " " + px + " " + py + " rectfill\n");
			if (doPercentage && zValues [_x][_y]/totalSum > 0.0005)
			{
				if (meColor > 0.3)
				{
					psDensityPlot * ("1 1 1 setrgbcolor\n");
				}
				else
				{
					psDensityPlot * ("0 0 0 setrgbcolor\n");
				}
				psDensityPlot * ("" + (plotOriginX+_x*px+secondaryShift) + " " + (plotOriginY+_y*py+secondaryShift) + " moveto ("+ Format(100*zValues [_x][_y]/totalSum,3,1) + ") show\n");
			}
		}
	}
	
	/* now do the coordinates */
	/* determine base scale */
	
	if (doPercentage)
	{
		psDensityPlot * _HYPSSetFont ("Times-Roman", plotDim[2]);
		psDensityPlot * "\n";
	}
	
	xscaler = determineCoordinateTicks (xyranges[0][0],xyranges[0][1]);
	_x	= ((xyranges[0][0]/xscaler)$1 + 1)*xscaler;
	px	= plotWidth/(xyranges[0][1]-xyranges[0][0]);
	psDensityPlot * ("0 0 0 setrgbcolor\n");
	while (_x < xyranges[0][1])
	{
		xStep = (plotOriginX + px*(_x-xyranges[0][0]));
		psDensityPlot * ("" +  xStep + " " + (2.5*plotDim[2]) + " (" + Format(_x,4,2) + ") centertext\n");  
		psDensityPlot * ("" +  xStep + " " + (plotOriginY+0.25*plotDim[2]) + " moveto 0 "+ (-0.25*plotDim[2]) +" rlineto stroke\n");  
		_x = _x + xscaler;
	}
	
	yscaler = determineCoordinateTicks (xyranges[1][0],xyranges[1][1]);
	_y	= ((xyranges[1][0]/yscaler)$1 + 1)*yscaler;
	py	= plotHeight/(xyranges[1][1]-xyranges[1][0]);
	while (_y < xyranges[1][1])
	{
		yStep = (plotOriginY + py*(_y-xyranges[1][0]));
		psDensityPlot * ("" +  (4*plotDim[2]) + " " + yStep + " (" + Format(_y,4,2) + ") righttext\n");  
		psDensityPlot * ("" +  plotOriginX    + " " + yStep + " moveto "+(0.25*plotDim[2]) +" 0 rlineto stroke\n");  
		_y = _y + yscaler;
	}
	
	psDensityPlot * ("0 0 0 setrgbcolor\n"+plotOriginX+" "+plotOriginY+" moveto\n"+
					 plotWidth + " " + plotHeight + " rlineto stroke\n");

	psDensityPlot * ("" + (plotOriginX+plotWidth/2) + " " + (0.5*plotDim[2]) +" (" + labels[1] + ") centertext\n");
	psDensityPlot * ("" + (plotOriginY+plotHeight/2) + " " + (-1.5*plotDim[2]) +" ("+ labels[2] + ") vcentertext\n");
	
	_x = Rows (circles);
	
	for (_y = 0; _y < _x; _y = _y + 1)
	{
		xStep = circles[_y][0]*px + plotOriginX;
		yStep = circles[_y][1]*py + plotOriginY;
		psDensityPlot * ("newpath " + xStep + " " + yStep + " " + circles[_y][2] + " 0 360 arc stroke\n");
	}
	
	psDensityPlot * "\nshowpage\n";
	psDensityPlot * 0;
	
	return psDensityPlot;
}

/*--------------------------------------------------------------------------------------------------------------------------*/


function generateHeatMap	 (data_matrix&, /* Nx4 matrix with x,y,z,plot/or not (0, >0) points to plot */
							  xyranges, 	/* 2x2 matrix {{x_bins, x_left, x_steps}{y_bins, y_bottom, y_steps} */
							  zrange,   	/* 1x2 matrix {{z_min, z_max}}, used for color scaling 
							  			   	z_min and z_max will be automatically adjusted to include
							  			   	at least the range of values from data_matrix*/
							  plotDim, 		/* 1x3 matrix {{width, height,font_size}} of the plot in points */
							  colors, 		/* 2x3 matrix {{R_base,G_base,B_base}{R_max,G_max,R_max
							  			 		the colors are linearly interpolated from base (min intensity)*/
							  labels  		/* 1x3 matrix of strings: plot-label, x-axis label, y-axis label*/
							  )
{
	
	psDensityPlot = ""; psDensityPlot*1024;
	
	plotHeight = Max (100, plotDim[1]);
	plotWidth  = Max (100, plotDim[0]);
	
	plotOriginX = 4.5*plotDim[2];
	plotOriginY = 3.5*plotDim[2];
	
	xBoxes		= xyranges[0][0];
	yBoxes		= xyranges[1][0];
	xStep		= xyranges[0][2];
	yStep		= xyranges[1][2];
	px			= plotWidth/xBoxes;
	py			= plotHeight/yBoxes;
	zMin		= 1e100;
	zMax		= -1e100;

	psDensityPlot * _HYPSPageHeader (plotWidth + 5*plotDim[2], plotHeight + 6*plotDim[2], "Density Plot");
	psDensityPlot * "\n";
	psDensityPlot * _HYPSSetFont ("Times-Roman", plotDim[2]);

	psDensityPlot * "\n";
	psDensityPlot * _HYPSTextCommands(0);
	psDensityPlot * _HYPSHatchRectanglePattern (45, Min(10,px$2), px, py, {{0.5,0.5,0.5}}, "hatch",plotOriginX,plotOriginY);
	
	psDensityPlot * "\n 1 setlinewidth 1 setlinecap 0 setlinejoin 0 0 0 setrgbcolor\n";
	
	zValues		= {xBoxes,yBoxes};
	doPlot		= {xBoxes,yBoxes};
	
	/* compute and condense the data matrix */
	
	
	_x = Rows (data_matrix);
	
	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		xSquare = _dataPoint $ yBoxes;
		ySquare = _dataPoint % yBoxes;
		plotMe  = data_matrix[_dataPoint][3];
		if		(plotMe)
		{
			zValues [xSquare][ySquare]  = data_matrix[_dataPoint][2]/plotMe;
			zMin						= Min(zMin,zValues [xSquare][ySquare]);
			zMax						= Max(zMax,zValues [xSquare][ySquare]);
			doPlot[xSquare][ySquare]    = 1;
		}
	}
	
	zMin = Min(zMin,zrange[0]);
	zMax = Max(zMax,zrange[1])-zMin;
	
	for (_x = 0; _x < xBoxes; _x = _x+1)
	{
		for (_y = 0; _y < yBoxes; _y = _y+1)
		{
			if (doPlot[_x][_y])
			{
				meColor = (zValues [_x][_y]-zMin)/zMax;
				meColor2 = colors[0][-1]*(1-meColor) + colors[1][-1]*meColor;
				psDensityPlot * (""+ meColor2[0] + " " + meColor2[1] + " " + meColor2[2] + " setrgbcolor\n");
				psDensityPlot * ("" + Format(plotOriginX+_x*px,20,10) + " " + Format(plotOriginY+_y*py,20,10) + " " + Format(px,20,10) + " " + Format(py,20,10) + " rectfill\n");
			}
			else
			{
				psDensityPlot * ("/Pattern setcolorspace hatch setcolor " + Format(plotOriginX+_x*px,20,10) + " " + Format(plotOriginY+_y*py,20,10) + " " + Format(px,20,10) + " " + Format(py,20,10) + " rectfill\n");	
			}
		}
	}
	
	/* now do the coordinates */
	/* determine base scale */
	
	if (doPercentage)
	{
		psDensityPlot * _HYPSSetFont ("Times-Roman", plotDim[2]);
		psDensityPlot * "\n";
	}
	
	x_max = xyranges[0][1]+xyranges[0][0]*xyranges[0][2];
	xscaler = determineCoordinateTicks (xyranges[0][1],x_max);
	_x	= ((xyranges[0][1]/xscaler)$1 + 1)*xscaler;
	px	= plotWidth/(x_max-xyranges[0][1]);
	psDensityPlot * ("0 0 0 setrgbcolor\n");
	while (_x < x_max)
	{
		xStep = (plotOriginX + px*(_x-xyranges[0][1]));
		psDensityPlot * ("" +  xStep + " " + (2.5*plotDim[2]) + " (" + Format(_x,4,2) + ") centertext\n");  
		psDensityPlot * ("" +  xStep + " " + (plotOriginY+0.25*plotDim[2]) + " moveto 0 "+ (-0.25*plotDim[2]) +" rlineto stroke\n");  
		_x = _x + xscaler;
	}
	
	y_max = xyranges[1][1]+xyranges[1][0]*xyranges[1][2];
	yscaler = determineCoordinateTicks (xyranges[1][1],y_max);
	_y	= ((xyranges[1][1]/yscaler)$1 + 1)*yscaler;
	py	= plotHeight/(y_max-xyranges[1][1]);
	while (_y < y_max)
	{
		yStep = (plotOriginY + py*(_y-xyranges[1][1]));
		psDensityPlot * ("" +  (4*plotDim[2]) + " " + yStep + " (" + Format(_y,4,2) + ") righttext\n");  
		psDensityPlot * ("" +  plotOriginX    + " " + yStep + " moveto "+(0.25*plotDim[2]) +" 0 rlineto stroke\n");  
		_y = _y + yscaler;
	}
	
	psDensityPlot * ("" + (plotOriginX+plotWidth/2) + " " + (0.5*plotDim[2]) +" (" + labels[1] + ") centertext\n");
	psDensityPlot * ("" + (plotOriginY+plotHeight/2) + " " + (-1.5*plotDim[2]) +" ("+ labels[2] + ") vcentertext\n");
		
	psDensityPlot * ("\n0 0 0 setrgbcolor " + plotOriginX + " " + plotOriginY + " " + plotWidth + " " + plotHeight + " rectstroke\n");
	
	_thermLabels = {{"0%",".2"}{"25%",".2"}{"50%",".2"}{"75%",".2"}{"100%",".2"}};
	_colorLabels = {5,6};
	for (_x = 0; _x < 5; _x = _x+1)
	{
		meColor = _x*0.25;
		meColor2 = colors[0][-1]*(1-meColor) + colors[1][-1]*meColor;
		_colorLabels[_x][0] = meColor2[0];
		_colorLabels[_x][1] = meColor2[1];
		_colorLabels[_x][2] = meColor2[2];
		if (Max(_colorLabels[_x][0],Max(_colorLabels[_x][1],_colorLabels[_x][2])) < 0.4)
		{
			_colorLabels[_x][3] = 1;
			_colorLabels[_x][4] = 1;
			_colorLabels[_x][5] = 1;		
		}
	}
	
	psDensityPlot * ("" + plotOriginX + " " + (plotHeight+plotOriginY+plotDim[2]/2) + " translate \n"); 
	psDensityPlot * (_HYPSLabeledBoxes (200,1.5*plotDim[2],plotDim[2],_thermLabels,_colorLabels))["PS"];
	psDensityPlot * "\nshowpage\n";
	psDensityPlot * 0;
	
	return psDensityPlot;
}