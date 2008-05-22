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
		px = Abs (dataLabels[_dataPoint]) * plotDim[2];
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

