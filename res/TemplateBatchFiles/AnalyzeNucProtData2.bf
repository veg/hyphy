_Genetic_Code = 0;
#include "simpleBootstrap.bf";

SetDialogPrompt ("Please specify a nucleotide or amino-acid data file:");

DataSet ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filteredData = CreateFilter (ds,1);

fprintf (stdout,"\n______________READ THE FOLLOWING DATA______________\n",ds);

SelectTemplateModel(filteredData);

_DO_TREE_REBALANCE_ = 1;

#include "queryTree.bf";

LikelihoodFunction lf = (filteredData,givenTree);

Optimize (res,lf);

params = res[1][1]+19;

fprintf (stdout, "\nLog(L)\t",res[1][0],"\nParameters:\t", params,"\nc-AIC Score = ", 
				  2(-res[1][0]+params*(filteredData.sites/(filteredData.sites-params-1))));
				  
branchLengths = BranchLength (givenTree,-1);

TL = 0;

for (k=Rows(branchLengths)*Columns(branchLengths)-1; k>=0; k=k-1)				  
{
	TL = TL + branchLengths[k];
}

fprintf (stdout, "\nTotal tree length: ", TL, "\n");

bName = TipName ( givenTree, 0);
ExecuteCommands ("GetInformation (rmx,givenTree."+bName+");");

SELECTED_CHART_DATA = {190,1};

k = 0;

for (h=0; h<20; h=h+1)
{
	for (v=h+1; v<20; v=v+1)
	{
		SELECTED_CHART_DATA[k] = rmx[h][v]/vectorOfFrequencies[v];
		k=k+1;
	}
}

SELECTED_CHART_DATA = SELECTED_CHART_DATA*(vectorOfFrequencies[9]/rmx[7][9]);


count = Rows (SELECTED_CHART_DATA);

sum  = 0;
sum2 = 0;
sum3 = 0;
sum4 = 0;

SELECTED_CHART_DATA = SELECTED_CHART_DATA%0;

data_min  = SELECTED_CHART_DATA [0];
data_max  = SELECTED_CHART_DATA [0];

if (count%2)
{
	median = SELECTED_CHART_DATA[count/2];
}
else
{
	counter = count/2-1;
	median = (SELECTED_CHART_DATA[counter]+SELECTED_CHART_DATA[counter+1])/2;
}

for (counter=0; counter < count; counter = counter+1)
{
	term = SELECTED_CHART_DATA [counter];
	sum  = sum+term;
	sum2 = sum2+term*term;
	sum3 = sum3+term^3;
	sum4 = sum4+term^4;
	if (term < data_min)
	{
		data_min = term;
	}
	else
	{
		if (term > data_max)
		{
			data_max = term;
		}
	}
}

counter = (sum2-sum*sum/count)/(count-1);

if (count > 3)
{
	sum4 = (-6*sum^4+12*count*sum^2*sum2-3(count-1)*count*sum2^2-4*count(count+1)*sum*sum3+count^2*(count+1)*sum4)/(count*(count-1)*(count-2)*(count-3));
	sum4 = sum4/counter^2+3;
}
else
{
	sum4 = 0;
}	

if (count > 2)
{
	sum3 = (2*sum^3-3*count*sum*sum2+count^2*sum3)/(count*(count-1)*(count-2))/Sqrt(counter^2);
}
else
{
	sum3 = 0;
}	

fprintf (stdout, "\nCount    = ", count,
				 "\nMean     = ", sum/count,
				 "\nMedian   = ", median,
				 "\nVariance = ", counter,
				 "\nStd.Dev  = ", Sqrt (counter),
				 "\nCOV      = ", Sqrt (counter)*count/sum,
				 "\nSum      = ", sum,
				 "\nSq. sum  = ", sum2,
				 "\nSkewness = ", sum3,
				 "\nKurtosis = ", sum4,
				 "\nMin      = ", data_min,
				 "\nMax      = ", data_max,"\n\n");


#include "categoryEcho.bf";
