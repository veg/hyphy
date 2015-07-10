//Built-in functions

function operators()
{
    //! (operator)

    //Undefined for Matrix and Trees
    a = 0.0;
    assert(!a==1.0,"!a==1.0 failed");

    a = 5.7;
    assert(!a==0.0,"!a==0.0 failed");

    a = "/bin/bash";
    assert(!a==1.0,"!a==1.0 failed");

    a = "/nonexistent/filepath";
    assert(!a==0.0, "!\"/nonexistent/filepath\"==0.0 failed");
    
    
    //!= (operator)
    assert(1 != 2 == 1,"1 != 2 == 1 failed");
    assert("Hyphy" != "HyPhy" == 1,"\"Hyphy\" != \"HyPhy\" == 1 failed");
    assert(1 != "1" == 1,"1 != \"1\" == 1 failed");
    assert("1" != 1 == 1,"\"1\" != 1 == 1 failed");

    //$ (operator)
    //Matrix
    a = {{1,2,3}}${{1,2,3}};
    assert(a[1]==4,"a[1]==4 failed");

    //Numeric
    assert(29.2$5.9==5,"29.2$5.9==5 failed");

    //String  
    assert(("ATATA"$"ATAT")[1]==3,"(\"ATATA\"$\"ATAT\")[1]==3 failed");
    assert(("ATATA"$"CCC")[1]==-1,"(\"ATATA\"$\"CCC\")[1]==-1 failed");
    assert(("ATATA"$"(TA)T")[1]==3,"(\"ATATA\"$\"(TA)T\")[1]==3 failed");

    //% (operator)

    //Matrix  
    a = {
        {    0.685433846825,    0.937009058412}
        {    0.757032230905,    0.530296328368}
        {     0.81098658657,    0.116591965341}
        {    0.485633417379,     0.28659384099}
        {    0.277806169884,    0.913600717418}
        {    0.153395332897,    0.101711391029}
        {    0.358844396043,   0.0366885005582}
        {    0.847414346376,    0.498019525664}
        {   0.0848120991338,    0.148932764807}
        {    0.848594787961,    0.940901103649}
        }; 

    b = a % 0;
    //{
    //{   0.0848120991338,    0.148932764807}
    //{    0.153395332897,    0.101711391029}
    //{    0.277806169884,    0.913600717418}
    //{    0.358844396043,   0.0366885005582}
    //{    0.485633417379,     0.28659384099}
    //{    0.685433846825,    0.937009058412}
    //{    0.757032230905,    0.530296328368}
    //{     0.81098658657,    0.116591965341}
    //{    0.847414346376,    0.498019525664}
    //{    0.848594787961,    0.940901103649}
    //}

    assert(b[0][1]==0.148932764807,"b[0][1]==0.148932764807 failed");

    //Numeric     
    assert(29%5==4,"29%5==4 failed");

    //String  
    assert("Hyphy"%"HYPHY"==1,"\"Hyphy\"%\"HYPHY\"==1 failed");

    //&& (operator)
    //Numeric     
    a = 29&&5;
    assert(a==1,"29&&5==1 failed");

    //String  
    a = "hyphy"&&1;
    assert(a=="HYPHY","a==\"HYPHY\" failed");

    //* (operator)

    //Associative List    
    //list is now {"key":"value","key2":"value2", "key3":"value3"}. Returns 3.
    a = {"key":"value"};
    b = {"key2":"value2", "key3":"value3"};
    a*b;

    assert(a["key3"]=="value3","a[\"key3\"]==\"value3\" failed");

    //Numeric     
    assert(5*2==10,"5*2==10 failed");

    //Matrix*Matrix   
    a = {{2,2,4,2}};
    some_matrix = {{1}{3}{5}{7}};
    c = a*some_matrix;
    assert(c[0][0]==42);    

    //Matrix*Number   
    assert(({{2,2,4,2}}*2)[2]==8,"({{2,2,4,2}}*2)[2]==8 failed");

    //String*String   
    a = "juxta";
    a*"position";
    assert(a=="juxtaposition","a==\"juxtaposition\" failed");
    
    //String*Matrix   
    assert(("tattarrattat"*{{"a","r","t"}})[0]==2,"(\"tattarrattat\"*{{\"a\",\"r\",\"t\"}})[0]==2 failed");

    //+ (operator)
    //Associative Array   
    list = {"key":"value"};
    list2 = {"key2":"value2", "key3":"value3"};
    a = list + list2;
    list3 = list["1"];
    assert(list3["key3"] == "value3","list3[\"key3\"] == \"value3\" failed");

    //Matrix  
    assert(({{1,2,3}} + {{3,2,1}})[0]==4,"{({1,2,3}} + {{3,2,1}})[0]==4 failed");

    //Numeric     
    assert((21+21)==42,"(21+21)==42 failed");

    //String  
    assert("Juxta" + "position"=="Juxtaposition","\"Juxta\" + \"position\"==\"Juxtaposition\" failed");

    //Topology/Tree   
    Topology T = ((a_node,f_node)N1,c_node,d_node);
    t = T + {"NAME":"e_node", "WHERE": "f_node", "PARENT": "g_node"};
    Topology expect = ((a_node,(f_node,e_node)g_node)N1,c_node,d_node);
    assert(t==0,"T==expect failed");

    //- (operator)

    //Associative List    
    list = {"key":"value", "key2":"value2"};
    list2 = {"key":"value"};
    list - list2;
    assert(list["key2"]=="value2","list[\"key2\"]==\"value2\" failed");

    //Matrix  
    assert(({{2,2,4,2}}-{{1,1,1,1}})[0] == 1,"({{2,2,4,2}}-{{1,1,1,1}})[0] == 1 failed");

    //Numeric     
    assert(5-3==2,"5-3==2 failed");

    /// (operator)
    //Numeric     
    assert(6/3==2,"6/3==2 failed");

    //String  
    assert("string"/"strin*"==1,"\"string\"/\"strin*\"==1 failed");

    //Greater Than

    //Matrix  
/*
 *    a = {{0,7,11,14}{7,0,6,9}{11,6,0,7}{14,9,7,0}}>1;
 *
 *    a1 = a>1;
 *    a2 = {
 *            {                 4,                 0,                 1}
 *            {                 4,                 0,                 1}
 *            {                 5,                 0,                 1}
 *            {                 5,            0.3125,                 1}
 *            {                 5,                 0,                 3}
 *            {                -1,                 0,                 6}
 *            {                 0,                 0,                 0}
 *            {                 0,                 0,                 0}
 *            {                 0,                 0,                 0}
 *            {                 0,                 0,                 0}
 *        };
 *
 *    assert(a1[0][1] == a2[0][1],"a1[0][1] == a2[0][1] failed");
 *
 */
    //Numeric     
    assert(5>4==1,"5>4==1 failed");

    //String  
    assert("Bears" > "Battlestar Galactica"==1,"\"Bears\" > \"Battlestar Galactica\"==1 failed");

    //Greater Than or Equal (operator)

    //Matrix  
/*
 *    a = {{0,7,11,14}{7,0,6,9}{11,6,0,7}{14,9,7,0}};
 *    a1 = a>=3;
 *
 *
 *    c = {
 *        {                 4,                 0,                 7,                 0,                -1}
 *        {                 0,                 0,                 0,                 0,                 0}
 *        {                 0,                 0,                 0,                 0,                 0}
 *        {                 0,                 0,                 0,                21,                 0}
 *        {                 0,                 0,                 0,                 0,                 0}
 *        {                 0,                 0,                 0,                 0,                 0}
 *        {                 0,                 0,                 0,                 0,                 0}
 *        {                 0,                 0,                 0,                 0,                 0}
 *        };
 *
 *    assert(a1[0][2]==c[0][2],"b[0][2],c[0][2] failed");
 *
 */
    //Numeric     
    assert(5>=4==1,"5>=4==1 failed");

    //String  
    assert("Bears" >= "Battlestar Galactica"==1,"\"Bears\" >= \"Battlestar Galactica\"==1 failed");

    //Less Than (operator)
    //Numeric     
    assert(4<5==1,"4<5==1 failed");

    //String  
    assert("Battlestar Galactica"<"Bears"==1,"\"Battlestar Galactica\"<\"Bears\"==1 failed");

    //Less Than or Equal (Operator)
    //Numeric     
    assert(4<=5==1,"4<=5==1 failed");

    //String  
    assert("Battlestar Galactica"<="Bears"==1,"\"Battlestar Galactica\"<=\"Bears\"==1 failed");

    //Topology    
    Topology T1 = ((a,b_node)N1,c,d,((g,h)N3,e,f)N2);
    Topology T2 = ((a,b_node)N1,c,d,(e,f)N2);
    assert(T2<=T1 == 1,"T2<=T1 == 1 failed");


    //Or (operator)

    //Numeric     
    assert(.5||1==1,".5||1==1 failed");

    //String  
    a = {
    {                 0}
    {                 1}
    {                 6}
    {                 7}
    };

    assert(("stringstring"||"st")[0][2]==a[0][2],"(\"stringstring\"||\"st\")[0][2]==a[0][2] failed");
        
    //^ (operator)
    //Matrix  
    m =
    {
    {                 1,                 2}
    {                 3,                 4}
    };

    assert(m^2>=(-6.73149-.01),"m^2==-6.73149 failed");

    //Numeric     
    assert(2^6>=64-1,"2^6==64 failed");

    return 1;
}

function tree_specific() 
{
    //BranchCount
    ACCEPT_ROOTED_TREES = 0;
    Tree T = ((a,b),(c,d));
    unrooted = BranchCount (T);
    assert(unrooted == 2,"unrooted == 2 failed);");
     
    ACCEPT_ROOTED_TREES = 1;
    Tree T = ((a,b),(c,d));
    rooted = BranchCount (T);
    assert(rooted == 1,"rooted == 1 failed);");

    //BranchLength
    Topology T = ((a:0.1,b:0.2):0.4,c:0.15,d:0.33);
    assert(BranchLength(T,1)==0.2,"BranchLength(T,1)==0.2 failed);");

    //BranchName
    Topology T = (((a:0.1,b:0.2)ab:0.4,e:0.1):0.2,c:0.15,d:0.33);
    assert(BranchName(T,1)=="Node1","BranchName(T,1)==\"Node1\" failed);");

    //TipCount
    Tree T1 = ((a,b)N1,c,d(N4),((g,h)N3,e,f)N2);
    assert(TipCount(T1)==8,"TipCount(T1)==8 failed);");

    //TipName
    Tree T1 = ((a,b)N1,c,d(N4),((g,h)N3,e,f)N2);
    assert(TipName(T1,1)=="b","TipName(T1,1)==\"b\" failed);");
}


function hbl_functions()
{
    //Absolute
    //Associative Array   
    list={"key":"value", "key2":"value2"};
    assert(Abs(list)==2,"Abs(list)==2 failed);");

    //Matrix (vector)     
    assert(Format(Abs({{1,2,3}}),0,5)=="3.74166","Abs({{1,2,3}})==3.74166 failed);");

    //Matrix (not a vector)   
    assert(Abs({{1,2}{3,-4}})==4,"Abs({{1,2}{3,-4}})==4 failed);");

    //Numeric     
    assert(Abs(-7.5)==7.5,"Abs(-7.5)==7.5 failed);");

    //String  
    assert(Abs("Cornholio")==9,"Abs(\"Cornholio\")==9 failed);");

    //Topology/Tree   
    Topology T = ((a,b)N1,c,d);
    m = Abs(T);

    //{{2,2,5, 5,5, -1}}
    assert(m[0][0]==2,"m[0][0]==2 failed);");

    //ArcTan
    assert(Format((4*Arctan(1)),0,5)=="3.14159","(4*Arctan(1))==3.14159 failed);");

    //Beta
    assert(Format(Beta(2,2),0,6)=="0.166667","Beta(2,2)==0.166667 failed);");

    //CChi2

    //Matrix  
    assert(Format(CChi2({{1,2}{3,0}},5),0,1) == "0.4","CChi2({{1,2}{3,0}},5) == 0.4 failed);");

    //Numeric     
    assert(Format(CChi2(1.44,1),0,6) == "0.769861","CChi2(1.44,1) == 0.769861 failed);");

    //CGammaDist
    assert(Format(CGammaDist(4,20,5),0,6) == "0.529743","CGammaDist(4,20,5) == 0.529743 failed);");

    //Columns
    assert( Columns({{1,2,3,4}{1,2,3,4}{1,2,3,4}}) == 4 ," Columns({{1,2,3,4}{1,2,3,4}{1,2,3,4}}) == 4  failed);");

    //Cos
    assert(Cos(0)==1,"Cos(0)==1 failed);");

    //Eigensystem
    a=Eigensystem({{19,3}{-2,26}});
    //assert(b[0]==20,"a[\"1\"][0]==20 failed);");

    //Erf
    assert(Format(Erf(.75),0,6)=="0.711156","Erf(.75)==0.711156 failed);");

    //Eval
    assert(Eval("3+3*13")==42,"Eval(\"3+3*13\")==42 failed);");

    //Exp
    assert(Format((Exp({{1,2}{2,1}}))[0][0],0,5)=="10.22671","(Exp({{1,2}{2,1}}))[0][0]==10.2267081822 failed);");
    assert(Format(Exp(1),0,5)=="2.71828","Exp(1)==2.71828 failed);");
    assert(Exp("1001111011000010")==6,"Exp(\"1001111011000010\")==6 failed);");

    //Format
    assert(Format(5,0,5)=="5.00000","Format(5,10,5)==5.00000 failed);");
    assert(Format("5",0,5)=="5.00000","Format(\"5\",10,5)==5.00000 failed);");

    Topology T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);
    assert(Format(T1,1,10) == "((a:-1,b:-1)N1:-1,c:-1,d:-1,((g:-1,h:-1)N3:-1,e:-1,f:-1)N2:-1)", "Line 368 failed");

    //Gamma
    assert(Format(Gamma(4),0,0)=="6","Gamma(4)==6 failed);");

    //GammaDist
    assert(Format(GammaDist(2,1,2),0,7)=="0.0366313","GammaDist(2,1,2)==0.0366313 failed);");

    //IBeta
    assert(Format(IBeta(.5,.5,5),0,5)=="0.98988","IBeta(.5,.5,5)==0.98988 failed);");

    //IGamma
    assert(Format(IGamma(3,4),0,5)=="0.76190","IGamma(3,4)==0.761897 failed);");

    //Inverse
    assert(Inverse("Dracula")=="alucarD","Inverse(\"Dracula\")==\"alucarD\"");
    assert((Inverse({{1,3,3}{1,4,3}{1,3,4}}))[0]==7,"(Inverse({{1,3,3}{1,4,3}{1,3,4}}))[0]==7 failed);");

    //Join
    list={"key":"value", "key2":"value2"};
    assert((Join(",", list))=="value,value2","(Join(\",\", list))==\"value,value2\" failed);");

    assert((Join(",",{{1,2,3}{4,5,6}{7,8,9}}))=="1,2,3,4,5,6,7,8,9","(Join(\",\",{{1,2,3}{4,5,6}{7,8,9}}))==\"1,2,3,4,5,6,7,8,9\" failed);");

    //LnGamma
    assert(Format(LnGamma(200),0,3)=="857.934","LnGamma(200)==857.934 failed);");

    //Log
    m = {
        {               1.2,               5.6,                 7}
        {                 3,                 4,                12}
        {                12,              3.23,                 8}
        };

    assert(Format((Log(m))[0][0],0,5)=="0.18232","(Log(m))[0][0]==0.182321556794 failed);");
    assert(Format(Log(2.71828183),0,0)=="1","Log(2.71828183)==1 failed);");
    assert(Log("hyphy")==109707827,"Log(\"hyphy\")==109707827 failed);");

    //LUDecompose
    assert((LUDecompose({{1.2,5.6,7}{3, 4,12}{12,3.23,8}}))[0][0]==12,"(LUDecompose({{1.2,5.6,7}{3, 4,12}{12,3.23,8}}))[0][0]==12 failed);");

    //LUSolve
    lu = LUDecompose({{1.2,5.6,7}{3, 4,12}{12,3.23,8}});
    b = {{1}{2}{3}};
    LUSolve(lu,b);

    //Max
    assert((Max({{1,2,3}{4,5,6}{7,8,9}},1))[0][0]==9,"Max{{1,2,3}{4,5,6}{7,8,9}},1)=={{9,8}} failed);");
    assert(Max({{1,2,3}{4,5,6}{7,8,9}},5)==9,"Max{{1,2,3}{4,5,6}{7,8,9}},5)==9 failed);");
    assert(Max(5,15)==15,"Max(5,15)==15 failed);");

    //Min
    assert((Min({{1,2,3}{4,5,6}{7,8,9}},1))[0][0]==1,"Min{{1,2,3}{4,5,6}{7,8,9}},1)=={{1,0}} failed);");
    assert(Min({{1,2,3}{4,5,6}{7,8,9}},5)==1,"Min{{1,2,3}{4,5,6}{7,8,9}},5)==1 failed);");
    assert(Min(5,15)==5,"Min(5,15)==5 failed);");

    //PSTreeString
    Tree T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);
    t = PSTreeString (T1,"STRING_SUPPLIED_LENGTHS",{{-1,-1}});

    //Random
    Random(3,15);
    Random({{1,2,3}{4,5,6}{7,8,9}},2);
    Random({{1,2}{3,4}},"LHS");

    mean = {{1,1,1}};
    cov = {{1,0}{0,1}};
    a = {"PDF":"Dirichlet","ARG0":cov};
    Random(mean,a);

    //RerootTree
    RerootTree("((a,b)N1,c,d)", "string");

    Tree T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);
    str = "N3";
    t = RerootTree(T1,str);

    //Rows
    assert(Rows({{1,2}{3,-4}})==2,"Rows({{1,2}{3,-4}})==2 failed);");

    //Simplex
    m = {
        {                 1,                 2,                 3,                 4,                 0,                 0,                 0}
        {                 0,                 3,                 2,                 1,                 1,                 0,                10}
        {                 0,                 2,                 5,                 3,                 0,                 1,                15}
        };

    Simplex(m);

    //Sin
    assert(Format(Sin(3.1415297/6),0,5)=="0.49999","Sin(3.1415297/6)==0.499991 failed);");

    //Sqrt
    assert(Format(Sqrt(2),0,5)=="1.41421","Sqrt(2)==1.41421 failed);");

    //Tan
    assert(Format(Tan(3.1415297/4),0,5)=="0.99997","Tan(3.1415297/4)==0.999969 failed);");

    //Time
    assert(Time(0),"Time(0) failed);");
    assert(Time(1),"Time(1) failed);");

    //Transpose
    m =
        {
            {                 1,                 2,                 3}
            {                 1,                 2,                 3}
        };

    Transpose(m);

    //Type
    assert(Type(Exp(1))=="Number","Type(Exp(1))==\"Number\" failed);");

    M = {{1,2}};
    assert(Type(M) == "Matrix","Type(M) == \"Matrix\" failed);");

    A = {"0":"A", "A":"0"};
    assert(Type(A)=="AssociativeList","Type(A)==\"AssociativeList\" failed);");

    Topology Top = (1,2,3);
    assert(Type(Top) == "Topology","Type(Top) == \"Topology\" failed);");

    Tree T = (1,2,3);
    assert(Type(T)=="Tree","Type(T)==\"Tree\" failed);");

    assert(Type("Hello world")=="String","Type(\"Hello world\")==\"String\" failed);");

    //ZCDF
    assert(Format(ZCDF(1),0,5)=="0.84134","ZCDF(1)==0.841345 failed);");

    return 1;
}

//HBL Commands
//AlignSequences
//Assert
//BGM
//Break
//Category
//ChoiceList
//ClearConstraints
//ConstructCategoryMatrix
//Continue
//CovarianceMatrix
//DataSetFilter
//DeleteObject
//Differentiate
//Do
//DoSQL
//ExecuteAFile
//ExecuteCommands
//Export
//Ffunction
//FindRoot
//For
//Fprintf
//Fscanf
//Function
//GetDataInfo
//GetInformation
//GetNeutralNull
//GetString
//GetURL
//HarvestFrequencies
//If
//Import
//Integrate
//LFCompute
//LikelihoodFunction
//LikelihoodFunction3
//LoadFunctionLibrary
//Model
//MolecularClock
//MPIReceive
//MPISend
//OpenDataPanel
//OpenWindow
//Optimize
//ReplicateConstraint
//RequireVersion
//Return
//SCFG
//SelectTemplateModel
//SetDialogPrompt
//SetParameter
//SimulateDataSet
//SpawnLikelihoodFunction
//StateCounter
//UseModel
//While

/* execution stub */

//fprintf    (stdout, "[Running COVERAGE TEST '", getTestName(), "']\n");
result1 = operators();
result2 = hbl_functions();
if (result1 && result2)
{
	fprintf (stdout, "[TEST PASSED]\n");
}
else
{
	fprintf (stdout, "[TEST FAILED]\n");
}

