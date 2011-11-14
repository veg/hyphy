steven = "Hi, Steven's test here\n";
fprintf (stdout, steven&&5);

//n=Abs(-7.5);
//fprintf (stdout, n);
//fprintf (stdout, "\n");
//s=Abs("Cornholio");
//fprintf (stdout, s);
//fprintf (stdout, "\n");
//v=Abs({{1,2,3}});
//fprintf (stdout, v);
//fprintf (stdout, "\n");
//m=Abs({{1,2}{3,-4}});
//fprintf (stdout, m);
//fprintf (stdout, "\n");
//list = {"key":"value", "key2":"value2"};
//l=Abs(list);
//fprintf (stdout, l);
//fprintf (stdout, "\n");
//
//Topology T = ((a,b)N1,c,d);
//t = Abs(T);
//fprintf (stdout, t);
//fprintf (stdout, "\n");


//Add
//a = 7;
//b = 8;
//c = a + b;
//fprintf (stdout, c);
//fprintf (stdout, "\n");
//
//s = "Juxta" + "position";
//fprintf (stdout, s);
//fprintf (stdout, "\n");
//
//v = {{1,2,3}} + {{3,2,1}};
//fprintf (stdout, v);
//fprintf (stdout, "\n");
//
//list = {"keyes":"value", "key2":"value2", "keynes":"kekeke", "five":"fold"} + {"key3":"value3", "kesdaynes":"kekeke" };
//fprintf (stdout, list);
//fprintf (stdout, "\n");
//
Topology T1 = ((a,b)N1,c,d);
////Topology T2 = ((e,f)N2,g,h);
//t = T1 + {"NAME":"e", "WHERE": "b", "PARENT":"f"};
//fprintf (stdout, T1);
//fprintf (stdout, "\n");

//Not
//Topology T1 = ((a,(e,(k,(m,(o,(q,(s,(u,(w,(y,(cc,(ee,ff)dd)z)x)v)t)r)p)n)l)f)b)N1,(c,(g,h)d)N2,(i,j)N3,(aa,bb)N4,(a1,(e1,(k1,(m1,(o1,(q1,(s1,(u1,(w1,(y1,(cc1,(ee1,ff1)dd1)z1)x1)v1)t1)r1)p1)n1)l1)f1)b1)N51);
//Topology T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);
//Topology T2 = ((a,b)N1,c,d,(e,f)N2);
//t = T1-T2;
//fprintf (stdout, t);
//fprintf (stdout, "\n");
//fprintf (stdout, T1);
//fprintf (stdout, "\n");
//fprintf (stdout, T2);
//fprintf (stdout, "\n");


//Tree T1 = ((a,b)N1,c,d,((g,h)N3,e,f)N2);
//t = PSTreeString (T1,"STRING_SUPPLIED_LENGTHS",{{-1,-1}});
//fprintf (stdout, t);
//fprintf (stdout, "\n");
//Topology T2 = ((a,b)N1,c,d,(e,f)N2);

Tree T1 = ((a,b)N1,(c,d)N4,((g,h)N3,e,f)N2);
t = T1^0;
fprintf (stdout, t);
fprintf (stdout, "\n");
