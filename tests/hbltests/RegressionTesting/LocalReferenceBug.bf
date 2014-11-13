lfunction foo (x) {
    fprintf (stdout, ^x, "\n");
}

lfunction bar () {
    y = "Hai!";
    foo (&y);
}

bar ();

x = "Canadian devil!";

fprintf (stdout, "*" , &x, " = ", x, " = ", *(&x), "\n");