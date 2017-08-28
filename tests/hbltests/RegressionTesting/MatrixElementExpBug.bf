a = {{0,1,2,3}};
ea = a ["Exp (_MATRIX_ELEMENT_VALUE_)"];

for (i = 0; i < Columns (a); i+=1) {
	assert (ea [i] == Exp (a[i]), "Exp(" + a[i] + ") != "+ ea[i]);
}
