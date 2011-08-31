function			checkMarginalReconstruction (expected_partials, expected_classes, lfID)
{
	ConstructCategoryMatrix (cmx, "`lfID`");
	cmx = Transpose(cmx);
	cmx = cmx["_MATRIX_ELEMENT_VALUE_*Exp(cmx.log_scale_multiplier)^cmx.site_scalers[_MATRIX_ELEMENT_ROW_]"];
	ConstructCategoryMatrix (cmx2, lfID, SHORT);
	ConstructCategoryMatrix (cmx3, lfID, WEIGHTS);
	matrixError1 = Max(expected_partials - cmx, 0);
	matrixError2 = Abs(expected_classes - cmx2);
	return {{matrixError1__, matrixError2__}};
}
