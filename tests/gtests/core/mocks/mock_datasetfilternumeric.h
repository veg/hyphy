class Mock_DataSetFilterNumeric : public _DataSetFilterNumeric {
 public:
  MOCK_METHOD0(IsNormalFilter,
      bool(void));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(GetDimension,
      long(bool));
  MOCK_METHOD3(CompareTwoSites,
      bool(unsigned long, unsigned long, unsigned long));
};
