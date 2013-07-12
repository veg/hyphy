class Mock_DataSetFilter : public _DataSetFilter {
 public:
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(toFileStr,
      void(FILE));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(FreeUpMemory,
      long(long));
  MOCK_METHOD0(IsNormalFilter,
      bool(void));
  MOCK_METHOD0(GetFullLengthSpecies,
      long(void));
  MOCK_METHOD0(GetSiteCount,
      long(void));
  MOCK_METHOD1(GetDimension,
      long(bool));
  MOCK_METHOD2(operator,
      _String&(unsigned long site, unsigned long pos));
  MOCK_METHOD4(RetrieveState,
      void(long, long, , ));
  MOCK_METHOD2(GetChar,
      char(unsigned long site, unsigned long pos));
  MOCK_METHOD3(CompareTwoSites,
      bool(unsigned long, unsigned long, unsigned long));
};
