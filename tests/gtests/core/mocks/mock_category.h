class Mock_CategoryVariable : public _CategoryVariable {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD0(IsGlobal,
      bool(void));
  MOCK_METHOD0(IsConstant,
      bool(void));
  MOCK_METHOD0(IsCategory,
      bool(void));
  MOCK_METHOD4(ScanForVariables,
      void(, , _AVLListX, long));
  MOCK_METHOD1(ScanForGVariables,
      void(_AVLList));
};
