class Mock_SimpleList : public _SimpleList {
 public:
  MOCK_METHOD1(=,
      operator _SimpleList(_SimpleList));
  MOCK_METHOD1(&,
      operator _SimpleList(_SimpleList));
  MOCK_METHOD1(<<,
      void operator(long));
  MOCK_METHOD1(>>,
      bool operator(long));
  MOCK_METHOD1(<<,
      void operator(_SimpleList));
  MOCK_METHOD2(BinaryFind,
      long(, long));
  MOCK_METHOD2(Compare,
      long(long, long));
  MOCK_METHOD2(Compare,
      long(BaseRef, long));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD1(DeleteList,
      void(const _SimpleList));
  MOCK_METHOD2(FilterRange,
      void(long, long));
  MOCK_CONST_METHOD2(Find,
      long(const, long));
  MOCK_CONST_METHOD3(FindStepping,
      long(const, const, const));
  MOCK_METHOD1(Initialize,
      void());
  MOCK_METHOD4(InsertElement,
      void(BaseRef, long, bool, bool));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD0(toStr,
      BaseRef(void));
};
