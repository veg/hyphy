class Mock_Polynomial : public _Polynomial {
 public:
  MOCK_METHOD4(Execute,
      _MathObject*(long, _MathObject, _MathObject, _hyExecutionContext));
  MOCK_METHOD0(makeDynamic,
      BaseObj*(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD1(Add,
      _MathObject*(_MathObject));
  MOCK_METHOD2(Plus,
      _MathObject*(, bool));
  MOCK_METHOD1(Sub,
      _MathObject*(_MathObject));
  MOCK_METHOD1(Raise,
      _MathObject*(_MathObject));
  MOCK_METHOD0(Minus,
      _MathObject*(void));
  MOCK_METHOD1(Mult,
      _MathObject*(_MathObject));
  MOCK_METHOD0(Compute,
      _MathObject*(void));
  MOCK_METHOD1(Equal,
      bool(_MathObject));
  MOCK_METHOD0(IsObjectEmpty,
      bool(void));
  MOCK_METHOD0(ObjectClass,
      unsigned long(void));
  MOCK_METHOD0(Value,
      _Parameter(void));
  MOCK_METHOD0(toStr,
      BaseObj*(void));
  MOCK_METHOD1(toFileStr,
      void(FILE));
  MOCK_METHOD4(ScanForVariables,
      void(_AVLList, bool, _AVLListX, long));
  MOCK_METHOD0(HasChanged,
      bool(void));
};
