class Mock_Operation : public _Operation {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseObj*(void));
  MOCK_METHOD1(StackDepth,
      void(long));
  MOCK_METHOD0(toStr,
      BaseObj*(void));
  MOCK_METHOD0(Initialize,
      void(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD1(IsAVariable,
      bool());
  MOCK_METHOD0(IsConstant,
      bool(void));
  MOCK_METHOD0(IsAFunctionCall,
      bool(void));
  MOCK_METHOD0(UserFunctionID,
      long(void));
  MOCK_METHOD0(GetAVariable,
      long(void));
  MOCK_METHOD1(SetAVariable,
      void(long d));
  MOCK_METHOD0(AssignmentVariable,
      bool(void));
  MOCK_METHOD0(HasChanged,
      bool(void));
  MOCK_METHOD1(SetTerms,
      void(long d));
  MOCK_METHOD0(GetANumber,
      _PMathObj(void));
  MOCK_METHOD1(SetNumber,
      void(_PMathObj d));
  MOCK_METHOD1(EqualOp,
      bool(_Operation));
};
