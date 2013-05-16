class Mock_CELInternals : public _CELInternals {
};

class Mock_ExecutionList : public _ExecutionList {
 public:
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
};
