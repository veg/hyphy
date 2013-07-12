class Mock_CString : public _CString {
 public:
  MOCK_METHOD1(<<,
      void operator(_String));
  MOCK_METHOD1(<<,
      void operator(char));
  MOCK_METHOD0(Finalize,
      void(void));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(FreeUpMemory,
      long(long));
  MOCK_METHOD1(Duplicate,
      void(BaseRef ref));
};
