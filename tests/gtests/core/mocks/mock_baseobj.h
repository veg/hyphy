class MockBaseObj : public BaseObj {
 public:
  MOCK_METHOD0(toStr,
      BaseObj*(void));
  MOCK_METHOD0(toErrStr,
      BaseObj*(void));
  MOCK_METHOD1(toFileStr,
      void(FILE));
  MOCK_METHOD0(makeDynamic,
      BaseObj*(void));
  MOCK_METHOD1(FreeUpMemory,
      long(long));
  MOCK_METHOD0(Initialize,
      void(void));
  MOCK_METHOD1(Duplicate,
      void(BaseObj *ref));
  MOCK_METHOD0(AddAReference,
      void(void));
  MOCK_METHOD0(RemoveAReference,
      void(void));
};
