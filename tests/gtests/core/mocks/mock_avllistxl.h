class Mock_AVLListXL : public _AVLListXL {
 public:
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(GetDataByKey,
      BaseRef(BaseRef));
  MOCK_METHOD3(InsertData,
      long(BaseRef, long, bool));
  MOCK_METHOD1(Clear,
      void());
  MOCK_METHOD1(DeleteXtra,
      void(long));
  MOCK_METHOD4(UpdateValue,
      long(, , , ));
};
