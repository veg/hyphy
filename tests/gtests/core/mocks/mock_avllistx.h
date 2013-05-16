class Mock_AVLListX : public _AVLListX {
 public:
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(Clear,
      void());
  MOCK_METHOD1(DeleteXtra,
      void(long));
  MOCK_METHOD1(PopulateFromList,
      void(_List));
  MOCK_METHOD3(InsertData,
      long(BaseRef, long, bool));
  MOCK_METHOD3(UpdateValue,
      long(BaseRef, long, long));
};
